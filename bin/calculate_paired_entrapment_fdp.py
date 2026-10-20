#!/usr/bin/env python3
"""Assess the paired novel-target/entrapment peptide list (r=1).

Ports the paired estimator from calculate_paired_entrapment_fdp.py and the
complete-tie target/decoy recalculation from recalculate_benchmark_global_peptidoform.py.
FDP = (N_E + N_E_only + 2*N_E_better_both)/(N_T + N_E), equation 4 of
10.1038/s41592-025-02719-x. The precursor cutoff does not establish novel-only FDR.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


def score(frame, rank, output):
    if not frame[rank].between(0, 1).all() or not frame.Decoy.isin([0, 1]).all():
        raise ValueError('Missing/invalid ranking values or decoy labels')
    keys = ['Precursor.Id', 'Decoy']
    if frame.groupby(keys)['Stripped.Sequence'].nunique().gt(1).any():
        raise ValueError('Inconsistent sequence mapping')
    entities = frame.sort_values(rank, kind='stable').drop_duplicates(keys).copy()
    curve = entities.groupby(rank, sort=True).Decoy.agg(['sum', 'count'])
    curve['D'] = curve['sum'].cumsum()
    curve['T'] = (curve['count'] - curve['sum']).cumsum()
    curve['FDR'] = (curve['D'] / curve['T'].replace(0, np.nan)).fillna(1).clip(upper=1)
    curve[output] = curve.FDR.iloc[::-1].cummin().iloc[::-1]
    assert curve[output].is_monotonic_increasing
    entities[output] = entities[rank].map(curve[output])
    lookup = entities.set_index(keys)[output]
    out = frame.join(lookup, on=keys, validate='many_to_one')
    assert len(out) == len(frame)
    return out, entities, curve.reset_index()


def count_fdp(pairs, cutoff):
    result = pairs.copy()
    target = result.target_q.le(cutoff)
    entrapment = result.entrapment_q.le(cutoff)
    alone = entrapment & ~target
    target_rank = result.get("target_rank", result.target_q)
    entrapment_rank = result.get("entrapment_rank", result.entrapment_q)
    wins = entrapment & target & entrapment_rank.lt(target_rank)
    ties = entrapment & target & entrapment_rank.eq(target_rank)
    if ties.any():
        raise ValueError("Accepted partners have tied ranking values; resolve score ordering before estimating FDP.")
    result["target_accepted"] = target
    result["entrapment_accepted"] = entrapment
    result["entrapment_only"] = alone
    result["entrapment_better_both_accepted"] = wins
    result["paired_numerator_contribution"] = (
        entrapment.astype(int) + alone.astype(int) + 2 * wins.astype(int)
    )
    nt, ne = int(target.sum()), int(entrapment.sum())
    numerator = int(result.paired_numerator_contribution.sum())
    summary = {
        "q_cutoff": cutoff,
        "target_peptides": nt,
        "entrapment_peptides": ne,
        "entrapment_only": int(alone.sum()),
        "entrapment_better_both_accepted": int(wins.sum()),
        "accepted_pair_ties": int(ties.sum()),
        "paired_numerator": numerator,
        "combined_discoveries": nt + ne,
        "paired_fdp": numerator / (nt + ne) if nt + ne else np.nan,
    }
    return result, summary


def calculate(report, pair_file, q_column='Global.Q.Value', cutoff=0.001,
              recalculate_q=False, rank_column=None):
    if q_column != 'Global.Q.Value':
        raise ValueError('Paired FDP requires the Global.Q.Value column')
    if not 0 < cutoff < 1:
        raise ValueError('Cutoff must be between 0 and 1')
    rank_column = rank_column or q_column
    columns = list(dict.fromkeys(['Run', 'Protein.Ids', 'Stripped.Sequence', 'Decoy',
                                  q_column, rank_column] + (['Precursor.Id'] if recalculate_q else [])))
    parquet = pq.ParquetFile(report)
    missing = set(columns) - set(parquet.schema_arrow.names)
    if missing:
        raise ValueError(f'Report is missing required columns: {sorted(missing)}')
    frames, rows, runs = [], 0, set()
    for batch in parquet.iter_batches(columns=columns):
        frame = batch.to_pandas()
        rows += len(frame)
        runs.update(frame.Run.unique())
        if not frame.Decoy.isin([0, 1]).all():
            raise ValueError('Missing/invalid decoy labels')
        if not recalculate_q:
            frame = frame.loc[frame.Decoy.eq(0) & frame['Protein.Ids'].str.contains('NOVEL_TARGET_', na=False)]
        frames.append(frame)
    frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=columns)
    selected_q = q_column
    if recalculate_q and not frame.empty:
        selected_q = 'Recalculated.Global.TD.Q.Value'
        frame, _, _ = score(frame, q_column, selected_q)
    observed = frame.loc[frame.Decoy.eq(0) & frame['Protein.Ids'].str.contains('NOVEL_TARGET_', na=False)].copy()
    if not observed[q_column].between(0, 1).all() or not np.isfinite(observed[rank_column]).all():
        raise ValueError('Invalid q-values or pairing rank scores')
    if observed['Protein.Ids'].str.contains(';').any():
        raise ValueError('Ambiguous novel/entrapment accessions')
    observed['accession'] = observed['Protein.Ids']
    observed_ids = set(observed.accession)
    base_ids = {value[:-9] if value.endswith('_p_target') else value for value in observed_ids}
    expected_ids = base_ids | {value + '_p_target' for value in base_ids}
    selected = []
    for chunk in pd.read_csv(pair_file, sep='\t', chunksize=50000):
        required = {'sequence', 'decoy', 'proteins', 'peptide_type'}
        if not required.issubset(chunk.columns):
            raise ValueError(f'Pair table requires columns: {sorted(required)}')
        chunk['accession'] = chunk.proteins.str.split('|').str[1]
        selected.append(chunk.loc[chunk.accession.isin(expected_ids)].copy())
    mapping = pd.concat(selected, ignore_index=True)
    if set(mapping.accession) != expected_ids or mapping.accession.duplicated().any():
        raise ValueError('Pair mapping is missing partners or contains duplicate IDs')
    if not mapping.decoy.eq('No').all():
        raise ValueError('Entrapments must be search targets, not decoys')
    mapping['pair_id'] = mapping.accession.str.replace(r'_p_target$', '', regex=True)
    mapping['kind'] = np.where(mapping.accession.str.endswith('_p_target'), 'entrapment', 'target')
    if not mapping.peptide_type.eq(mapping.kind.map({'target': 'target', 'entrapment': 'p_target'})).all():
        raise ValueError('Pair roles disagree with accession suffixes')
    if mapping.sequence.str.replace('I', 'L', regex=False).duplicated().any():
        raise ValueError('Nonunique peptide identities in evaluated pairs')
    joined = observed.merge(mapping[['accession', 'sequence']], on='accession', validate='many_to_one')
    if len(joined) != len(observed) or not joined['Stripped.Sequence'].eq(joined.sequence).all():
        raise ValueError('Report sequences disagree with the pair mapping')
    mapping['q'] = mapping.accession.map(observed.groupby('accession')[selected_q].min())
    mapping['rank'] = mapping.accession.map(observed.groupby('accession')[rank_column].min())
    if mapping.empty:
        pairs = pd.DataFrame(columns=['pair_id', 'target_sequence', 'entrapment_sequence',
                                     'target_q', 'entrapment_q', 'target_rank', 'entrapment_rank'])
    else:
        pairs = mapping.pivot(index='pair_id', columns='kind', values=['sequence', 'q', 'rank'])
        pairs.columns = [f'{kind}_{value}' for value, kind in pairs.columns]
        pairs = pairs.reset_index()
        if not pairs.target_sequence.map(sorted).eq(pairs.entrapment_sequence.map(sorted)).all():
            raise ValueError('Paired sequences do not preserve amino-acid composition')
    result, summary = count_fdp(pairs, cutoff)
    summary.update(
        input_report=str(report.resolve()), pair_file=str(pair_file.resolve()),
        q_column=selected_q, ranking_column=q_column, pair_ranking_column=rank_column,
        recalculated_q_values=recalculate_q, input_rows=rows, input_runs=len(runs),
        evaluated_nondecoy_rows=len(observed), validated_observed_peptides=len(observed_ids),
        pairs_with_observed_member=len(pairs), entrapment_ratio=1,
        aggregation='minimum q-value per unmodified peptide across charges, modifications and runs',
        scope='novel targets plus paired shuffled entrapments; excludes reference and decoy entries',
        method_doi='10.1038/s41592-025-02719-x',
        equation='(N_E + N_E_only + 2*N_E_better_both_accepted)/(N_T + N_E)',
        target_decoy_method=('best score per (Precursor.Id, Decoy); complete ascending score ties; '
                             'cumulative D/T capped at 1; reverse cumulative minimum; all decoys retained'
                             if recalculate_q else 'use supplied q-values'),
    )
    return result, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--pairs', type=Path, required=True)
    parser.add_argument('--output-prefix', type=Path, required=True)
    parser.add_argument('--q-column', choices=['Global.Q.Value'], default='Global.Q.Value')
    parser.add_argument('--rank-column', help='Lower-is-better score for ordering accepted partners')
    parser.add_argument('--cutoff', type=float, default=0.001)
    parser.add_argument('--recalculate-q-values', action='store_true')
    args = parser.parse_args()
    paths = [Path(str(args.output_prefix) + suffix)
             for suffix in ('_summary.tsv', '_method.json')]
    for path in paths:
        if path.exists():
            raise FileExistsError(f'Refusing to overwrite {path}')
        if not path.parent.is_dir():
            raise FileNotFoundError(path.parent)
    _, summary = calculate(args.report, args.pairs, args.q_column, args.cutoff,
                           args.recalculate_q_values, args.rank_column)
    pd.DataFrame([summary]).to_csv(paths[0], sep='\t', index=False, mode='x')
    with paths[1].open('x') as handle:
        json.dump({key: (None if isinstance(value, float) and np.isnan(value) else value)
                   for key, value in summary.items()}, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(f"Paired FDP = {summary['paired_numerator']}/{summary['combined_discoveries']} "
          f"= {summary['paired_fdp']:.6%}")


if __name__ == '__main__':
    main()
