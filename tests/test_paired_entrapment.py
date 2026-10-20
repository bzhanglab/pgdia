"""Small-data checks for paired database construction and FDP assessment."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
import build_paired_entrapment_db as builder
import calculate_paired_entrapment_fdp as assessment


class PairedEntrapmentTest(unittest.TestCase):
    def setUp(self):
        self.directory = Path(tempfile.mkdtemp(prefix='pgdia-paired-test-'))
        self.reference = self.directory / 'reference.fa'
        self.novel = self.directory / 'novel.fa'
        self.reference.write_text('>sp|REF|REF GN=GENE\nACDEFGHIK\n')
        self.novel.write_text('>var_one\nPEPTIDER\n>var_two\nMNPQSTVK\n'
                              '>var_repeat\nAAAAAAAK\n>reference_IL\nACDEFGHLK\n')
        self.prefix = self.directory / 'sample'
        self.manifest = builder.build(self.reference, self.novel, 'sample', self.prefix)
        self.pairs = pd.read_csv(self.directory / 'sample.pairs.tsv', sep='\t')
        self.pairs['accession'] = self.pairs.proteins.str.split('|').str[1]
        self.report = self.directory / 'report.parquet'

    def report_rows(self):
        rows = []
        # Pair 0: both accepted, entrapment ranks better. Pair 1: entrapment only.
        for index, row in self.pairs.iterrows():
            q = [.004, .003, .04, .002][index]
            rows.append({'Run': 'run1', 'Protein.Ids': row.accession,
                         'Stripped.Sequence': row.sequence, 'Decoy': 0,
                         'Global.Q.Value': q, 'Precursor.Id': row.sequence + '2'})
        return rows

    def test_database_constraints_and_reproducibility(self):
        self.assertEqual(self.manifest['nr_targets'], 2)
        self.assertEqual(self.manifest['entrapments'], 2)
        self.assertEqual(self.manifest['excluded_nr_targets'], 1)
        targets = self.pairs[self.pairs.peptide_type.eq('target')]
        entrapments = self.pairs[self.pairs.peptide_type.eq('p_target')]
        for target, entrapment in zip(targets.sequence, entrapments.sequence):
            self.assertEqual(sorted(target), sorted(entrapment))
            self.assertEqual(target[-1], entrapment[-1])
        self.assertFalse(builder.database_substring_hits(set(entrapments.sequence),
                                                         [self.reference, self.novel]))
        self.assertFalse(self.pairs.sequence.str.replace('I', 'L', regex=False).duplicated().any())
        self.assertNotIn('ACDEFGHLK', set(targets.sequence))
        builder.build(self.reference, self.novel, 'sample', self.directory / 'repeat')
        self.assertEqual((self.directory / 'sample.fasta').read_bytes(),
                         (self.directory / 'repeat.fasta').read_bytes())
        with self.assertRaises(FileExistsError):
            builder.build(self.reference, self.novel, 'sample', self.prefix)

    def test_full_protein_substring_exclusion(self):
        # The candidate occurs inside a longer protein, away from cleavage sites.
        self.reference.write_text('>ref\nAAAPEDITPERGGGK\n')
        self.assertEqual(builder.database_substring_hits({'PEDITPER'}, [self.reference]), {'PEDLTPER'})

    def test_complete_ranking_ties_and_precursor_deduplication(self):
        frame = pd.DataFrame({'Precursor.Id': list('abcde') + ['a'],
                              'Decoy': [0, 0, 1, 0, 0, 0],
                              'Stripped.Sequence': list('ABCDE') + ['A'],
                              'Global.Q.Value': [.01, .02, .02, .03, .04, .05]})
        scored, entities, curve = assessment.score(frame, 'Global.Q.Value', 'q')
        np.testing.assert_allclose(entities.q, [0, .25, .25, .25, .25])
        self.assertEqual(curve.D.tolist(), [0, 1, 1, 1])
        self.assertEqual(curve['T'].tolist(), [1, 2, 3, 4])
        self.assertEqual(scored.q.iloc[0], scored.q.iloc[-1])

    def test_paired_numerator_and_missing_partner(self):
        pd.DataFrame(self.report_rows()).to_parquet(self.report, index=False)
        audit, summary = assessment.calculate(self.report, self.directory / 'sample.pairs.tsv', cutoff=.01)
        self.assertEqual(summary['target_peptides'], 1)
        self.assertEqual(summary['entrapment_peptides'], 2)
        self.assertEqual(summary['entrapment_only'], 1)
        self.assertEqual(summary['entrapment_better_both_accepted'], 1)
        self.assertEqual(summary['paired_numerator'], 5)
        self.assertAlmostEqual(summary['paired_fdp'], 5 / 3)
        self.assertEqual(audit.paired_numerator_contribution.sum(), 5)
        # A missing report partner must remain an undiscovered member of the pair.
        pd.DataFrame(self.report_rows()[1:2]).to_parquet(self.directory / 'missing.parquet', index=False)
        _, missing = assessment.calculate(self.directory / 'missing.parquet',
                                          self.directory / 'sample.pairs.tsv', cutoff=.01)
        self.assertEqual(missing['paired_numerator'], 2)
        self.assertEqual(missing['combined_discoveries'], 1)

    def test_reference_and_decoy_retention_in_recalculation(self):
        rows = self.report_rows()
        rows += [dict(rows[0], **{'Protein.Ids': 'T_reference', 'Stripped.Sequence': 'ACDEFGHIK',
                                  'Precursor.Id': 'reference', 'Global.Q.Value': .001}),
                 dict(rows[0], **{'Protein.Ids': 'decoy', 'Decoy': 1,
                                  'Precursor.Id': 'decoy', 'Global.Q.Value': .003})]
        pd.DataFrame(rows).to_parquet(self.report, index=False)
        _, summary = assessment.calculate(self.report, self.directory / 'sample.pairs.tsv',
                                          cutoff=.01, recalculate_q=True)
        self.assertEqual(summary['target_peptides'], 0)
        self.assertEqual(summary['entrapment_peptides'], 1)
        self.assertEqual(summary['paired_numerator'], 2)
        self.assertEqual(summary['input_rows'], 6)

    def test_rejects_tied_accepted_partners(self):
        rows = self.report_rows()
        rows[0]['Global.Q.Value'] = rows[1]['Global.Q.Value']
        pd.DataFrame(rows).to_parquet(self.report, index=False)
        with self.assertRaisesRegex(ValueError, 'tied ranking'):
            assessment.calculate(self.report, self.directory / 'sample.pairs.tsv', cutoff=.01)

    def test_rejects_ambiguous_ids_and_sequence_mismatch(self):
        rows = self.report_rows()
        rows[0]['Protein.Ids'] = 'REF;' + rows[0]['Protein.Ids']
        pd.DataFrame(rows).to_parquet(self.report, index=False)
        with self.assertRaisesRegex(ValueError, 'Ambiguous'):
            assessment.calculate(self.report, self.directory / 'sample.pairs.tsv')
        rows = self.report_rows()
        rows[0]['Stripped.Sequence'] = 'INVALID'
        pd.DataFrame(rows).to_parquet(self.directory / 'mismatch.parquet', index=False)
        with self.assertRaisesRegex(ValueError, 'sequences disagree'):
            assessment.calculate(self.directory / 'mismatch.parquet', self.directory / 'sample.pairs.tsv')

    def test_empty_discoveries_and_cli_outputs(self):
        rows = [dict(self.report_rows()[0], **{'Protein.Ids': 'T_reference'})]
        pd.DataFrame(rows).to_parquet(self.report, index=False)
        script = Path(assessment.__file__)
        result = subprocess.run([sys.executable, str(script), '--report', str(self.report),
                                 '--pairs', str(self.directory / 'sample.pairs.tsv'),
                                 '--output-prefix', str(self.directory / 'empty')],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        method = json.loads((self.directory / 'empty_method.json').read_text())
        self.assertIsNone(method['paired_fdp'])
        self.assertEqual(method['combined_discoveries'], 0)
        self.assertEqual(sorted(path.name for path in self.directory.glob('empty*')),
                         ['empty_method.json', 'empty_summary.tsv'])


if __name__ == '__main__':
    unittest.main()
