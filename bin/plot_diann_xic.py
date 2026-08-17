#!/usr/bin/env python3
"""Plot precursor transition XICs exported by DIA-NN with ``--xic``."""

from __future__ import annotations

import argparse
import math
import re
import textwrap
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd
import pyarrow.parquet as pq


REQUIRED_COLUMNS = {"pr", "feature", "info", "rt", "value"}
FRAGMENT_PATTERN = re.compile(r"^([by])(\d+)(?:\^\d+)?$", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot transition-level extracted-ion chromatograms from the "
            "*.xic.parquet files written by DIA-NN with --xic."
        )
    )
    parser.add_argument("--xic-dir", required=True, type=Path,
                        help="Directory containing DIA-NN *.xic.parquet files.")
    parser.add_argument("--xic-glob", default="*.xic.parquet",
                        help="Input filename glob within --xic-dir (default: *.xic.parquet).")
    parser.add_argument("--precursor",
                        help="Exact DIA-NN precursor identifier from the pr column.")
    parser.add_argument("--peptide",
                        help="Unmodified peptide sequence; use with --charge if --precursor is omitted.")
    parser.add_argument("--charge", type=int,
                        help="Precursor charge; required with --peptide when --precursor is omitted.")
    parser.add_argument(
        "--variant-position", type=int,
        help=("Optional 1-based variant position within --peptide. Fragment ions "
              "containing this residue are marked with an asterisk."),
    )
    parser.add_argument("--transition", action="append", default=[],
                        help="Plot only this feature (for example y7^1); may be supplied repeatedly.")
    parser.add_argument(
        "--top-fragments", type=int, default=0,
        help="Keep the N fragment ions with greatest integrated intensity (default: all).",
    )
    parser.add_argument("--include-ms1", action="store_true",
                        help="Include the DIA-NN ms1 trace in addition to fragment-ion traces.")
    parser.add_argument(
        "--normalize", action="store_true",
        help="Normalize each run to its maximum plotted intensity; exported data remain raw.",
    )
    parser.add_argument("--format", choices=("png", "pdf", "both"), default="both",
                        help="Plot format (default: both).")
    parser.add_argument("--output-dir", type=Path, default=Path("diann_xic_plots"),
                        help="Output directory (default: diann_xic_plots).")
    args = parser.parse_args()

    if args.precursor is None:
        if args.peptide is None or args.charge is None:
            parser.error("provide --precursor or both --peptide and --charge")
        args.precursor = f"{args.peptide}{args.charge}"
    if args.peptide is None:
        args.peptide = infer_sequence(args.precursor)
    if args.charge is not None and not 1 <= args.charge <= 9:
        parser.error("--charge must be between 1 and 9")
    if args.variant_position is not None and not 1 <= args.variant_position <= len(args.peptide):
        parser.error("--variant-position must fall within the unmodified peptide sequence")
    if args.top_fragments < 0:
        parser.error("--top-fragments cannot be negative")
    return args


def infer_sequence(precursor: str) -> str:
    without_charge = re.sub(r"\d+$", "", precursor)
    without_modifications = re.sub(r"\([^)]*\)|\[[^]]*]", "", without_charge)
    return "".join(re.findall(r"[A-Z]", without_modifications.upper()))


def run_name(path: Path) -> str:
    suffix = ".xic.parquet"
    return path.name[: -len(suffix)] if path.name.endswith(suffix) else path.stem


def read_precursor(path: Path, precursor: str) -> pd.DataFrame:
    parquet = pq.ParquetFile(path)
    columns = set(parquet.schema_arrow.names)
    missing = REQUIRED_COLUMNS - columns
    if missing:
        raise ValueError(f"{path} is missing DIA-NN XIC columns: {', '.join(sorted(missing))}")
    table = pq.read_table(path, columns=sorted(REQUIRED_COLUMNS),
                          filters=[("pr", "=", precursor)])
    frame = table.to_pandas()
    if not frame.empty:
        frame.insert(0, "run", run_name(path))
    return frame


def fragment_contains_site(feature: str, peptide_length: int,
                           position: Optional[int]) -> bool:
    if position is None:
        return False
    match = FRAGMENT_PATTERN.match(str(feature))
    if match is None:
        return False
    ion_type, ion_number_text = match.groups()
    ion_number = int(ion_number_text)
    if ion_type.lower() == "b":
        return ion_number >= position
    return ion_number >= peptide_length - position + 1


def select_transitions(frame: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    selected = frame.loc[frame["feature"].astype(str).str.lower() != "index"].copy()
    if not args.include_ms1:
        selected = selected.loc[selected["feature"].astype(str).str.lower() != "ms1"]
    if args.transition:
        selected = selected.loc[selected["feature"].astype(str).isin(args.transition)]
    if args.top_fragments:
        fragments = selected.loc[selected["feature"].astype(str).str.lower() != "ms1"]
        keep = fragments.groupby("feature", observed=True)["value"].sum().nlargest(
            args.top_fragments
        ).index.tolist()
        if args.include_ms1:
            keep.append("ms1")
        selected = selected.loc[selected["feature"].isin(keep)]
    if selected.empty:
        raise ValueError("no transition rows remain after applying the requested filters")
    selected["contains_variant_site"] = selected["feature"].map(
        lambda feature: fragment_contains_site(feature, len(args.peptide),
                                                args.variant_position)
    )
    return selected.sort_values(["run", "feature", "rt"]).reset_index(drop=True)


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")


def transition_label(feature: str, contains_site: bool) -> str:
    return f"{feature} *" if contains_site else str(feature)


def wrapped_run_name(run: object) -> str:
    return "\n".join(textwrap.wrap(str(run), width=48))


def plot_run(axis: plt.Axes, frame: pd.DataFrame, colors: dict[str, object],
             normalize: bool) -> float:
    scale = float(frame["value"].max()) if normalize else 1.0
    if scale <= 0:
        scale = 1.0
    for feature, transition in frame.groupby("feature", sort=True, observed=True):
        transition = transition.sort_values("rt")
        contains_site = bool(transition["contains_variant_site"].iloc[0])
        axis.plot(
            transition["rt"], transition["value"] / scale,
            color=colors[str(feature)], linewidth=2.0 if contains_site else 1.2,
            label=transition_label(str(feature), contains_site),
        )
    summed = frame.groupby("rt", observed=True)["value"].sum()
    apex_rt = float(summed.idxmax())
    axis.axvline(apex_rt, color="0.35", linestyle="--", linewidth=0.9, alpha=0.65)
    axis.set_xlabel("Retention time (min)")
    axis.set_ylabel("Relative intensity" if normalize else "Intensity")
    return apex_rt


def save_figure(figure: plt.Figure, stem: Path, output_format: str) -> None:
    formats = ("png", "pdf") if output_format == "both" else (output_format,)
    for extension in formats:
        figure.savefig(Path(f"{stem}.{extension}"), dpi=300, bbox_inches="tight")


def plot_xics(frame: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    features = sorted(frame["feature"].astype(str).unique())
    palette = plt.get_cmap("tab20")
    colors = {feature: palette(index % 20) for index, feature in enumerate(features)}
    precursor_name = safe_name(args.precursor)
    summaries: list[dict[str, object]] = []

    for run, run_frame in frame.groupby("run", sort=True, observed=True):
        figure, axis = plt.subplots(figsize=(7.2, 4.6))
        apex_rt = plot_run(axis, run_frame, colors, args.normalize)
        axis.set_title(f"{args.precursor}\n{wrapped_run_name(run)}", fontsize=10)
        axis.legend(title="Transition", fontsize=8, ncols=2)
        figure.tight_layout()
        save_figure(figure,
                    args.output_dir / f"xic_{precursor_name}_{safe_name(str(run))}",
                    args.format)
        plt.close(figure)
        summaries.append({
            "run": run, "precursor": args.precursor, "apex_rt_min": apex_rt,
            "transition_count": run_frame["feature"].nunique(),
            "point_count": len(run_frame),
        })

    runs = list(frame["run"].drop_duplicates())
    column_count = min(3, len(runs))
    row_count = math.ceil(len(runs) / column_count)
    figure, axes = plt.subplots(row_count, column_count,
                                figsize=(6.1 * column_count, 4.2 * row_count),
                                squeeze=False)
    for axis, run in zip(axes.flat, runs):
        run_frame = frame.loc[frame["run"] == run]
        plot_run(axis, run_frame, colors, args.normalize)
        axis.set_title(wrapped_run_name(run), fontsize=9)
        axis.legend(fontsize=7, ncols=2)
    for axis in axes.flat[len(runs):]:
        axis.set_visible(False)
    figure.suptitle(f"DIA-NN transition XICs | {args.precursor}", y=1.01)
    figure.tight_layout()
    save_figure(figure, args.output_dir / f"xic_{precursor_name}_all_runs", args.format)
    plt.close(figure)
    return pd.DataFrame(summaries)


def main() -> None:
    args = parse_args()
    input_files = sorted(args.xic_dir.glob(args.xic_glob))
    if not input_files:
        raise SystemExit(f"no files matching {args.xic_glob!r} were found in {args.xic_dir}")
    frames = [read_precursor(path, args.precursor) for path in input_files]
    frames = [frame for frame in frames if not frame.empty]
    if not frames:
        raise SystemExit(
            f"precursor {args.precursor!r} was not found in {len(input_files)} XIC file(s)"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    selected = select_transitions(pd.concat(frames, ignore_index=True), args)
    summaries = plot_xics(selected, args)
    precursor_name = safe_name(args.precursor)
    selected.rename(columns={"pr": "precursor", "feature": "transition",
                             "rt": "rt_min", "value": "intensity"}).to_csv(
        args.output_dir / f"xic_{precursor_name}_transition_data.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    summaries.to_csv(args.output_dir / f"xic_{precursor_name}_summary.tsv",
                     sep="\t", index=False)
    print(
        f"Wrote XIC plots and filtered transition data for {args.precursor} "
        f"across {len(summaries)} run(s) to {args.output_dir}"
    )


if __name__ == "__main__":
    main()
