#!/usr/bin/env python3
"""Generate a small deterministic DIA-MS mzML file for the PG-DIA example.

The output is synthetic and is intended only for workflow smoke testing. It is
not intended for benchmarking identification sensitivity or quantitative
accuracy.
"""

from __future__ import annotations

import argparse
import base64
import math
import random
import struct
import zlib
from pathlib import Path
from xml.etree import ElementTree as ET


NS = "http://psi.hupo.org/ms/mzml"
XSI = "http://www.w3.org/2001/XMLSchema-instance"
ET.register_namespace("", NS)
ET.register_namespace("xsi", XSI)

PROTON = 1.007276466621
WATER = 18.010564684
CARBAMIDOMETHYL = 57.021463735
AMINO_ACID_MASS = {
    "A": 71.037113805,
    "R": 156.101111050,
    "N": 114.042927470,
    "D": 115.026943065,
    "C": 103.009184505,
    "E": 129.042593135,
    "Q": 128.058577540,
    "G": 57.021463735,
    "H": 137.058911875,
    "I": 113.084063975,
    "L": 113.084063975,
    "K": 128.094963015,
    "M": 131.040484645,
    "F": 147.068413945,
    "P": 97.052763875,
    "S": 87.032028435,
    "T": 101.047678505,
    "W": 186.079312980,
    "Y": 163.063328575,
    "V": 99.068413945,
}

# Sequence, precursor charge, retention-time apex (minutes), peak height.
PEPTIDES = (
    ("VLGHENR", 2, 3.4, 850_000.0),
    ("PEPTIDER", 2, 4.6, 1_100_000.0),
    ("SDFGHLEK", 2, 5.7, 760_000.0),
    ("VAVVQAEK", 2, 6.8, 930_000.0),
    ("FQGQVLDSLR", 2, 8.0, 1_250_000.0),
    ("LSEPAELTDAVK", 2, 9.3, 980_000.0),
    ("LGEYGFQNALIVR", 2, 10.5, 1_300_000.0),
    ("AGLCQTFVYGGCR", 2, 11.7, 880_000.0),
    ("TMTDIDGDGTITTK", 2, 13.0, 1_020_000.0),
    ("YLQEIYNSNNQK", 2, 14.2, 1_180_000.0),
)

DIA_WINDOWS = ((400.0, 500.0), (500.0, 600.0), (600.0, 700.0), (700.0, 800.0))


def qname(tag: str) -> str:
    return f"{{{NS}}}{tag}"


def add_cv(parent: ET.Element, accession: str, name: str, value: object | None = None,
           cv_ref: str = "MS", unit: tuple[str, str, str] | None = None) -> ET.Element:
    attributes = {"cvRef": cv_ref, "accession": accession, "name": name}
    if value is not None:
        attributes["value"] = str(value)
    if unit is not None:
        unit_cv_ref, unit_accession, unit_name = unit
        attributes.update(
            unitCvRef=unit_cv_ref,
            unitAccession=unit_accession,
            unitName=unit_name,
        )
    return ET.SubElement(parent, qname("cvParam"), attributes)


def peptide_mass(sequence: str) -> float:
    return WATER + sum(
        AMINO_ACID_MASS[aa] + (CARBAMIDOMETHYL if aa == "C" else 0.0)
        for aa in sequence
    )


def precursor_mz(sequence: str, charge: int) -> float:
    return (peptide_mass(sequence) + charge * PROTON) / charge


def fragment_ions(sequence: str) -> list[float]:
    residue_masses = [
        AMINO_ACID_MASS[aa] + (CARBAMIDOMETHYL if aa == "C" else 0.0)
        for aa in sequence
    ]
    ions: list[float] = []
    for cut in range(1, len(sequence)):
        b_ion = sum(residue_masses[:cut]) + PROTON
        y_ion = sum(residue_masses[cut:]) + WATER + PROTON
        if 300.0 <= b_ion <= 1500.0:
            ions.append(b_ion)
        if 300.0 <= y_ion <= 1500.0:
            ions.append(y_ion)
    return sorted(set(round(mz, 6) for mz in ions))


def encode(values: list[float]) -> str:
    packed = struct.pack(f"<{len(values)}d", *values)
    return base64.b64encode(zlib.compress(packed)).decode("ascii")


def add_binary_array(parent: ET.Element, values: list[float], array_type: str) -> None:
    encoded = encode(values)
    binary_array = ET.SubElement(
        parent,
        qname("binaryDataArray"),
        {"encodedLength": str(len(encoded))},
    )
    add_cv(binary_array, "MS:1000523", "64-bit float")
    add_cv(binary_array, "MS:1000574", "zlib compression")
    if array_type == "mz":
        add_cv(
            binary_array,
            "MS:1000514",
            "m/z array",
            unit=("MS", "MS:1000040", "m/z"),
        )
    else:
        add_cv(
            binary_array,
            "MS:1000515",
            "intensity array",
            unit=("MS", "MS:1000131", "number of detector counts"),
        )
    ET.SubElement(binary_array, qname("binary")).text = encoded


def add_spectrum(
    spectrum_list: ET.Element,
    index: int,
    scan_number: int,
    rt: float,
    ms_level: int,
    mz_values: list[float],
    intensities: list[float],
    preceding_ms1_id: str | None = None,
    isolation_window: tuple[float, float] | None = None,
) -> None:
    spectrum_id = f"controllerType=0 controllerNumber=1 scan={scan_number}"
    spectrum = ET.SubElement(
        spectrum_list,
        qname("spectrum"),
        {
            "index": str(index),
            "id": spectrum_id,
            "defaultArrayLength": str(len(mz_values)),
        },
    )
    add_cv(spectrum, "MS:1000511", "ms level", ms_level)
    add_cv(spectrum, "MS:1000127", "centroid spectrum")
    add_cv(spectrum, "MS:1000130", "positive scan")
    add_cv(spectrum, "MS:1000285", "total ion current", f"{sum(intensities):.4f}")
    if intensities:
        apex_index = max(range(len(intensities)), key=intensities.__getitem__)
        add_cv(spectrum, "MS:1000504", "base peak m/z", f"{mz_values[apex_index]:.6f}")
        add_cv(spectrum, "MS:1000505", "base peak intensity", f"{intensities[apex_index]:.4f}")
        add_cv(spectrum, "MS:1000528", "lowest observed m/z", f"{min(mz_values):.6f}")
        add_cv(spectrum, "MS:1000527", "highest observed m/z", f"{max(mz_values):.6f}")

    scan_list = ET.SubElement(spectrum, qname("scanList"), {"count": "1"})
    add_cv(scan_list, "MS:1000795", "no combination")
    scan = ET.SubElement(scan_list, qname("scan"), {"instrumentConfigurationRef": "IC1"})
    add_cv(
        scan,
        "MS:1000016",
        "scan start time",
        f"{rt:.4f}",
        unit=("UO", "UO:0000031", "minute"),
    )
    scan_window_list = ET.SubElement(scan, qname("scanWindowList"), {"count": "1"})
    scan_window = ET.SubElement(scan_window_list, qname("scanWindow"))
    lower_limit, upper_limit = (350.0, 1200.0) if ms_level == 1 else (300.0, 1500.0)
    add_cv(scan_window, "MS:1000501", "scan window lower limit", lower_limit,
           unit=("MS", "MS:1000040", "m/z"))
    add_cv(scan_window, "MS:1000500", "scan window upper limit", upper_limit,
           unit=("MS", "MS:1000040", "m/z"))

    if ms_level == 2 and isolation_window is not None and preceding_ms1_id is not None:
        lower, upper = isolation_window
        center = (lower + upper) / 2.0
        precursor_list = ET.SubElement(spectrum, qname("precursorList"), {"count": "1"})
        precursor = ET.SubElement(
            precursor_list,
            qname("precursor"),
            {"spectrumRef": preceding_ms1_id},
        )
        isolation = ET.SubElement(precursor, qname("isolationWindow"))
        add_cv(isolation, "MS:1000827", "isolation window target m/z", center,
               unit=("MS", "MS:1000040", "m/z"))
        add_cv(isolation, "MS:1000828", "isolation window lower offset", center - lower,
               unit=("MS", "MS:1000040", "m/z"))
        add_cv(isolation, "MS:1000829", "isolation window upper offset", upper - center,
               unit=("MS", "MS:1000040", "m/z"))
        activation = ET.SubElement(precursor, qname("activation"))
        add_cv(activation, "MS:1000422", "beam-type collision-induced dissociation")
        add_cv(activation, "MS:1000045", "collision energy", 30,
               unit=("UO", "UO:0000266", "electronvolt"))

    array_list = ET.SubElement(spectrum, qname("binaryDataArrayList"), {"count": "2"})
    add_binary_array(array_list, mz_values, "mz")
    add_binary_array(array_list, intensities, "intensity")


def make_ms1(rt: float, rng: random.Random) -> tuple[list[float], list[float]]:
    peaks: list[tuple[float, float]] = []
    for sequence, charge, apex, height in PEPTIDES:
        mz = precursor_mz(sequence, charge)
        if not 350.0 <= mz <= 1200.0:
            continue
        profile = math.exp(-0.5 * ((rt - apex) / 0.34) ** 2)
        intensity = height * profile
        if intensity > 50.0:
            peaks.extend(
                [
                    (mz, intensity),
                    (mz + 1.003355 / charge, intensity * 0.46),
                    (mz + 2.006710 / charge, intensity * 0.18),
                ]
            )
    for _ in range(20):
        peaks.append((rng.uniform(350.0, 1200.0), rng.uniform(40.0, 350.0)))
    peaks.sort()
    return [peak[0] for peak in peaks], [peak[1] for peak in peaks]


def make_ms2(rt: float, window: tuple[float, float], rng: random.Random) -> tuple[list[float], list[float]]:
    lower, upper = window
    peaks: list[tuple[float, float]] = []
    for sequence, charge, apex, height in PEPTIDES:
        mz = precursor_mz(sequence, charge)
        if not lower <= mz < upper:
            continue
        profile = math.exp(-0.5 * ((rt - apex) / 0.34) ** 2)
        fragments = fragment_ions(sequence)
        for rank, fragment_mz in enumerate(fragments):
            relative = 0.75 / (1.0 + 0.10 * rank)
            intensity = height * profile * relative
            if intensity > 25.0:
                peaks.append((fragment_mz, intensity))
    for _ in range(24):
        peaks.append((rng.uniform(300.0, 1500.0), rng.uniform(25.0, 220.0)))
    peaks.sort()
    return [peak[0] for peak in peaks], [peak[1] for peak in peaks]


def build_mzml(output_path: Path) -> None:
    root = ET.Element(
        qname("mzML"),
        {
            "id": "PGDIA_SIMULATED_DIA",
            "version": "1.1.0",
            f"{{{XSI}}}schemaLocation": (
                "http://psi.hupo.org/ms/mzml "
                "http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd"
            ),
        },
    )
    cv_list = ET.SubElement(root, qname("cvList"), {"count": "2"})
    ET.SubElement(cv_list, qname("cv"), {
        "id": "MS", "fullName": "Proteomics Standards Initiative Mass Spectrometry Ontology",
        "version": "4.1.0", "URI": "https://github.com/HUPO-PSI/psi-ms-CV/raw/master/psi-ms.obo",
    })
    ET.SubElement(cv_list, qname("cv"), {
        "id": "UO", "fullName": "Unit Ontology", "version": "2023-05-25",
        "URI": "http://purl.obolibrary.org/obo/uo.obo",
    })

    file_description = ET.SubElement(root, qname("fileDescription"))
    file_content = ET.SubElement(file_description, qname("fileContent"))
    add_cv(file_content, "MS:1000579", "MS1 spectrum")
    add_cv(file_content, "MS:1000580", "MSn spectrum")
    add_cv(file_content, "MS:1000127", "centroid spectrum")
    add_cv(file_content, "MS:1003215", "data-independent acquisition")
    ET.SubElement(file_description, qname("sourceFileList"), {"count": "0"})

    software_list = ET.SubElement(root, qname("softwareList"), {"count": "1"})
    software = ET.SubElement(software_list, qname("software"), {"id": "PGDIA_SIM", "version": "1.0"})
    add_cv(software, "MS:1000799", "custom unreleased software tool")

    instrument_list = ET.SubElement(root, qname("instrumentConfigurationList"), {"count": "1"})
    instrument = ET.SubElement(instrument_list, qname("instrumentConfiguration"), {"id": "IC1"})
    component_list = ET.SubElement(instrument, qname("componentList"), {"count": "3"})
    source = ET.SubElement(component_list, qname("source"), {"order": "1"})
    add_cv(source, "MS:1000073", "electrospray ionization")
    analyzer = ET.SubElement(component_list, qname("analyzer"), {"order": "2"})
    add_cv(analyzer, "MS:1000484", "orbitrap")
    detector = ET.SubElement(component_list, qname("detector"), {"order": "3"})
    add_cv(detector, "MS:1000624", "inductive detector")

    data_processing_list = ET.SubElement(root, qname("dataProcessingList"), {"count": "1"})
    data_processing = ET.SubElement(data_processing_list, qname("dataProcessing"), {"id": "DP1"})
    method = ET.SubElement(data_processing, qname("processingMethod"), {
        "order": "0", "softwareRef": "PGDIA_SIM",
    })
    add_cv(method, "MS:1000544", "Conversion to mzML")

    run = ET.SubElement(root, qname("run"), {
        "id": "simulated_dia_run", "defaultInstrumentConfigurationRef": "IC1",
    })
    cycle_times = [2.0 + 0.2 * index for index in range(71)]
    spectrum_count = len(cycle_times) * (1 + len(DIA_WINDOWS))
    spectrum_list = ET.SubElement(run, qname("spectrumList"), {
        "count": str(spectrum_count), "defaultDataProcessingRef": "DP1",
    })

    rng = random.Random(314159)
    spectrum_index = 0
    scan_number = 1
    for rt in cycle_times:
        mz_values, intensities = make_ms1(rt, rng)
        ms1_id = f"controllerType=0 controllerNumber=1 scan={scan_number}"
        add_spectrum(
            spectrum_list, spectrum_index, scan_number, rt, 1, mz_values, intensities,
        )
        spectrum_index += 1
        scan_number += 1
        for window in DIA_WINDOWS:
            mz_values, intensities = make_ms2(rt + 0.02, window, rng)
            add_spectrum(
                spectrum_list,
                spectrum_index,
                scan_number,
                rt + 0.02,
                2,
                mz_values,
                intensities,
                preceding_ms1_id=ms1_id,
                isolation_window=window,
            )
            spectrum_index += 1
            scan_number += 1

    ET.SubElement(run, qname("chromatogramList"), {
        "count": "0", "defaultDataProcessingRef": "DP1",
    })

    ET.indent(root, space="  ")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(root).write(output_path, encoding="utf-8", xml_declaration=True)


def write_ground_truth(output_dir: Path) -> None:
    with (output_dir / "simulated_dia_ground_truth.tsv").open("w", encoding="utf-8") as handle:
        handle.write("peptide\tcharge\tprecursor_mz\trt_apex_min\tpeak_height\n")
        for sequence, charge, apex, height in PEPTIDES:
            handle.write(
                f"{sequence}\t{charge}\t{precursor_mz(sequence, charge):.6f}\t"
                f"{apex:.2f}\t{height:.0f}\n"
            )

    with (output_dir / "simulated_proteins.fasta").open("w", encoding="utf-8") as handle:
        for index, (sequence, _, _, _) in enumerate(PEPTIDES, start=1):
            handle.write(
                f">SIM_P{index:03d}|SIM_TX{index:03d}|SIM_GENE{index:03d} "
                f"synthetic protein containing {sequence}\n"
            )
            handle.write(f"M{sequence}\n")


def parse_args() -> argparse.Namespace:
    default_output = Path(__file__).resolve().parent / "data" / "simulated_dia.mzML"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=default_output,
        help=f"Output mzML path (default: {default_output})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_path = args.output.resolve()
    build_mzml(output_path)
    write_ground_truth(output_path.parent)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
