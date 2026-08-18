#!/usr/bin/env python3
"""Compare ExoCross cross sections with benchmark reference files.

The first supported benchmark type is absorption/Voigt-super with a
temperature-list and a pressure-list.  The implementation intentionally uses
only the Python standard library so it can be copied into ExoCross test jobs.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass, field
from itertools import zip_longest
from pathlib import Path
from typing import Iterable, Sequence


SUPPORTED_PROFILE = "voigt-super"


class BenchmarkError(Exception):
    """Raised for an invalid benchmark configuration or data file."""


@dataclass(frozen=True)
class BenchmarkSpec:
    input_path: Path
    task: str
    profile: str
    output: str
    temperatures: tuple[float, ...]
    pressures: tuple[float, ...]
    resolving_power: float | None


@dataclass
class FileResult:
    pressure: float
    reference: str
    candidate: str
    rows: int = 0
    values: int = 0
    max_absolute_error: float = 0.0
    max_relative_error: float = 0.0
    max_grid_absolute_error: float = 0.0
    max_grid_relative_error: float = 0.0
    mismatches: int = 0
    examples: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return self.mismatches == 0


def _strip_comment(line: str) -> str:
    """Remove common ExoCross comment forms from one input line."""
    for marker in ("!", "#", "("):
        if marker in line:
            line = line.split(marker, 1)[0]
    return line.strip()


def _number(text: str, *, context: str) -> float:
    try:
        value = float(text.replace("D", "E").replace("d", "e"))
    except ValueError as exc:
        raise BenchmarkError(f"invalid number {text!r} in {context}") from exc
    if not math.isfinite(value):
        raise BenchmarkError(f"non-finite number {text!r} in {context}")
    return value


def parse_input(path: Path) -> BenchmarkSpec:
    """Read the subset of an ExoCross input needed to select this benchmark."""
    path = path.resolve()
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise BenchmarkError(f"cannot read input file {path}: {exc}") from exc

    task = ""
    profile = ""
    output = ""
    temperatures: list[float] = []
    pressures: list[float] = []
    scalar_temperature: float | None = None
    scalar_pressure: float | None = None
    resolving_power: float | None = None
    active_list: str | None = None

    for line_number, raw_line in enumerate(lines, 1):
        line = _strip_comment(raw_line)
        if not line:
            continue
        fields = line.split()
        keyword = fields[0].lower()

        if active_list is not None:
            if keyword == "end":
                active_list = None
                continue
            value = _number(fields[0], context=f"{active_list}, line {line_number}")
            if active_list == "temperature-list":
                temperatures.append(value)
            else:
                pressures.append(value)
            continue

        if keyword in {"temperature-list", "pressure-list"}:
            active_list = keyword
        elif keyword in {"absorption", "emission"}:
            task = keyword
        elif keyword == SUPPORTED_PROFILE:
            profile = keyword
        elif keyword == "output":
            if len(fields) < 2:
                raise BenchmarkError(f"missing output name on line {line_number}")
            output = fields[1]
        elif keyword in {"temperature", "temp"}:
            if len(fields) < 2:
                raise BenchmarkError(f"missing temperature on line {line_number}")
            scalar_temperature = _number(fields[1], context=f"temperature, line {line_number}")
        elif keyword == "pressure":
            if len(fields) < 2:
                raise BenchmarkError(f"missing pressure on line {line_number}")
            scalar_pressure = _number(fields[1], context=f"pressure, line {line_number}")
        elif keyword == "r":
            if len(fields) < 2:
                raise BenchmarkError(f"missing resolving power on line {line_number}")
            resolving_power = _number(fields[1], context=f"R, line {line_number}")

    if active_list is not None:
        raise BenchmarkError(f"unterminated {active_list} block")
    if not temperatures and scalar_temperature is not None:
        temperatures.append(scalar_temperature)
    if not pressures and scalar_pressure is not None:
        pressures.append(scalar_pressure)

    missing = []
    if not output:
        missing.append("output")
    if not temperatures:
        missing.append("temperature-list (or temperature)")
    if not pressures:
        missing.append("pressure-list (or pressure)")
    if missing:
        raise BenchmarkError("missing required input: " + ", ".join(missing))
    if task != "absorption" or profile != SUPPORTED_PROFILE:
        shown_task = task or "(default/unspecified)"
        shown_profile = profile or "(unrecognised)"
        raise BenchmarkError(
            "unsupported ExoCross task/profile: "
            f"{shown_task}/{shown_profile}; currently supported: "
            f"absorption/{SUPPORTED_PROFILE}"
        )
    if any(value <= 0.0 for value in temperatures):
        raise BenchmarkError("temperatures must be positive")
    if any(value < 0.0 for value in pressures):
        raise BenchmarkError("pressures must be non-negative")

    return BenchmarkSpec(
        input_path=path,
        task=task,
        profile=profile,
        output=output,
        temperatures=tuple(temperatures),
        pressures=tuple(pressures),
        resolving_power=resolving_power,
    )


def fortran_pressure(value: float) -> str:
    """Format pressure like ExoCross's ES16.8-style filename suffix."""
    if not math.isfinite(value):
        raise BenchmarkError("pressure used in a filename must be finite")
    mantissa, exponent = f"{value:.8E}".split("E")
    return f"{mantissa}E{int(exponent):+04d}"


def output_name(spec: BenchmarkSpec, pressure: float) -> str:
    return f"{spec.output}__{fortran_pressure(pressure)}.xsec"


def significant_figure_rtol(digits: int) -> float:
    if digits < 1:
        raise BenchmarkError("significant figures must be at least 1")
    return 0.5 * 10.0 ** (1 - digits)


def _parse_row(line: str, path: Path, row_number: int) -> list[float]:
    values = []
    for column, token in enumerate(line.split(), 1):
        try:
            value = float(token.replace("D", "E").replace("d", "e"))
        except ValueError as exc:
            raise BenchmarkError(
                f"{path}, row {row_number}, column {column}: invalid number {token!r}"
            ) from exc
        if not math.isfinite(value):
            raise BenchmarkError(
                f"{path}, row {row_number}, column {column}: non-finite value {token!r}"
            )
        values.append(value)
    return values


def _close(
    expected: float,
    actual: float,
    *,
    rtol: float,
    atol: float,
    relative_floor: float,
) -> tuple[bool, float, float]:
    difference = abs(actual - expected)
    scale = max(abs(expected), abs(actual), relative_floor)
    relative = difference / scale
    return difference <= atol + rtol * scale, difference, relative


def _add_example(result: FileResult, message: str, limit: int) -> None:
    result.mismatches += 1
    if len(result.examples) < limit:
        result.examples.append(message)


def compare_xsec(
    reference: Path,
    candidate: Path,
    *,
    pressure: float,
    temperatures: Sequence[float],
    rtol: float,
    grid_rtol: float,
    atol: float,
    relative_floor: float,
    max_examples: int,
) -> FileResult:
    result = FileResult(
        pressure=pressure,
        reference=str(reference),
        candidate=str(candidate),
    )
    expected_columns = len(temperatures) + 1
    try:
        reference_stream = reference.open("r", encoding="utf-8")
    except OSError as exc:
        raise BenchmarkError(f"cannot read reference file {reference}: {exc}") from exc
    try:
        candidate_stream = candidate.open("r", encoding="utf-8")
    except OSError as exc:
        reference_stream.close()
        raise BenchmarkError(f"cannot read candidate file {candidate}: {exc}") from exc

    with reference_stream, candidate_stream:
        pairs: Iterable[tuple[str | None, str | None]] = zip_longest(
            reference_stream, candidate_stream
        )
        for row_number, (reference_line, candidate_line) in enumerate(pairs, 1):
            if reference_line is None:
                _add_example(result, f"row {row_number}: extra candidate row", max_examples)
                for _ in candidate_stream:
                    result.mismatches += 1
                break
            if candidate_line is None:
                _add_example(result, f"row {row_number}: candidate row is missing", max_examples)
                for _ in reference_stream:
                    result.mismatches += 1
                break

            result.rows += 1
            expected = _parse_row(reference_line, reference, row_number)
            actual = _parse_row(candidate_line, candidate, row_number)
            if len(expected) != expected_columns:
                raise BenchmarkError(
                    f"{reference}, row {row_number}: expected {expected_columns} columns "
                    f"from {len(temperatures)} temperatures, found {len(expected)}"
                )
            if len(actual) != expected_columns:
                _add_example(
                    result,
                    f"row {row_number}: candidate has {len(actual)} columns; "
                    f"expected {expected_columns}",
                    max_examples,
                )
                continue

            grid_ok, grid_abs, grid_rel = _close(
                expected[0],
                actual[0],
                rtol=grid_rtol,
                atol=0.0,
                relative_floor=1.0e-300,
            )
            result.max_grid_absolute_error = max(result.max_grid_absolute_error, grid_abs)
            result.max_grid_relative_error = max(result.max_grid_relative_error, grid_rel)
            if not grid_ok:
                _add_example(
                    result,
                    f"row {row_number}: grid expected {expected[0]:.9E}, "
                    f"got {actual[0]:.9E} (relative error {grid_rel:.3E})",
                    max_examples,
                )

            for column, temperature in enumerate(temperatures, 1):
                result.values += 1
                expected_value = expected[column]
                actual_value = actual[column]
                ok, absolute, relative = _close(
                    expected_value,
                    actual_value,
                    rtol=rtol,
                    atol=atol,
                    relative_floor=relative_floor,
                )
                result.max_absolute_error = max(result.max_absolute_error, absolute)
                result.max_relative_error = max(result.max_relative_error, relative)
                if expected_value < 0.0 or actual_value < 0.0:
                    ok = False
                if not ok:
                    _add_example(
                        result,
                        f"row {row_number}, T={temperature:g} K: "
                        f"expected {expected_value:.9E}, got {actual_value:.9E} "
                        f"(relative error {relative:.3E})",
                        max_examples,
                    )
    return result


def compare_benchmark(
    spec: BenchmarkSpec,
    reference_dir: Path,
    candidate_dir: Path,
    *,
    significant_figures: int,
    grid_significant_figures: int,
    absolute_tolerance: float,
    relative_floor: float,
    max_examples: int,
) -> list[FileResult]:
    if absolute_tolerance < 0.0:
        raise BenchmarkError("absolute tolerance cannot be negative")
    if relative_floor <= 0.0:
        raise BenchmarkError("relative floor must be positive")
    if max_examples < 0:
        raise BenchmarkError("max examples cannot be negative")

    reference_dir = reference_dir.resolve()
    candidate_dir = candidate_dir.resolve()
    if reference_dir == candidate_dir:
        raise BenchmarkError(
            "reference and candidate directories are the same; provide independently "
            "generated results with --candidate-dir"
        )
    rtol = significant_figure_rtol(significant_figures)
    grid_rtol = significant_figure_rtol(grid_significant_figures)
    results = []
    for pressure in spec.pressures:
        filename = output_name(spec, pressure)
        results.append(
            compare_xsec(
                reference_dir / filename,
                candidate_dir / filename,
                pressure=pressure,
                temperatures=spec.temperatures,
                rtol=rtol,
                grid_rtol=grid_rtol,
                atol=absolute_tolerance,
                relative_floor=relative_floor,
                max_examples=max_examples,
            )
        )
    return results


def _display_number(value: float) -> str:
    return f"{value:g}"


def _print_report(
    spec: BenchmarkSpec,
    results: Sequence[FileResult],
    *,
    significant_figures: int,
    grid_significant_figures: int,
    absolute_tolerance: float,
    relative_floor: float,
) -> None:
    print(f"Input: {spec.input_path}")
    print(f"Benchmark: {spec.task}/{spec.profile}")
    print("Temperatures (K): " + ", ".join(map(_display_number, spec.temperatures)))
    print("Pressures (bar): " + ", ".join(map(_display_number, spec.pressures)))
    if spec.resolving_power is not None:
        print(f"Resolving power: {spec.resolving_power:g}")
    print(
        f"Cross-section tolerance: {significant_figures} significant figures "
        f"(rtol={significant_figure_rtol(significant_figures):.3E}, "
        f"atol={absolute_tolerance:.3E}, relative floor={relative_floor:.3E})"
    )
    print(
        f"Grid tolerance: {grid_significant_figures} significant figures "
        f"(rtol={significant_figure_rtol(grid_significant_figures):.3E})"
    )
    print()

    for result in results:
        status = "PASS" if result.passed else "FAIL"
        print(
            f"{status}  P={result.pressure:g} bar  rows={result.rows:,}  "
            f"values={result.values:,}  max_rel={result.max_relative_error:.3E}  "
            f"mismatches={result.mismatches:,}"
        )
        for example in result.examples:
            print(f"      {example}")

    total_rows = sum(result.rows for result in results)
    total_values = sum(result.values for result in results)
    total_mismatches = sum(result.mismatches for result in results)
    status = "PASS" if total_mismatches == 0 else "FAIL"
    print()
    print(
        f"{status}: {len(results)} pressure file(s), {total_rows:,} rows, "
        f"{total_values:,} cross-section values, {total_mismatches:,} mismatch(es)"
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Read an ExoCross input file, select its benchmark, and compare "
            "candidate .xsec files with reference .xsec files."
        )
    )
    parser.add_argument("input", type=Path, help="ExoCross .inp file")
    parser.add_argument(
        "--candidate-dir",
        type=Path,
        required=True,
        help="directory containing newly generated ExoCross output files",
    )
    parser.add_argument(
        "--reference-dir",
        type=Path,
        help="directory containing reference outputs (default: input directory)",
    )
    parser.add_argument(
        "--significant-figures",
        type=int,
        default=6,
        help="required significant figures for cross sections (default: 6)",
    )
    parser.add_argument(
        "--grid-significant-figures",
        type=int,
        default=9,
        help="required significant figures for wavenumbers (default: 9)",
    )
    parser.add_argument(
        "--absolute-tolerance",
        type=float,
        default=0.0,
        help="optional absolute tolerance for cross sections (default: 0)",
    )
    parser.add_argument(
        "--relative-floor",
        type=float,
        default=1.0e-300,
        help="minimum scale used in relative error near zero (default: 1e-300)",
    )
    parser.add_argument(
        "--max-examples",
        type=int,
        default=10,
        help="maximum mismatch examples printed per file (default: 10)",
    )
    parser.add_argument("--json", action="store_true", help="print machine-readable JSON")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        spec = parse_input(args.input)
        reference_dir = args.reference_dir or spec.input_path.parent
        results = compare_benchmark(
            spec,
            reference_dir,
            args.candidate_dir,
            significant_figures=args.significant_figures,
            grid_significant_figures=args.grid_significant_figures,
            absolute_tolerance=args.absolute_tolerance,
            relative_floor=args.relative_floor,
            max_examples=args.max_examples,
        )
    except BenchmarkError as exc:
        print(f"benchmark error: {exc}", file=sys.stderr)
        return 2

    passed = all(result.passed for result in results)
    if args.json:
        report = {
            "passed": passed,
            "spec": {
                **asdict(spec),
                "input_path": str(spec.input_path),
            },
            "tolerance": {
                "significant_figures": args.significant_figures,
                "relative_tolerance": significant_figure_rtol(args.significant_figures),
                "grid_significant_figures": args.grid_significant_figures,
                "grid_relative_tolerance": significant_figure_rtol(
                    args.grid_significant_figures
                ),
                "absolute_tolerance": args.absolute_tolerance,
                "relative_floor": args.relative_floor,
            },
            "files": [asdict(result) | {"passed": result.passed} for result in results],
        }
        print(json.dumps(report, indent=2))
    else:
        _print_report(
            spec,
            results,
            significant_figures=args.significant_figures,
            grid_significant_figures=args.grid_significant_figures,
            absolute_tolerance=args.absolute_tolerance,
            relative_floor=args.relative_floor,
        )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
