#!/usr/bin/env python3
"""Convert one generic LAMMPS data file to a VMD-readable XYZ file.
Usage:  python data2xyz_generic.py path/to/lammps.data
"""

import argparse
from pathlib import Path

# ========================== USER SETTINGS ==========================
# By default, X, Y, Z coordinates are the LAST THREE columns of an atom row.
# Set these to 1-based positive column numbers to select different columns.
# Negative numbers count from the end: -1 is the last column.
X_COLUMN = -3
Y_COLUMN = -2
Z_COLUMN = -1

# Generic XYZ needs an atom name in column 1. Without any chemical element
# information, use a fixed label ("A") for all atoms, which VMD can display.
ATOM_LABEL = "A"
# ===================================================================


def column_value(parts, column):
    """Fetch a token with 1-based positive or Python-like negative indexing."""
    if column == 0:
        raise ValueError("Column numbering cannot be zero")
    index = column - 1 if column > 0 else column
    return parts[index]


def extract_atoms(data_file):
    """Read only the first uninterrupted block of rows after 'Atoms'."""
    rows = []
    found_atoms = False
    started = False

    with open(data_file, encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not found_atoms:
                if line.split("#", 1)[0].strip().split()[:1] == ["Atoms"]:
                    found_atoms = True
                continue

            # Blank lines immediately after 'Atoms' are ignored; the NEXT
            # blank line after coordinates start terminates the Atoms section.
            if not line:
                if started:
                    break
                continue

            # A row may contain at most ONE contiguous alphabetic field
            # (such as an atom label). Further text is allowed after '#',
            # but '#' must follow at least TWO numeric columns.
            body, separator, comment = line.partition("#")
            parts = body.split()
            if not parts:
                raise ValueError(f"Line {line_number}: no atom data before '#' comment")
            if separator:
                numeric_before_comment = sum(is_number(part) for part in parts)
                if numeric_before_comment < 2:
                    raise ValueError(
                        f"Line {line_number}: '#' comment must follow at least two numeric columns"
                    )
            alphabetic_fields = sum(any(ch.isalpha() for ch in part) and not is_number(part) for part in parts)
            if alphabetic_fields > 1:
                raise ValueError(f"Line {line_number}: more than one alphabetic field")

            try:
                atom_id = int(parts[0])
                x = float(column_value(parts, X_COLUMN))
                y = float(column_value(parts, Y_COLUMN))
                z = float(column_value(parts, Z_COLUMN))
            except (ValueError, IndexError) as exc:
                raise ValueError(f"Line {line_number}: invalid atom ID or coordinate columns") from exc
            rows.append((atom_id, x, y, z))
            started = True

    if not found_atoms:
        raise ValueError("No 'Atoms' section found")
    if not rows:
        raise ValueError("No atom coordinate rows found after 'Atoms'")
    if len({row[0] for row in rows}) != len(rows):
        raise ValueError("Duplicate atom IDs found")
    rows.sort(key=lambda row: row[0])
    return rows


def is_number(token):
    try:
        float(token)
        return True
    except ValueError:
        return False


def main():
    parser = argparse.ArgumentParser(description="Convert generic LAMMPS DATA coordinates to XYZ")
    parser.add_argument("data_file", help="Path to input LAMMPS .data file")
    args = parser.parse_args()
    source = Path(args.data_file)
    if not source.is_file():
        parser.error(f"Input file not found: {source}")
    if not ATOM_LABEL or any(ch.isspace() for ch in ATOM_LABEL):
        parser.error("ATOM_LABEL must be one nonempty token")

    try:
        rows = extract_atoms(source)
        destination = source.with_suffix(".xyz")
        with destination.open("w", encoding="utf-8") as handle:
            handle.write(f"{len(rows)}\nFrame 0\n")
            for _, x, y, z in rows:
                handle.write(f"{ATOM_LABEL} {x:.15g} {y:.15g} {z:.15g}\n")
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Error: {exc}\n")
    print(f"Wrote {destination} ({len(rows)} atoms)")


if __name__ == "__main__":
    main()
