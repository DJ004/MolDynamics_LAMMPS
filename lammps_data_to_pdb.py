#!/usr/bin/env python3
"""
Convert CHARMM-GUI generated LAMMPS atom_style full DATA coordinates to PDB using CHARMM PSF labels.

Usage:
    python lammpsdata_to_pdb_v4.py frame.data topology.psf
    python lammpsdata_to_pdb_v4.py frame.data topology.psf -o frame.pdb

Requirements: LAMMPS atom IDs and PSF atom indices must refer to the same atoms.
Coordinates are copied without scaling or unwrapping (assumed Angstroms).
"""

import argparse
from pathlib import Path

CHAIN_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
SECTIONS = {
    "Masses", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals",
    "Impropers", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
    "Dihedral Coeffs", "Improper Coeffs", "PairIJ Coeffs",
    "BondBond Coeffs", "BondAngle Coeffs", "MiddleBondTorsion Coeffs",
    "EndBondTorsion Coeffs", "AngleTorsion Coeffs", "AngleAngleTorsion Coeffs",
    "BondBond13 Coeffs", "AngleAngle Coeffs", "Ellipsoids", "Lines",
    "Triangles", "Bodies", "Molecules",
}


def read_psf(path):
    atoms = {}
    with open(path, encoding="utf-8") as file:
        for line in file:
            if "!NATOM" in line:
                count = int(line.split()[0])
                break
        else:
            raise ValueError("PSF does not contain a !NATOM section")
        for _ in range(count):
            parts = next(file).split()
            if len(parts) < 8:
                raise ValueError(f"Malformed PSF atom record: {parts}")
            atom_id = int(parts[0])
            if atom_id in atoms:
                raise ValueError(f"Duplicate PSF atom index: {atom_id}")
            # Standard PSF: index segid resid resname atomname atomtype charge mass
            atoms[atom_id] = (parts[1], parts[2], parts[3], parts[4])
    return atoms


def read_lammps_full(path):
    atoms = {}
    in_atoms = False
    with open(path, encoding="utf-8") as file:
        for line in file:
            content = line.split("#", 1)[0].strip()
            if not content:
                continue
            if not in_atoms:
                if content == "Atoms":
                    in_atoms = True
                continue
            if content in SECTIONS:
                break
            parts = content.split()
            if len(parts) < 7:
                raise ValueError(f"Bad Atoms record: {line.rstrip()}")
            try:
                atom_id = int(parts[0])
                molecule_id = int(parts[1])
                atom_type = int(parts[2])
                x, y, z = map(float, parts[4:7])
            except ValueError as exc:
                raise ValueError(f"Bad Atoms record: {line.rstrip()}") from exc
            if atom_id in atoms:
                raise ValueError(f"Duplicate LAMMPS atom ID: {atom_id}")
            atoms[atom_id] = (molecule_id, atom_type, x, y, z)
    if not in_atoms or not atoms:
        raise ValueError("No atoms found in the LAMMPS Atoms section")
    return atoms


def chain_ids(psf_atoms):
    segments = list(dict.fromkeys(item[0] for item in psf_atoms.values()))
    if len(segments) > len(CHAIN_ALPHABET):
        raise ValueError("More than 62 PSF segments; PDB chain IDs are insufficient")
    return {segment: CHAIN_ALPHABET[i] for i, segment in enumerate(segments)}


def write_pdb(path, lammps_atoms, psf_atoms):
    missing = set(lammps_atoms) - set(psf_atoms)
    extra = set(psf_atoms) - set(lammps_atoms)
    if missing or extra:
        raise ValueError(
            f"PSF/LAMMPS atom IDs do not match: {len(missing)} missing from PSF, "
            f"{len(extra)} missing from DATA. Examples: "
            f"missing={sorted(missing)[:5]}, extra={sorted(extra)[:5]}. "
            "Use a PSF built for exactly this atom-ID mapping."
        )
    chains = chain_ids(psf_atoms)
    overflow_residues = 0
    with open(path, "w", encoding="ascii") as output:
        for atom_id in sorted(lammps_atoms):
            _, _, x, y, z = lammps_atoms[atom_id]
            seg, resid, resname, atomname = psf_atoms[atom_id]
            # Standard PDB atom and residue names fit in 4 and 3 columns.
            # Keep valid atom names such as OH2; truncate TIP3 to TIP.
            atomname = atomname[:4]
            resname = resname[:3]
            try:
                residue_number = int(resid)
            except ValueError as exc:
                raise ValueError(f"Non-integer PSF residue ID {resid!r}; PDB output requires a numeric ID") from exc
            if residue_number > 9999 or residue_number < 0:
                overflow_residues += 1
            # Repeats after 99,999 atoms.
            serial = (atom_id - 1) % 99999 + 1
            pdb_resid = residue_number % 10000
            # PDB fixed columns 13-16: atom name; 18-20: resname;
            # 22: chain; 23-26: resseq; 31-54: XYZ.
            name_field = f" {atomname:<3}" if len(atomname) < 4 else atomname
            record = (
                f"ATOM  {serial:5d} {name_field:4s} {resname:>3s} "
                f"{chains[seg]}{pdb_resid:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{1.00:6.2f}{0.00:6.2f}          "
                f"  "  # Unknown element: PSF atom type is NOT an element
                f"  "  # Blank charge
            )
            if len(record) != 80 or any(len(f"{v:8.3f}") > 8 for v in (x, y, z)):
                raise ValueError(
                    f"Atom {atom_id}: coordinate cannot fit in PDB 8.3 fields; "
                    "use an extended coordinate format instead."
                )
            output.write(record + "\n")
        output.write("END\n")
    return chains, overflow_residues


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_file", type=Path, help="LAMMPS DATA file (atom_style full)")
    parser.add_argument("psf_file", type=Path, help="Matching CHARMM PSF topology")
    parser.add_argument("-o", "--output", type=Path, help="Output PDB (default: DATA stem + .pdb)")
    args = parser.parse_args()
    output = args.output or Path(args.data_file.stem + ".pdb")
    topology = read_psf(args.psf_file)
    atoms = read_lammps_full(args.data_file)
    chains, overflow_residues = write_pdb(output, atoms, topology)
    print(f"Wrote {len(atoms):,} atoms to {output}")
    print("PSF segment -> PDB chain:", ", ".join(f"{s}:{c}" for s, c in chains.items()))
    if len(atoms) > 99999 or max(atoms) > 99999:
        print("WARNING: PDB serial numbers repeat after 99,999")
    if overflow_residues:
        print(f"WARNING: {overflow_residues:,} atoms have residue IDs outside 0-9999; IDs wrapped.")
    print("NOTE: Atom IDs are sorted, XYZ coordinates are unchanged, and element fields are blank.")
    print("NOTE: PSF atom order/IDs must correspond exactly to LAMMPS atom IDs.")


if __name__ == "__main__":
    main()
