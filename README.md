## MolDynamics_LAMMPS

### `lammps_data_to_pdb.py`
- Successful conversion does not guarantee that the PSF atom indices match the LAMMPS atom IDs. If they do not, the coordinates can look correct while residue and atom identities are assigned incorrectly.
- Recommend loading the original PSF as the topology and the LAMMPS-derived coordinates separately in VMD (as below). This avoids the PDB naming and numbering limitations and preserves the full molecular topology. 
  - VMD TkConsole:
    - mol new `PSF_FILE` type psf waitfor all
    - mol addfile `PDB_FILE` type pdb waitfor all
  - VMD GUI:
    - File → New Molecule → Load `PSF_FILE`
    - Select the same molecule in the Molecule File Browser.
    - Click Browse → select `PDB_FILE` → Load

### `lammpstrj_to_xyz.py`
- Converts a generiz LAMMPS trajectory (.lammpstrj extension) into an XYZ trajectory.
- Requires MDAnalysis to be installed in the env.

### `lammps_data_to_xyz.py`
- Converts a generic LAMMPS data file to VMD-readable XYZ coordinate file.
  - Searches for line containing "Atoms", then ignores the following empty line and extracts coordinates from the following line. Details are provided as comments in the script.
- Assumes by default that the last 3 columns of coordinate-containing rows are X, Y, Z coordinates. User can set custom columns for X, Y, Z in the user settings of the script.
- Assumes all atom names as "A" - this is the first column of the XYZ file.