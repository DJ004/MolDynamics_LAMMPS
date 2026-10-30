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
