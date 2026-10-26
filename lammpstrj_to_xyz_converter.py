# Need MDA installed in conda env with python=3.12

import MDAnalysis as mda

# 1. Load using .lammpstrj file 
# (format="LAMMPSDUMP" tells Python how to read it despite the extension)
u = mda.Universe("top.data", "traj.lammpstrj", format="LAMMPSDUMP")

# 2. Select atom range ("all" to convert full lammpstrj to xyz)
atom_selection = "all"
subset = u.select_atoms(atom_selection)

# 3. Save selection into .xyz file
with mda.Writer("traj.xyz", subset.n_atoms) as W:
    for ts in u.trajectory:
        W.write(subset)
