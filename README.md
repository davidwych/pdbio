# pdbio
A minimalist package for reading, manipulating, and writing PDB files

Easilly incorporated in to a conda environemnt by adding to `lib/python3.#/site_packages/`

Main class is the `PDBFile` class, which can be instantialized with:
`>>> system = PDBFile(ifilename=input_pdb_file.pdb)`

Crystallographic information is stored in attribute: `system.crystinfo`
- e.g. unit cell box length `a` accessible through `system.crystinfo.a`

Contents of the PDB file, including `TER` records, stored in array `system.contents`
Atoms stored in array `system.atoms`: e.g. one can loop through all atoms in the system with:
```
all_glycines = []
for atom in system.atoms:
    if atom.resname == "GLY":
        all_glycines.append(atom)
```

PDBFile objects can be concatenated:
`>>> system1.combine(system2)`

Write new PDBFiles with:
`>>> system.write("output_pdb_file.pdb")`

Contains an `Atom` class which stores the elements of each `ATOM` and `HETATM` line
- index
- atomtype
- resname
- chainid
- resid
- x, y, and z
- occupancy
- bfactor
- element

Contains a `CrystInfo` class which stores the info from the `CRYST1` record
- Unit cell box lengths a, b, and c
- Unit cell angles alpha, beta, and gamma
- space group
- z value
