# Observations on SMILES parsing and 3D generation

1. The `--gen3d` option is handled by `OpGen3D::Do` in `src/ops/gen3d.cpp`.
2. The core logic for building 3D structures from 0D/2D resides in `OBBuilder` (`src/builder.cpp`) and distance geometry (`src/distgeom.cpp`).
3. SMILES strings are read by `OBSmilesParser::ParseSmiles` in `src/formats/smilesformat.cpp` or via `src/formats/smiley.h` depending on the exact format used (`smi` vs `smiley`).
4. The generation of 3D coordinates involves several steps:
   a. Perceiving stereo if missing (`StereoFrom0D`).
   b. Setup `OBGen3DStereoHelper`.
   c. Use `OBBuilder::Build` to generate initial 3D structure.
   d. If `OBBuilder` fails or if distance geometry is requested (like with `--gen3d dg`), fallback to `OBDistanceGeometry::GetGeometry`.
   e. Perform Force Field (FF) cleanup with MMFF94 or UFF.
   f. Various combinations of rotor search and conjugate gradients are applied depending on speed requested (1-5).
   g. Final step involves `stereoHelper.Check` to ensure correctness.

5. The XYZ file generation logic is handled by `src/formats/xyzformat.cpp`, where `WriteMolecule` writes out the number of atoms, an optional title, and then each atom's symbol and its x, y, and z coordinates.

6. The generation process heavily relies on Open Babel specific internal classes (`OBMol`, `OBAtom`, `OBBond`, `OBBuilder`, `OBForceField`, `OBDistanceGeometry`), making an exact 1:1 rewrite in Go non-trivial without a substantial cheminformatics library or writing CGo bindings to Open Babel itself.

### Relevant Test Cases for SMILES string parsing
From `test/smilestest.cpp`:

1. `C[C@H](O)N`
2. `Cl[C@@](CCl)(I)Br`
3. `Cl/C=C/F`
4. `F[Po@SP1](Cl)(Br)I`
5. `F[Po@SP2](Br)(Cl)I`
6. `F[Po@SP3](Cl)(I)Br`
7. `CCC[C@@H](O)CC\C=C\C=C\C#CC#C\C=C\CO`
8. `OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1`
9. `OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2`
10. `CC(=O)OCCC(/C)=C\C[C@H](C(C)=C)CCC=C`
11. `CC[C@H](O1)CC[C@@]12CCCO2`
12. `CN1CCC[C@H]1c2cccnc2`
13. `CC(C)[C@@]12C[C@@H]1[C@@H](C)C(=O)C2`
14. `CC(C)[C@H]1CC[C@]([C@@H]2[C@@H]1C=C(COC2=O)C(=O)O)(CCl)O`
15. `C(CS[14CH2][14C@@H]1[14C@H]([14C@H]([14CH](O1)O)O)O)[C@@H](C(=O)O)N`
16. `CCC[C@@H]1C[C@H](N(C1)C)C(=O)NC([C@@H]2[C@@H]([C@@H]([C@H]([C@H](O2)SC)OP(=O)(O)O)O)O)C(C)Cl`

Other tests that might be relevant for 3D coordinate generation and XYZ format:
- `test/testdistgeom.py` tests distance geometry coordinate generation via `-osdf --gen3d dg`.
- `test/testsym.py` runs `obabel -ismi -oxyz --gen3d` and tests symmetry.
