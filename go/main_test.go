package main

import (
	"strings"
	"testing"
)

// Tests derived from `test/smilestest.cpp` in the OpenBabel repository.
// These SMILES strings involve various stereochemistry scenarios, chiral centers, and complex rings.
// Using genericSmilesCanonicalTest strings and explicit test cases:
func TestSMILESGen3D(t *testing.T) {
	smilesStrings := []string{
		"C[C@H](O)N",
		"Cl[C@@](CCl)(I)Br",
		"Cl/C=C/F",
		"F[Po@SP1](Cl)(Br)I",
		"F[Po@SP2](Br)(Cl)I",
		"F[Po@SP3](Cl)(I)Br",
		"CCC[C@@H](O)CC\\C=C\\C=C\\C#CC#C\\C=C\\CO",
		"OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1",
		"OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2",
		"CC(=O)OCCC(/C)=C\\C[C@H](C(C)=C)CCC=C",
		"CC[C@H](O1)CC[C@@]12CCCO2",
		"CN1CCC[C@H]1c2cccnc2",
		"CC(C)[C@@]12C[C@@H]1[C@@H](C)C(=O)C2",
		"CC(C)[C@H]1CC[C@]([C@@H]2[C@@H]1C=C(COC2=O)C(=O)O)(CCl)O",
		"C(CS[14CH2][14C@@H]1[14C@H]([14C@H]([14CH](O1)O)O)O)[C@@H](C(=O)O)N",
		"CCC[C@@H]1C[C@H](N(C1)C)C(=O)NC([C@@H]2[C@@H]([C@@H]([C@H]([C@H](O2)SC)OP(=O)(O)O)O)O)C(C)Cl",
	}

	for _, smiles := range smilesStrings {
		t.Run("SMILES: "+smiles, func(t *testing.T) {
			output, err := Generate3DXYZFromSMILES(smiles)
			if err != nil {
				t.Fatalf("Failed to generate 3D coordinates for SMILES %q: %v", smiles, err)
			}

			// Verify output has basic XYZ format characteristics
			lines := strings.Split(strings.TrimSpace(output), "\n")
			if len(lines) < 3 {
				t.Errorf("Expected at least 3 lines for an XYZ file (number of atoms, comment, and atoms), got %d lines", len(lines))
			}

			// Further validation of XYZ format could be done here.
		})
	}
}
