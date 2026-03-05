package main

import (
	"strings"
	"testing"
)

// Tests simplified for native Go SMILES pipeline.
// These SMILES strings focus on non-stereo chains, rings, branches and common atom types.
func TestSMILESGen3D(t *testing.T) {
	smilesStrings := []string{
		"C",                  // Methane
		"CC",                 // Ethane
		"C=C",                // Ethene
		"C#C",                // Ethyne
		"CCO",                // Ethanol
		"CCN(C)C",            // Dimethylethylamine
		"C1CCCCC1",           // Cyclohexane
		"c1ccccc1",           // Benzene
		"C1=CC=CC=C1",        // Benzene (Kekule)
		"c1ccncc1",           // Pyridine
		"C1COCCO1",           // 1,4-Dioxane
		"CCC(C)CC(O)C",       // Branched chain
		"C1CC2CCC1C2",        // Norbornane
		"CC(=O)O",            // Acetic acid
		"N#CC#N",             // Dicyanogen
		"C12C3C4C1C5C4C3C25", // Cubane (simplified)
		"ClC(F)(Br)I",        // Halogens
		"[nH]1cccc1",         // Pyrrole (explicit H)
		"[O-]C(=O)C",         // Acetate (explicit charge)
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
