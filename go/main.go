package main

import (
	"fmt"
	"os"

	"smiles3d/builder"
	"smiles3d/smiles"
	"smiles3d/uff"
	"smiles3d/xyz"
)

// Generate3DXYZFromSMILES takes a SMILES string and generates a 3D XYZ representation
// using the native Go implementation.
func Generate3DXYZFromSMILES(smilesStr string) (string, error) {
	// Parse SMILES
	mol, err := smiles.ParseSMILES(smilesStr)
	if err != nil {
		return "", fmt.Errorf("failed to parse SMILES: %v", err)
	}

	// Build Initial 3D Geometry
	builder.Build3D(mol)

	// Optimize Geometry using UFF
	// Using 200 iterations and a conservative step size
	uff.Optimize(mol, 200, 0.05)

	// Write to XYZ format
	return xyz.WriteXYZ(mol), nil
}

func main() {
	if len(os.Args) < 2 {
		fmt.Println("Usage: go run main.go <SMILES string>")
		os.Exit(1)
	}
	smilesStr := os.Args[1]

	fmt.Printf("Generating 3D coordinates for SMILES: %s\n", smilesStr)

	xyzOutput, err := Generate3DXYZFromSMILES(smilesStr)
	if err != nil {
		fmt.Printf("Error: %v\n", err)
		os.Exit(1)
	}

	fmt.Println("Result:")
	fmt.Println(xyzOutput)

	err = os.WriteFile("output.xyz", []byte(xyzOutput), 0644)
	if err != nil {
		fmt.Printf("Error writing to file: %v\n", err)
		os.Exit(1)
	}
	fmt.Println("Wrote 3D representation to output.xyz")
}
