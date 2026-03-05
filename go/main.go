package main

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
)

// Generate3DXYZFromSMILES takes a SMILES string and generates a 3D XYZ representation
// using the Open Babel executable.
// It reproduces the behavior of `obabel -ismi -oxyz --gen3d`
func Generate3DXYZFromSMILES(smiles string) (string, error) {
	// Prepare the command
	// We use echo to pipe the SMILES string into obabel
	// Or we can pass it directly. obabel allows format -:<SMILES>
	// However, standard way is passing as a file or stdin

	cmd := exec.Command("obabel", "-ismi", "-oxyz", "--gen3d")

	// Provide the SMILES string via stdin
	cmd.Stdin = bytes.NewBufferString(smiles)

	var outb, errb bytes.Buffer
	cmd.Stdout = &outb
	cmd.Stderr = &errb

	err := cmd.Run()
	if err != nil {
		return "", fmt.Errorf("obabel failed: %v, stderr: %s", err, errb.String())
	}

	return outb.String(), nil
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
