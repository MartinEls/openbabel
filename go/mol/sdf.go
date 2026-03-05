package mol

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// ReadSDF reads a multi-molecule SDF file and returns a slice of *Mol.
func ReadSDF(filename string) ([]*Mol, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var mols []*Mol
	scanner := bufio.NewScanner(file)

	for {
		m, done, err := readSDFMol(scanner)
		if err != nil {
			return nil, err
		}
		if m != nil {
			mols = append(mols, m)
		}
		if done {
			break
		}
	}

	return mols, nil
}

func readSDFMol(scanner *bufio.Scanner) (*Mol, bool, error) {
	// Read title
	if !scanner.Scan() {
		return nil, true, nil // EOF
	}
	title := strings.TrimSpace(scanner.Text())
	_ = title

	// Skip 2 header lines
	for i := 0; i < 2; i++ {
		if !scanner.Scan() {
			return nil, true, nil
		}
	}

	// Read counts line
	if !scanner.Scan() {
		return nil, true, nil
	}
	countsLine := scanner.Text()
	if len(countsLine) < 6 {
		return nil, false, fmt.Errorf("invalid counts line: %s", countsLine)
	}

	numAtomsStr := strings.TrimSpace(countsLine[0:3])
	numBondsStr := strings.TrimSpace(countsLine[3:6])

	numAtoms, err := strconv.Atoi(numAtomsStr)
	if err != nil {
		return nil, false, fmt.Errorf("invalid atom count: %s", numAtomsStr)
	}

	numBonds, err := strconv.Atoi(numBondsStr)
	if err != nil {
		return nil, false, fmt.Errorf("invalid bond count: %s", numBondsStr)
	}

	m := NewMol()

	// Read atoms
	for i := 0; i < numAtoms; i++ {
		if !scanner.Scan() {
			return nil, false, fmt.Errorf("unexpected EOF reading atoms")
		}
		line := scanner.Text()
		if len(line) < 39 {
			return nil, false, fmt.Errorf("invalid atom block line: %s", line)
		}

		x, _ := strconv.ParseFloat(strings.TrimSpace(line[0:10]), 64)
		y, _ := strconv.ParseFloat(strings.TrimSpace(line[10:20]), 64)
		z, _ := strconv.ParseFloat(strings.TrimSpace(line[20:30]), 64)
		symbol := strings.TrimSpace(line[31:34])

		atom := m.AddAtom(GetAtomicNum(symbol))
		atom.X = x
		atom.Y = y
		atom.Z = z
	}

	// Read bonds
	for i := 0; i < numBonds; i++ {
		if !scanner.Scan() {
			return nil, false, fmt.Errorf("unexpected EOF reading bonds")
		}
		line := scanner.Text()
		if len(line) < 9 {
			return nil, false, fmt.Errorf("invalid bond block line: %s", line)
		}

		idx1, _ := strconv.Atoi(strings.TrimSpace(line[0:3]))
		idx2, _ := strconv.Atoi(strings.TrimSpace(line[3:6]))
		order, _ := strconv.Atoi(strings.TrimSpace(line[6:9]))
		isAromatic := order == 4
		if isAromatic {
			order = 1
		}

		m.AddBond(idx1, idx2, order, isAromatic)
	}

	// Read until $$$$
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "$$$$" {
			break
		}
	}

	return m, false, nil
}
