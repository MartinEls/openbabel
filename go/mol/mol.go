package mol

import (
	"fmt"
)

// Atom represents a chemical atom.
type Atom struct {
	Idx          int // Unique ID within the molecule
	AtomicNum    int // Z (e.g., 6 for Carbon)
	FormalCharge int
	IsAromatic   bool
	ImplicitH    int
	X, Y, Z      float64 // 3D Coordinates
}

func (a *Atom) Symbol() string {
	sym := GetSymbol(a.AtomicNum)
	if sym == "Xx" {
		return fmt.Sprintf("X%d", a.AtomicNum)
	}
	return sym
}

// Bond represents a chemical bond.
type Bond struct {
	Idx        int
	BeginAtom  *Atom
	EndAtom    *Atom
	Order      int // 1=Single, 2=Double, 3=Triple, 4=Aromatic
	IsAromatic bool
}

// Mol represents a chemical molecule.
type Mol struct {
	Atoms []*Atom
	Bonds []*Bond
}

// NewMol creates a new empty molecule.
func NewMol() *Mol {
	return &Mol{
		Atoms: make([]*Atom, 0),
		Bonds: make([]*Bond, 0),
	}
}

// AddAtom adds an atom to the molecule and returns it.
func (m *Mol) AddAtom(atomicNum int) *Atom {
	atom := &Atom{
		Idx:       len(m.Atoms) + 1,
		AtomicNum: atomicNum,
	}
	m.Atoms = append(m.Atoms, atom)
	return atom
}

// AddBond adds a bond between two atoms in the molecule.
func (m *Mol) AddBond(beginAtomIdx, endAtomIdx int, order int, isAromatic bool) *Bond {
	if beginAtomIdx < 1 || beginAtomIdx > len(m.Atoms) || endAtomIdx < 1 || endAtomIdx > len(m.Atoms) {
		return nil
	}

	// OpenBabel handles aromatic bonds generally as order 1 (or special flag).
	// We'll use Order 4 or IsAromatic flag.
	if isAromatic {
		order = 4
	}

	bond := &Bond{
		Idx:        len(m.Bonds) + 1,
		BeginAtom:  m.Atoms[beginAtomIdx-1],
		EndAtom:    m.Atoms[endAtomIdx-1],
		Order:      order,
		IsAromatic: isAromatic,
	}
	m.Bonds = append(m.Bonds, bond)
	return bond
}

// GetAtom returns an atom by its index (1-based).
func (m *Mol) GetAtom(idx int) *Atom {
	if idx < 1 || idx > len(m.Atoms) {
		return nil
	}
	return m.Atoms[idx-1]
}

// GetNeighbors returns all atoms bonded to the given atom.
func (m *Mol) GetNeighbors(a *Atom) []*Atom {
	var neighbors []*Atom
	for _, b := range m.Bonds {
		if b.BeginAtom.Idx == a.Idx {
			neighbors = append(neighbors, b.EndAtom)
		} else if b.EndAtom.Idx == a.Idx {
			neighbors = append(neighbors, b.BeginAtom)
		}
	}
	return neighbors
}

// NumAtoms returns the number of atoms.
func (m *Mol) NumAtoms() int {
	return len(m.Atoms)
}

// NumBonds returns the number of bonds.
func (m *Mol) NumBonds() int {
	return len(m.Bonds)
}
