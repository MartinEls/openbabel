package smiles

import (
	"fmt"
	"strconv"
	"strings"

	"smiles3d/mol"
)

// ParseSMILES parses a subset of SMILES strings.
// Supports: elements (C,N,O,S,P,F,Cl,Br,I), aromaticity (lowercase c,n,o,s,p),
// rings (1-9, %10-%99), branches (), and basic bond types (-,=,#,:).
// Ignores stereochemical symbols (@, @@, /, \).
func ParseSMILES(s string) (*mol.Mol, error) {
	m := mol.NewMol()

	// Clean stereocenters and other unneeded markers.
	// Since stereochemistry does not affect basic atom/bond graph connectivity,
	// we just ignore / \ @.
	s = strings.ReplaceAll(s, "@", "")
	s = strings.ReplaceAll(s, "/", "")
	s = strings.ReplaceAll(s, "\\", "")

	// A simple stack-based parser.
	type stackEntry struct {
		atomIdx int
		idx     int // pos in string
	}
	var stack []stackEntry

	// Ring lookup map for closing rings: ringNumber -> atomIdx
	rings := make(map[int]int)

	var lastAtom *mol.Atom

	var i int
	n := len(s)

	for i < n {
		c := s[i]

		switch c {
		case '(': // Branch start
			if lastAtom == nil {
				return nil, fmt.Errorf("branch start without previous atom at pos %d", i)
			}
			stack = append(stack, stackEntry{lastAtom.Idx, i})
			i++
		case ')': // Branch end
			if len(stack) == 0 {
				return nil, fmt.Errorf("unmatched branch end at pos %d", i)
			}
			top := stack[len(stack)-1]
			stack = stack[:len(stack)-1]
			lastAtom = m.GetAtom(top.atomIdx)
			i++
		case '-', '=', '#', ':': // Explicit bonds
			// Just remember the bond type.
			bondOrder := 1
			isAromatic := false
			if c == '=' {
				bondOrder = 2
			} else if c == '#' {
				bondOrder = 3
			} else if c == ':' {
				isAromatic = true
				bondOrder = 4
			}
			i++

			// Get next atom and attach it with the specified bond type.
			if i >= n {
				return nil, fmt.Errorf("unexpected end of string after bond type at pos %d", i-1)
			}

			nextAtom, nextI, err := readAtom(s, i, m)
			if err != nil {
				return nil, fmt.Errorf("error reading atom after bond at pos %d: %v", i, err)
			}
			if lastAtom != nil {
				m.AddBond(lastAtom.Idx, nextAtom.Idx, bondOrder, isAromatic)
			}
			lastAtom = nextAtom
			i = nextI
		case '%': // Two digit ring number
			if i+2 >= n {
				return nil, fmt.Errorf("unexpected end of string parsing ring at pos %d", i)
			}
			ringNum, err := strconv.Atoi(string(s[i+1 : i+3]))
			if err != nil {
				return nil, fmt.Errorf("invalid ring number at pos %d", i)
			}
			i += 3
			handleRing(m, lastAtom, ringNum, rings)
		default:
			// One digit ring
			if c >= '1' && c <= '9' {
				ringNum := int(c - '0')
				i++
				handleRing(m, lastAtom, ringNum, rings)
			} else {
				// Regular atom (e.g. C, c, [nH])
				atom, nextI, err := readAtom(s, i, m)
				if err != nil {
					return nil, fmt.Errorf("error reading atom at pos %d: %v", i, err)
				}
				if lastAtom != nil {
					// Implicit bond logic: single by default, or aromatic if both are aromatic.
					bondOrder := 1
					isAromatic := false
					if lastAtom.IsAromatic && atom.IsAromatic {
						isAromatic = true
						bondOrder = 4
					}
					m.AddBond(lastAtom.Idx, atom.Idx, bondOrder, isAromatic)
				}
				lastAtom = atom
				i = nextI
			}
		}
	}

	AddHydrogens(m)

	return m, nil
}

func handleRing(m *mol.Mol, currentAtom *mol.Atom, ringNum int, rings map[int]int) {
	if currentAtom == nil {
		return
	}

	if startAtomIdx, exists := rings[ringNum]; exists {
		// Close ring
		startAtom := m.GetAtom(startAtomIdx)

		bondOrder := 1
		isAromatic := false
		if currentAtom.IsAromatic && startAtom.IsAromatic {
			isAromatic = true
			bondOrder = 4
		}
		m.AddBond(currentAtom.Idx, startAtom.Idx, bondOrder, isAromatic)
		delete(rings, ringNum)
	} else {
		// Open ring
		rings[ringNum] = currentAtom.Idx
	}
}

func readAtom(s string, startIdx int, m *mol.Mol) (*mol.Atom, int, error) {
	i := startIdx
	c := s[i]

	// Handle explicit atom definitions like [nH], [14C], [O-]
	if c == '[' {
		// Simple parsing for [...] block. We only care about the symbol and charge.
		end := strings.Index(s[i:], "]")
		if end == -1 {
			return nil, i, fmt.Errorf("unmatched bracket at pos %d", i)
		}

		content := s[i+1 : i+end]
		i += end + 1 // skip to next char after ]

		// Strip leading numbers (isotopes)
		contentStart := 0
		for contentStart < len(content) && content[contentStart] >= '0' && content[contentStart] <= '9' {
			contentStart++
		}
		content = content[contentStart:]

		// Very basic parsing
		// Look for charge
		charge := 0
		if strings.HasSuffix(content, "+") {
			charge = 1
			content = content[:len(content)-1]
		} else if strings.HasSuffix(content, "-") {
			charge = -1
			content = content[:len(content)-1]
		}

		// Look for explicit H
		explicitH := 0
		if idx := strings.Index(content, "H"); idx != -1 {
			if idx+1 < len(content) && content[idx+1] >= '0' && content[idx+1] <= '9' {
				explicitH = int(content[idx+1] - '0')
				content = content[:idx] // strip H#
			} else {
				explicitH = 1
				content = content[:idx] // strip H
			}
		}

		sym := content
		atomicNum, isAromatic := symbolToAtomicNum(sym)

		atom := m.AddAtom(atomicNum)
		atom.IsAromatic = isAromatic
		atom.FormalCharge = charge
		atom.ImplicitH = explicitH // Not fully implicit, but we store it here for simplicity

		return atom, i, nil
	}

	// Normal 1 or 2 character symbol
	var sym string
	if i+1 < len(s) && isLower(s[i+1]) && s[i+1] != 'c' && s[i+1] != 'n' && s[i+1] != 'o' && s[i+1] != 's' && s[i+1] != 'p' {
		// Multi-char atom like Cl, Br
		sym = s[i : i+2]
		i += 2
	} else {
		sym = s[i : i+1]
		i += 1
	}

	atomicNum, isAromatic := symbolToAtomicNum(sym)
	if atomicNum == 0 {
		return nil, startIdx, fmt.Errorf("unknown symbol %q at pos %d", sym, startIdx)
	}

	atom := m.AddAtom(atomicNum)
	atom.IsAromatic = isAromatic
	return atom, i, nil
}

func isLower(c byte) bool {
	return c >= 'a' && c <= 'z'
}

func symbolToAtomicNum(sym string) (int, bool) {
	isArom := false
	if sym == strings.ToLower(sym) && len(sym) == 1 {
		isArom = true
	}
	symLower := strings.ToLower(sym)

	switch symLower {
	case "h":
		return 1, false
	case "c":
		return 6, isArom
	case "n":
		return 7, isArom
	case "o":
		return 8, isArom
	case "f":
		return 9, false
	case "p":
		return 15, isArom
	case "s":
		return 16, isArom
	case "cl":
		return 17, false
	case "br":
		return 35, false
	case "i":
		return 53, false
	default:
		return 0, false
	}
}

// AddHydrogens calculates implicit valences and adds explicit hydrogen atoms to the Mol.
func AddHydrogens(m *mol.Mol) {
	targetValence := map[int]int{
		6:  4, // Carbon
		7:  3, // Nitrogen
		8:  2, // Oxygen
		9:  1, // Fluorine
		15: 3, // Phosphorus
		16: 2, // Sulfur
		17: 1, // Chlorine
		35: 1, // Bromine
		53: 1, // Iodine
	}

	// First pass: Calculate sum of bond orders for each atom
	// We iterate through all atoms and sum their connected bond orders.
	currentBonds := make(map[int]int)

	for _, b := range m.Bonds {
		order := b.Order
		if b.IsAromatic {
			// Aromatic bonds average to ~1.5. For valence math, we typically count it as 1.5,
			// but a simpler integer approach for CHNO works:
			// an aromatic carbon usually has 3 neighbors (2 aromatic bonds, 1 single).
			// Let's just count aromatic as order 1.5 rounded somehow, or just handle it differently.
			order = 1
			// To properly add H to aromatic systems without full Kekulization,
			// a simple rule is: Aromatic C normally has 1 H. Aromatic N with no charge/H has 0 H.
		}
		currentBonds[b.BeginAtom.Idx] += order
		currentBonds[b.EndAtom.Idx] += order
	}

	// Make a snapshot of atoms to avoid adding hydrogens to newly added hydrogens.
	numAtoms := m.NumAtoms()
	for i := 1; i <= numAtoms; i++ {
		atom := m.GetAtom(i)
		if atom.AtomicNum == 1 {
			continue // skip hydrogens
		}

		if atom.ImplicitH > 0 {
			// Explicitly defined in SMILES (e.g. [nH])
			for h := 0; h < atom.ImplicitH; h++ {
				hAtom := m.AddAtom(1)
				m.AddBond(atom.Idx, hAtom.Idx, 1, false)
			}
			continue
		}

		// If explicit H was not defined via bracket, infer from valence
		target, ok := targetValence[atom.AtomicNum]
		if !ok {
			continue // Unrecognized atom type for H addition
		}

		val := currentBonds[atom.Idx]

		// Adjust target for formal charge
		if atom.AtomicNum == 6 && atom.FormalCharge != 0 {
			target -= 1
		} // Carbocations/carbanions
		if atom.AtomicNum == 7 && atom.FormalCharge == 1 {
			target = 4
		} // Ammonium
		if atom.AtomicNum == 8 && atom.FormalCharge == -1 {
			target = 1
		} // Alkoxide
		if atom.AtomicNum == 8 && atom.FormalCharge == 1 {
			target = 3
		} // Oxonium

		var numHToAdd int
		if atom.IsAromatic {
			// Simple heuristic for aromatic rings:
			// e.g., Benzene C has 2 aromatic neighbors. Total 'bonds' counted = 2.
			// True valence needed = 4. 2 aromatic bonds roughly equal 3 single bonds.
			// So numH = 4 - 3 = 1.

			// Count neighbors
			neighbors := len(m.GetNeighbors(atom))
			if atom.AtomicNum == 6 {
				numHToAdd = 3 - neighbors // aromatic C bonded to 2 other atoms in ring needs 1 H.
			} else if atom.AtomicNum == 7 {
				numHToAdd = 2 - neighbors // aromatic N bonded to 2 atoms in ring (like pyridine) needs 0 H.
			}
		} else {
			numHToAdd = target - val
		}

		if numHToAdd < 0 {
			numHToAdd = 0
		}

		for h := 0; h < numHToAdd; h++ {
			hAtom := m.AddAtom(1)
			m.AddBond(atom.Idx, hAtom.Idx, 1, false)
		}
	}
}
