package uff

import (
	"smiles3d/mol"
)

// assignUFFTypes returns a map from Atom Index to UFF parameter string type.
func assignUFFTypes(m *mol.Mol) map[int]string {
	types := make(map[int]string)

	// Create adjacency for bond orders
	neighbors := make(map[int][]*mol.Bond)
	for _, b := range m.Bonds {
		neighbors[b.BeginAtom.Idx] = append(neighbors[b.BeginAtom.Idx], b)
		neighbors[b.EndAtom.Idx] = append(neighbors[b.EndAtom.Idx], b)
	}

	for _, a := range m.Atoms {
		typ := assignTypeForAtom(a, neighbors[a.Idx])
		types[a.Idx] = typ
	}

	return types
}

func assignTypeForAtom(a *mol.Atom, bonds []*mol.Bond) string {
	// Calculate sum of bond orders
	// (Note: in SDF or OpenBabel, aromatic bonds might have order 1 with a flag,
	// or order 4. In our builder, we map 4 -> 1 but set IsAromatic.)
	sumOrder := 0
	hasDouble := false
	hasTriple := false
	hasAromatic := a.IsAromatic

	for _, b := range bonds {
		order := b.Order
		if b.IsAromatic {
			hasAromatic = true
		}
		if order == 2 {
			hasDouble = true
		} else if order == 3 {
			hasTriple = true
		}
		sumOrder += order
	}

	sym := a.Symbol()

	switch sym {
	case "H":
		return "H_"
	case "C":
		if hasAromatic {
			return "C_R"
		}
		if hasTriple || sumOrder >= 4 && hasDouble && len(bonds) == 2 {
			return "C_1"
		}
		if hasDouble || sumOrder == 3 && len(bonds) == 3 {
			return "C_2"
		}
		return "C_3"
	case "N":
		if hasAromatic {
			return "N_R"
		}
		if hasTriple {
			return "N_1"
		}
		if hasDouble {
			return "N_2"
		}
		return "N_3"
	case "O":
		if hasAromatic {
			return "O_R"
		}
		if hasDouble {
			return "O_2"
		}
		return "O_3"
	case "S":
		if hasAromatic {
			return "S_R"
		}
		if hasDouble {
			return "S_2"
		}
		return "S_3+2" // Default generic S
	case "P":
		return "P_3+3"
	case "F":
		return "F_"
	case "Cl":
		return "Cl"
	case "Br":
		return "Br"
	case "I":
		return "I_"
	default:
		// Attempt to fallback to generic name if available in params
		if _, ok := params[sym]; ok {
			return sym
		}
		// Some atoms have generic types like Element_
		if _, ok := params[sym+"_"]; ok {
			return sym + "_"
		}
		// Some like Element3+something
		// We'll just return Du (Dummy) if not found, to avoid crash.
		return "Du"
	}
}
