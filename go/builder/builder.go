package builder

import (
	"math"
	"math/rand"
	"time"

	"smiles3d/mol"
)

// Build3D assigns initial 3D coordinates to a molecule.
// We use a simple rule-based random walk to place atoms, followed by a naive spring layout
// to untangle them slightly before handing off to the UFF optimizer.
func Build3D(m *mol.Mol) {
	rand.Seed(time.Now().UnixNano())

	// Step 1: Assign initial coordinates
	// If it's a single atom, put it at the origin.
	if m.NumAtoms() == 0 {
		return
	}

	m.GetAtom(1).X = 0.0
	m.GetAtom(1).Y = 0.0
	m.GetAtom(1).Z = 0.0

	// We place atoms randomly within a sphere to avoid exact overlaps.
	// Standard bond length is ~1.5 Angstroms.
	for i := 2; i <= m.NumAtoms(); i++ {
		a := m.GetAtom(i)

		// Find a bonded neighbor that already has coordinates
		var placedNeighbor *mol.Atom
		for _, b := range m.Bonds {
			if b.BeginAtom.Idx == a.Idx && b.EndAtom.Idx < a.Idx {
				placedNeighbor = b.EndAtom
				break
			} else if b.EndAtom.Idx == a.Idx && b.BeginAtom.Idx < a.Idx {
				placedNeighbor = b.BeginAtom
				break
			}
		}

		if placedNeighbor != nil {
			// Place at a random direction, ~1.5A away from neighbor
			theta := rand.Float64() * 2 * math.Pi
			phi := math.Acos(2*rand.Float64() - 1)
			r := 1.5 // typical bond length

			a.X = placedNeighbor.X + r*math.Sin(phi)*math.Cos(theta)
			a.Y = placedNeighbor.Y + r*math.Sin(phi)*math.Sin(theta)
			a.Z = placedNeighbor.Z + r*math.Cos(phi)
		} else {
			// Disconnected fragment (shouldn't happen in valid SMILES, but handle just in case)
			a.X = rand.Float64() * 5.0
			a.Y = rand.Float64() * 5.0
			a.Z = rand.Float64() * 5.0
		}
	}

	// Step 2: A quick 3D untangling step (simple spring model)
	// This helps UFF avoid getting stuck in terrible local minima.
	untangle(m, 50)
}

func untangle(m *mol.Mol, steps int) {
	kSpring := 1.0 // bond strength
	kRepel := 0.5  // repulsion strength
	idealDist := 1.5
	damping := 0.1

	numAtoms := m.NumAtoms()
	if numAtoms < 2 {
		return
	}

	velocities := make([][3]float64, numAtoms+1)

	for step := 0; step < steps; step++ {
		forces := make([][3]float64, numAtoms+1)

		// Repulsive forces between ALL atoms (prevent overlap)
		for i := 1; i <= numAtoms; i++ {
			for j := i + 1; j <= numAtoms; j++ {
				a1 := m.GetAtom(i)
				a2 := m.GetAtom(j)

				dx := a1.X - a2.X
				dy := a1.Y - a2.Y
				dz := a1.Z - a2.Z
				dist2 := dx*dx + dy*dy + dz*dz
				dist := math.Sqrt(dist2)

				if dist > 0.01 && dist < 5.0 { // only repel nearby atoms
					f := kRepel / (dist * dist)
					fx := f * (dx / dist)
					fy := f * (dy / dist)
					fz := f * (dz / dist)

					forces[i][0] += fx
					forces[i][1] += fy
					forces[i][2] += fz
					forces[j][0] -= fx
					forces[j][1] -= fy
					forces[j][2] -= fz
				}
			}
		}

		// Attractive spring forces for bonded atoms
		for _, b := range m.Bonds {
			a1 := b.BeginAtom
			a2 := b.EndAtom

			dx := a2.X - a1.X
			dy := a2.Y - a1.Y
			dz := a2.Z - a1.Z
			dist := math.Sqrt(dx*dx + dy*dy + dz*dz)

			if dist > 0.01 {
				// Hooke's Law: F = k(x - x0)
				f := kSpring * (dist - idealDist)
				fx := f * (dx / dist)
				fy := f * (dy / dist)
				fz := f * (dz / dist)

				forces[a1.Idx][0] += fx
				forces[a1.Idx][1] += fy
				forces[a1.Idx][2] += fz
				forces[a2.Idx][0] -= fx
				forces[a2.Idx][1] -= fy
				forces[a2.Idx][2] -= fz
			}
		}

		// Apply forces to update positions
		for i := 1; i <= numAtoms; i++ {
			a := m.GetAtom(i)

			velocities[i][0] = (velocities[i][0] + forces[i][0]) * damping
			velocities[i][1] = (velocities[i][1] + forces[i][1]) * damping
			velocities[i][2] = (velocities[i][2] + forces[i][2]) * damping

			a.X += velocities[i][0]
			a.Y += velocities[i][1]
			a.Z += velocities[i][2]
		}
	}
}
