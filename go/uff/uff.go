package uff

import (
	"math"

	"smiles3d/mol"
)

// CalculateEnergy calculates the current UFF energy of the molecule (in kcal/mol).
func CalculateEnergy(m *mol.Mol) float64 {
	n := m.NumAtoms()
	if n <= 1 {
		return 0.0
	}

	atomTypes := assignUFFTypes(m)
	neighbors := make(map[int][]*mol.Atom)
	for i := 1; i <= n; i++ {
		neighbors[i] = m.GetNeighbors(m.GetAtom(i))
	}

	totalEnergy := 0.0

	// 1. Bond Stretch Energy
	for _, b := range m.Bonds {
		a1 := b.BeginAtom
		a2 := b.EndAtom
		t1 := atomTypes[a1.Idx]
		t2 := atomTypes[a2.Idx]

		p1, ok1 := params[t1]
		p2, ok2 := params[t2]
		if !ok1 || !ok2 {
			continue
		}

		bo := float64(b.Order)
		if b.IsAromatic {
			bo = 1.5
		}

		r0 := p1.BondRadius + p2.BondRadius
		kb := 664.12 * (p1.Zeff * p2.Zeff) / math.Pow(r0, 3) * bo

		dx := a2.X - a1.X
		dy := a2.Y - a1.Y
		dz := a2.Z - a1.Z
		r := math.Sqrt(dx*dx + dy*dy + dz*dz)

		// Harmonic oscillator energy: E = 0.5 * kb * (r - r0)^2
		totalEnergy += 0.5 * kb * (r - r0) * (r - r0)
	}

	// 2. Angle Bend Energy
	for i := 1; i <= n; i++ {
		center := m.GetAtom(i)
		tCenter := atomTypes[center.Idx]
		pCenter, okC := params[tCenter]
		if !okC {
			continue
		}

		nbrs := neighbors[center.Idx]
		if len(nbrs) < 2 {
			continue
		}

		theta0 := pCenter.Angle
		ka := 100.0 // Simplified fixed force constant for bending

		for j := 0; j < len(nbrs); j++ {
			for k := j + 1; k < len(nbrs); k++ {
				a1 := nbrs[j]
				a2 := nbrs[k]

				v1x := a1.X - center.X
				v1y := a1.Y - center.Y
				v1z := a1.Z - center.Z
				r1 := math.Sqrt(v1x*v1x + v1y*v1y + v1z*v1z)

				v2x := a2.X - center.X
				v2y := a2.Y - center.Y
				v2z := a2.Z - center.Z
				r2 := math.Sqrt(v2x*v2x + v2y*v2y + v2z*v2z)

				if r1 < 0.001 || r2 < 0.001 {
					continue
				}

				dot := v1x*v2x + v1y*v2y + v1z*v2z
				cosTheta := dot / (r1 * r2)
				if cosTheta > 1.0 {
					cosTheta = 1.0
				}
				if cosTheta < -1.0 {
					cosTheta = -1.0
				}

				theta := math.Acos(cosTheta)

				// E = 0.5 * ka * (theta - theta0)^2
				totalEnergy += 0.5 * ka * (theta - theta0) * (theta - theta0)
			}
		}
	}

	// 3. Van der Waals Energy
	for i := 1; i <= n; i++ {
		for j := i + 1; j <= n; j++ {
			bonded := false
			for _, a := range neighbors[i] {
				if a.Idx == j {
					bonded = true
					break
				}
			}

			is13 := false
			if !bonded {
				for _, a := range neighbors[i] {
					for _, b := range neighbors[j] {
						if a.Idx == b.Idx {
							is13 = true
							break
						}
					}
				}
			}

			if bonded || is13 {
				continue
			}

			a1 := m.GetAtom(i)
			a2 := m.GetAtom(j)
			t1 := atomTypes[i]
			t2 := atomTypes[j]
			p1, ok1 := params[t1]
			p2, ok2 := params[t2]
			if !ok1 || !ok2 {
				continue
			}

			dx := a1.X - a2.X
			dy := a1.Y - a2.Y
			dz := a1.Z - a2.Z
			r := math.Sqrt(dx*dx + dy*dy + dz*dz)
			if r < 0.1 {
				continue
			}

			xij := math.Sqrt(p1.VdwRadius * p2.VdwRadius)
			dij := math.Sqrt(p1.VdwDepth * p2.VdwDepth)

			xr := xij / r
			xr6 := xr * xr * xr * xr * xr * xr
			xr12 := xr6 * xr6

			// E = D * ( (x/r)^12 - 2 * (x/r)^6 )
			totalEnergy += dij * (xr12 - 2.0*xr6)
		}
	}

	return totalEnergy
}

// Optimize takes a molecule and runs a simplified steepest descent UFF optimization.
func Optimize(m *mol.Mol, maxIter int, stepSize float64) {
	n := m.NumAtoms()
	if n <= 1 {
		return
	}

	// Assign UFF string types
	atomTypes := assignUFFTypes(m)

	// Prepare connectivity (for angles and exclusion)
	neighbors := make(map[int][]*mol.Atom)
	for i := 1; i <= n; i++ {
		neighbors[i] = m.GetNeighbors(m.GetAtom(i))
	}

	for iter := 0; iter < maxIter; iter++ {
		forces := make([][3]float64, n+1)

		// 1. Bond Stretch Energy (Harmonic)
		for _, b := range m.Bonds {
			a1 := b.BeginAtom
			a2 := b.EndAtom
			t1 := atomTypes[a1.Idx]
			t2 := atomTypes[a2.Idx]

			p1, ok1 := params[t1]
			p2, ok2 := params[t2]

			if !ok1 || !ok2 {
				continue
			}

			// BO = bond order
			bo := float64(b.Order)
			if b.IsAromatic {
				bo = 1.5
			}

			// Ideal bond length (r0) = ri + rj + bond-order correction (ignored for simplicity)
			r0 := p1.BondRadius + p2.BondRadius

			// Simple Pauling electronegativity bond order correction
			// r0 -= 0.1332 * math.Log(bo)
			// But for simplicity, we just use r0 directly for a robust initial geometry

			// Force constant (k_b) roughly 600-700 for single bonds, scaled by BO
			kb := 664.12 * (p1.Zeff * p2.Zeff) / math.Pow(r0, 3) * bo

			dx := a2.X - a1.X
			dy := a2.Y - a1.Y
			dz := a2.Z - a1.Z
			r := math.Sqrt(dx*dx + dy*dy + dz*dz)

			if r < 0.001 {
				continue
			} // avoid singularity

			// F = k_b * (r - r0)
			fMag := kb * (r - r0)

			fx := fMag * dx / r
			fy := fMag * dy / r
			fz := fMag * dz / r

			forces[a1.Idx][0] += fx
			forces[a1.Idx][1] += fy
			forces[a1.Idx][2] += fz
			forces[a2.Idx][0] -= fx
			forces[a2.Idx][1] -= fy
			forces[a2.Idx][2] -= fz
		}

		// 2. Angle Bend Energy (Simplified)
		for i := 1; i <= n; i++ {
			center := m.GetAtom(i)
			tCenter := atomTypes[center.Idx]
			pCenter, okC := params[tCenter]
			if !okC {
				continue
			}

			nbrs := neighbors[center.Idx]
			if len(nbrs) < 2 {
				continue
			}

			theta0 := pCenter.Angle

			// Force constant
			ka := 100.0 // Simplified fixed force constant for bending

			for j := 0; j < len(nbrs); j++ {
				for k := j + 1; k < len(nbrs); k++ {
					a1 := nbrs[j]
					a2 := nbrs[k]

					v1x := a1.X - center.X
					v1y := a1.Y - center.Y
					v1z := a1.Z - center.Z
					r1 := math.Sqrt(v1x*v1x + v1y*v1y + v1z*v1z)

					v2x := a2.X - center.X
					v2y := a2.Y - center.Y
					v2z := a2.Z - center.Z
					r2 := math.Sqrt(v2x*v2x + v2y*v2y + v2z*v2z)

					if r1 < 0.001 || r2 < 0.001 {
						continue
					}

					dot := v1x*v2x + v1y*v2y + v1z*v2z
					cosTheta := dot / (r1 * r2)

					// Clamp to [-1, 1]
					if cosTheta > 1.0 {
						cosTheta = 1.0
					}
					if cosTheta < -1.0 {
						cosTheta = -1.0
					}

					theta := math.Acos(cosTheta)

					// E = 0.5 * ka * (theta - theta0)^2
					// F = ka * (theta - theta0)
					fMag := ka * (theta - theta0)

					// Calculate derivatives using simplified central difference for forces
					delta := 0.001

					// derivative wrt a1.X
					v1xNew := v1x + delta
					r1New := math.Sqrt(v1xNew*v1xNew + v1y*v1y + v1z*v1z)
					cosNew := (v1xNew*v2x + v1y*v2y + v1z*v2z) / (r1New * r2)
					if cosNew > 1.0 {
						cosNew = 1.0
					}
					if cosNew < -1.0 {
						cosNew = -1.0
					}
					thetaNew := math.Acos(cosNew)
					da1x := (thetaNew - theta) / delta

					v1yNew := v1y + delta
					r1NewY := math.Sqrt(v1x*v1x + v1yNew*v1yNew + v1z*v1z)
					cosNewY := (v1x*v2x + v1yNew*v2y + v1z*v2z) / (r1NewY * r2)
					if cosNewY > 1.0 {
						cosNewY = 1.0
					}
					if cosNewY < -1.0 {
						cosNewY = -1.0
					}
					thetaNewY := math.Acos(cosNewY)
					da1y := (thetaNewY - theta) / delta

					v1zNew := v1z + delta
					r1NewZ := math.Sqrt(v1x*v1x + v1y*v1y + v1zNew*v1zNew)
					cosNewZ := (v1x*v2x + v1y*v2y + v1zNew*v2z) / (r1NewZ * r2)
					if cosNewZ > 1.0 {
						cosNewZ = 1.0
					}
					if cosNewZ < -1.0 {
						cosNewZ = -1.0
					}
					thetaNewZ := math.Acos(cosNewZ)
					da1z := (thetaNewZ - theta) / delta

					forces[a1.Idx][0] -= fMag * da1x
					forces[a1.Idx][1] -= fMag * da1y
					forces[a1.Idx][2] -= fMag * da1z

					// Force on a2
					v2xNew := v2x + delta
					r2New := math.Sqrt(v2xNew*v2xNew + v2y*v2y + v2z*v2z)
					cosNew2 := (v1x*v2xNew + v1y*v2y + v1z*v2z) / (r1 * r2New)
					if cosNew2 > 1.0 {
						cosNew2 = 1.0
					}
					if cosNew2 < -1.0 {
						cosNew2 = -1.0
					}
					thetaNew2 := math.Acos(cosNew2)
					da2x := (thetaNew2 - theta) / delta

					v2yNew := v2y + delta
					r2NewY := math.Sqrt(v2x*v2x + v2yNew*v2yNew + v2z*v2z)
					cosNew2Y := (v1x*v2x + v1y*v2yNew + v1z*v2z) / (r1 * r2NewY)
					if cosNew2Y > 1.0 {
						cosNew2Y = 1.0
					}
					if cosNew2Y < -1.0 {
						cosNew2Y = -1.0
					}
					thetaNew2Y := math.Acos(cosNew2Y)
					da2y := (thetaNew2Y - theta) / delta

					v2zNew := v2z + delta
					r2NewZ := math.Sqrt(v2x*v2x + v2y*v2y + v2zNew*v2zNew)
					cosNew2Z := (v1x*v2x + v1y*v2y + v1z*v2zNew) / (r1 * r2NewZ)
					if cosNew2Z > 1.0 {
						cosNew2Z = 1.0
					}
					if cosNew2Z < -1.0 {
						cosNew2Z = -1.0
					}
					thetaNew2Z := math.Acos(cosNew2Z)
					da2z := (thetaNew2Z - theta) / delta

					forces[a2.Idx][0] -= fMag * da2x
					forces[a2.Idx][1] -= fMag * da2y
					forces[a2.Idx][2] -= fMag * da2z

					// Force on center (conservation of momentum: sum F = 0)
					forces[center.Idx][0] += fMag * (da1x + da2x)
					forces[center.Idx][1] += fMag * (da1y + da2y)
					forces[center.Idx][2] += fMag * (da1z + da2z)
				}
			}
		}

		// 3. Van der Waals (Lennard-Jones)
		// We compute VdW between all non-bonded atoms (1-4 distance and beyond)
		for i := 1; i <= n; i++ {
			for j := i + 1; j <= n; j++ {
				// Are they bonded?
				bonded := false
				for _, a := range neighbors[i] {
					if a.Idx == j {
						bonded = true
						break
					}
				}

				// Are they 1-3 bonded (sharing a common neighbor)?
				is13 := false
				if !bonded {
					for _, a := range neighbors[i] {
						for _, b := range neighbors[j] {
							if a.Idx == b.Idx {
								is13 = true
								break
							}
						}
					}
				}

				if bonded || is13 {
					continue
				}

				a1 := m.GetAtom(i)
				a2 := m.GetAtom(j)
				t1 := atomTypes[i]
				t2 := atomTypes[j]
				p1, ok1 := params[t1]
				p2, ok2 := params[t2]
				if !ok1 || !ok2 {
					continue
				}

				dx := a1.X - a2.X
				dy := a1.Y - a2.Y
				dz := a1.Z - a2.Z
				r2 := dx*dx + dy*dy + dz*dz
				r := math.Sqrt(r2)

				if r < 0.1 {
					continue
				} // avoid explosion

				xij := math.Sqrt(p1.VdwRadius * p2.VdwRadius)
				dij := math.Sqrt(p1.VdwDepth * p2.VdwDepth)

				// Lennard-Jones Force: F = -24 * D * (2*(X/r)^12 - (X/r)^6) / r
				// Force is negative for attraction, positive for repulsion
				xr := xij / r
				xr6 := xr * xr * xr * xr * xr * xr
				xr12 := xr6 * xr6

				fMag := 24.0 * dij * (2.0*xr12 - xr6) / r

				fx := fMag * dx / r
				fy := fMag * dy / r
				fz := fMag * dz / r

				forces[i][0] += fx
				forces[i][1] += fy
				forces[i][2] += fz
				forces[j][0] -= fx
				forces[j][1] -= fy
				forces[j][2] -= fz
			}
		}

		// Update positions (Steepest Descent)
		maxForce := 0.0
		for i := 1; i <= n; i++ {
			a := m.GetAtom(i)
			fx, fy, fz := forces[i][0], forces[i][1], forces[i][2]

			fTotal := math.Sqrt(fx*fx + fy*fy + fz*fz)
			if fTotal > maxForce {
				maxForce = fTotal
			}

			// Clamp force vector to avoid exploding in highly strained initial geometries
			clamp := 5.0
			if fTotal > clamp {
				fx = fx * clamp / fTotal
				fy = fy * clamp / fTotal
				fz = fz * clamp / fTotal
			}

			a.X -= stepSize * fx
			a.Y -= stepSize * fy
			a.Z -= stepSize * fz
		}

		// Optional convergence criteria could be added here
		// if maxForce < 0.01 { break }
	}
}
