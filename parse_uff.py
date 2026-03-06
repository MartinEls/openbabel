import sys

def main():
    print("package uff")
    print("")
    print("import \"math\"")
    print("")
    print("// UFFParams defines the parameters for a specific UFF atom type.")
    print("type UFFParams struct {")
    print("\tBondRadius float64 // R1")
    print("\tAngle      float64 // Theta0 (radians)")
    print("\tVdwRadius  float64 // x1")
    print("\tVdwDepth   float64 // D1")
    print("\tZeff       float64 // Z*")
    print("}")
    print("")
    print("var params = map[string]UFFParams{")

    with open("data/UFF.prm", "r") as f:
        for line in f:
            if line.startswith("param "):
                parts = line.split()
                if len(parts) >= 7:
                    atom_type = parts[1]
                    try:
                        r1 = float(parts[2])
                        theta0_deg = float(parts[3])
                        x1 = float(parts[4])
                        d1 = float(parts[5])
                        zeta = float(parts[6]) # not used in go/uff/uff.go but keep in struct? Wait, zeff is parts[7]
                        zeff = float(parts[7])
                        print(f'\t"{atom_type}": {{BondRadius: {r1}, Angle: {theta0_deg} * math.Pi / 180.0, VdwRadius: {x1}, VdwDepth: {d1}, Zeff: {zeff}}},')
                    except ValueError:
                        pass
    print("}")

if __name__ == "__main__":
    main()
