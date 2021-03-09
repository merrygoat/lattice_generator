from sys import argv, exit

def main():
    filetype = 1  # if 0 write lammps output, if 1 write standard xyz output

    if filetype == 0:
        outputfile = open("startconfiguration.xyz", "w")
    else:
        outputfile = open("xyzstartconfiguration.xyz", "w")

    if len(argv) != 5:
        print("Incorrect syntax. Use: arguments lattice_size density, lattice_type (0 for Al2Cu 1 for cubic) and particle_ratio (2 for 2:1, 3 for 3:1 etc). e.g: ./lammpslatticegenerator.py 7 1.5 0 2")
        print("Note that ratio will only be used for cubic lattices. Al2Cu lattices are always 2:1 ratio.")
        exit()

    cubic_box_size = int(argv[1])
    density = float(argv[2])

    if int(argv[3]) == 0:  # if setting up Al2Cu lattice
        unitcelldimensions = [6.067, 6.067, 4.877]
        numberofparticles = 12 * cubic_box_size ** 3
    else:        # if using cubic lattice
        unitcelldimensions = [1, 1, 1]
        numberofparticles = cubic_box_size ** 3

    volume = numberofparticles / density
    sidelenratio = unitcelldimensions[0] / unitcelldimensions[2]
    sidez = (volume / (sidelenratio * sidelenratio)) ** (1.0 / 3)
    sidexy = sidez * sidelenratio
    unitcellsidexy = sidexy / cubic_box_size
    unitcellsidez = sidez / cubic_box_size

    print("unit cell side x + y = " + str(unitcellsidexy))
    print("unit cell side z = " + str(unitcellsidez))

    particlecounter = 1

    al2cucoordarray = [[0.65810, 0.15810, 0.50000], [0.15810, 0.65810, 0.00000], [0.84190, 0.34190, 0.00000], [0.34190, 0.84190, 0.50000], [0.15810, 0.34190, 0.50000], [0.65810, 0.84190, 0.00000], [.84190, 0.65810, 0.50000], [0.34190, 0.15810, 0.00000], [0.50000, 0.50000, 0.75000], [0.00000, 0.00000, 0.25000], [0.50000, 0.50000, 0.25000], [0.00000, 0.00000, 0.75000]]

    if filetype == 1:
        outputfile.write(str(numberofparticles) + "\n")

    outputfile.write("Auto generated 11A crystal lattice for LAMMPS\n")

    if filetype == 0:
        outputfile.write(f"\n{numberofparticles} atoms\n\n")
        outputfile.write("2 atom types\n\n")
        outputfile.write(f"{-0.5 * sidexy:.7f} {sidexy * 0.5:.7f} xlo xhi\n")
        outputfile.write(f"{-0.5 * sidexy:.7f} {sidexy * 0.5:.7f} ylo yhi\n")
        outputfile.write(f"{-0.5 * sidez:.7f} {sidez * 0.5:.7f} zlo zhi\n\n")
        outputfile.write("Masses\n\n1 1\n2 1\n\nAtoms\n\n")

    if int(argv[3]) == 0:   # Al2Cu lattice
        for x in range(cubic_box_size):
            for y in range(cubic_box_size):
                for z in range(cubic_box_size):
                    for unit_cell_counter in range(12):
                        lammpsprefix = f"{particlecounter}\t"
                        if unit_cell_counter > 7:
                            particle_type = 2
                        else:
                            particle_type = 1
                        x_pos = al2cucoordarray[unit_cell_counter][0] * unitcellsidexy + x * unitcellsidexy - 0.5 * sidexy
                        y_pos = al2cucoordarray[unit_cell_counter][1] * unitcellsidexy + y * unitcellsidexy - 0.5 * sidexy
                        z_pos = al2cucoordarray[unit_cell_counter][2] * unitcellsidez + z* unitcellsidez - 0.5 * sidez
                        if filetype == 0:
                            outputfile.write(lammpsprefix)
                        outputfile.write(f"{particle_type}\t{x_pos:.7f}\t{y_pos:.7f}\t{z_pos:.7f}\n")
                        particlecounter += 1
    else:   # Cubic lattice
        for x in range(cubic_box_size):
            for y in range(cubic_box_size):
                for z in range(cubic_box_size):
                    if particlecounter % (int(argv[4]) + 1) == 0:
                        particle_type = 2
                    else:
                        particle_type = 1
                    x_pos = x * unitcellsidexy - 0.5 * sidexy
                    y_pos = y * unitcellsidexy - 0.5 * sidexy
                    z_pos = z * unitcellsidez - 0.5 * sidez
                    lammpsprefix = f"{particlecounter}\t"
                    lineout = f"{particle_type}\t{x_pos:.7f}\t{y_pos:.7f}\t{z_pos:.7f}\n"
                    if filetype == 0:
                        outputfile.write(lammpsprefix + lineout)
                    else:
                        outputfile.write(lineout)
                    particlecounter += 1
    outputfile.close()

main()