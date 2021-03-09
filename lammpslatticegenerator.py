import argparse


def main(args: dict):

    if args["output_type"] == "lammps":
        output_file = open("out.lmp", "w")
    else:
        output_file = open("out.xyz", "w")

    if args["type"] == "al2cu":
        unitcelldimensions = [6.067, 6.067, 4.877]
        num_particles = 12 * args["lattice_size"] ** 3
    else:
        # cubic lattice
        unitcelldimensions = [1, 1, 1]
        num_particles = args["lattice_size"] ** 3

    volume = num_particles / args["density"]
    sidelenratio = unitcelldimensions[0] / unitcelldimensions[2]
    sidez = (volume / (sidelenratio * sidelenratio)) ** (1.0 / 3)
    sidexy = sidez * sidelenratio
    unitcellsidexy = sidexy / args["lattice_size"]
    unitcellsidez = sidez / args["lattice_size"]

    print("unit cell side x + y = " + str(unitcellsidexy))
    print("unit cell side z = " + str(unitcellsidez))

    particlecounter = 1

    al2cucoordarray = [[0.65810, 0.15810, 0.50000],
                       [0.15810, 0.65810, 0.00000],
                       [0.84190, 0.34190, 0.00000],
                       [0.34190, 0.84190, 0.50000],
                       [0.15810, 0.34190, 0.50000],
                       [0.65810, 0.84190, 0.00000],
                       [.84190, 0.65810, 0.50000],
                       [0.34190, 0.15810, 0.00000],
                       [0.50000, 0.50000, 0.75000],
                       [0.00000, 0.00000, 0.25000],
                       [0.50000, 0.50000, 0.25000],
                       [0.00000, 0.00000, 0.75000]]

    if args["output_type"] == "xyz":
        output_file.write(str(num_particles) + "\n")

    output_file.write("Auto generated 11A crystal lattice for LAMMPS\n")

    if args["output_type"] == "lammps":
        output_file.write(f"\n{num_particles} atoms\n\n")
        output_file.write("2 atom types\n\n")
        output_file.write(f"{-0.5 * sidexy:.7f} {sidexy * 0.5:.7f} xlo xhi\n")
        output_file.write(f"{-0.5 * sidexy:.7f} {sidexy * 0.5:.7f} ylo yhi\n")
        output_file.write(f"{-0.5 * sidez:.7f} {sidez * 0.5:.7f} zlo zhi\n\n")
        output_file.write("Masses\n\n1 1\n2 1\n\nAtoms\n\n")

    if args["type"] == "al2cu":
        for x in range(args["lattice_size"]):
            for y in range(args["lattice_size"]):
                for z in range(args["lattice_size"]):
                    for unit_cell_counter in range(12):
                        lammpsprefix = f"{particlecounter}\t"
                        if unit_cell_counter > 7:
                            particle_type = 2
                        else:
                            particle_type = 1
                        x_pos = al2cucoordarray[unit_cell_counter][0] * unitcellsidexy + x * unitcellsidexy - 0.5 * sidexy
                        y_pos = al2cucoordarray[unit_cell_counter][1] * unitcellsidexy + y * unitcellsidexy - 0.5 * sidexy
                        z_pos = al2cucoordarray[unit_cell_counter][2] * unitcellsidez + z* unitcellsidez - 0.5 * sidez
                        if args["output_type"] == "lammps":
                            output_file.write(lammpsprefix)
                        output_file.write(f"{particle_type}\t{x_pos:.7f}\t{y_pos:.7f}\t{z_pos:.7f}\n")
                        particlecounter += 1
    else:   # Cubic lattice
        for x in range(args["lattice_size"]):
            for y in range(args["lattice_size"]):
                for z in range(args["lattice_size"]):
                    if particlecounter % (args["ratio"] + 1) == 0:
                        particle_type = 2
                    else:
                        particle_type = 1
                    x_pos = x * unitcellsidexy - 0.5 * sidexy
                    y_pos = y * unitcellsidexy - 0.5 * sidexy
                    z_pos = z * unitcellsidez - 0.5 * sidez
                    lammpsprefix = f"{particlecounter}\t"
                    lineout = f"{particle_type}\t{x_pos:.7f}\t{y_pos:.7f}\t{z_pos:.7f}\n"
                    if args["output_type"] == "lammps":
                        output_file.write(lammpsprefix + lineout)
                    else:
                        output_file.write(lineout)
                    particlecounter += 1
    output_file.close()


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser()

    PARSER.add_argument("-o", "--output_type", help="The format of the output, 'xyz' or 'lammps'.",
                        default="xyz")
    PARSER.add_argument("-s", "--lattice_size", help="The number of unit cells in each axis.",
                        default=5, type=int)
    PARSER.add_argument("-d", "--density", help="The number density of the lattice.", default=1,
                        type=float)
    PARSER.add_argument("-t", "--type", help="The lattice type, 'cubic' or 'al2cu'",
                        default="al2cu")
    PARSER.add_argument("-r", "--ratio", help="Ratio of A to B particles. This argument is ignored"
                                              "for Al2Cu lattices.", default=4, type=int)
    ARGS = vars(PARSER.parse_args())

    main(ARGS)
