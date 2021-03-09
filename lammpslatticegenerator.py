import argparse
from typing import List, TextIO


def main(args: dict):

    if args["type"] == "al2cu":
        unit_cell_dimensions = [6.067, 6.067, 4.877]
        num_particles = 12 * args["lattice_size"] ** 3
    else:
        # cubic lattice
        unit_cell_dimensions = [1, 1, 1]
        num_particles = args["lattice_size"] ** 3

    volume = num_particles / args["density"]
    side_len_ratio = unit_cell_dimensions[0] / unit_cell_dimensions[2]
    side_z = (volume / (side_len_ratio * side_len_ratio)) ** (1.0 / 3)
    side_xy = side_z * side_len_ratio
    unit_cell_side_xy = side_xy / args["lattice_size"]
    unit_cell_side_z = side_z / args["lattice_size"]

    print(f"unit cell side x + y = {unit_cell_side_xy}")
    print(f"unit cell side z = {unit_cell_side_z}")

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
    if args["output_type"] == "lammps":
        file_name = "out.lmp"
    else:
        file_name = "out.xyz"

    with open(file_name, 'w') as output_file:
        write_header(args["output_type"], num_particles, output_file, side_xy, side_z)

        if args["type"] == "al2cu":
            write_al2cu(al2cucoordarray, args, output_file, side_xy, side_z, unit_cell_side_xy,
                        unit_cell_side_z)
        else:   # Cubic lattice
            write_cubic(args, output_file, side_xy, side_z, unit_cell_side_xy, unit_cell_side_z)


def write_cubic(args: dict, output_file: TextIO, side_xy: float, side_z: float,
                cell_side_xy: float, cell_side_z: float):
    particle_counter = 1
    for x in range(args["lattice_size"]):
        for y in range(args["lattice_size"]):
            for z in range(args["lattice_size"]):
                if args["output_type"] == "lammps":
                    output_file.write(f"{particle_counter}\t")

                if particle_counter % (args["ratio"] + 1) == 0:
                    particle_type = 2
                else:
                    particle_type = 1

                x_pos = x * cell_side_xy - 0.5 * side_xy
                y_pos = y * cell_side_xy - 0.5 * side_xy
                z_pos = z * cell_side_z - 0.5 * side_z
                output_file.write(f"{particle_type}\t{x_pos:.7f}\t{y_pos:.7f}\t{z_pos:.7f}\n")
                particle_counter += 1


def write_al2cu(al2cu_coords: List[List[int]], args: dict, output_file: TextIO, side_xy: float,
                side_z: float, unitcellside_xy: float, unitcellside_z: float):
    particle_counter = 1
    for x in range(args["lattice_size"]):
        for y in range(args["lattice_size"]):
            for z in range(args["lattice_size"]):
                for unit_cell_counter in range(12):
                    if args["output_type"] == "lammps":
                        output_file.write(f"{particle_counter}\t")

                    if unit_cell_counter > 7:
                        particle_type = 2
                    else:
                        particle_type = 1

                    x_pos = al2cu_coords[unit_cell_counter][0] * unitcellside_xy + x * unitcellside_xy - 0.5 * side_xy
                    y_pos = al2cu_coords[unit_cell_counter][1] * unitcellside_xy + y * unitcellside_xy - 0.5 * side_xy
                    z_pos = al2cu_coords[unit_cell_counter][2] * unitcellside_z + z * unitcellside_z - 0.5 * side_z
                    output_file.write(f"{particle_type}\t{x_pos:.7f}\t{y_pos:.7f}\t{z_pos:.7f}\n")
                    particle_counter += 1


def write_header(output_type: str, num_particles: int, output_file: TextIO, side_xy: float,
                 side_z: float):
    if output_type == "xyz":
        write_xyz_header(num_particles, output_file)
    if output_type == "lammps":
        write_lammps_header(num_particles, output_file, side_xy, side_z)


def write_lammps_header(num_particles: int, output_file: TextIO, side_xy: float, side_z: float):
    output_file.write(f"\n{num_particles} atoms\n\n")
    output_file.write("2 atom types\n\n")
    output_file.write(f"{-0.5 * side_xy:.7f} {side_xy * 0.5:.7f} xlo xhi\n")
    output_file.write(f"{-0.5 * side_xy:.7f} {side_xy * 0.5:.7f} ylo yhi\n")
    output_file.write(f"{-0.5 * side_z:.7f} {side_z * 0.5:.7f} zlo zhi\n\n")
    output_file.write("Masses\n\n1 1\n2 1\n\nAtoms\n\n")
    output_file.write("Auto generated crystal lattice for LAMMPS\n")


def write_xyz_header(num_particles: int, output_file: TextIO):
    output_file.write(f"{num_particles}\n")
    output_file.write("Auto generated crystal lattice\n")


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
