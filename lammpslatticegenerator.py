from sys import argv, exit

filetype = 0  # if 0 write lammps output, if 1 write standard xyz output

if filetype == 0:
    outputfile = open("startconfiguration.xyz", "w")
else:
    outputfile = open("xyzstartconfiguration.xyz", "w")

if len(argv) != 5:
    print("Incorrect syntax. Use: arguments lattice_size density, lattice_type (0 for Al2Cu 1 for cubic) and particle_ratio (2 for 2:1, 3 for 3:1 etc). e.g: ./lammpslatticegenerator.py 7 1.5 0 2")
    print("Note that ratio will only be used for cubic lattices. Al2Cu lattices are always 2:1 ratio.")
    exit()

cubicboxsize = int(argv[1])
density = float(argv[2])

if int(argv[3]) == 0:  # if setting up Al2Cu lattice
    unitcelldimensions = [6.067, 6.067, 4.877]
    numberofparticles = 12 * cubicboxsize ** 3
else:        # if using cubic lattice
    unitcelldimensions = [1, 1, 1]
    numberofparticles = cubicboxsize ** 3

volume = numberofparticles / density
sidelenratio = unitcelldimensions[0] / unitcelldimensions[2]
sidez = (volume / (sidelenratio * sidelenratio)) ** (1.0 / 3)
sidexy = sidez * sidelenratio
unitcellsidexy = sidexy / cubicboxsize
unitcellsidez = sidez / cubicboxsize
lammpsprefix = ""

print("unit cell side x + y = " + str(unitcellsidexy))
print("unit cell side z = " + str(unitcellsidez))

i = 1
j = 1
k = 1
l = 1
particlecounter = 1
particletype = 1
lineout = ""

al2cucoordarray = [[0.65810, 0.15810, 0.50000], [0.15810, 0.65810, 0.00000], [0.84190, 0.34190, 0.00000], [0.34190, 0.84190, 0.50000], [0.15810, 0.34190, 0.50000], [0.65810, 0.84190, 0.00000], [.84190, 0.65810, 0.50000], [0.34190, 0.15810, 0.00000], [0.50000, 0.50000, 0.75000], [0.00000, 0.00000, 0.25000], [0.50000, 0.50000, 0.25000], [0.00000, 0.00000, 0.75000]]

if filetype == 1:
    outputfile.write(str(numberofparticles) + "\n")

outputfile.write("Auto generated 11A crystal lattice for LAMMPS\n")

if filetype == 0:
    outputfile.write("\n" + str(numberofparticles) + " atoms\n\n")
    outputfile.write("2 atom types\n\n")
    outputfile.write('%.7f' % (-0.5 * sidexy) + " " + '%.7f' % (sidexy * 0.5) + " xlo xhi\n")
    outputfile.write('%.7f' % (-0.5 * sidexy) + " " + '%.7f' % (sidexy * 0.5) + " ylo yhi\n")
    outputfile.write('%.7f' % (-0.5 * sidez) + " " + '%.7f' % (sidez * 0.5) + " zlo zhi\n\n")
    outputfile.write("Masses\n\n1 1\n2 1\n\nAtoms\n\n")

if int(argv[3]) == 0:   # Al2Cu lattice
    for i in range(0, cubicboxsize):
        for j in range(0, cubicboxsize):
            for k in range(0, cubicboxsize):
                for l in range(0, 12):
                    lammpsprefix = str(particlecounter) + "\t"
                    if l > 7:
                        particletype = 2
                    else:
                        particletype = 1
                    lineout = str(particletype) + "\t" + '%.7f' % (al2cucoordarray[l][0] * unitcellsidexy + i * unitcellsidexy - 0.5 * sidexy) + "\t" + '%.7f' % (al2cucoordarray[l][1] * unitcellsidexy + j * unitcellsidexy - 0.5 * sidexy) + "\t" + '%.7f' % (al2cucoordarray[l][2] * unitcellsidez + k*unitcellsidez - 0.5 * sidez) + "\n"
                    if filetype == 0:
                        outputfile.write(lammpsprefix + lineout)
                    else:
                        outputfile.write(lineout)
                    particlecounter += 1
else:   # Cubic lattice
    for i in range(0, cubicboxsize):
        for j in range(0, cubicboxsize):
            for k in range(0, cubicboxsize):
                if particlecounter % (int(argv[4]) + 1) == 0:
                    particletype = 2
                else:
                    particletype = 1
                lammpsprefix = str(particlecounter) + "\t"
                lineout = str(particletype) + "\t" + '%.7f' % (i * unitcellsidexy - 0.5 * sidexy) + "\t" + '%.7f' % (j * unitcellsidexy - 0.5 * sidexy) + "\t" + '%.7f' % (k * unitcellsidez - 0.5 * sidez) + "\n"
                if filetype == 0:
                    outputfile.write(lammpsprefix + lineout)
                else:
                    outputfile.write(lineout)
                particlecounter += 1

outputfile.close()