from math import pow
from sys import argv
import sys

"""The program lammpslatticegenerator takes 3 arguments: <cubicboxsize> <density> <filelammps>
	cubicboxsize: 	integer, number of unit cells per direction 
	density: 		float, density of the box
	filelammps: 	boolean, whether the file type is lmp or xyz
The program creates a lmp or xyz file with the initial configuration of the crystal
"""

	### IMPORT OF THE COORDINATES FILE ###

if len(argv) == 4:
	filelammps=argv[3]
elif len(argv) == 3:
	filelammps=True			#default value
else:						#in that case, the number of arguments is not good, so an error is raised
    print("Incorrect syntax. Use: arguments 'lattice size' and density. e.g: ./lammpslatticegenerator.py 7 1.5 True "
		  "or ./lammpslatticegenerator.py 7 1.5")
    sys.exit()

#opening of the file
if filelammps:
	print "Writing LMP"
	outputfile = open("startconfiguration.lmp","w")
else:
	print "Writing XYZ"
	outputfile = open("xyzstartconfiguration.xyz","w")

	### SETTINGS ###

unitcelldimensions 	= [6.067, 6.067, 4.877]							#centered cubic cell
cubicboxsize 		= int(argv[1])									#size of the box = 1st arg given
numberofparticles 	= 12 * cubicboxsize**3							#number of particles = number per box * number of cells
density 			= float(argv[2])								#density of the box = 2nd arg given
volume 				= numberofparticles/density						#volume
sidelenratio 		= unitcelldimensions[0]/unitcelldimensions[2]	#aspect ratio of the cell
sidez 				= (volume/(sidelenratio*sidelenratio))**(1.0/3)	#side along z of the box
sidexy 				= sidez * sidelenratio							#side along x and y of the box
unitcellsidez 		= sidez/cubicboxsize							#redefinition of the unit cell's size along Z
unitcellsidexy 		= sidexy/cubicboxsize 							#					"						x and y


#definition of the coordinates of the sites in the unit cell
coordarray = [[0.65810, 0.15810, 0.50000],
			  [0.15810, 0.65810, 0.00000],
			  [0.84190, 0.34190, 0.00000],
			  [0.34190, 0.84190, 0.50000],
			  [0.15810, 0.34190, 0.50000],
			  [0.65810, 0.84190, 0.00000],
			  [0.84190, 0.65810, 0.50000],
			  [0.34190, 0.15810, 0.00000],
			  [0.50000, 0.50000, 0.75000],
			  [0.00000, 0.00000, 0.25000],
			  [0.50000, 0.50000, 0.25000],
			  [0.00000, 0.00000, 0.75000]]


	### PREAMBULE OF THE FILE ###

outputfile.write("Auto generated 11A crystal lattice for LAMMPS\n")

if filelammps:
	outputfile.write(
					 "\n" + str(numberofparticles) + " atoms\n\n"
					 "2 atom types\n\n"
					 '%.7f' % (-0.5 * sidexy) + " " + '%.7f' % (sidexy * 0.5) + " xlo xhi\n"
					 '%.7f' % (-0.5 * sidexy) + " " + '%.7f' % (sidexy * 0.5) + " ylo yhi\n"
					 '%.7f' % (-0.5 * sidez)  + " " + '%.7f' % (sidez  * 0.5) + " zlo zhi\n\n"
	                 "Masses\n\n1 1\n2 1\n\nAtoms\n\n"
					)
else:
	outputfile.write(str(numberofparticles) + "particles\n")


	### POSITIONNING OF THE ATOMS ###

#counters and preliminary variables
particlecounter = 1
number=0
lineout = ""
lammpsprefix = ""
listing=range(cubicboxsize)

#loop over each position of the unit cells in the 3d box
for (i,j,k) in [(i,j,k) for i in listing for j in listing for k in listing]:
	#loop over each site of the cell
	for l in range (12):
		#distinction between particle types
		if l > 7:
			number=2
		else:
			number=1

		#computation of the position of the particle
		lammpsprefix = str(particlecounter) + "\t"
		lineout = "{}".format(number) + "\t" \
				  "%.7f" % ( (coordarray[l][0] + i) * unitcellsidexy - 0.5 * sidexy) + "\t" \
				  "%.7f" % ( (coordarray[l][1] + j) * unitcellsidexy - 0.5 * sidexy) + "\t" \
				  "%.7f" % ( (coordarray[l][2] + k) * unitcellsidez  - 0.5 * sidez ) + "\n"

		#writing of the position of the particle
		if filelammps:
			outputfile.write(lammpsprefix + lineout)
		else:
			outputfile.write(lineout)

		#one particle added
		particlecounter += 1

outputfile.close()

#number of particles displayed
print "The number of particles is", particlecounter-1
