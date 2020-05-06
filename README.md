# Introduction

This c program can be used to evolve mildly-hierarchical three body systems under secular approximation. It can evolve a three body system in which a central object of mass (m=1) is being orbited by a test particle(mass=0) and a perturber of mass m2. While the orbit of the perturber is constant, orbital elements of the test particle are evolved using a double averaged hamiltonian. 

# Compile

The main program file is called  "gda_three_body.c". It can be compiled using:

		gcc gda_three_body.c -lm -lgsl -lgslcblas -o gda3bd
                
This program uses GNU Scientific Library(GSL). You need to install it to compile and run the program. See https://www.gnu.org/software/gsl/ for details. 

# Run

The compiled binary can be run using:

		./gda3bd input_file 1 output_file integration_time
                
You may try “./gda3bd IC1.txt 1 outputfile.txt 100” for a simple example.
		
# Input file:		

The first argument is the input file which contains the initial orbital parameters of the test particle and the perturber(see below). The second argument is the name of the output file to which the program will write it's output. The third argument is the end time(in years). 

This program expects the input file to provide initial conditions for the test particle in the following format:

---
HEADER
        
a2	e2	i2	omega2	Omega2	m2
        
<HEADER>
        
<a1>	<e1>	<i1>	<omega1>	<Omega1>

---
The first and the third lines are headers which are used to improve readability and are not parsed by the program. The second line contains the orbital elements of the perturber. The first argument is the semi-major axis(in AU) followed by e2, i2, argument of pericenter, longitude of ascending node and mass(in solar masses)(e2 is the eccentricity and i2 is the inclination of the perturber).  They can be separated by spaces or tabs. All angles should be in radians. Line 4 contains the orbital elements of the test particles. The first argument is the semi-major axis(in AU) followed by e1, i1, argument of pericenter and longitude of ascending node (e1 is the eccentricity and i1 is the inclination  of the test-particle). Please see the sample input file "IC1.txt". 


# Output file:

The format for the output file is given below:

---
<HEADER>
        
<a1>	<a2>	<e2>	<i2>	<omega2>	<Omega2>	<m2>
        
<HEADER>
        
<time>,<e>,<i>,<omega>,<Omega>,<H>,<dplus>,<dminus>
        
		:
                
		:
                
<time>,<e>,<i>,<omega>,<Omega>,<H>,<dplus>,<dminus>
        
---
The first and the third lines are headers. The second line contains the configuration of the system. It contains semi-major axis of the test-particle, semi-major axis, eccentricity, inclination, argument of pericenter, longitude of ascending node and the mass of the pertuber(in solar masses). From Line 4 till the end of the file, output of the simulation is printed. Each line contains time, eccentricity, inclination, argument of pericenter and longitude of ascending node of the test particle as well as the value of the hamiltonian and the closest distance between the test-particle's orbit and perturber's orbit.  


We provide a python file "plot_traj.py" to plot trajectories from the output file. To use it you need to install python3 and matplotlib. Use the following command to run the python file:

python3 plot_traj.py <output filename>
        
Plots will be stored in the file named "trajectory.pdf"



Following are the brief descriptions of various files in the directory:

1.gda_three_body.c

  The main file which contains all the code to calculate hamiltonian and it's derivatives.
  
2.gHenon.c

  This file contains code to implement henon's method. This is used to resolve orbital crossings.
  
3.gSL2cquad.c

  This file contains code to calculate double integral used in hamiltonian calculation. It uses GSL's CQUAD method.
  
4.gJmatrix.c

  This file contains code to calculate jacobian of the Hamiltonian.
  
5.gStack3.c

  This file contains code to implement timestep mechanism during crossings.
  
6.plot_traj.py

  Python file to plot the trajectories.

