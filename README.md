# Introduction

This c program can be used to evolve mildly-hierarchical three body systems under secular approximation. It can evolve a three body system in which a central object of mass (m=1) is being orbited by a test particle(mass=0) and a perturber of mass m2. While the orbit of the perturber is constant, orbital elements of the test particle are evolved using a double averaged hamiltonian. Please see Bhaskar et al. in prep for more details about the integration methods and the reliability of the double average method. This work is supported by NASA grant: 80NSSC20K0641.

# Compile

The main program file is called  "gda_three_body.c". It can be compiled using:

		gcc gda_three_body.c -lm -lgsl -lgslcblas -o gda3bd
                
This program uses GNU Scientific Library(GSL). You need to install it to compile and run the program. See https://www.gnu.org/software/gsl/ for details. 

# Run

The compiled binary can be run using:

		./gda3bd input_file 1 output_file integration_time
                
You may try “./gda3bd IC1.txt 1 outputfile.txt 100” for a simple example.
		
# Input file		

The input file contains the initial orbital parameters of the test particle and the perturber in the following format:

---
HEADER
        
a2	e2	i2	omega2	Omega2	m2
        
HEADER
        
a1	e1	i1	omega1	Omega1

---
The first and the third lines are headers which are used to improve readability and are not parsed by the program. The second line contains the orbital elements of the perturber. The first argument is the semi-major axis (in AU) followed by the eccentricity, inclination, argument of pericenter, longitude of ascending node and mass (in solar masses). They can be separated by spaces or tabs. All angles should be in radians. Line 4 contains the orbital elements of the test particles. Please see file "IC1.txt" for an example. 


# Output file

The format for the output file is given below:

---
HEADER
        
a1	a2	e2	i2	omega2	Omega2	m2
        
HEADER
        
time,e,i,omega,Omega,H,dmin
        
		:
                
		:
                
time,e,i,omega,Omega,H,dmin
        
---
The first and the third lines are headers. The second line contains the configuration of the system. It includes the semi-major axis of the test-particle, as well as the semi-major axis, eccentricity, inclination, argument of pericenter, longitude of ascending node and the mass of the pertuber. The units are the same as those included in the input file. Output of the simulation is printed from line 4 till the end of the file. Each line contains the time, eccentricity, inclination, argument of pericenter and longitude of ascending node of the test particle as well as the value of the hamiltonian and the closest distance between the orbits of the test-particle and the perturber.  

## Plotting output
We provide a python file "plot_traj.py" to plot results of the simulation. To use it you need to install python3 and matplotlib. Use the following command to run the python file:

python3 plot_traj.py <output filename>
        
Plots will be stored in the file named "trajectory.pdf"



# Descriptions of Other Functions

1. gda_three_body.c
  &nbsp;&nbsp;&nbsp;&nbsp;The main file which contains all the code to calculate hamiltonian and it's derivatives.
  
2. gHenon.c
  &nbsp;&nbsp;&nbsp;&nbsp;This file contains code to implement henon's method. This is used to resolve orbital crossings.
  
3. gSL2cquad.c
  &nbsp;&nbsp;&nbsp;&nbsp;This file contains code to calculate double integral used in hamiltonian calculation. It uses GSL's CQUAD method.
  
4. gJmatrix.c
  &nbsp;&nbsp;&nbsp;&nbsp;This file contains code to calculate jacobian of the Hamiltonian.
  
5. gStack3.c
  &nbsp;&nbsp;&nbsp;&nbsp;This file contains code to implement timestep mechanism during crossings.
  
6. plot_traj.py
  &nbsp;&nbsp;&nbsp;&nbsp;Python file to plot the trajectories.

# Acknowledgement
When you use this code or parts of this code, we would greatly appreciate a citation to Bhaskar et al. in prep. This work is supported by NASA grant: 80NSSC20K0641.
