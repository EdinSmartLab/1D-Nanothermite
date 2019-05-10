# 1D Heat Conduction Code

This repository contains the Python code to solve the 1D Heat conduction equations. Based on 2D conduction code architechture.

# Current state:
-solve 1D heat conduction with uniform heat generation (Explicit)

-Combustion source term from Kim (Explicit)

-two species model with Darcy's law, species tracker (gas and solid only); TESTING PHASE

-can be run from command prompt and must be run in parallel

-can restart a simulation using variable data from previous run

-post-processing script outputs Temperature, reaction progress, reaction rate and species data

# Run code from cmd:
cd [directory]

mpiexec -n [proc] python main.py [input file name] [Output directory]

where:

[proc]-number of processors used; must be an even number

[input file name]-name of input file including extension in name (.txt files have been tested); based on current directory

[Output directory]-directory to save data files to; based on current directory; will create if non-existent

# Post-processing data:
python Post.py [Output directory]

where:

[Output directory]-directory where data files are stored
