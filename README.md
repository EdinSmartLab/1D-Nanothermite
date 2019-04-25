# 1D Heat Conduction Code

This repository contains the Python code to solve the 1D Heat conduction equations. Based on 2D conduction code architechture.

# Current state:
-solve 1D Planar heat conduction with uniform heat generation (Explicit)

-Combustion source term from Kim (Explicit)

-can be run from command prompt and must be run in parallel

-can restart a simulation using Temperature data from previous run

-post-processing script outputs Temperature, reaction progress and reaction rate contours

# Run code from cmd:
cd [directory]

mpiexec -n [proc] python main.py [input file name] [Output directory]

where:

[proc]-number of processors used; must be an even number

[input file name]-name of input file including extension in name (.txt files have been tested); based on current directory

[Output directory]-directory to save data files to; based on current directory; will create if non-existent

# Post-processing data:
python Post.py [Output directory] [1D graphs]

where:

[Output directory]-directory where data files are stored
[1D graphs] indicates whether 1D graphs should be output (1 or 0); default is 0