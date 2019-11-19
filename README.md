# 1D Heat Conduction Code

This repository contains the Python code to solve the 1D Heat conduction equations with a new model for nano-thermite combustion.

# Solver Details
-Vertex-centred, finite volume method, Explicit

-2nd order central differences for diffusion flux; harmonic or linear interpolation at control surfaces

-Solve heat conduction equations (Heat model) or nano-thermite model (Species model)

-Customizable specific heat capacity based on reaction progress (Arrhenius source term), temperature or a constant


-Customizable thermal conductivity models and calculation methods


-Boundary Conditions: Constant, Flux or convective


-Run from command prompt; must be in parallel with equal nodes in each process

-Can restart a simulation using variable data from previous run

# Heat Model
-Uniform heat generation or exponential (Arrhenius) source term options

# Species Model (In Development)
-Two species tracker of solid material with a fraction of it converted to a gas phase


-Exponential (Arrhenius) source term for heat generation and species conversion rate


-Darcy's law implemented in convective terms with permeability calculated via Carmen-Kozeny

# Run code from **cmd**
cd [directory]

mpiexec -n [proc] python main.py [input file name] [Output directory]

where:

[directory]-is the directory where python files for code are located

[proc]-number of processors used; must be an even number

[input file name]-name and relative path of input file including extension in name (.txt files have been tested)

[Output directory]-relative path to directory to save data files to; will create if non-existent

# Post-processing data
-post-processing script outputs Temperature, reaction progress, reaction rate and species distribution

python Post.py [input file name]

where:

[input file name]-relative path to Post.py input file
