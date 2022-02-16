"# GridWizard" 
GridWizard is a python class that create a 2D-grid that can be automatically refined around an obstacle. This grid is adapted for 2D diffusion solving, for example solving the thermodynamic equilibrium in the medium surrounding the obstacle.

The repository contains the GridWizard.py script which define the python class and an example of use of this class in a small thermodynamic solving case, with the output as a pdf figure of the result. The thermodynamic case is a wing foiler profile NACA0015 at 200 Celsius degrees embedded in a 50 degrees fluid (something like air). The simulation calculate the thermodynamic equilibrium when considering only the thermic diffusion of the fluid. The aim of this simulation is to show the behaviour of the grid for such a diffusion case.

![view](https://user-images.githubusercontent.com/57091666/154318943-51918562-ae3d-47c7-b66c-8b87590580c2.png)
