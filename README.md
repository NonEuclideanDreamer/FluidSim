# FluidSim
2d- fluid simulation with velocity field

Solving Navier Stokes numerically, figured out the workflow with the help of this paper by Jos Stam: https://pages.cs.wisc.edu/~chaol/data/cs777/stam-stable_fluids.pdf

This is a Java Project consisting of 3 classes. It produces a series of images that can be put together to a video.

The main class is FluidFlow.java, where one can do modifications to change the output.

-In the class head on can change image size, and various parameters. 
  The boolean "fourier" should remain false, I didn't get it to work
  The boolean "viscosity" can be made true, but I only implemented a one-step numerical solution, so visc can't be set too high lest it explodes
  The boolean "obst" specifies wether we have an obstacle with boundary conditions

-In the constructor one can specify the initial state of the colors

-In the main method, one can specify the type of obstacle

-In the method "update" one can specify the force pattern disturbing the fluid


The other two classes handle Vector Fields and Scalar Fields respectively.

If you want to play with it and have questions, feel free to ask, for example in my discord: https://discord.gg/Q2CJfeMJn
