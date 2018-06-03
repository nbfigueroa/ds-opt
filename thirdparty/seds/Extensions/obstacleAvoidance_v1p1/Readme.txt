Obstacle avoidance Package: version 1.1 issued on March 26, 2012

This packages contains the obstacle avoidance algorithm presented in the
following paper:

     S.M. Khansari Zadeh and A. Billard, "A Dynamical System Approach to
     Realtime Obstacle Avoidance", Autonomous Robots, 2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Copyright (c) 2011 S. Mohammad Khansari-Zadeh, LASA Lab, EPFL,    %%%
%%%          CH-1015 Lausanne, Switzerland, http://lasa.epfl.ch         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The presented source codes contains two MATLAB functions:
'Tutorial_Obstacle_Format.m' and 'Tutorial_Obstacle_Avoidance.m', and two
subdirectories 'lib_obstacle_avoidance' and 'Example_DSs'.

- Tutorial_Obstacle_Format.m: A tutorial on how to define obstacles in the
  format that is consistent with the obstacle avoidance library.

- Tutorial_Obstacle_Avoidance.m: A tutorial including a set of examples on
  how to use the obstacle avoidance library.

- The folder 'lib_obstacle_avoidance' contains three functions:
    1) 'obs_modulation_ellipsoid.m' is the main file and includes the implementation
        of the proposed obstacle avoidance approach.
    2) 'obs_draw_ellipsoid.m' is used to illustrate the obstacle shape.
    3) 'Simulation.m' is a simple simulation function, and is used to
        simulate robot motions. This function is taken from the SEDS library.

- The folder 'Example_DSs' contains a set of DS that are used in the
tutorials.

When running the tutorials, it is assumed that your current directory is the
'sc_obs_AR' directory. Otherwise, you should manually add both the
'lib_obstacle_avoidance' and 'Example_DSs' directories to the MATLAB path.

    ***This source code was written and tested with MATLAB 2011b.***
