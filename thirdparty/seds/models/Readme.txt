recorded_motions: contains a library of 24 handwriting motions recorded from 
	          Tablet-PC. Demonstrations are saved as '.mat' file and 
		  contains two variables:

            o dt: time steps which is 0.02 second for all demonstrations
            o demos: A variable containing all demonstrations (only
              trajectories). The variable 'demos' follows the following
              format:
              - demos{n}: d x T^n matrix representing the d dimensional
                          trajectories. T^n is the number of datapoint in
                          this demonstration (1 < n < N)

SEDS_models: cotains the same models as the folder 'recorded_motions' plus learnt
             models of demonstrations using SEDS library.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Copyright (c) 2010 S. Mohammad Khansari-Zadeh, LASA Lab, EPFL,   %%%
%%%          CH-1015 Lausanne, Switzerland, http://lasa.epfl.ch         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This library of motion is free for non-commercial academic use. The library 
must not be modified or distributed without prior permission of the authors.
Please acknowledge the authors in any academic publications that have made 
use of this library or part of it. Please use this BibTex for reference:

   S. M. Khansari Zadeh and A. Billard, "Imitation learning of Globally 
   Stable Non-Linear Point-to-Point Robot Motions using Nonlinear
   Programming", in Proceeding of the 2010 IEEE/RSJ International
   Conference on Intelligent Robots and Systems (IROS 2010), Taipei,
   Taiwan, October 2010

To get latest upadate of the software please visit
                          http://lasa.epfl.ch/khansari

Please send your feedbacks or questions to:
                          mohammad.khansari_at_epfl.ch
