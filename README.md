# ds-opt
Toolbox including optimization techniques for estimation of Globally Asymptotically Stable Dynamical Systems focused on (1) Linear Parameter Varying formulation with GMM-based mixing function and different Lyapunov candidate functions as proposed in [1,2]. 
For comparison purposes, this toolbox also includes implementations and demo scripts for DS learning with SEDS [3] and the diffeomorphic matching approach [4].

<p align="center">
  <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Lshape_lpvO3.png" width="220">
<img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Ashape_lpvO3.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Sshape_lpvO3.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Ashape_lpvO3.png" width="220">
</>
  
### Installation Instructions
This package needs the **physically-consisent** GMM fitting proposed in [1] and implemented in [phys-gmm](https://github.com/nbfigueroa/phys-gmm.git). If you do not already have this package, you can download it as a submodule. After cloning this repo, one must initialize/download the submodule with the following commands:
```
cd ~./ds_opt
git submodule init
git submodule update
```
In case you want to update the submodule to its latest version, you can do so with the following command:
```
git submodule update --remote
```
### Instructions and Content


##### Running the demo scripts
There are three important demo scripts:

### Example Datasets
These examples + more datasets are provided in ```
./datasets``` folder.


**References**     
> [1] Figueroa, N. and Billard, A. (2018) A Physically-Consistent Bayesian Non-Parametric Mixture Model for Dynamical System Learning. In Proceedings of the 2nd Conference on Robot Learning (CoRL). Accepted.   
> [2] Mirrazavi Salehian, S. S. (2018) Compliant control of Uni/ Multi- robotic arms with dynamical systems. PhD Thesis.   
> [3] Khansari Zadeh, S. M. and Billard, A. (2011) Learning Stable Non-Linear Dynamical Systems with Gaussian Mixture Models. IEEE Transaction on Robotics, vol. 27, num 5, p. 943-957.    
> [4] N. Perrin and P. Schlehuber-Caissier, “Fast diffeomorphic matching to learn globally asymptotically stable nonlinear dynamical systems,” Systems & Control Letters (2016).

**Contact**: [Nadia Figueroa](http://lasa.epfl.ch/people/member.php?SCIPER=238387) (nadia.figueroafernandez AT epfl dot ch)
