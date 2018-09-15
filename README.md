# ds-opt
Toolbox including optimization techniques for estimation of Globally Asymptotically Stable Dynamical Systems focused on (1) Linear Parameter Varying formulation with GMM-based mixing function and different Lyapunov candidate functions as proposed in [1]. For comparison purposes, this toolbox also includes implementations and demo scripts for DS learning with SEDS [2] and the diffeomorphic matching approach [3].

<p align="center">
  <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Lshape_lpvO3.png" width="220">
<img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Ashape_lpvO3.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Sshape_lpvO3.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Ashape_lpvO3.png" width="220">
</>
  
### Installation Instructions
This package needs the **physically-consisent** GMM (PC-GMM) fitting proposed in [1] and implemented in [phys-gmm](https://github.com/nbfigueroa/phys-gmm.git). If you do not already have this package, you can download it as a submodule. After cloning this repo, one must initialize/download the submodule with the following commands:
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
.. Comments here.. introduce approach..

##### Running the demo scripts
There are three important demo scripts:

### Example Datasets
These examples + more datasets are provided in ```
./datasets``` folder. Following we show some **notably challenging motions** that cannot be accurately encoded with SEDS [3] **(1st Column)** or a PC-GMM-based LPV-DS [1] with a simple Quadradtic Lyapunov Function (QLF) **(2nd Column)**, yet can be correctly encoded with the proposed PC-GMM-based LPV-DS with a parametrized QLF (P-QLF) [1] **(3rd Column)** yielding comparable **(in some cases BETTER)** results than the global diffeomorphic matching approach **(4th Column)**, which is the state-of-the-art approach known to outperform SEDS.

-  **2D S-shape Dataset**  
<p align="center">
  <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Sshape_seds.png" width="220">
<img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Sshape_lpvO1.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Sshape_lpvO3.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Sshape_diff.png" width="215">
</>

-  **2D Messy Snake Dataset**  
<p align="center">
  <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Messy-snake_seds.png" width="215">
<img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Messy-snake_lpvO1.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Messy-snake_lpvO3.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Messy-snake_diff.png" width="225">
</>

-  **2D Multi-Behavior (Single Target) Dataset**  
<p align="center">
  <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Multi_seds.png" width="215">
<img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Multi_lpv01.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Multi_lpv03.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Multi_diff.png" width="220">
</>  
  
-  **2D SharpC-shape from LASA Handwriting Dataset**  
<p align="center">
  <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/CSharp_seds.png" width="220">
<img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/CSharp_lpvO1.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/CSharp_lpvO3.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/CSharp_diff.png" width="220">
</>

-  **2D N-shape from LASA Handwriting Dataset**  
<p align="center">
  <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Nshape_seds.png" width="220">
<img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Nshape_lpvO1.png" width="225"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Nshape_lpvO3.png" width="215"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Nshape_diff.png" width="220">
</>


-  **2D Hee-shape from LASA Handwriting Dataset**  
<p align="center">
  <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Hee_seds.png" width="220">
<img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Hee_lpvO1.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Hee_lpvO3.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/Hee_diff.png" width="220">
</>
  
-  **2D Snake-shape from LASA Handwriting Dataset**  
<p align="center">
  <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/LASASnake_seds.png" width="220">
<img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/LASASnake_lpvO1.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/LASASnake_lpvO3.png" width="220"><img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/LASASnake_diff.png" width="220">
</>
  


**References**     
> [1] Figueroa, N. and Billard, A. (2018) A Physically-Consistent Bayesian Non-Parametric Mixture Model for Dynamical System Learning. In Proceedings of the 2nd Conference on Robot Learning (CoRL). Accepted.     
> [2] Khansari Zadeh, S. M. and Billard, A. (2011) Learning Stable Non-Linear Dynamical Systems with Gaussian Mixture Models. IEEE Transaction on Robotics, vol. 27, num 5, p. 943-957.    
> [3] N. Perrin and P. Schlehuber-Caissier, “Fast diffeomorphic matching to learn globally asymptotically stable nonlinear dynamical systems,” Systems & Control Letters (2016).

**Contact**: [Nadia Figueroa](http://lasa.epfl.ch/people/member.php?SCIPER=238387) (nadia.figueroafernandez AT epfl dot ch)
