# Tube-Certified Trajectory Tracking  With Robust Control Contraction Metrics (CCM)


## Codes for the paper 

P. Zhao, A. Lakshmanan, K. Ackerman, A. Gahlawat, M. Pavone, and N. Hovakimyan, “[Tube-certified trajectory tracking for nonlinear systems with robust control contraction metrics](https://arxiv.org/abs/2109.04453),” arXiv preprint arXiv:2109.04453, 2021

## Dependencies
- [YALMIP](https://yalmip.github.io/) + [Mosek](https://yalmip.github.io/solver/mosek/) solver for search of CCM or robust CCM (RCCM) using sum of squares (SOS) relaxation.
- [OptimTraj](https://github.com/MatthewPeterKelly/OptimTraj) for planning trajectories.
- [OPTI](https://inverseproblem.co.nz/OPTI/) + Matlab `fmincon` slover for solving the nonlinear programming (NLP) problem associated with geodesic computation. 


## Usage 
- Run `main_synthesis.m` for desigining the CCM/RCCM controller.
- Run `main_simulate.m` for simulating the designed controller. 
- Accelerate the geodesic computation by using generated C codes
  - Use `generate_c_ccm.m` to generate C codes for a controller right after designing it. 
  - Copy the `.mexw64` files from a specific folder storing the pre-designed controller, e.g., `ccm_0.8_plim_0.33pi` to the `pvtol` folder for simulation. 

