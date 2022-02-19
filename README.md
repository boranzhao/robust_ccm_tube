# Tube-Certified Trajectory Tracking With Robust Control Contraction Metrics (CCM)


## Codes for the paper 

P. Zhao, A. Lakshmanan, K. Ackerman, A. Gahlawat, M. Pavone, and N. Hovakimyan, “[Tube-certified trajectory tracking for nonlinear systems with robust control contraction metrics](https://arxiv.org/abs/2109.04453),” IEEE Robotics and Automation Letter, 2022. arXiv:2109.04453.

## Dependencies
- [YALMIP](https://yalmip.github.io/) + [Mosek](https://yalmip.github.io/solver/mosek/) solver for search of CCM or robust CCM (RCCM) using sum of squares (SOS) relaxation.
- [SPOTLESS](https://github.com/spot-toolbox/spotless) also needed for CCM/RCCM synthesis for the 3D quadrotor example
- [OptimTraj](https://github.com/MatthewPeterKelly/OptimTraj) for planning trajectories.
- [OPTI](https://inverseproblem.co.nz/OPTI/) + Matlab `fmincon` slover for solving the nonlinear programming (NLP) problem associated with geodesic computation. 


## Usage 
For each example, 
- Run `main.m` under the `metric` folder to design the CCM/RCCM controller.
- Run `main.m` under the `sim` folder to simulate the behavior of the system under a CCM/RCCM controller. 
- Accelerate the geodesic and control law computation by generating C codes with `generate_code_for_geodesic_cal.m`.

## Suggestions
SOS programming seems not reliable for CCM/RCCM synthesis for high dimensional systems (e.g., a 3D quadrotor). If you struggle a lot in getting a CCM/RCCM for your system, we suggest to try gridding the state space and solving an LMI problem instead of solving a SOS problem, as we did for the quadrotor example.


