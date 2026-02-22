%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation of planning and control of a 3-D quadrotor with (robust) CCM

%  Author: Pan Zhao, University of Alabama
%  pan.zhao@ua.edu; boranzhao9@gmail.com
  
%  Last update: Feb 22, 2026.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
addpath quad_aux_files quad_plotting quad_traj_opt create_env;
addpath('../metric');
%% 
controller_name = 'ccm';            % {'open-loop','ccm','rccm','rccm-p'} 
print_file = 0; 

%% simulation settings 
sim_config.tight_input_bnd = 0;     % whether to tighten the input bounds in planning a nominal trajectory
sim_config.include_dist = 1;        % include the disturbance  
sim_config.save_sim_rst = 0;        % whether to save simulation results
sim_config.replan_nom_traj = 0;     % whether to replan a trajectory
sim_config.dt_sim = 1/250;          % sampling interval for controller 
use_generated_code = 0;             % whether to use the generated codes for simulations: using generated codes can accelerate by at least ten fold

%% wind disturbance settings
w_max = 1;                         % maximum amplitude of wind disturbance
T_w = 10;                          % period of disturbance
dist_config.sim_config.include_dist = sim_config.include_dist;
dist_config.w_max = w_max;
dist_config.T_w = T_w;
dist_config.wind_ang = pi/4;
dist_config.gen_dist= @(t) w_max*(0.8+0.2*sin(2*pi/dist_config.T_w.*t));

%% load dynamics, controller and setup simulation environment
% ---------------------- load plant and controller ------------------------
if strcmp(controller_name, 'ccm')
    file_controller = 'ccm_2.8.mat';  % 3.4   
elseif strcmp(controller_name, 'rccm') 
    file_controller = 'rccm_4.4.mat';
elseif strcmp(controller_name, 'rccm-p') 
    file_controller = 'rccm_3.4_pos.mat';
else    
    error('Please choose a controller');
end
load(file_controller);
controller.dW_dxi_fcn = @(i,x) (i==7)*dW_dT(x)+(i==8)*dW_dphi(x)+(i==9)*dW_dtheta(x);

%set to 1 to generate new obs
new_obs = 0;
load_quad_config;

%% Load solver/trajectory generator
%set to 1 to  generate new path and trajectory (set 1 if new_obs = 1)
gen_new_traj = new_obs;
if strcmp(controller_name, 'ccm')
    file_quad_traj = 'quad_traj_opt/quad_traj_0.58.mat';     
elseif strcmp(controller_name, 'rccm-p') 
    file_quad_traj = 'quad_traj_opt/quad_traj_0.32.mat';  
else
    error('Please select the right controller');
end
load_solvers;

controller.MP_t = MP_t;
controller.MP_state = MP_state;
controller.MP_ctrl = MP_ctrl;
controller.x_nom_fcn = @(t) interp1(MP_t,MP_state,t);
controller.u_nom_fcn = @(t) interp1(MP_t,MP_ctrl,t);

controller.use_generated_code = use_generated_code;

%% Formulate the NLP problem for geodesic computation
D = 3;      % degree of the polynomial
N = 6;      % stopping index for the CGL (Chebyshev-Gauss-Lobatto) nodes: #notes N+1
geodesic = setup_geo_opt_quad(plant.nc,D,N,controller);
geodesic.nlprob = []; % setting geodesic.nlprob to [] indicates that the geodesic will be approximated by a straight line 

controller.geodesic = geodesic;
controller.w_nom = 0;  % nominal value for disturbances

% simulate
t_vec = 0:sim_config.dt_sim:controller.MP_t(end);
T_steps = length(t_vec)-1;

% note that 
% x: pos,vel, roll, pitch, yaw, om; u: T, moments
% xc: pos, vel, thrust, roll, pitch; 
% uc: thrust_dot, roll_dot, pitch_dot,yaw_dot

x = test_state;
xc = x2xc(x,plant.g);

xTraj = zeros(plant.n,T_steps+1);
xcTraj = zeros(plant.nc,T_steps+1);
ucTraj = zeros(plant.nuc,T_steps);
xcnomTraj = xcTraj;
ucnomTraj = ucTraj;
energyTraj = zeros(1,T_steps);

xTraj(:,1) = x;
xcTraj(:,1) = xc;

plant.n = plant.nc; % for use in the function "ccm_law".
for i=1:T_steps
    t = t_vec(i);   
    if (mod(i,200)==0)
        fprintf('t = %.2f s\n',t);
    end
    xc = x2xc(x,xc(7)); % xc(7) is thrust;
    
    % get the nominal state and input
    xc_nom = controller.x_nom_fcn(t); xc_nom = xc_nom(1:plant.nc)';
    uc_nom = controller.u_nom_fcn(t); uc_nom = uc_nom(1:plant.nuc)';   
    [uc,energy] = ccm_law(t,xc_nom,uc_nom,xc,plant,controller);
    
    % record
    xcnomTraj(:,i) = xc_nom;
    ucTraj(:,i) = uc; 
    ucnomTraj(:,i) = uc_nom;
    energyTraj(i) = energy;
        
    % propagate with zero-hold for control inputs
    euler_dot_des = [uc(2);uc(3);0]; % set desired yaw_dot to constant zero
    [d_t,d_state] = ode113(@(t,x) quad_dyn(t,x,xc(7),euler_dot_des,plant,dist_config),[t_vec(i) t_vec(i+1)],x); %,ode_options
    
    % update and record;
    x = d_state(end,:)';
    xc(7) = xc(7) + uc(1)*sim_config.dt_sim; 
    
    xTraj(:,i+1) = x;
    xcTraj(:,i+1) = xc;    
end
xc_nom = controller.x_nom_fcn(t_vec(i+1));xc_nom = xc_nom(1:plant.nc)';
uc_nom = controller.u_nom_fcn(t_vec(i+1));uc_nom = uc_nom(1:plant.nuc)';   
xcnomTraj(:,i+1) = xc_nom;
[~,energy(i+1)] = ccm_law(t_vec(i+1),xc_nom,uc_nom,xc,plant,controller);

%% plot the result
plot_trajs;

%% save_data
if sim_config.save_sim_rst == 1
    if controller.type == CtrlDesignOpts.ccm    
        file_name = 'sim_ccm';  
    else
        if  size(plant.C,1)>=6
            file_name = 'sim_rccm';
        else
            file_name = 'sim_rccm_p';
        end
    end
    file_name = [file_name '_lam_' num2str(controller.lambda,2)];
    if dist_config.sim_config.include_dist == 1
        file_name = [file_name '_w_dist_' num2str(w_max)];
    end    
    file_name  = [file_name '.mat'];
    save(file_name,'t_vec','xTraj','xcTraj','ucTraj','xcnomTraj','ucnomTraj','energyTraj','dist_config','sim_config','plant','controller');
end

%%
function dx_dot = quad_dyn(t,x,thrust,euler_dot_des,plant,dist_config)
global kp_omg
%get moments
omg_des = R_omg(x(7:9))*euler_dot_des;
omg = x(10:12);
M = kp_omg*(omg_des - omg);

%assemble
u = [thrust;M];

%propagate
wdist = dist_config.gen_dist(t);
dx_dot = plant.f_sim(x) + plant.B_sim(x)*u + ...
   (dist_config.sim_config.include_dist == 1)*plant.Bw_sim*...
   (wdist*[sin(dist_config.wind_ang);-cos(dist_config.wind_ang);0]);
end

function x = xc2x(xc,yaw,omg)
% convert the (9) states for CCM control to (12) full states
x = [xc(1:6);xc(8:9);yaw;omg];
end

function xc = x2xc(x,thrust)
% convert the (9) states for CCM control to (12) full states
xc = [x(1:6);thrust;x(7:8)];
end