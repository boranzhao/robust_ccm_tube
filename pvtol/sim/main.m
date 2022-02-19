%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation of planning and control of a planar quadrotor with (robust) CCM

%  Author: Pan Zhao, UIUC, Advanced Controls Research Lab,
%  panzhao2@illinois.edu
%  Codes for the paper:
%  P. Zhao, et al. Tube-certified trajectory tracking for nonliner systems
%  with robust control contraction metrics. IEEE Robotics and
%  Automation Letters, 2022. 
%  Last update: Feb 18, 2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
addpath plot_pvtol trajGen_pvtol trajGen_pvtol/nominal_trajs

%% simulation settings
controller_type = 'rccm';            % {'open-loop','ccm','rccm','rccm-p'} 
% ---------------------- load plant and controller ------------------------
rmpath ccm_0.8_plim_0.33pi rccm_1.4_wmax_1_plim_0.33pi rccm_1.2_wmax_1_plim_0.33pi_pos
if strcmp(controller_type, 'ccm')
    file_controller = 'ccm_0.8_plim_0.33pi.mat'; 
    addpath ccm_0.8_plim_0.33pi
elseif strcmp(controller_type, 'rccm') 
    file_controller = 'rccm_1.4_wmax_1_plim_0.33pi.mat';
    addpath rccm_1.4_wmax_1_plim_0.33pi
elseif strcmp(controller_type, 'rccm-p') 
    file_controller = 'rccm_1.2_wmax_1_plim_0.33pi_pos.mat';
    addpath rccm_1.2_wmax_1_plim_0.33pi_pos
else    
    error('Please choose a controller');
end
load(file_controller);
%  -----------------------simulation settings -----------------------------
sim_config.include_tube = 1;        % whether to tighten the state bounds in planning a nominal trajectory
sim_config.tight_input_bnd = 1;     % whether to tighten the input bounds in planning a nominal trajectory
sim_config.include_obs = 0;         % whether to include the obstacles
sim_config.include_dist = 0;        % include the disturbance  
sim_config.save_sim_rst = 0;        % whether to save simulation results
sim_config.replan_nom_traj = 0;     % whether to replan a trajectory

use_generated_code = 1;             % whether to use the generated codes for simulations: using generated codes can accelerate by at least one fold

n = 6; nu = 2;
x0 = zeros(6,1);                    % initial state
xF = [10 10 0 0 0 0]';              % final state
duration = 13;                      % (estimated) time 
umax = 3*plant.m*plant.g;           % control limit
% ----- bounds for input and states for using OptimTraj to plan trajs.-----
u_bnd = [0 0; umax umax]';
x_bnd = [-inf -inf -state_set.p_lim -state_set.vx_lim, -state_set.vz_lim, -state_set.pd_lim;
          inf  inf  state_set.p_lim  state_set.vx_lim   state_set.vz_lim   state_set.pd_lim]';

% --------------------- wind disturbance settings -------------------------
w_max = 1;                          % maximum amplitude of wind disturbance
T_w = 10;                           % period of disturbance
dist_config.sim_config.include_dist = sim_config.include_dist;
dist_config.w_max = w_max;
dist_config.T_w = T_w;
dist_config.gen_dist= @(t) w_max*(0.8+0.2*sin(2*pi/dist_config.T_w.*t));

if ~isfield(controller,'tube_gain_u')
    controller.tube_gain_u = controller.u_bnd/state_set.w_lim;
end
tube_u = controller.tube_gain_u*w_max;

if w_max >1
    error('RCCM controllers were not designed for such larger disturbance!');
end
if sim_config.tight_input_bnd == 1
    u_bnd = u_bnd +[0 0; -tube_u -tube_u]';    
end

%% Plan or load a nominal trajecotory 
gen_nominal_traj
controller.x_nom_fcn = x_nom_fcn;
controller.u_nom_fcn = u_nom_fcn;
controller.w_nom = 0;  % nominal value for disturbances

%% Formulate the NLP problem for geodesic computation
controller.use_generated_code = use_generated_code;
lambda = controller.lambda;
%  problem setting for geodesic computation
D = 2;      % degree of the polynomial
N = D+6;    % stopping index for the CGL (Chebyshev-Gauss-Lobatto) nodes: #nodes N+1

% --------------- re-generate the code is necessary after change of ------
% controller or geodesic optimization settings: remember to re-generate the
% m-file functions, e.g., dW_dphi, dW_dvx, etc., first.------------------  
if controller.use_generated_code 
    load geodesic_setting_for_codegen.mat
    if D ~= geodesic_setting_for_codegen.D ||  N ~= geodesic_setting_for_codegen.N
       error('The geodesic setting here does not match the one used for code generation.'); 
    end
    answer = questdlg('Are the generated codes for this particular scenario?','Question for using C-code in simulation','Yes','No','No');
    switch answer 
        case 'Yes'            
        case 'No'
            error('You cannot continue without including the generated codes for this scenario!');
    end
end
controller.geodesic = set_opt_prob_for_geodesic_computation(n,D,N,controller);

%% simulate
dist0 = norm(x0); 
% --------for additionally outputing control inputs and Reim. energy-------
% compute the initial Riemann energy function value
ue = ccm_law(0,x0,zeros(2,1),x0,plant,controller);
xu0 = [x0;controller.u_nom_fcn(0);ue(end)]; % state, input, energy
Klp = 500*[1 1 1]';
tic;
% ode23 is fastest, followed by ode45
[times,xuTraj] = ode23(@(t,xu) pvtol_dyn(t,xu,Klp,plant,controller,dist_config),[0 duration],xu0); %,ode_opts)
toc;
xuTraj = xuTraj';
xTraj = xuTraj(1:n,:);
uTraj = xuTraj(n+1:end-1,:);
energyTraj = xuTraj(end,:); % Riem. Energy

%% plot the result
plot_trajs;

%% Save data
if sim_config.save_sim_rst == 1
    if controller.type == CtrlDesignOpts.ccm    
        file_name = 'sim_ccm';  
    else
        if  size(plant.C,1)>=6
            file_name = 'sim_rccm';
        else
            file_name = 'sim_rccm_pos';
        end
    end
    file_name = [file_name '_lam_' num2str(controller.lambda,2)];
    if dist_config.sim_config.include_dist == 1
        file_name = [file_name '_w_dist_' num2str(w_max)];
    end
    
    if sim_config.include_obs == 1   
        file_name  = [file_name '_w_obs.mat'];
    else
        file_name  = [file_name '.mat'];
    end
    save(file_name,'times','xTraj','uTraj','xnomTraj','unomTraj','energyTraj','dist_config','sim_config','plant','controller');
end

%% PVTOL dynamics
function dxdt = pvtol_dyn(t,xu,Klp,plant,controller,dist_config)
n = plant.n;
x = xu(1:n);
ue_filter = xu(n+1:end);
% tic;
xnom = controller.x_nom_fcn(t);
unom = controller.u_nom_fcn(t);
[u,Erem] = ccm_law(t,xnom,unom,x,plant,controller);
% toc;
wt = dist_config.gen_dist(t);
dxdt = [plant.f_fcn(x); -Klp.*ue_filter]+[plant.B_fcn(x)*u;  Klp.*[u;Erem]];
if dist_config.sim_config.include_dist == 1
   dxdt(1:n,:) = dxdt(1:n,:) + plant.Bw_fcn(x).*wt;
end
end

