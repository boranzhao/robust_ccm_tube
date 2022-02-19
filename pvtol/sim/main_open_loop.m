% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% Simulate planning and control of a planar quadrotor with (robust) CCM

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

clc;clear;

controller_type = 'open-loop';  % {'open-loop','ccm','rccm','rccm-p'} 
% ---------------------- load plant and controller ------------------------

file_controller = 'rccm_1.2_wmax_1_plim_0.33pi_pos.mat';           
load(file_controller);
controller.type = CtrlDesignOpts.open_loop;

%  -----------------------simulation settings -----------------------------
sim_config.include_tube = 0;    % whether to tighten the state bounds in planning a nominal trajectory
sim_config.tight_input_bnd = 0;    % whether to tighten the input bounds in planning a nominal trajectory
sim_config.include_obs = 1;        % whether to include the obstacles
sim_config.include_dist = 1;       % include the disturbance  
sim_config.save_sim_rst = 1;      % whether to save simulation results
sim_config.replan_nom_traj = 0;   % whether to replan a trajectory

use_generated_code = 1; % whether to use the generated codes for simulations: using generated codes can accelerate by at least one fold

n = 6; nu = 2;
x0 = zeros(6,1);            % initial state
xF = [10 10 0 0 0 0]';      % final state
duration = 13;
umax = 3*plant.m*plant.g;
u_bnd = [0 0; umax umax]';
x_bnd = [-inf -inf -state_set.p_lim -state_set.vx_lim, -state_set.vz_lim, -state_set.pd_lim;
          inf  inf  state_set.p_lim  state_set.vx_lim   state_set.vz_lim   state_set.pd_lim]';

% wind disturbance settings
w_max = 0.7;              % maximum amplitude of disturbance
T_w = 10;                % period of disturbance
dist_config.sim_config.include_dist = sim_config.include_dist;
dist_config.w_max = w_max;
dist_config.T_w = T_w;
dist_config.gen_dist= @(t) w_max*(0.8+0.2*sin(2*pi/dist_config.T_w.*t));

if sim_config.tight_input_bnd == 1
    u_bnd = u_bnd +[0 0; -tube_u -tube_u]';    
end
%% Plan or load a nominal trajecotory 
file_traj = 'nomTraj_w_obs_open_loop';    

if sim_config.replan_nom_traj == 1
    trajGen_config.x0 = x0;
    trajGen_config.xF = xF;
    trajGen_config.x_bnd = x_bnd;
    trajGen_config.u_bnd = u_bnd;
    trajGen_config.include_obs = sim_config.include_obs;
    trajGen_config.include_tube = sim_config.include_tube;
    trajGen_config.duration = duration;
    % ------------------------ Specify the obstacles--------------
    %     obs = [4.5 2.5 1;
    %            4 6 1;
    %            10 7.5 1;
    %            7.5 4.5 1;
    %            6.5 8.5 1];
    % obs = [3 7.5 0.8;
    %        3 5 0.8;
    %        5.5 3.5 0.8;
    %        5.5 7.5 0.8;
    %        7.5 5.5 0.8];
    %    
    obs = [3 5 0.8;           
               5.5 7 0.8;
               6 4 0.8;
               11 5 0.8];

    figure(1);clf;hold on;
    visualize_obs(obs);
    xlim([0 12]);
    ylim([0 12]);
    trajGen_config.obs = obs;
    soln = plan_traj_pvtol(plant,controller,trajGen_config);

    tF = soln.grid.time(end); trajGen_config.tF = tF;
    save(file_traj,'trajGen_config','soln');
else
    load(file_traj);
end
duration = trajGen_config.tF;   % modify the duration

% show the planned traj
x_nom_fcn = soln.interp.state;
u_nom_fcn = soln.interp.control;
times = 0:0.05:duration;
simuLen = length(times);
xnomTraj = zeros(n,simuLen);
unomTraj = zeros(nu,simuLen);
for t =1:simuLen
    xnomTraj(:,t) = x_nom_fcn(times(t));
    unomTraj(:,t) = u_nom_fcn(times(t));
end
figure(1);clf
hold on;
plot(xnomTraj(1,:),xnomTraj(2,:),'linewidth',1);
if trajGen_config.include_obs == 1
    visualize_obs(trajGen_config.obs);
end    
sim_config.trajGen_config = trajGen_config;

controller.x_nom_fcn = x_nom_fcn;
controller.u_nom_fcn = u_nom_fcn;
controller.w_nom = 0;  % nominal value for disturbances

%% simulate
dist0 = norm(x0); 
% -------------------------------------------------------
tic;
% ode23 is fastest, followed by ode45
[times,xuTraj] = ode23(@(t,xu) pvtol_dyn(t,xu,plant,controller,dist_config),[0 duration],x0); %,ode_opts)
toc;
xuTraj = xuTraj';
xTraj = xuTraj(1:n,:);

%% plot the result
n = plant.n; nu = plant.nu;
% ccm_sim_data = 'ccm_sim_w_dist_1.mat';
% rccm_sim_data = 'ccm_sim_w_dist_1.mat';
%-------------------------------------------------
% nominal trajectory
simuLen = length(times);
xnomTraj = zeros(n,simuLen);
unomTraj = zeros(nu,simuLen);
for t =1:simuLen
    xnomTraj(:,t) = x_nom_fcn(times(t));
    unomTraj(:,t) = u_nom_fcn(times(t));
end
addText = 'Disturbed';

close all;
figure(1); clf;
subplot(2,2,1); hold on;
plot(times,xnomTraj(1,:),'b--',times,xnomTraj(2,:),'r--');
plot(times,xTraj(1,:),'b-',times,xTraj(2,:),'r-');
xlabel('Time (s)')
ylabel('x & z (m)')
legend('x_{nom}', 'z_{nom}',['x: ' addText],['z: ' addText]);

subplot(2,2,2); hold on;
plot(times,xnomTraj(3,:)*180/pi,'--');
plot(times,xTraj(3,:)*180/pi);
xlabel('Time (s)')
ylabel('\phi (deg)')
legend('Nominal',addText);

figure(2);clf;
hold on;
h1 = plot(xnomTraj(1,:),xnomTraj(2,:),'b--');
h2 = plot(xTraj(1,:),xTraj(2,:),'r-');
scatter(xnomTraj(1,1),xnomTraj(1,1),'bo')
scatter(xnomTraj(1,end),xnomTraj(1,end),'b*')
scatter(xTraj(1,end),xTraj(2,end),'r*')
xlabel('x (m)')
ylabel('z (m)')
legend('Nominal', addText);


if sim_config.save_sim_rst == 1    
    file_name = 'sim_open_loop_w_dist_0.5_w_obs.mat';  
    save(file_name,'times','xTraj','xnomTraj','unomTraj','dist_config','sim_config','plant','controller');
end

return;


function dxdt = pvtol_dyn(t,xu,plant,controller,dist_config)
n = plant.n;
x = xu(1:n);
% tic;
u_nom = controller.u_nom_fcn(t);
% toc;
u = u_nom; 
wt = dist_config.gen_dist(t);
dxdt = plant.f_fcn(x)+plant.B_fcn(x)*u;
if dist_config.sim_config.include_dist == 1
   dxdt(1:n,:) = dxdt(1:n,:) + plant.Bw_fcn(x).*wt;
end

end

