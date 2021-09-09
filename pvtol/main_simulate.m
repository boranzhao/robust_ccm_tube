%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation of planning and control of a planar quadrotor with (robust) CCM

%  Author: Pan Zhao, UIUC, Advanced Controls Research Lab,
%  panzhao2@illinois.edu
%  Codes for the paper:
%  P. Zhao, et al. Tube-certified trajectory tracking for nonliner systems
%  with robust control contraction metrics. Submitted to IEEE Robotics and
%  Automation Letters, 2021. 
%  Last update: Sep 9, 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
controller_type = 'ccm';            % {'open-loop','ccm','rccm','rccm-p'} 
% ---------------------- load plant and controller ------------------------
if strcmp(controller_type, 'ccm')
    file_controller = 'ccm_0.8_plim_0.33pi.mat';           
elseif strcmp(controller_type, 'rccm') 
    file_controller = 'rccm_1.4_wmax_1_plim_0.33pi.mat';
elseif strcmp(controller_type, 'rccm-p') 
    file_controller = 'rccm_1.2_wmax_1_plim_0.33pi_pos.mat';
else    
    error('Please choose a controller');
end
load(file_controller);
%  -----------------------simulation settings -----------------------------
sim_config.include_tube = 1;        % whether to tighten the state bounds in planning a nominal trajectory
sim_config.tight_input_bnd = 1;     % whether to tighten the input bounds in planning a nominal trajectory
sim_config.include_obs = 1;         % whether to include the obstacles
sim_config.include_dist = 1;        % include the disturbance  
sim_config.save_sim_rst = 0;        % whether to save simulation results
sim_config.replan_nom_traj = 1;     % whether to replan a trajectory

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
file_traj = 'nomTraj';
if sim_config.include_obs == 0
    file_traj = ['nomTraj' '.mat']; 
elseif sim_config.include_obs == 1
    file_traj = ['nomTraj_w_obs_' file_controller];    
end
if sim_config.replan_nom_traj == 1
    trajGen_config.x0 = x0;
    trajGen_config.xF = xF;
    trajGen_config.x_bnd = x_bnd;
    trajGen_config.u_bnd = u_bnd;
    trajGen_config.include_obs = sim_config.include_obs;
    trajGen_config.include_tube = sim_config.include_tube;
    trajGen_config.duration = duration;
    % ------------------------ Specify the obstacles-----------------------
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
duration = trajGen_config.tF;   % modify the duration according to the computed trajectory

% --------------------- show the planned traj -----------------------------
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

%% Formulate the NLP problem for geodesic computation
controller.use_generated_code = use_generated_code;
lambda = controller.lambda;
%  problem setting for geodesic computation
D = 2;      % degree of the polynomial
N = D+6;    % stopping index for the CGL (Chebyshev-Gauss-Lobatto) nodes: #notes N+1

% optimization variables: chebyshev coefficients for geodesics
% [c_10, c_11,..., c_1D, ..., c_n0, c_n1,..., c_nD]

% --------------------obtain chebyshev pseudospectral numerics-------------
% --------------------(to be used for computing the integral)--------------
[s,w_cheby] = clencurt(N); % t is 1 by N+1, with values lying b/t 0 and 1.
% evaluate the value of the CCM
W = zeros(n,n,N+1);

% compute Cheby basis at all points
[T, Tdot] = compute_cheby(N,D,s); % Both T and T_dot are D+1 by N+1
% for equality constraints
Aeq = [kron(eye(n),T(:,1)'); kron(eye(n),ones(1,D+1))];
Aeq = sparse(Aeq);

% --------- formulate and solve the NLP problem using OPTI --------------
ndec = n*(D+1);
if controller.use_generated_code == 1
    costf = @(c) RiemannEnergy1_mex(c,n,D,N,T,Tdot,w_cheby);
    grad = @(c) energyGradient1_mex(c,n,D,N,T,Tdot,w_cheby);     
else    
    costf = @(c) RiemannEnergy(c,n,D,N,T,Tdot,w_cheby,controller.W_fcn);
    grad = @(c) energyGradient(c,n,D,N,T,Tdot,w_cheby,controller.W_fcn,controller.dW_dxi_fcn); 
end

geodesic.D = D; geodesic.N = N; geodesic.ndec = ndec;
geodesic.T = T; geodesic.Tdot = Tdot;
geodesic.Aeq = Aeq; 
geodesic.costf = costf;
geodesic.grad = grad;
geodesic.w_cheby = w_cheby;

beq = zeros(2*n,1);c0 = zeros(n*(D+1),1);
% add some bounds to mitigate numerical issues when computing the geodesic
% lb = -20*ones(size(c0));
% ub = 20*ones(size(c0));
lb = -20*ones(n,D+1);
ub = 20*ones(n,D+1);
lb(3:4,:) = -5*ones(2,D+1);  % phi and vx
ub(3:4,:) = 5*ones(2,D+1);   % phi and vx
lb = lb';lb = lb(:);
ub = ub';ub= ub(:);

% i = 3;
% lb((i-1)*(D+1)+(1:D+1),1)= -10*ones(D+1,1);
% ub((i-1)*(D+1)+(1:D+1),1)= 10*ones(D+1,1);
% i = 4;
% lb((i-1)*(D+1)+(1:D+1),1)= -10*ones(D+1,1);
% ub((i-1)*(D+1)+(1:D+1),1)= 10*ones(D+1,1);

% --------------- re-generate the code is necessary after change of ------
% controller or geodesic optimization settings: remember to re-generate the
% m-file functions, e.g., dW_dphi, dW_dvx, etc., first.------------------  
if controller.use_generated_code 
    answer = questdlg('Are the generated codes for this particular scenario?','Question for using C-code in simulation','Yes','No','No');
    switch answer 
        case 'Yes'
        case 'No'
            error('You cannot continue without including the generated codes for this scenario!');
    end
end

%     [copt1,Erem,exitflag,info] = solve(Opt,c0);

% ---------- for using ipopt solver ---------------------------------------
% opts_opti = optiset('solver','ipopt','maxiter', 500,'display','iter');  %,,'derivCheck','on'
% Opt = opti('fun',geodesic.costf,'grad',geodesic.grad,'eq',geodesic.Aeq,beq,...
%             'bounds',lb,ub,'ndec',geodesic.ndec,'x0',c0,'options',geodesic.opts_opti);
% geodesic.nlprob = convIpopt(Opt.prob,geodesic.opts_opti); 
% -------------------------------------------------------------------------

% ----------for using matlab fmincon solver--------------------------------
opts_opti = optiset('solver','matlab','maxiter',500,'tolrfun',1e-4,'tolafun',1e-4,'display','off','derivCheck','off'); 
Opt = opti('fun',costf,'grad',grad,'eq',Aeq,beq,'bounds',lb,ub,'ndec',ndec,'x0',c0,'options',opts_opti);
geodesic.nlprob = convMatlab(Opt.prob,opts_opti); 
geodesic.nlprob.options = ...
optimoptions(@fmincon,'Display','off','HessianApproximation','lbfgs',...
'MaxIterations',opts_opti.maxiter,'SpecifyObjectiveGradient',true,'CheckGradients',false,...
'OptimalityTolerance',opts_opti.tolrfun,'FunctionTolerance',opts_opti.tolrfun,'FunValCheck','on','StepTolerance',1.0e-8);
% -------------------------------------------------------------------------

geodesic.opts_opti = opts_opti;
controller.x_nom_fcn = x_nom_fcn;
controller.u_nom_fcn = u_nom_fcn;
controller.geodesic = geodesic;
controller.w_nom = 0;  % nominal value for disturbances

% simulate
dist0 = norm(x0); 
% -------------------------------------------------------------------------
% [times,xTraj] = ode23s(@(t,x) plant.f_fcn(x)+plant.B*ccm_law(t,x,plant,controller),[0 duration],x0); 
% toc;
% return;
% ode_opts = odeset('MaxStep',5e-1);
% --------for additionally outputing control inputs and Reim. energy-------
% compute the initial Riemann energy function value
ue = ccm_law(0,x0,plant,controller);
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
if controller.type == CtrlDesignOpts.ccm
    addText = 'CCM';
else
    addText = 'RCCM';
end
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

subplot(2,2,3);hold on;
plot(times,unomTraj(1,:),'b--',times,unomTraj(2,:),'r--');
plot(times,uTraj(1,:),'b-',times,uTraj(2,:),'r-');
xlabel('Time (s)');
ylabel('u (N)')
legend('u_{nom,1}', 'u_{nom,2}',['u_1: ' addText],['u_2: ' addText]);
% plot(times,distTraj);
% xlabel('Time (s)');
% ylabel('$\|x-x^\star\|/\|x_0-x^\star_0\|$','interpreter','latex');

subplot(2,2,4)
plot(times, energyTraj);
xlabel('Time (s)');
ylabel('Riemann energy')
legend(addText)

figure(2);clf;
hold on;
h1 = plot(xnomTraj(1,:),xnomTraj(2,:),'b--');
h2 = plot(xTraj(1,:),xTraj(2,:),'r-');
scatter(xnomTraj(1,1),xnomTraj(1,1),'bo')
scatter(xnomTraj(1,end),xnomTraj(1,end),'b*')
scatter(xTraj(1,end),xTraj(2,end),'r*')
xlabel('z (m)')
ylabel('y (m)')
legend('Nominal', addText);


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


function dxdt = pvtol_dyn(t,xu,Klp,plant,controller,dist_config)
n = plant.n;
x = xu(1:n);
ue_filter = xu(n+1:end);
% tic;
ue = ccm_law(t,x,plant,controller);
% toc;
u = ue(1:end-1); 
wt = dist_config.gen_dist(t);
dxdt = [plant.f_fcn(x); -Klp.*ue_filter]+[plant.B_fcn(x)*u;  Klp.*ue];
if dist_config.sim_config.include_dist == 1
   dxdt(1:n,:) = dxdt(1:n,:) + plant.Bw_fcn(x).*wt;
end

end

