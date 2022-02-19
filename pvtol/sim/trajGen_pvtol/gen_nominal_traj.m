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