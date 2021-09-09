%% plot the result
close all;
print_file = 1;     % whether to save the fig and print it to pdf. 
line_wd = 1.5;      %
plot_rst_w_obs = 1; % whether to plot the results with obstacles
w_max = 1;          % amplitude of disturbance

ccm_sim_data = 'sim_ccm_lam_0.8_w_dist_1_w_obs.mat';
rccm_sim_data = 'sim_rccm_lam_1.4_w_dist_1_w_obs.mat';
rccm_pos_sim_data = 'sim_rccm_pos_lam_1.2_w_dist_1_w_obs.mat';
ccm_nomTraj_data = 'nomTraj_w_obs_ccm_0.8_plim_0.33pi.mat';
rccm_nomTraj_data = 'nomTraj_w_obs_rccm_1.4_wmax_1_plim_0.33pi.mat';
rccm_pos_nomTraj_data = 'nomTraj_w_obs_rccm_1.2_wmax_1_plim_0.33pi_pos.mat';
ol_sim_data = 'sim_open_loop_w_dist_0.5_w_obs.mat';

n = 6; nu =2;

% load data

% CCM
load(ccm_sim_data);
times_ccm = times;
xTraj_ccm = xTraj;
uTraj_ccm = uTraj;
energyTraj_ccm = energyTraj;
controller_ccm = controller;
controller_ccm.tube_x = controller_ccm.tube_gain_x*w_max;
controller_ccm.tube_z = controller_ccm.tube_gain_z*w_max;
load(ccm_nomTraj_data);
nomTrajSoln = soln; 
simuLen = length(times);
xnomTraj_ccm = zeros(n,simuLen);
unomTraj_ccm = zeros(nu,simuLen);
for t =1:simuLen
    xnomTraj_ccm(:,t) = nomTrajSoln.interp.state(times(t));
    unomTraj_ccm(:,t) = nomTrajSoln.interp.control(times(t));
end
times_ccm = times;


if sim_config.include_obs == 1
    fig_name_com = '_w_obs';
    obs = sim_config.trajGen_config.obs;
else
    fig_name_com = '';
end

% RCCM
load(rccm_sim_data);
times_rccm = times;
xTraj_rccm = xTraj;
uTraj_rccm = uTraj;
energyTraj_rccm = energyTraj;
% tube_gain_rccm = tube_gain;
% tube_rccm = tube_gain*w_max*artificial_gain; 
controller_rccm = controller;
controller_rccm.tube_xz = controller_rccm.tube_gain_xz*w_max;
controller_rccm.tube_x = controller_rccm.tube_xz;
controller_rccm.tube_z = controller_rccm.tube_xz; 
load(rccm_nomTraj_data);
nomTrajSoln = soln; 
simuLen = length(times);
xnomTraj_rccm = zeros(n,simuLen);
unomTraj_rccm = zeros(nu,simuLen);
for t =1:simuLen
    xnomTraj_rccm(:,t) = nomTrajSoln.interp.state(times(t));
    unomTraj_rccm(:,t) = nomTrajSoln.interp.control(times(t));
end
times_rccm = times;

% open loop
load(ol_sim_data);
times_ol = times;
xTraj_ol = xTraj;
uTraj_ol = uTraj;
% tube_gain_ol = tube_gain;
% tube_ol = tube_gain*w_max*artificial_gain;
controller_ol = controller;
load('nomTraj_w_obs_open_loop');
xnomTraj_ol =  xnomTraj;
unomTraj_ol = unomTraj;
n = plant.n; nu = plant.nu;

%% ----------- Plot the trajectory -----------------

% color testing
color = {'k','b',[0 0.5 0],'r',[0.8 0.9 0.9741],[0.8 0.98 0.9],[1 0.8 0.8]};
% figure;hold on;
% plot(times',xnomTraj(1,:),'color',color{5},'Linewidth',5);
% plot(times',xnomTraj(1,:)+0.5,'color',color{6},'Linewidth',5);
% plot(times',xnomTraj(1,:)-0.5,'color',color{7},'Linewidth',5);
 
linestyle = {':','--','-.','-'};
% close all;

figure(1);clf;
hold on;
% add the tube using fill functions

% CCM tube
scale = 20;
Len = length(times_ccm);
for i=1:Len/scale
    center = [xnomTraj_ccm(1,i*scale-1),xnomTraj_ccm(2,i*scale-1)];
    h2_tube = ellipse(controller_ccm.tube_x, controller_ccm.tube_z,0,center(1),center(2),color{5},[],1);
end
% RCCM tube
Len = length(times_rccm);
for i=1:Len/scale
    center = [xnomTraj_rccm(1,i*scale-1),xnomTraj_rccm(2,i*scale-1)];
   h3_tube= ellipse(controller_rccm.tube_xz,controller_rccm.tube_xz,0,center(1),center(2),color{6},[],1);
end
h2 = plot(xTraj_ccm(1,:),xTraj_ccm(2,:),[color{2} linestyle{4}],'Linewidth',line_wd);
% open loop

h1 = plot(xTraj_ol(1,1:162),xTraj_ol(2,1:162),[color{4} linestyle{4}],'Linewidth',line_wd);


% nominal trajectories
h3 = plot(xTraj_rccm(1,:),xTraj_rccm(2,:),'-','color',color{3},'Linewidth',line_wd);
h0 = plot(xnomTraj_ccm(1,:),xnomTraj_ccm(2,:),'k--','Linewidth',line_wd);
plot(xnomTraj_ccm(1,1),xnomTraj_ccm(1,1), 'ro','MarkerSize',10);
plot(xnomTraj_ccm(1,end),xnomTraj_ccm(1,end),'rh','MarkerSize',10);

plot(xnomTraj_rccm(1,:),xnomTraj_rccm(2,:),'k--','Linewidth',line_wd);

plot(xnomTraj_ol(1,:),xnomTraj_ol(2,:),'k--','Linewidth',line_wd);

visualize_obs(obs); 

xlabel('$p_x$ (m)','interpreter','latex');
ylabel('$p_z$ (m)','interpreter','latex')
plot(xTraj_ol(1,159),xTraj_ol(2,159),'rp','MarkerSize',10,'MarkerFaceColor','r'); 
legend([h0 h1,h2,h3,h2_tube,h3_tube],{'Planned traj.','Actual traj.: OL','Actual traj.: CCM','Actual traj.: RCCM'},'Location','northwest');
ylim(ylim-[2 0]);
annotation('arrow',[0.15 0.25],[0.55 0.55],'Linewidth',2);
text(-4.8,5.5,'Wind','Fontsize',15,'color','r');
text(-2,-1,'Start','Fontsize',15);
text(9,12,'Goal','Fontsize',15);
axis off
% add the axes
shift = 0.65;
annotation('arrow',[0.1+shift 0.17+shift],[0.15 0.15],'Linewidth',2); 
annotation('arrow',[0.1+shift 0.1+shift],[0.15 0.25],'Linewidth',2);
text(-4+shift*37,-8+shift*2,'X', 'Fontsize',15);text(-7+shift*37,-4.5+shift*2,'Z', 'Fontsize',15);

goodplot([6 5]);
fig_name = 'traj_w_open_loop';
return;
if print_file == 1
    print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
end


