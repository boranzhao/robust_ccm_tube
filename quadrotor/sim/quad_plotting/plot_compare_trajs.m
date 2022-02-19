% compare the performance of different controllers under different scenarios. 
%-------------------------------------------------
close all;

%% plot the result
print_file = 0;     % whether to save the fig and print it to pdf. 
line_wd = 1.5;      %

ccm_sim_data = 'sim_ccm_lam_3.4_w_dist_1.mat';
rccm_p_sim_data  = {'sim_rccm_lam_3.4_w_dist_1.mat', ...
                    'sim_rccm_lam_3.4_w_dist_2.mat',...
                    'sim_rccm_lam_3.4_w_dist_5.mat',...
                    'sim_rccm_lam_3.4_w_dist_6.mat'}; 
                
rccm_p_sim_data  = {'sim_rccm_lam_3.4_w_dist_1.mat'}; 
     
num_data = length(rccm_p_sim_data);
% load data
XcTraj = cell(num_data+1,1);
XTraj = cell(num_data+1,1);
UcTraj = cell(num_data+1,1);
EnergyTraj = cell(num_data+1,1);
W_max = nan(num_data+1,1);

i = num_data+1;
load(ccm_sim_data);
XcTraj{i} = xcTraj;
XTraj{i} = xTraj;
UcTraj{i} = ucTraj;
EnergyTraj{i} = energyTraj;
W_max(i) = dist_config.w_max;
t_vec_ccm = t_vec;
xcnomTraj_ccm = xcnomTraj;
for i = 1:num_data
    load(rccm_p_sim_data{i});
    XcTraj{i} = xcTraj;
    XTraj{i} = xTraj;
    UcTraj{i} = ucTraj;
    EnergyTraj{i} = energyTraj;
    W_max(i) = dist_config.w_max;
end

T_steps = length(t_vec);

n = plant.nc; nu = plant.nuc;
color = {'r','b','k',[0 0.5 0],[0.8 0.9 0.9741],[0.8 0.98 0.9],[1 0.8 0.8],[0.7, 0.7 0.9]};
line = {'-','--','-.','-','--'};

%% trajectory
% figure(1);clf;
% plot3dObstacles(goal_box,'g',0.2);
% plot_all_obstacles();hold on;
% scale = 20;
% tube_xyz = controller.tube_gain_xyz*dist_config.w_max;
% for i=1:T_steps/scale
%    center = [xcnomTraj(1,i*scale-1),-xcnomTraj(2,i*scale-1),-xcnomTraj(3,i*scale-1)];
%    [X,Y,Z] = ellipsoid(center(1),center(2),center(3),tube_xyz,tube_xyz,tube_xyz);
%    surf(X,Y,Z,'FaceColor',color{8},'FaceAlpha',0.2,'EdgeColor','none'); %'FaceLighting','flat'
% end
% h2 = plot3(xcTraj(1,:),-xcTraj(2,:),-xcTraj(3,:),'r-','Linewidth',1.5);
% h1 = plot3(xcnomTraj(1,:),-xcnomTraj(2,:),-xcnomTraj(3,:),'k--','Linewidth',1.5);hold on;
% 
% scatter3(xcnomTraj(1,1),-xcnomTraj(2,1),-xcnomTraj(3,1),'ko')
% scatter3(xcnomTraj(1,end),-xcnomTraj(2,end),-xcnomTraj(3,end),'ko')
% scatter3(xcTraj(1,end),-xcTraj(2,end),-xcTraj(3,end),'r*')
% xlabel('X (m)')
% ylabel('Y (m)')
% zlabel('Z (m)')
% view(52.5,30);
% legend([h1 h2],{'nominal', 'actual'},'Location','northoutside');
% goodplot([6 6]);
% fig_name = 'traj_3d';
% if print_file == 1
%     savefig([fig_name '.fig']);
%     print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
%     view(2);
%     legend off
%     savefig('traj_2d_10.fig');
%     print('traj_2d_10.pdf', '-painters', '-dpdf', '-r150');
% end

%% trajectory error
figure();hold on;
% plot(t_vec,xcTraj(1,:)-xcnomTraj(1,:),'k-',t_vec,xcTraj(2,:)-xcnomTraj(2,:),'b-',t_vec,xcTraj(3,:)-xcnomTraj(3,:),'r-','Linewidth',1.5); hold on;
err = zeros(num_data,T_steps);
legend_txt0 = cell(num_data,1);
for i = 1:num_data 
   err_tmp = XcTraj{i}(1:3,:)-xcnomTraj(1:3,:); err(i,:) = sqrt(sum(err_tmp.^2,1));
    plot(t_vec,err(i,:),'color',color{i},'LineStyle',line{i},'Linewidth',line_wd);
    legend_txt0{i} = sprintf('$\\bar w= %d$',W_max(i));
end
i = i+ 1;
err_tmp = XcTraj{i}(1:3,:)-xcnomTraj_ccm(1:3,:); err_ccm = sqrt(sum(err_tmp.^2,1));
plot(t_vec_ccm,err_ccm,'color',color{i},'LineStyle',line{i},'Linewidth',line_wd);
legend_txt0{i} = sprintf('CCM: $\\bar w= %d$',W_max(i));

plot(t_vec,controller.tube_gain_xyz*1*ones(1,T_steps),'g--','Linewidth',1);
xlabel('Time (s)')
ylabel('$||p-p^\star||$','interpreter','latex');
legend(legend_txt0,'interpreter','latex','Location','best','Orientation','Horizontal');% #1\right\rVert $||[e_x;e_y;e_z]||$
goodplot([6 2.5]);
if print_file == 1
%     savefig('rccm_p_traj_err_diffW.fig');
    print('rccm_p_traj_err_diffW.pdf','-painters', '-dpdf', '-r150');
end
%% Rotation angles
figure(); clf;
legend_txt = cell(num_data+1,1); legend_txt{1} = 'nominal';
subplot(3,1,1)
plot(t_vec,xcnomTraj(8,:)*180/pi,'k:','Linewidth',1.5);hold on;
for i = 1:num_data     
    plot(t_vec,XcTraj{i}(8,:)*180/pi,'color',color{i},'LineStyle',line{i},'Linewidth',1.5);
    legend_txt{i+1} = sprintf('$\\bar w= %d$',W_max(i));
end
% plot(t_vec,xcnomTraj(8,:)*180/pi + controller.tube_gain_rp*1*180/pi*ones(1,T_steps),'g--','Linewidth',1);
% plot(t_vec,xcnomTraj(8,:)*180/pi - controller.tube_gain_rp*1*180/pi*ones(1,T_steps),'g--','Linewidth',1);
ylim(ylim+[0 10]);
ylabel('$\phi$ (deg)','interpreter','latex');
legend(legend_txt,'interpreter','latex','location','north','orientation','horizontal');%
goodplot([6 6]);
subplot(3,1,2)
plot(t_vec,xcnomTraj(9,:)*180/pi,'k:','Linewidth',1.5);hold on;
for i = 1:num_data
    plot(t_vec,XcTraj{i}(9,:)*180/pi,'color',color{i},'LineStyle',line{i},'Linewidth',1.5);
end
% plot(t_vec,xcnomTraj(9,:)*180/pi+controller.tube_gain_rp*1*180/pi*ones(1,T_steps),'g--','Linewidth',1);
% plot(t_vec,xcnomTraj(9,:)*180/pi-controller.tube_gain_rp*1*180/pi*ones(1,T_steps),'g--','Linewidth',1);
ylim(ylim+[0 10]);
% xlabel('Time (s)')
ylabel('$\theta$ (deg)','interpreter','latex');
goodplot([6 6]);
subplot(3,1,3); hold on;
for i = 1:num_data 
    err_tmp = XcTraj{i}(8:9,:)-xcnomTraj(8:9,:); err(i,:) = sqrt(sum(err_tmp.^2,1))*180/pi;
    plot(t_vec,err(i,:),'color',color{i},'LineStyle',line{i},'Linewidth',line_wd);
end
plot(t_vec,controller.tube_gain_rp*1*180/pi*ones(1,T_steps),'g--','Linewidth',1);
ylim(ylim+[-10 10]);
xlabel('Time (s)')
ylabel('$|\!|[\phi\!,\theta]\!- \![\phi^\star\!\!,\theta^\star]|\!|\!$ (deg)','interpreter','latex');
legend(legend_txt0{1:num_data},'interpreter','latex','location','north','orientation','horizontal');%
goodplot([6 6]);
if print_file == 1
%     savefig('rccm_p_rot_angle_diffW.fig');
    print('rccm_p_rot_angle_diffW.pdf','-painters', '-dpdf', '-r150');
end

%% Control effort
figure()
% subplot(3,1,1);
h  = cell(num_data+1,1);
for i = 1:num_data
    h{i+1} = plot(t_vec(1:end-1), UcTraj{i}(1,:),'color',color{i},'LineStyle',line{i},'linewidth',1.5);hold on
end
h{1} = plot(t_vec(1:end-1), ucnomTraj(1,:),'k:','linewidth',1.5); 

xlabel('Time (s)');
ylabel('$\dot{\tau}$ (m/s$^3$)','interpreter','latex');
legend([h{:}],legend_txt,'interpreter','latex','location','best','orientation','horizontal');%'numcolumns',2
goodplot([6 6]);

% subplot(3,1,2);hold on
% for i = 1:num_data
%     plot(t_vec(1:end-1), UcTraj{i}(2,:),'color',color{i},'LineStyle',line{i},'linewidth',1.5);hold on
% end
% plot(t_vec(1:end-1), ucnomTraj(2,:),'k--','linewidth',1.5); 
% % xlabel('Time (s)');
% % legend('nominal','actual');
% ylabel('$\dot{\phi}$ (rad/s)','interpreter','latex');
% grid on
% % set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
% goodplot([6 6]);
% subplot(3,1,3); hold on
% for i = 1:num_data
%     plot(t_vec(1:end-1), UcTraj{i}(3,:),'color',color{i},'LineStyle',line{i},'linewidth',1.5);hold on
% end
% plot(t_vec(1:end-1), ucnomTraj(3,:),'k--','linewidth',1.5);
% 
% xlabel('Time (s)');
% ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex');
% % grid on
% % set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
% goodplot([6 6],[0.1 0.01]);
goodplot([6 2.5],[0.1 0.01]);
if print_file == 1
%     savefig('rccmp_p_control_diffW.fig');
    print('rccm_p_control_diffW.pdf','-painters', '-dpdf', '-r150');
end
%% Geodesic Energy
figure()
for i = 1:num_data
    plot(t_vec(1:end-1),EnergyTraj{i},'color',color{i},'LineStyle',line{i},'linewidth',1.5); hold on
end
% grid on
xlabel('Time (s)'); ylabel('Energy');
legend(legend_txt0,'interpreter','latex','orientation','horizontal');% #1\right\rVert $||[e_x;e_y;e_z]||$
goodplot([6 2.5],[0.1 0.01])
if print_file == 1
%     savefig('rccm_p_energy_diffW.fig');
    print('rccm_p_energy_diffW.pdf','-painters', '-dpdf', '-r150');
end
%% Movie
plot_quad_movie2(MP_state(:,1),MP_state(:,2),MP_state(:,3),MP_t(1:end-1),xTraj',round(0.1/sim_config.dt),controller,goal,...
                tree_obstacles,tower_obstacles,obstacles_infl,World_dim);

