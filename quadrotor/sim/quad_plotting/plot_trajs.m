n = plant.nc; nu = plant.nuc;
color = {'k','b',[0 0.5 0],'r',[0.8 0.9 0.9741],[0.8 0.98 0.9],[1 0.8 0.8],[0.7, 0.7 1]};
%-------------------------------------------------
if controller.type == CtrlDesignOpts.ccm
    addtext = 'ccm';
else
    addtext = 'rccm';
end

%% trajectory
figure(1);clf;
plot3dObstacles(goal_box,'g',0.2);
plot_all_obstacles();hold on;
scale = 10;
tube_xyz = controller.tube_gain_xyz*dist_config.w_max;
for i=1:length(t_vec)/scale
   center = [xcnomTraj(1,i*scale-1),-xcnomTraj(2,i*scale-1),-xcnomTraj(3,i*scale-1)];
   if i== 1 || norm(center-center_prev) >0.03
       [X,Y,Z] = ellipsoid(center(1),center(2),center(3),tube_xyz,tube_xyz,tube_xyz);
       surf(X,Y,Z,'FaceColor',color{8},'FaceAlpha',0.1,'EdgeColor','none'); %'FaceLighting','flat'
   else
       aa = 1;
   end
   center_prev = center;
end
h2 = plot3(xcTraj(1,:),-xcTraj(2,:),-xcTraj(3,:),'r-','Linewidth',0.8);
h1 = plot3(xcnomTraj(1,:),-xcnomTraj(2,:),-xcnomTraj(3,:),'k--','Linewidth',1.5);hold on;


scatter3(xcnomTraj(1,1),-xcnomTraj(2,1),-xcnomTraj(3,1),'ko')
scatter3(xcnomTraj(1,end),-xcnomTraj(2,end),-xcnomTraj(3,end),'ko')
scatter3(xcTraj(1,end),-xcTraj(2,end),-xcTraj(3,end),'r*')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
view(52.5,30);
legend([h1 h2],{'nominal', 'actual'},'Location','northoutside');
goodplot([6 6]);
set(gca,'Box','off');
figname = [addtext '_traj_3d'];
if print_file == 1
%     savefig([fig_name '.fig']);
%     print(figname, '-painters', '-dpdf', '-r150');
    export_fig(figname,'-pdf','-q101');

    h2.ZData = h2.ZData +0.7;
    h1.ZData = h1.ZData +0.7;
    view(2);
    legend off
%     savefig('traj_2d_10.fig');
    figname = [addtext '_traj_2d'];
%     print(figname, '-painters', '-dpdf', '-r150');
    export_fig(figname,'-pdf','-q101');
end

%% trajectory error
figure();hold on;
plot(t_vec,xcTraj(1,:)-xcnomTraj(1,:),'k-',t_vec,xcTraj(2,:)-xcnomTraj(2,:),'b-',t_vec,xcTraj(3,:)-xcnomTraj(3,:),'r-','Linewidth',1.5); hold on;
err = xcTraj(1:3,:)-xcnomTraj(1:3,:); err = sqrt(sum(err.^2,1));
plot(t_vec,err,'m--','Linewidth',2);
plot(t_vec,controller.tube_gain_xyz*dist_config.w_max*ones(1,length(t_vec)),'g--','Linewidth',1);
xlabel('Time (s)')
ylabel('Tracking error (m)');
legend('$p_x-p_x^\star$','$p_z-p_z^\star$','$p_z-p_z^\star$','$||p-p^\star||$','interpreter','latex','Location','best');% #1\right\rVert $||[e_x;e_y;e_z]||$
goodplot([6 2.5]);
figname = [addtext '_traj_err'];
if print_file == 1
%     savefig('traj_err_10.fig');
    print(figname,'-dpdf', '-r300');
%     export_fig(figname,'-pdf', '-q101');
end
%% Rotation angles
figure(); clf;
subplot(2,1,1)
plot(t_vec,xcnomTraj(8,:)*180/pi,'--','Linewidth',1.5);hold on;
plot(t_vec,xcTraj(8,:)*180/pi,'r-','Linewidth',1.5);
ylabel('$\phi$ (deg)','interpreter','latex');
legend('nominal','actual','location','north','orientation','horizontal');
goodplot([6 5]);
subplot(2,1,2)
plot(t_vec,xcnomTraj(9,:)*180/pi,'--','Linewidth',1.5);hold on;
plot(t_vec,xcTraj(9,:)*180/pi,'r-','Linewidth',1.5);
xlabel('Time (s)')
ylabel('$\theta$ (deg)','interpreter','latex');
goodplot([6 5]);
if print_file == 1
    figname = [addtext '_traj_err'];
%     savefig('rot_angle.fig');
    print(figname,'-painters', '-dpdf', '-r150');
end

%% Control effort
figure()
subplot(3,1,1);
h2 = plot(t_vec(1:end-1), ucTraj(1,:),'r-','linewidth',1.5);hold on
h1 = plot(t_vec(1:end-1), ucnomTraj(1,:),'k--','linewidth',1.5); 

% xlabel('Time (s)');
ylabel('$\dot{\tau}$ (m/s$^3$)','interpreter','latex');
legend([h1 h2],{'nominal','actual'},'location','best','orientation','horizontal');
goodplot([6 6]);

subplot(3,1,2);hold on
plot(t_vec(1:end-1), ucTraj(2,:),'r-','linewidth',1.5);
plot(t_vec(1:end-1), ucnomTraj(2,:),'k--','linewidth',1.5); 
% xlabel('Time (s)');
% legend('nominal','actual');
ylabel('$\dot{\phi}$ (rad/s)','interpreter','latex');
grid on
% set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
goodplot([6 6]);
subplot(3,1,3); hold on
plot(t_vec(1:end-1), ucTraj(3,:),'r-','linewidth',1.5);
plot(t_vec(1:end-1), ucnomTraj(3,:),'k--','linewidth',1.5);

xlabel('Time (s)');
ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex');
% grid on
% set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
goodplot([6 6],[0.1 0.01]);
if print_file == 1
    figname = [addtext '_control'];
%     savefig('control.fig');
    print(figname,'-painters', '-dpdf', '-r150');
end
%% Geodesic Energy
figure()
plot(t_vec(1:end-1),energyTraj,'b-','linewidth',1.5); hold on
% grid on
xlabel('Time (s)'); ylabel('Energy');
goodplot([6 2.5],[0.1 0.01])
if print_file == 1
    figname = [addtext '_energy'];
%     savefig('energy.fig');
    print(figname,'-painters', '-dpdf', '-r150');
end
%% Movie
%plot_quad_movie2(xcnomTraj',t_vec,xTraj',round(0.1/sim_config.dt_sim),controller,goal,...
%                tree_obstacles,tower_obstacles,obstacles_infl,World_dim,tube_xyz);


