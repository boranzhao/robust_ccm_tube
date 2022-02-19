print_file = 0;
line_wd = 0.8;
color = {'k','b',[0 0.5 0],'r',[0.8 0.9 0.9741],[0.8 0.98 0.9],[1 0.8 0.8]};
c_g = 'g';[0 0.5 0];
% [lambda  ccm:all_states ccm:xyz rccm:all_states rccm:xyz
% rccm-p:all_states rccm-p:xyz];
tubes = [0.4    6.256   5.287   2.10    2.06    44.87   0.57; 
         0.6    5.646   3.709   1.81    1.757   39.29   0.485;
         0.8    5.83    2.665   1.70    1.561   35.07   0.441;
         1.0    6.263   2.234   1.643   1.437   32.09   0.413;
         1.2    6.808   1.599   1.612   1.356   29.66   0.395;
         1.4    7.42    1.81    1.598   1.307   27.85   0.382;
         1.6    8.086   1.840   1.595   1.286   26.69   0.376;
         1.8    8.78    1.538   1.60    1.309   25.69   0.37;
         2.0    9.503   1.238   1.61    1.30    24.96   0.361;
         2.2    10.243  1.368   1.628    1.263   24.56   0.354;
         2.4    10.997  1.128    1.652   1.262   24.116  0.350;
         2.6    11.763  0.896    1.666   1.263   23.644  0.3403;
         2.8    12.538  0.768    1.717   1.250   22.89   0.33;
         3.0    13.321  0.807    1.76    1.229   22.011  0.322;
         3.2    14.110  0.615    1.809   1.194   21.13   0.318;
         3.4    14.905  0.580    1.867   1.165   20.349  0.316;
         3.6    15.703  0.765    1.932   1.138   19.677  0.318];

lambdas = tubes(:,1);
ccm_tubes = tubes(:,2:3);
rccm_tubes = tubes(:,4:5);
rccm_p_tubes = tubes(:,6:7);
close all;

figure;
subplot(2,1,1); 
semilogy(lambdas,ccm_tubes(:,1),'bo--','linewidth',line_wd);hold on;
semilogy(lambdas,rccm_tubes(:,1),'*-.','color',c_g,'linewidth',line_wd);
semilogy(lambdas,rccm_p_tubes(:,1),'rs-','linewidth',line_wd);
ylabel('$\alpha_x$', 'interpreter','latex');
legend('CCM','RCCM','RCCM-P','Orientation','Horizontal');
% ylim(ylim+[-0.2 100]);
% yticks([10 100]);
grid on;
goodplot([6 4.5]);

subplot(2,1,2);
semilogy(lambdas,ccm_tubes(:,2),'bo--','linewidth',line_wd); hold on;
semilogy(lambdas,rccm_tubes(:,2),'*-.','color',c_g,'linewidth',line_wd);
semilogy(lambdas,rccm_p_tubes(:,2),'rs-','linewidth',line_wd);
ylabel('$\alpha_{p}$', 'interpreter','latex');
xlabel('$\lambda$ ($\hat \lambda$)', 'interpreter','latex');
ylim([0.05 11]);
yticks([0.1 10]);
grid on;
goodplot([6 4.5]);

% subplot(3,1,3); 
% semilogy(lambdas,ccm_tubes(:,3),'bo--','linewidth',line_wd);hold on;
% semilogy(lambdas,rccm_tubes(:,3),'*-.','color',c_g,'linewidth',line_wd);
% semilogy(lambdas,rccm_p_tubes(:,3),'rs-','linewidth',line_wd);
% ylabel('$\alpha_u \bar w$', 'interpreter','latex');
% xlabel('$\lambda$ ($\hat \lambda$)','interpreter','latex');
% ylim([0.05 1100]);
% yticks([0.1 10 100]);
% grid on;
% goodplot([6 5]);
fig_name = 'lambda_vs_tube_size';
if print_file == 1
    savefig([fig_name '.fig']);
    print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
end
return 
   
%% load data and plot the tubes

% CCM
load(ccm_data);
controller_ccm  = controller;
controller_ccm.tube_x = controller_ccm.tube_gain_x*w_max;
controller_ccm.tube_z = controller_ccm.tube_gain_z*w_max;
controller_ccm.tube_vx = controller_ccm.tube_gain_vx*w_max;
controller_ccm.tube_vz = controller_ccm.tube_gain_vz*w_max;
controller_ccm.tube_phi = controller_ccm.tube_gain_phi*w_max;
controller_ccm.tube_phidot = controller_ccm.tube_gain_phidot*w_max;

% RCCM
load(rccm_data);
controller_rccm = controller;
controller_rccm.tube_xz = controller_rccm.tube_gain_xz*w_max;
controller_rccm.tube_x = controller_rccm.tube_xz;
controller_rccm.tube_z = controller_rccm.tube_xz; 

% RCCM (pos)
load(rccm_pos_data);
controller_rccm_pos = controller;
controller_rccm_pos.tube_xz = controller_rccm_pos.tube_gain_xz*w_max;
controller_rccm_pos.tube_x = controller_rccm_pos.tube_xz;
controller_rccm_pos.tube_z = controller_rccm_pos.tube_xz;
n = plant.n; nu = plant.nu;

% Plot the tubes 
close all
figId = 1;
figure(figId);
subplot(2,2,1);hold on;
ellipse(controller_ccm.tube_x, controller_ccm.tube_z,0,0,0,color{5},[],1);
ellipse(controller_rccm.tube_xz, controller_rccm.tube_xz,0,0,0,'g',[],1);
ellipse(controller_rccm_pos.tube_xz, controller_rccm_pos.tube_xz,0,0,0,'r',[],1);
xlabel('$p_x$ (m)','interpreter','latex');
ylabel('$p_z$ (m)','interpreter','latex');
legend('CCM','RCCM','RCCM-P','Orientation','Horizontal');
axis equal
% ylim([-5 5]);
% subplot(2,2,2);hold on;
% ellipse(d_bar_ccm*sqrt(W_bar_ccm(4,4)), d_bar_ccm*sqrt(W_bar_ccm(5,5)),0,0,0);
% ellipse(controller_rccm.tube_xz, controller_rccm.tube_xz,0,0,0,'m');
% xlabel('$v_x$ (m/s)','interpreter','latex');
% ylabel('$v_z$ (m/s)','interpreter','latex');
% legend('CCM','RCCM');
% axis equal
% 
% subplot(2,2,3); hold on;
% ellipse(d_bar_ccm*sqrt(W_bar_ccm(3,3)), d_bar_ccm*sqrt(W_bar_ccm(6,6)),0,0,0);
% ellipse(controller_rccm.tube_xz, controller_rccm.tube_xz,0,0,0,'m');
% legend('CCM','RCCM');
% xlabel('$\phi$ (rad/s)','interpreter','latex');
% ylabel('$\dot{\phi}$ (rad)','interpreter','latex');
% axis equal
goodplot([6 5]);
subplot(2,2,2);hold on;
ellipse(controller_ccm.u_bnd, controller_ccm.u_bnd,0,0,0,color{5},[],1);
ellipse(controller_rccm.u_bnd, controller_rccm.u_bnd,0,0,0,'g',[],1);
ellipse(controller_rccm_pos.u_bnd, controller_rccm_pos.u_bnd,0,0,0,'r',[],1);
xlabel('$u_1$ (N)','interpreter','latex');
ylabel('$u_2$ (N)','interpreter','latex');
% legend('CCM','RCCM','RCCM-P');
axis equal
% ylim([-5 5]);
goodplot([6 5])
subplot(2,2,3);hold on;
ellipse(controller_ccm.tube_vx, controller_ccm.tube_vz,0,0,0,color{5},[],1);
ellipse(controller_rccm_pos.tube_vxvz, controller_rccm_pos.tube_vxvz,0,0,0,'r',[],1);
ellipse(controller_rccm.tube_vxvz, controller_rccm.tube_vxvz,0,0,0,'g',[],1);
xlabel('$v_x$ (m/s)','interpreter','latex');
ylabel('$v_z$ (m/s)','interpreter','latex');
% legend('CCM','RCCM','RCCM-P','Orientation','Horizontal');
axis equal
goodplot([6 5])
subplot(2,2,4);hold on;
ellipse(controller_ccm.tube_phi, controller_ccm.tube_phidot,0,0,0,color{5},[],1);
ellipse(controller_rccm_pos.tube_phi, controller_rccm_pos.tube_phidot,0,0,0,'r',[],1);
ellipse(controller_rccm.tube_phi, controller_rccm.tube_phidot,0,0,0,'g',[],1);
xlabel('$\phi$ (rad)','interpreter','latex');
ylabel('$\dot{\phi}$ (rad/s)','interpreter','latex');
% legend('CCM','RCCM','RCCM-P','Orientation','Horizontal');
% axis equal
ylim([-15 15]);
xlim([-5 5]);
goodplot([6 5])

fig_name = 'tube_4';
if print_file == 1
    savefig([fig_name '.fig']);
    print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
end


    
        
