lambdas = [0.2 0.4 0.6 0.8 1 1.2 1.4 1.6]';

ccm_data = 'ccm_0.8_plim_0.33pi.mat';
rccm_data = 'rccm_1.4_wmax_1_plim_0.33pi.mat';
rccm_pos_data = 'rccm_1.2_wmax_1_plim_0.33pi_pos.mat';

w_max = 1;  % amplitude of disturbance

print_file = 0;
line_wd = 0.8;
color = {'k','b',[0 0.5 0],'r',[0.8 0.9 0.9741],[0.8 0.98 0.9],[1 0.8 0.8]};
c_g = 'g';[0 0.5 0];
% 1st column: all states; 2nd column: xz; 3nd column: u
ccm_tubes = [35.44 10.07 15.84; 
        20.8 5.97 8.5;
        16.2 4.77 5.24;
        14.07 4.26 4.65;
        13.25 5.44 9.53;
        14.95 5.15 26.14;
        24 9.53 105.67;
        41.26 24.95 516.81];
rccm_tubes = [3.38 3.30 3.27;
    2.25 2.13 2.18;
    1.92 1.77 1.84;
    1.77 1.57 1.71;
    1.71 1.48 1.66;
    1.69 1.43 1.66;
    1.705 1.41 1.65;
    1.75 1.41 1.68];
rccm_p_tubes = [46.66 1.255 1.254;
    26.15 0.787 0.788;
    20.27 0.653 0.653;
    17.18 0.60 0.60;
    15.72 0.57 0.57;
    14.76 0.563 0.56;
    13.78 0.570 0.571;
    13.11 0.589 0.588];

close all;
figure;
subplot(3,1,1); 
semilogy(lambdas,ccm_tubes(:,1),'bo--','linewidth',line_wd);hold on;
semilogy(lambdas,rccm_tubes(:,1),'*-.','color',c_g,'linewidth',line_wd);
semilogy(lambdas,rccm_p_tubes(:,1),'rs-','linewidth',line_wd);
ylabel('$\alpha_x \bar w$', 'interpreter','latex');
legend('CCM','RCCM','RCCM-P','Orientation','Horizontal');
ylim(ylim+[-0.2 100]);
yticks([10 100]);
grid on;
goodplot([6 5]);

subplot(3,1,2);
semilogy(lambdas,ccm_tubes(:,2),'bo--','linewidth',line_wd); hold on;
semilogy(lambdas,rccm_tubes(:,2),'*-.','color',c_g,'linewidth',line_wd);
semilogy(lambdas,rccm_p_tubes(:,2),'rs-','linewidth',line_wd);
ylabel('$\alpha_{p_x,p_z}\bar w$', 'interpreter','latex');
ylim([0.05 110]);
yticks([0.1 10 100]);
grid on;
goodplot([6 5]);

subplot(3,1,3); 
semilogy(lambdas,ccm_tubes(:,3),'bo--','linewidth',line_wd);hold on;
semilogy(lambdas,rccm_tubes(:,3),'*-.','color',c_g,'linewidth',line_wd);
semilogy(lambdas,rccm_p_tubes(:,3),'rs-','linewidth',line_wd);
ylabel('$\alpha_u \bar w$', 'interpreter','latex');
xlabel('$\lambda$ ($\hat \lambda$)','interpreter','latex');
ylim([0.05 1100]);
yticks([0.1 10 100]);
grid on;
goodplot([6 5]);
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


    
        
