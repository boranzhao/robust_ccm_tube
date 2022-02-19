
if (gen_new_traj)    
    %% Setup FMT* to get geometric path (NWU frame)    
    %connection radius
    plan_radius = 2;
    
    %Generate samples and determine nearest neighbors
    N_nodes = 500;
    nodes = generateHaltonSamples(3,N_nodes);
    nodes = nodes.*repmat(World_dim,N_nodes,1);
    %get nearest neighbor graph
    disp('Computing FMT NN graph');
    [N_nn,C_nn] = compute_NN(nodes,plan_radius);
    disp('Done');
    
    vis_nn = 0;
    if (vis_nn)
        figure(fig)
        for i = 1:N_nodes
            NN_i = find(N_nn(i,:));
            for j = 1:length(NN_i)
                plot3(nodes([i;NN_i(j)],1),nodes([i;NN_i(j)],2),nodes([i;NN_i(j)],3),'bo-','linewidth',1.5);
            end
        end
    end    
    %get geometric path
    [FMT_time, geom_path] =  FMTStar_Geom3D(N_nodes,nodes,N_nn,C_nn,...
        plan_radius,obstacles_coll,pos_init',goal,fig);
    
    %% Run Richter smoothing algorithm (NWU frame)
    [pos_coef,vel_coef,acc_coef,jer_coef,T_seg,coeff] = min_snap(geom_path,obstacles_coll,fig);
    save(file_quad_traj,'pos_coef','vel_coef','acc_coef','jer_coef','T_seg','coeff');
    
else
    load(file_quad_traj);
end

%% Extract nominal state and control trajectories
dt = 0.001;
[MP_t,MP_state,MP_ctrl, Thrust_nom] = generate_quad_traj_10(dt,mq,plant.g,[pos_init(1);
                                                                 -pos_init(2);
                                                                 -pos_init(3)],file_quad_traj);

%% Plot
figure(fig)
plot3(MP_state(:,1),-MP_state(:,2),-MP_state(:,3),'k-','linewidth',2.0);
% for i = 1:(0.05/(MP_t(2)-MP_t(1))):length(MP_t)
%     proj_Ellipse([1:3],M_ccm_pos_infl,1,[MP_state(i,1);-MP_state(i,2);-MP_state(i,3)],30,'b');
% end
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
view(52.5,30);

figure()
subplot(2,1,1)
plot(MP_t,MP_state(:,8:9)*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Attitude [deg]');
legend('Roll','Pitch');

subplot(2,1,2)
plot(MP_t,MP_ctrl(:,2:4)*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Angular rate [deg/s]');
hl = legend('$\dot{\phi}$','$\dot{\theta}$','$\dot{\psi}$');
set(hl,'interpreter','latex');

figure()
subplot(2,1,1)
plot(MP_t, MP_state(:,4:6),'linewidth',2);
grid on
xlabel('Time [s]'); legend('v_x','v_y','v_z');

subplot(2,1,2)
plot(MP_t, MP_ctrl(:,1),'linewidth',2);
grid on
xlabel('Time [s]'); legend('Thrust derivative');


