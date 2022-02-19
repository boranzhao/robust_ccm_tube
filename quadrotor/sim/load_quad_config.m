% Source: https://github.com/stanfordASL/RobustMP

%% Constants

Jx = 0.03; Jy =  0.04; Jz =  0.1;
global mq Jq;

Jq = diag([Jx;Jy;Jz]);
mq = 0.9574;

%% Dynamics

plant.n = 12; plant.nu = 4;
plant.nc= 9; plant.nuc = 3;
%dynamics for CCM controller
%pos,vel,thrust,roll,pitch

%full-state dynamics
%pos,vel, roll, pitch, yaw, omg
plant.f_sim = @(x) [x(4:6);
          [0;0;plant.g];
          R_eul(x(7:9))*x(10:12);
          -Jq\cross(x(10:12),Jq*x(10:12))];
      
plant.B_sim =  @(x) [zeros(3,4);
          -rot_matrix(x(7),x(8),x(9))*[0;0;1],zeros(3);
          zeros(3,4);
          zeros(3,1), Jq\eye(3)];  

plant.Bw_sim = [zeros(3);
       eye(3);
       zeros(6,3)];

%% Setup lower-level controller
%Angular rate P controller gain
global kp_omg;
kp_omg = 6;


%% Setup Sim environment
close all;

%Set world space dimensions
World_dim = [8,8,8];

%set initial state
pos_init = [0;0;1]; %in NWU frame
test_state = [pos_init(1);-pos_init(2);-pos_init(3); %NED frame
              zeros(9,1)];
thrust_init = plant.g;          
xc_init = [test_state(1:6);
           thrust_init;
           test_state(7:9)];
          
%Define obstacle environment (if not already done) in NWU frame
if (new_obs)
    define_obstacles;
    save('create_env/Obstacles.mat','tree_obstacles','tower_obstacles',...
          'obstacles','obstacles_infl','obstacles_coll','n_obstacles');
else
    load('create_env/Obstacles.mat');
    fig = figure();
    plot_all_obstacles();
end  

%define goal
goal_box = [World_dim-[1,1,1];
            World_dim];
goal_V =  corner2Hull(goal_box,1);
goal = Polyhedron('V',goal_V);
figure(fig)
plot3dObstacles(goal_box,'b',0.2);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
