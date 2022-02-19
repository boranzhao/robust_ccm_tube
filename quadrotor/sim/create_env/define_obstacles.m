%% setting
% tube size associated with the controller 
tube_infl = controller.tube_gain_xyz*dist_config.w_max;

% define numbers and positions of tress and towers
halton_pos = generateHaltonSamples(2,100);
numTrees = 8; numTowers = 6;
% randomly generating the location of the trees and towers
% tree_pos_xy = halton_pos(randperm(100,numTrees),:);
% tower_pos_xy = halton_pos(randperm(100,numTowers),:);
% tree_pos_xy = tree_pos_xy.*repmat(World_dim(1:2),numTrees,1);
% tower_pos_xy = tower_pos_xy.*repmat(World_dim(1:2)-tower_width,numTowers,1);

% manually generate the location of obstacles
tree_pos_xy = [... 
    1.0625    4.5679;
    2.0500    0.8877;
    2.6250    7.5062;
    4.3750    4.2469;    
    6.1875    0.6255;
    6.5625    7.0123;
    6.6250    3.8395;
    7.8125    6.5514];

tower_pos_xy = [...
        0.0789    6.8519;
    0.358    1.0207;
    2.6966    0.8548;
    2.7952    5.5760;

    4.5891    5.1126;
    6.6495    4.3719];

%% create obstacles
% create tree obstacles (ll/ur format)
pine_rad = 1/4;
pine_height = 3;
tree_rad = 1.5/4;
tree_height = 3;
numTop = 6;
numDisc = 4;

tree_bounds = [];
tree_obstacles = [];

for i = 1:numTrees
    if rand(1) > 0.5
        tree_obstacles = [tree_obstacles;
                         createPineTreeObs([tree_pos_xy(i,:) 0], pine_rad, pine_height, numTop, numDisc)];
        tree_bounds = [tree_bounds;
                     tree_pos_xy(i,1)-pine_rad,tree_pos_xy(i,2)-pine_rad,0;
                     tree_pos_xy(i,1)+pine_rad,tree_pos_xy(i,2)+pine_rad,pine_height];
    else
        tree_obstacles = [tree_obstacles;
                         createTreeObs([tree_pos_xy(i,:) 0],tree_rad, tree_height, numTop, numDisc)];
        tree_bounds = [tree_bounds;
                     tree_pos_xy(i,1)-tree_rad,tree_pos_xy(i,2)-tree_rad,0;
                     tree_pos_xy(i,1)+tree_rad,tree_pos_xy(i,2)+tree_rad,tree_height];
    end
end

% Create towers (ll/ur format)
tower_width = (9/4-0.1)*(World_dim(1)/15); tower_height = World_dim(3);
tower_obstacles = [];
for i = 1:numTowers
    tower = createBoxObs([tower_pos_xy(i,1),tower_pos_xy(i,2),0],[tower_width,tower_width,tower_height]);
    tower_obstacles = [tower_obstacles; tower];
end

%% Inflate
%tightest bounding box form
obstacles = [tree_bounds; tower_obstacles];
n_obstacles = size(obstacles,1)/2;
size_infl = [0.0, 0.0, 0.0];

%add inflation by size of quad
obstacles_infl = obstacles  - kron(repmat([size_infl(1:2),0],n_obstacles,1),[1;0]) +...
                              kron(repmat(size_infl,n_obstacles,1),[0;1]);                       

%add inflation by size of tube bound
obstacles_coll = obstacles_infl  - kron(repmat([tube_infl(1:2),0],n_obstacles,1),[1;0]) +...
                                   kron(repmat(tube_infl,n_obstacles,1),[0;1]);
% Plot 
% fig = figure();
% plot_all_obstacles();
% plot3dObstacles(obstacles_coll,'r',0.0,1);
% view(2)

fig = figure();
plot_all_obstacles();
plot3dObstacles(obstacles_coll,'r',0.0,1);
view(52.5,30)
