function plot_quad_movie2(xnomTraj,solve_t,xTraj,~,controller,goal,...
                         tree_obstacles,tower_obstacles,obstacles_infl,World_dim,tube_xyz)
color = {'k','b',[0 0.5 0],'r',[0.8 0.9 0.9741],[0.8 0.98 0.9],[1 0.8 0.8],[0.7, 0.7 1]};

close all;

fig = figure();
% create the environment
plot3dObstacles(tree_obstacles,'g',0.6); 
plot3dObstacles(tower_obstacles,[0.5,0.5,0.5],1.0);
% plot3dObstacles(obstacles_infl,'r',0.01,1);
plot3dObstacles([0 0 0; World_dim],'k',0);
patch([0,World_dim(1),World_dim(1),0],...
      [0, 0, World_dim(2), World_dim(2)],...
      [0, 0,  0, 0],'FaceColor',[0,0.5,0],'FaceAlpha',0.3);
axis equal; axis tight 
view(3);
hold on
% plot the tubes
if controller.type == CtrlDesignOpts.rccm
    scale = 25; tube_interval = 0.08; vel_cst= 1.3;
elseif controller.type == CtrlDesignOpts.ccm
    scale = 50; tube_interval = 0.15; vel_cst= 1.3;
end
numtubes = 0;
for i=1:length(solve_t)/scale
   center = [xnomTraj(i*scale-1,1),-xnomTraj(i*scale-1,2),-xnomTraj(i*scale-1,3)];
   if i== 1 ||  norm(center-center_prev) > tube_interval || (norm(xnomTraj(i*scale-1,4:6))<= vel_cst && mod(i,10)==0 )
       [Xsurf,Ysurf,Zsurf] = ellipsoid(center(1),center(2),center(3),tube_xyz,tube_xyz,tube_xyz); 
       surf(Xsurf,Ysurf,Zsurf,'FaceColor',color{8},'FaceAlpha',0.1,'EdgeColor','none'); %'FaceLighting','flat'
       numtubes = numtubes +1;
   else
       aa = 1;
   end
   center_prev = center;
end
numtubes
plot3(xTraj(:,1),-xTraj(:,2),-xTraj(:,3),'r-','linewidth',1);hold on
plot3(xnomTraj(:,1),-xnomTraj(:,2),-xnomTraj(:,3),'k--','linewidth',1); 

goal.plot('color','green','alpha',0.2); 
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); 

%% FPV view:
fpv_view = 1;
set(gca,'CameraViewAngleMode','manual')
camproj('perspective')
camva(60)

%%
dt = solve_t(2)-solve_t(1);
t_step = round((1/20)/dt);
t_step_b = round(0.3/dt);
T_b = 2; %T_b second lookahead for bound
n_b = 1+round((T_b/0.3));

i = 1;
pos = xTraj(i,1:3)';
att = xTraj(i,7:9);

%% create objects

% hold on
% if (~fpv_view)
    h_b = [];
%     for j = 1:n_b
%         j_idx = min(i + (j-1)*t_step_b, length(solve_t));
%         h_b(j) = proj_Ellipse([1:3],M,1,[x_des(j_idx);-y_des(j_idx);-z_des(j_idx)],30,'b',0);
%     end
% end
[h_quad,~] = plot_quad_obj(pos,att,fig);
rot = rot_matrix(att(1),att(2),atan2(xTraj(i,5),xTraj(i,4)));
hold off

% view(3);
view(52.5,30);

Flip_M = [1,0,0;0,-1,0;0,0,-1];

%set view
if (fpv_view)
    if controller.type == CtrlDesignOpts.rccm
        pos_tail = [-1;0.1;-0.5]; 
    elseif controller.type == CtrlDesignOpts.ccm
        pos_tail = [-1;0.1;-0.5]; 
    end
    pos_tail0 = pos_tail;
    pos_view = [0.1;0;0.0];
     c_pos = (Flip_M*(pos+rot*pos_tail))';
    campos(c_pos);
%         c_target =(Flip_M*pos)'+(Flip_M*rot*pos_view)';
    c_target =(Flip_M*(pos+rot*pos_view))';
    camtarget(c_target);
else
    %iso-view
    campos((Flip_M*pos)'+[-1.5,-2.5,1]);
    camtarget((Flip_M*pos)');
end
drawnow;

%% Setup Movie
record_vid = 1;
if controller.type == CtrlDesignOpts.rccm
    video_name = 'Quad_sim_rccm';
else
    video_name = 'Quad_sim_ccm';
end
if (record_vid)
    writerObj = VideoWriter(video_name);
    writerObj.FrameRate = 1/(t_step*dt);
    writerObj.Quality = 100;
    open(writerObj);
    set(gcf, 'renderer', 'zbuffer')
end

%% Record
%record
view(52.5,30);camva(9.2)
if (record_vid)
    thisFrame = getframe(gcf);
    for i = 1:50
        writeVideo(writerObj, thisFrame);
    end
end
camva(60)

keyboard;
for i = 1:t_step:length(solve_t)
    pos = xTraj(i,1:3)';
    att = xTraj(i,7:9);
    
    
    delete(h_quad); 
    hold on
    %update quad and bound
%     if (~fpv_view)
%         delete(h_b);
%         h_b = [];
%         for j = 1:n_b
%             j_idx = min(i + (j-1)*t_step_b, length(solve_t));
%             h_b(j) = proj_Ellipse([1:3],M,1,[x_des(j_idx);-y_des(j_idx);-z_des(j_idx)],30,'b',0);
%         end
%     end
    [h_quad,~] = plot_quad_obj(pos,att,fig);
    rot = rot_matrix(att(1),att(2),atan2(xTraj(i,5),xTraj(i,4)));    
    hold off

    %update view
    if controller.type == CtrlDesignOpts.rccm
        if i>1050 && i<=1300
            pos_tail = [-1;-0.1;-0.3];
        elseif i>1300 && i<=1600
            pos_tail = [-1;0;-0.3];
        end
    elseif controller.type == CtrlDesignOpts.ccm
        if i>1400 && i<=2094
            pos_tail = [-1;-0.1;-0.2];
        elseif i>2094 
            pos_tail = [-1;0;-0.3];
        end
    end
        
    if (fpv_view)
        c_pos = (Flip_M*(pos+rot*pos_tail))';
        campos(c_pos);
%         c_target =(Flip_M*pos)'+(Flip_M*rot*pos_view)';
        c_target = (Flip_M*(pos+rot*pos_view))';
        camtarget(c_target);
    else
        %iso-view
        campos((Flip_M*pos)'+[-1.5,-2.5,1]);
        camtarget((Flip_M*pos)');
    end
    drawnow;
    
    %record
    if (record_vid)
        thisFrame = getframe(gcf);
        writeVideo(writerObj, thisFrame);
    end
end 
% get an iso view for the last frame
if (record_vid)
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); 
    view(52.5,30);camva(9.2);
    thisFrame = getframe(gcf);
    for i = 1:20
        writeVideo(writerObj, thisFrame); 
    end
    
    pause(0.5);
    close(writerObj);
end
end