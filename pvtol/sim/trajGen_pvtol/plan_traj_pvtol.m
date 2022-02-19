function [soln] =plan_traj_pvtol(plant,controller,trajGen_config)
x0 = trajGen_config.x0;
xF = trajGen_config.xF;
duration = trajGen_config.duration;
u_bnd = trajGen_config.u_bnd;
x_bnd = trajGen_config.x_bnd;

% initial guess of the solution: important to find a feasible solution    
initGuess.time = [0,duration/2 duration];
initGuess.state = [x0, (x0+xF)/2 xF];
initGuess.control = plant.g*plant.m*ones(2,3)/2;
if trajGen_config.include_obs == 1
    % ----- formulate the state constraints for collision avoidance 
     if trajGen_config.include_tube == 1
        tube_xz = controller.tube_gain_xz*w_max;
    else
        tube_xz = 0;
    end
    
    stateCst = @(x)  cst_avoid_obstacle(x,trajGen_config.obs,tube_xz,plant.nu);
    
    if ~((controller.type == CtrlDesignOpts.rccm && controller.opt_pos_dev ==1)) % indicates 'rccm-p' 
%         initGuess.time = [0, 3  6  10 15];
% %         initGuess.state = [x0, [-3 5 zeros(1,4)]'  xF];
%         initGuess.state =[x0, [-3 3 zeros(1,4)]' [-3 6 zeros(1,4)]' [1 10 zeros(1,4)]'  xF];
%         initGuess.control = plant.g*plant.m*ones(2,5)/2;
        
        if controller.type == CtrlDesignOpts.ccm
            initGuess.time = [0,  8  15];
            initGuess.state = [x0, [15 1 zeros(1,4)]'  xF];
            initGuess.control = plant.g*plant.m*ones(2,3)/2;
        elseif controller.type == CtrlDesignOpts.rccm
            initGuess.time = [0,  7  10 12];
            initGuess.state = [x0, [8.5 3 zeros(1,4)]' [8.5 7 zeros(1,4)]'  xF];
            initGuess.control = plant.g*plant.m*ones(2,4)/2;
        end
        
    end
    soln = trajOpt_pvtol(plant,x0,xF,duration,u_bnd,x_bnd,initGuess,stateCst);
else
    soln = trajOpt_pvtol(plant,x0,xF,duration,u_bnd,x_bnd,initGuess);
end    
fprintf(1,'cost = %.2f\n',soln.info.objVal);
