%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Synthesis of a (robust) CCM for control of a 3D quadrotor

%  Author: Pan Zhao, UIUC, Advanced Controls Research Lab,
%  panzhao2@illinois.edu
  
%  Last update: Feb 18, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %clc; 
yalmip('clear')
addpath('metrics_search');
addpath('../../../rccm/metric_search_offline/');
addpath('../../../rccm/control_law_online/');
addpath('../../../utilities/');
%% general settings 
save_rsts = 1;                                  % whether to save the results to a file
controller.opt_pos_dev = 1;                     % whether to focus on mitigate the effect of disturbance on the position states

%% load plant
load_plant;

%% settings for searching CCM & RCCM 
controller.type = CtrlDesignOpts.rccm;          % ccm
controller.ccm_law_form = ...
    CtrlDesignOpts.ccm_min_norm;                % ccm_min_norm, ccm_integration
controller.ccm_mat_ineq_form = ...
    CtrlDesignOpts.ccm_mat_ineq_use_B_perp;     % For ccm approach: ccm_mat_ineq_use_rho, or ccm_mat_ineq_use_B_perp.

lambda = 0.8;                                  % a line search or bisection search should be performed for lambda
consider_state_set = 1;                         % whether to consider a compact set for the states when formulating the constraints

if controller.type ==  CtrlDesignOpts.rccm
    controller.ccm_law_form = ...
        CtrlDesignOpts.ccm_integration;                % ccm_min_norm, ccm_integration
    w_lower_bound = .01;
elseif controller.type ==  CtrlDesignOpts.ccm
    w_lower_bound = 1;
end
w_states_index = [7 8 9];
w_order = 3; y_order = 3;
w_states = x(w_states_index);
[w_poly, w_poly_mat] = monomials(w_states,0:w_order);  % monomials of T, phi, theta, 

w_poly= w_poly(w_poly_mat(:,1)<=2);%only keep monomials whose degree w.r.t T(thrust) <=2
w_poly_mat = w_poly_mat(w_poly_mat(:,1)<=2,:);

% if controller.type ==  CtrlDesignOpts.rccm
%     [w_poly, w_poly_mat] = monomials([w_states;u],0:w_order);  % monomials of T, phi, theta, and u
%     II = find(w_poly_mat(:,1)<=1 & sum(w_poly_mat(:,length(w_states)+1:end),2)<=1);
%     w_poly= w_poly(II);%only keep monomials whose degree w.r.t T(thrust) <=2, w.r.t u<=1
%     w_poly_mat = w_poly_mat(II,:);
% end
n_monos_W = length(w_poly); 

% ------------state constraints imposed when searching CCM/RCCM-----------
phi_lim = pi/3;  
theta_lim = pi/3; 
T_lim = [0.5 2]*plant.g;    
Tdot_lim = 5*plant.g;
phidot_lim = pi;
thetadot_lim = pi; 
state_set.wdist_lim = 1;                        % w (disturbance)
state_set.ccm_states =  w_states;               % states invovled in the CCM condition
state_set.T_lim = T_lim;
state_set.phi_lim = phi_lim;                
state_set.theta_lim = theta_lim;
state_set.Tdot_lim = Tdot_lim;
state_set.phidot_lim = phidot_lim;
state_set.thetadot_lim = thetadot_lim;
state_set.w_states = w_states;
% -------------------------------------------------------------------------
     
dw_poly_dx = diff(w_poly,w_states);
if controller.type == CtrlDesignOpts.ccm 
    dw_poly_dt = dw_poly_dx*(f(w_states_index)+plant.B(w_states_index,:)*u);
%     dw_poly_dt = dw_poly_dx*(f(w_states_index));
elseif controller.type == CtrlDesignOpts.rccm
    % dW_poly_dt could depend on w
    dw_poly_dt = dw_poly_dx*(f(w_states_index)+plant.B(w_states_index,:)*u+Bw(w_states_index,:)*wdist);  
%     dw_poly_dt = dw_poly_dx*(f(w_states_index));  
end

% create the functions
w_poly_fcn = mss2fnc(w_poly,x,randn(length(x),2));
dw_poly_dt_fcn = mss2fnc(dw_poly_dt,[x;u],randn(length([x;u]),2));

if controller.type == CtrlDesignOpts.rccm && ...
        controller.ccm_mat_ineq_form == CtrlDesignOpts.ccm_mat_ineq_use_rho && n_monos_W>1
    error('Please use the CCM condition involving B_perp; otherwise, the condition will depends on u!');
end

% following [Singh et al 2019, IJRR], W is selected to be
% [W1  W2;
%  W2' W3];
prog = spotsosprog;

if controller.type  == CtrlDesignOpts.ccm 
    %---------------- W1 is independent of any states ---------------------
    W1_list = cell(n_monos_W,1);
    W2_list = cell(n_monos_W,1);
    W3_list = cell(n_monos_W,1);
    W_list  = cell(n_monos_W,1);
    [prog, W1_list{1}] = prog.newSym(6);
    [prog, W2_list{1}] = prog.newFree(6,3);
    [prog, W3_list{1}] = prog.newSym(3);
    W_list{1} = [W1_list{1}, W2_list{1};
                W2_list{1}', W3_list{1}];
    for i = 2:length(w_poly) 
        W1_list{i} = zeros(6);
        [prog, W2_list{i}] = prog.newFree(6,3);
        [prog, W3_list{i}] = prog.newSym(3);    
        W_list{i} = [W1_list{i}, W2_list{i};
                     W2_list{i}', W3_list{i}];
    end
    % ------------------------------------------------------------------
elseif controller.type  == CtrlDesignOpts.rccm  
    % ---------------- W1 is dependent on some states -------------------
    W_list  = cell(n_monos_W,1);
    for i = 1:length(w_poly)
        [prog, W_list{i}] = prog.newSym(n);
    end
    % ------------------------------------------------------------------
end
% create the function handle
W_exec = 'W_eval = @(ml)';
for i = 1:length(w_poly)
    if i<length(w_poly)
        W_exec = strcat(W_exec,sprintf('W_list{%d}*ml(%d) +',i,i));
    else
        W_exec = strcat(W_exec,sprintf('W_list{%d}*ml(%d);',i,i));
    end
end
eval(W_exec);    

controller.w_lower_bound = w_lower_bound;
state_set.consider_state_set = consider_state_set;
state_set.w_states = w_states;
state_set.w_states_index = w_states_index;

%% Search CCM/RCCM
tic;
if controller.type == CtrlDesignOpts.ccm  % search for CCM
    controller.lambda = lambda; 
    [cond_num_W,W_bar,w_lower,w_upper,W_coef,w_poly] = ccm_synthesis_quad(prog,plant,controller,W_list,W_eval,w_poly,w_poly_fcn,dw_poly_dt_fcn,lambda,state_set); 
    controller.w_upper = w_upper;
    controller.w_lower = w_lower;    
    controller.W_bar = W_bar;
elseif controller.type == CtrlDesignOpts.rccm % search for robust CCM minimizing PPG
    controller.lambda = lambda;  
    [y_poly, y_poly_mat] = monomials([w_states],0:y_order); 
    II = find(y_poly_mat(:,1)<=2 & sum(y_poly_mat(:,length(w_states)+1:end),2)<=1);
    y_poly= y_poly(II);%only keep monomials whose degree w.r.t u<=1
    y_poly_mat = y_poly_mat(II,:);
    n_monos_y = length(y_poly);
    Y_list  = cell(n_monos_y,1); 
    for i = 1:length(y_poly)
        [prog, Y_list{i}] = prog.newFree(nu,n);    
    end   
   
    Y_exec = 'Y_eval = @(ml)';
    for i = 1:length(y_poly)
        if i<length(y_poly)
            Y_exec = strcat(Y_exec,sprintf('Y_list{%d}*ml(%d) +',i,i));
        else
            Y_exec = strcat(Y_exec,sprintf('Y_list{%d}*ml(%d);',i,i));
        end
    end
   
    % create a few function handles
    y_poly_fcn = mss2fnc(y_poly,x,randn(length(x),2));  
    
    % create the function handle
    eval(W_exec); eval(Y_exec);     
    
    tic;
    [alpha_opt,mu_opt,W_coef,Y_coef,w_poly,y_poly] = rccm_synthesis_quad(prog,plant,controller,W_list,W_eval,Y_list,Y_eval,w_poly,w_poly_fcn,dw_poly_dt_fcn,y_poly,y_poly_fcn,lambda,state_set);
    toc;
       
    lambda_opt = lambda;
%     fprintf('The smallest gamma is %.3f obtained at lambda= %.3e\n',alpha_opt,lambda_opt);       
    fprintf('RCCM, lambda = %.2f, tube gain = %.4f\n',lambda_opt,alpha_opt);
    
end
toc;
%% extract functions -------------------------
extract_funcs;

controller.W_fcn = W_fcn;
controller.dW_dxi_fcn = dW_dxi_fcn;
controller.dW_dt_fcn = dW_dt_fcn;
if controller.type == CtrlDesignOpts.rccm
    tube_gain = alpha_opt; 
    controller.mu = mu_opt;    
    controller.Y_fcn = Y_fcn;
    controller.lambda = lambda_opt;
    controller.alpha = alpha_opt;
end

%% check CCM conditions, compute tubes 
compute_tubes;

%% save data
if save_rsts == 1
    if controller.type == CtrlDesignOpts.ccm 
        if n_monos_W == 1
            file_name = ['ccm_const_' num2str(lambda0) '.mat'];
        else
            file_name = ['ccm_' num2str(lambda0) '.mat'];
        end        
    elseif controller.type == CtrlDesignOpts.rccm
        if  controller.opt_pos_dev == 0
            file_name = ['rccm_' num2str(lambda0,3) '.mat'];
        else
            file_name = ['rccm_' num2str(lambda0,3) '_pos.mat'];
        end
    end
    save(file_name,'plant','controller','state_set');
end
%% generate the c codes for accelerating geodesic computation
% To use the generated codes, copy .mex and .mat files to the sim folder

% parameters used in the pseudospectral method for geodesic computation
geodesic_setting_for_codegen.D = 3; 
geodesic_setting_for_codegen.N = 6;
answer = questdlg('Do you want to generate the C codes for accelerating geodesic computation used for determining the control law?','Question for code generation','Yes','No','No');
switch answer 
    case 'Yes'
        generate_code_for_geodesic_cal(plant.n,plant.nu,plant.nw,geodesic_setting_for_codegen);        
        save('geodesic_setting_for_codegen.mat','geodesic_setting_for_codegen');
    case 'No'
end
