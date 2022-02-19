%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Synthesis of (robust) CCM for control of a planar VTOL

%  Author: Pan Zhao, UIUC, Advanced Controls Research Lab,
%  panzhao2@illinois.edu
%  Codes for the paper:
%  P. Zhao, et al. Tube-certified trajectory tracking for nonliner systems
%  with robust control contraction metrics. Submitted to IEEE Robotics and
%  Automation Letters, 2022. 
%  Last update: Feb 16, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
yalmip('clear');
%% add paths
addpath('../../../rccm/metric_search_offline/');
addpath('../../../rccm/control_law_online/');
addpath('../../../utilities/');
%% general settings 
save_rsts = 1;                              % whether to save the results to a file 
controller.opt_pos_dev = 1;                 % whether to focus on position states only when optimizing the tube size

%% load plant
load_plant;  

%% settings for searching CCM & RCCM
controller.type = CtrlDesignOpts.rccm;      % lqr,ccm, rccm
controller.ccm_law_form = ...
    CtrlDesignOpts.ccm_min_norm;            % ccm_min_norm, ccm_integration
controller.ccm_mat_ineq_form = ...
    CtrlDesignOpts.ccm_mat_ineq_use_B_perp; % For ccm approach: ccm_mat_ineq_use_rho, or ccm_mat_ineq_use_B_perp . 
lambda = 1.2;                               % a line search or bisection search should be performed for lambda
consider_state_set = 1;                     % whether to consider a compact set for the states when formulating the constraints
biSearch_lambda = 0;                        % 1 for conducting a bisection search for lambda in designing the robust CCM controller 
maxNum_bs = 10;                             % maximum number of trials for bisection-search 
lambdaSearchRange = [0.1 10^1];             % search range for lambda
W_lower_bound = 1e-2;                       % lower bound for the dual metric W
Wstates_index = [3 4];

% ------------state constraints imposed when searching CCM/RCCM-----------
p_lim = pi/3;  % phi
pd_lim = pi/3; % phi_dot 
vx_lim = 2;    % vx
vz_lim = 1;    % vz
w_lim = 1;     % w (disturbance)
state_set.box_lim = [p_lim^2-x(3)^2; vx_lim^2-x(4)^2; pd_lim^2-x(6)^2;  vz_lim^2-x(5)^2]*0.001;
state_set.num_consts_4_W_states = 2;        % #constraints from box_lim that involve states on which the metric W depends
state_set.other_lim_states = [x(6);x(5)]; 
state_set.lagrange_deg_W = 4;               % degree of the lafor the bounds of W
state_set.lagrange_deg_ccm = 4;             % for ccm condition
state_set.p_lim = p_lim;
state_set.pd_lim = pd_lim;
state_set.vx_lim = vx_lim;
state_set.vz_lim = vz_lim;
state_set.w_lim = w_lim;
% -------------------------------------------------------------------------

W_states = x(Wstates_index);
v_W = monolist(W_states,4);         % monomials of W_states up to degree 4
n_monos_W = length(v_W);
dv_W_dx = jacobian(v_W,W_states);

W_coef = sdpvar(n,n,n_monos_W);     
W = zeros(n);
for i=1:n_monos_W
    W = W+ W_coef(:,:,i)*v_W(i);
end
if controller.type == CtrlDesignOpts.ccm
    dv_W_dt = dv_W_dx*f(Wstates_index);
elseif controller.type == CtrlDesignOpts.rccm
    % dv_W_dt could depend on w
    dv_W_dt = dv_W_dx*(f(Wstates_index)+Bw(Wstates_index)*w);  
end

dW_dt = zeros(n);
for i=1:n_monos_W
    dW_dt = dW_dt+ W_coef(:,:,i)*dv_W_dt(i); 
end

if controller.ccm_mat_ineq_form == CtrlDesignOpts.ccm_mat_ineq_use_rho 
    [rho,c_rho,v_rho] = polynomial(W_states,2);
    controller.rho = rho;
    controller.c_rho = c_rho; 
    controller.v_rho = v_rho;
    n_monos_rho = length(v_rho);
end

controller.W_lower_bound = W_lower_bound;
state_set.consider_state_set = consider_state_set;
state_set.W_states = W_states;
state_set.W_states_index = Wstates_index;

if controller.type == CtrlDesignOpts.rccm
   state_set.box_lim = [state_set.box_lim; w_lim^2-w^2]; 
end

%% Search CCM/RCCM 
if controller.type == CtrlDesignOpts.ccm  % search for CCM
    % YALMIP (SOS) and Mosek are used to solve the SOS problem
    controller.lambda = lambda;%     
    paras_W = W_coef(:);    
    [cond_num_W,w_upper,w_lower,W_bar,max_res] = ccm(plant,controller,W,dW_dt,paras_W,lambda,state_set);

    controller.w_upper = w_upper;
    controller.w_lower = w_lower;    
    controller.W_bar = W_bar;
elseif controller.type == CtrlDesignOpts.rccm % search for robust CCM minimizing Linf gain
    controller.lambda = lambda;       
    Y_coef = sdpvar(nu,n,n_monos_W);
    Y = zeros(nu,n);
    for i=1:n_monos_W
        Y = Y+ Y_coef(:,:,i)*v_W(i);
    end 
    vars_WY = [W_coef(:); Y_coef(:)];
    if biSearch_lambda    
        % implement a bisection search to search for the best lambda
        Result_bs = ones(maxNum_bs+2,2)*inf; % 1st column for lambda, second column for Linf gain bounds
        Result_bs([1 2],1) = lambdaSearchRange';
        for i = 1:maxNum_bs+2
            if i <=2
                lambda = Result_bs(i,1);           
            else
                % sort the lambda value in ascending order 
                [~,index] = sort(Result_bs(:,1));
                Result_bs = Result_bs(index,:);
                % sort the Linf gain bounds in ascending order 
                [~,index] = sort(Result_bs(:,2));
                % find the lambda value corresponding to the smallest and
                % second smallest bounds and compute their mean
                lambda = sum(Result_bs(index(1:2),1))/2;
                Result_bs(i,1) = lambda;
            end
            [alpha_opt,mu_opt,max_residual] = rccm(plant,controller,W,dW_dt,Y,vars_WY,lambda,x1,state_set); 
            Result_bs(i,2) = alpha_opt;
        end
        Result_bs
        [alpha_opt, index] = min(Result_bs(:,2));
        lambda_opt = Result_bs(index,1);
    else
       [alpha_opt,mu_opt] = rccm(plant,controller,W,dW_dt,Y,vars_WY,lambda,state_set); 
       lambda_opt = lambda;
    end     
    fprintf('RCCM, lambda = %.2f, tube gain = %.4f\n',lambda_opt,alpha_opt);
    tube_gain = alpha_opt; 

% ------------------------------ extract Y_fcn-----------------------------
    Y_coef = value(Y_coef);  
    x = x_store; % must ensure that v_W and s contain "x" instead of "x_store"

    Y_fcn = zeros(nu,n);
    for i=1:n_monos_W
        Y_fcn = Y_fcn+ Y_coef(:,:,i)*v_W(i);
    end
    Y_fcn = clean(Y_fcn,1e-10);
    s = sdisplay(Y_fcn);
    syms x [n 1]
    syms Y_fcn [nu nu]
    for i=1:nu
        for j=1:n
            Y_fcn(i,j) = eval(s{i,j});
        end
    end
    matlabFunction(Y_fcn,'File','Y_fcn1','Vars',{x});
    Y_fcn = matlabFunction(Y_fcn,'Vars',{x});
    
    controller.mu = mu_opt;    
    controller.Y_fcn = Y_fcn;
    controller.lambda = lambda_opt;
    controller.alpha = alpha_opt;
end

% ------------------------- extract W_fcn & dW_fcn ------------------------
W_coef = value(W_coef);  
W_coef(abs(W_coef)<=1e-10) = 0;
x = x_store; % must ensure that v_W and s contain "x" instead of "x_store"
W_fcn = zeros(n);
for i=1:n_monos_W
    W_fcn = W_fcn+ W_coef(:,:,i)*v_W(i);
end
W_fcn = clean(W_fcn, 1e-10);
dv_W_dx = clean(dv_W_dx, 1e-10);
s = sdisplay(W_fcn);
s2 = sdisplay(dv_W_dx);
syms x [n 1]
syms W_fcn [n n]
for i=1:n
    for j=1:n
        W_fcn(i,j) = eval(s{i,j});
    end
end
%     W_fcn = arrayfun(@(si) eval(si), s);
matlabFunction(W_fcn,'File','W_fcn1','Vars',{x});
W_fcn = matlabFunction(W_fcn,'Vars',{x});
% W_fcn = @(x) W_fcn1(x(3),x(4));


[n1,n2]= size(dv_W_dx);
syms dv_W_dx_sym [n1 n2]

for i=1:n1
    for j=1:n2
        dv_W_dx_sym(i,j) = eval(s2{i,j});
    end
end    
dW_dphi = zeros(n);
dW_dvx = zeros(n);
for i=1:n_monos_W
    dW_dphi = dW_dphi + W_coef(:,:,i)*dv_W_dx_sym(i,1); 
    dW_dvx = dW_dvx + W_coef(:,:,i)*dv_W_dx_sym(i,2);
end

matlabFunction(dW_dphi,'File','dW_dphi','Vars',{x});
matlabFunction(dW_dvx,'File','dW_dvx','Vars',{x});
dW_dphi = matlabFunction(dW_dphi,'Vars',{x});
dW_dvx = matlabFunction(dW_dvx,'Vars',{x});
dW_dxi_fcn = @(i,x) (i==3)*dW_dphi(x)+(i==4)*dW_dvx(x);

if controller.type == CtrlDesignOpts.ccm
    dW_dt_fcn = @(x) dW_dphi(x)*f_phi_fcn(x) + dW_dvx(x)*f_vx_fcn(x); 
elseif controller.type == CtrlDesignOpts.rccm
    dW_dt_fcn = @(x,w) dW_dphi(x)*(f_phi_fcn(x)+Bw_phi_fcn(x)*w) + ...
                       dW_dvx(x)*(f_vx_fcn(x)+Bw_vx_fcn(x)*w); 
    %%%%%%%%%%%for testing with (Chebyshev poly.) approximated fcns. %%%%%%
    dW_dt_approx_fcn = @(x,w) dW_dphi(x)*(f_phi_fcn(x)+Bw_phi_fcn(x)*w) + ...
                       dW_dvx(x)*(f_vx_approx_fcn(x)+Bw_vx_approx_fcn(x)*w);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
controller.W_fcn = W_fcn;
controller.dW_dxi_fcn = dW_dxi_fcn;
controller.dW_dt_fcn = dW_dt_fcn;

%% check CCM conditions, compute tubes
compute_tubes;

%% save data
if save_rsts == 1
    if controller.type == CtrlDesignOpts.ccm
        file_name = ['ccm_' num2str(lambda0) '_plim_' num2str(p_lim/pi,2) 'pi.mat'];
    elseif controller.type == CtrlDesignOpts.rccm
        if  size(plant.C,1)>=6
            file_name = ['rccm_' num2str(lambda0,3) '_wmax_' num2str(w_lim,1) '_plim_' num2str(p_lim/pi,2) 'pi.mat'];
        else
            file_name = ['rccm_' num2str(lambda0,3) '_wmax_' num2str(w_lim,1) '_plim_' num2str(p_lim/pi,2) 'pi_pos.mat'];
        end
    end
    save(file_name,'plant','controller','state_set');
end
%% generate the c codes for accelerating geodesic computation
% To use the generated codes, copy .mex and .mat files to the sim folder

% parameters used in the pseudospectral method for geodesic computation
geodesic_setting_for_codegen.D = 2; 
geodesic_setting_for_codegen.N = 8;
answer = questdlg('Do you want to generate the C codes for accelerating geodesic computation used for determining the control law?','Question for code generation','Yes','No','No');
switch answer 
    case 'Yes'
        generate_code_for_geodesic_cal(plant.n,plant.nu,plant.nw,geodesic_setting_for_codegen);        
        save('geodesic_setting_for_codegen.mat','geodesic_setting_for_codegen');
    case 'No'
end
