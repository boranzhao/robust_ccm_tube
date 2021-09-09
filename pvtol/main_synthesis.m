%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Synthesis of (robust) CCM for control of a planar quadrotor

%  Author: Pan Zhao, UIUC, Advanced Controls Research Lab,
%  panzhao2@illinois.edu
%  Codes for the paper:
%  P. Zhao, et al. Tube-certified trajectory tracking for nonliner systems
%  with robust control contraction metrics. Submitted to IEEE Robotics and
%  Automation Letters, 2021. 
%  Last update: Sep 9, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
%% 
yalmip('clear');
controller.type = CtrlDesignOpts.rccm;      % lqr,ccm, rccm minimizing the PPG
controller.ccm_law_form = ...
    CtrlDesignOpts.ccm_min_norm;            % ccm_min_norm, ccm_integration
controller.ccm_mat_ineq_form = ...
    CtrlDesignOpts.ccm_mat_ineq_use_B_perp; % For ccm approach: ccm_mat_ineq_use_rho, or ccm_mat_ineq_use_B_perp . 
controller.opt_pos_dev = 1;                 % whether to focus on mitigate the effect of disturbance on the position states
lambda = 1.2;                               % a line search or bisection search should be performed for lambda
consider_state_set = 1;                     % whether to consider a compact set for the states when formulating the constraints
W_lower_bound = 1e-2;
biSearch_lambda = 0;                        % 1 for conducting a bisection search for lambda in designing the robust CCM controller 
maxNum_bs = 10;                             % maximum number of trials for bisection-search 
lambdaSearchRange = [0.1 10^1];             % search range for lambda
save_rsts = 1;                              % whether to save the results to a file

%----------------- settings for searching CCM & RCCM ----------------------
n = 6;nw = 1; nu =2;
x = sdpvar(n,1); x_store = x;
W_states_index = [3 4];
if controller.opt_pos_dev == 0
      C= [eye(n); zeros(nu,n)]; D =  [zeros(n,nu); 1*eye(nu)];  % consider both states and inputs
elseif controller.opt_pos_dev == 1
    C = [eye(2) zeros(2,4);zeros(nu,n)]; 
    D = [zeros(2,nu); 1*eye(nu)];                               % consider only position states and inputs
end
% ------------------ Stateconstraints for metric synthesis ----------------
p_lim = pi/3;  % phi
pd_lim = pi/3; % phi_dot 
vx_lim = 2;    % vx
vz_lim = 1;    % vz
w_lim = 1;     % w (disturbance)
state_set.box_lim = [p_lim^2-x(3)^2; vx_lim^2-x(4)^2; pd_lim^2-x(6)^2;  vz_lim^2-x(5)^2]*0.001;
state_set.num_consts_4_W_states = 2;        % # constraints from box_lim that involve states on which the metric W depends
state_set.other_lim_states = [x(6);x(5)]; 
state_set.lagrange_deg_W = 4;               % for the bounds of W
state_set.lagrange_deg_ccm = 4;             % for ccm condition
state_set.p_lim = p_lim;
state_set.pd_lim = pd_lim;
state_set.vx_lim = vx_lim;
state_set.vz_lim = vz_lim;
state_set.w_lim = w_lim;
% -------------------------------------------------------------------------

% ------------------ load system parameters -------------------------------
load_system_parameters; 

% approximating sin_x/cos_x with Chebshev polynomials
sinx = @(x) 0.9101*(x./(pi/3)) - 0.04466*(4*(x./(pi/3)).^3 - 3*(x./(pi/3))); % 0.8799 (close to 0.9101*3/pi, -0.03915
cosx = @(x) 0.7441 -0.2499*(2*(x./(pi/3)).^2 -1);                            % 0.7652, -0.2299

% ------------------------- system dynamics -------------------------------
w = sdpvar(nw,1); % disturbance
sin_p = sinx(x(3)); cos_p = cosx(x(3));
f = [x(4)*cos_p - x(5)*sin_p;    %px
    x(4)*sin_p + x(5)*cos_p;     %pz
    x(6);                        %phi
    x(6)*x(5)-plant.g*sin_p;         %vx
    -x(6)*x(4)-plant.g*cos_p;        %vz
    0];                          %phi_dot
% f written as a function handle: can also work when x has multiple columns
f_fcn = @(x) [x(4,:).*cos(x(3,:)) - x(5,:).*sin(x(3,:));    %px
            x(4,:).*sin(x(3,:)) + x(5,:).*cos(x(3,:));     %pz
            x(6,:);                        %phi
            x(6,:).*x(5,:)-plant.g*sin(x(3,:));         %vx
            -x(6,:).*x(4,:)-plant.g*cos(x(3,:));        %vz
            zeros(1,size(x,2))];          
        
B = [zeros(4,2); 1/plant.m 1/plant.m; plant.l/plant.J -plant.l/plant.J]; 
B_perp = [eye(4); zeros(2,4)];
Bw = [zeros(1,3),cosx(x(3)),-sinx(x(3)),0]'; 
Bw_fcn = @(x)[zeros(1,3),cos(x(3)),-sin(x(3)),0]';
df_dx = jacobian(f,x);
dBw_dx = jacobian(Bw,x);
A = df_dx + dBw_dx*w;
nz = size(C,1);

%%%%%%%%%% for testing with(Chebyshev polynomials) approximated fcns%%%%%%%
% using Chebyshev polynomial approximation.
Bw_approx_fcn = @(x)[zeros(1,3),cosx(x(3)),-sinx(x(3)),0]';
% approximated f_approx, df_dx_fcn
x = x_store;
s = sdisplay(f);
s2 = sdisplay(df_dx);
s3 = sdisplay(dBw_dx);
syms x [n 1]
syms f_approx_fcn [n 1]
syms df_dx_approx_fcn [n n]
syms dBw_dx_approx_fcn [n n]
for i=1:n    
    f_approx_fcn(i,1) = eval(s{i});    
    for j=1:n
        df_dx_approx_fcn(i,j) = eval(s2{i,j});
        dBw_dx_approx_fcn(i,j) =eval(s3{i,j}); 
    end
end
f_approx_fcn = matlabFunction(f_approx_fcn,'Vars',{x});
df_dx_approx_fcn = matlabFunction(df_dx_approx_fcn,'Vars',{x});
dBw_dx_approx_fcn = matlabFunction(dBw_dx_approx_fcn,'Vars',{x});
x = x_store;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plant.sinx = sinx;
plant.cosx = cosx;
plant.df_dx = df_dx;
plant.f_fcn = f_fcn;
plant.A = A;
plant.B = B;
plant.B_fcn = @(x) B;
plant.dynamics = @(x,u) f_fcn(x)+ B*u;
plant.w = w;

plant.B_perp = B_perp;
plant.Bw = Bw;
plant.Bw_fcn = Bw_fcn;
plant.C = C;
plant.D = D;
plant.n = n; plant.nu=nu; plant.nw = nw; plant.nz = nz; 

W_states = x(W_states_index);
f_phi_fcn = @(x) x(6);
f_vx_fcn = @(x) x(6)*x(5)-plant.g*sin(x(3));
Bw_phi_fcn = @(x) 0;
Bw_vx_fcn = @(x) cos(x(3));

%%%%%%%% for testing using (Chebyshev polynomials) approximated fcns%%%%%%%
f_vx_approx_fcn = @(x) x(6)*x(5)-plant.g*sinx(x(3));
Bw_vx_approx_fcn = @(x) cosx(x(3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_W = monolist(W_states,4);         % monomials of phi and vx up to degree 4
dv_W_dx = jacobian(v_W,W_states);
n_monos_W = length(v_W);
W_coef = sdpvar(n,n,n_monos_W);
W = zeros(n);
for i=1:n_monos_W
    W = W+ W_coef(:,:,i)*v_W(i);
end

if controller.type == CtrlDesignOpts.ccm
    dv_W_dt = dv_W_dx*f(W_states_index);
elseif controller.type == CtrlDesignOpts.rccm
    % dv_W_dt could depend on w
    dv_W_dt = dv_W_dx*(f(W_states_index)+Bw(W_states_index)*w);  
end

dW_dt = zeros(n);
for i=1:n_monos_W
    dW_dt = dW_dt+ W_coef(:,:,i)*dv_W_dt(i); 
end
[rho,c_rho,v_rho] = polynomial(W_states,2);
controller.rho = rho;
controller.c_rho = c_rho; 
controller.v_rho = v_rho;
n_monos_rho = length(v_rho);

controller.W_lower_bound = W_lower_bound;
state_set.consider_state_set = consider_state_set;
state_set.W_states = W_states;
state_set.W_states_index = W_states_index;

if controller.type == CtrlDesignOpts.rccm
   state_set.box_lim = [state_set.box_lim; w_lim^2-w^2]; 
end

%% Design CCM/RCCM controllers
if controller.type == CtrlDesignOpts.ccm  % search for CCM
    % YALMIP (SOS) and Mosek are used to solve the SOS problem
    % W0 = sdpvar(n,n);
    controller.lambda = lambda;%     
%     residual = res
    paras_W = W_coef(:);    
    [cond_num_W,w_upper,w_lower,W_bar,max_res] = ccm(plant,controller,W,dW_dt,paras_W,lambda,state_set);
        

    % ---------using a constant W matrix for testing or debugging----------
    % ---------e.g., whether the optimized geodesic is a straight line.---- 
    % W_fcn = @(x) Pinv; 
    % dW_fcn = {@(x) zeros(n), @(x) zeros(n), @(x) zeros(n)};
    % ---------------------------------------------------------------------
    
    controller.w_upper = w_upper;
    controller.w_lower = w_lower;    
    controller.W_bar = W_bar;
elseif controller.type == CtrlDesignOpts.rccm % search for robust CCM minimizing PPG
    controller.lambda = lambda;       
    Y_coef = sdpvar(nu,n,n_monos_W);
    Y = zeros(nu,n);
    for i=1:n_monos_W
        Y = Y+ Y_coef(:,:,i)*v_W(i);
    end 
    vars_WY = [W_coef(:); Y_coef(:)];
    if biSearch_lambda    
        % implement a bisection search for the optimal lambda to minimize gamma 
        Result_bs = ones(maxNum_bs+2,2)*inf; % 1st column for mu, second column for Gam
        Result_bs([1 2],1) = lambdaSearchRange';
        for i = 1:maxNum_bs+2
            if i <=2
                lambda = Result_bs(i,1);           
            else
                % sort the mu value in ascending order 
                [~,index] = sort(Result_bs(:,1));
                Result_bs = Result_bs(index,:);
                % sort the gamma value in ascending order 
                [~,index] = sort(Result_bs(:,2));
                % find the mu value corresponding to the smallest two gamma values and take their average as the mu value for next trial
                lambda = sum(Result_bs(index(1:2),1))/2;
                Result_bs(i,1) = lambda;
            end
            [alpha_opt,mu_opt,max_residual] = rccm_ppg(plant,controller,W,dW_dt,Y,vars_WY,lambda,x1,state_set); 
            Result_bs(i,2) = alpha_opt;
        end
        Result_bs
        [alpha_opt, index] = min(Result_bs(:,2));
        lambda_opt = Result_bs(index,1);
    else
       [alpha_opt,mu_opt] = rccm_ppg(plant,controller,W,dW_dt,Y,vars_WY,lambda,state_set); 
       lambda_opt = lambda;
    end
%     fprintf('The smallest gamma is %.3f obtained at lambda= %.3e\n',alpha_opt,lambda_opt);       
    fprintf('RCCM, lambda = %.2f, tube gain = %.4f\n',lambda_opt,alpha_opt);
    tube_gain = alpha_opt; % may include states and inputs

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
%     dW_dx_fcn = (i==3)*dW_fcn1 + (i==4)*dW_fcn2;
%     s = sdisplay(dW_dx_fcn);
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




%% --------------check CCM conditions, compute tubes on--------------------
%--------------- states and feedback control effort and save data----------
lambda0 = lambda;
if isfield(plant,'df_dx')
    plant= rmfield(plant,{'df_dx','Bw','w'});
end
if isfield(controller,'rho')
    controller= rmfield(controller,{'rho','c_rho','v_rho'});
end
if isfield(state_set,'box_lim')
    state_set = rmfield(state_set,{'box_lim','other_lim_states','W_states'});
end
plant.state_set = state_set;

disp('Checking CCM conditions ...');

lambda = 0.998*lambda; % slightly reduce lambda to avoid infeasible problems due to numerical issues    
ctrl_N = 10;
p_range = linspace(-p_lim, p_lim, ctrl_N);
vx_range = linspace(-vx_lim, vx_lim, ctrl_N);
vz_range = linspace(-vz_lim, vz_lim, 8);
pd_range = linspace(-pd_lim, pd_lim, 8);

df_dx_fcn = @(x) [0,0,-x(4)*sin(x(3))-x(5)*cos(x(3)),cos(x(3)),-sin(x(3)),0; 
               0,0, x(4)*cos(x(3))-x(5)*sin(x(3)),sin(x(3)), cos(x(3)),0;
               zeros(1,5),1;
               0,0,-plant.g*cos(x(3)),0,x(6),x(5);
               0,0, plant.g*sin(x(3)),-x(6),0,-x(4);
               zeros(1,6)];
           
dBw_dx_fcn = @(x) [zeros(3,6); zeros(2) [-sin(x(3)); -cos(x(3))] zeros(2,3);
                     zeros(1,6)];
A_fcn = @(x,w)  df_dx_fcn(x) + dBw_dx_fcn(x)*w;   
A_approx_fcn = @(x,w)  df_dx_approx_fcn(x) + dBw_dx_approx_fcn(x)*w;   
controller.df_dx_fcn = df_dx_fcn;
controller.A_fcn = A_fcn;
controller.A_approx_fcn = A_approx_fcn;
if controller.type == CtrlDesignOpts.ccm  
    % --------------- following the approach in the following paper -------
    % S.  Singh, et al. Robust  feedback  motion  planning  via
    % contraction  theory. IJRR, 2019.
    delta_u = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);
    eig_CCM = zeros(ctrl_N, ctrl_N, ctrl_N, ctrl_N);
    eig_W = zeros(ctrl_N,ctrl_N,2);
    sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);
    sigma_Bw = zeros(ctrl_N,1); % Note that Bw depends on only x(3)

    for i = 1:length(p_range)
        for j = 1:length(vx_range)
            x = [randn(2,1);p_range(i);vx_range(j);0;0];                
            W = W_fcn(x); % note that W depends only on phi and vx
            M = W\eye(n);
            eig_W(i,j,1) = min(eig(W));
            eig_W(i,j,2) = max(eig(W));
            for k = 1:length(vz_range)
                for l = 1:length(pd_range)
                    x = [randn(2,1);p_range(i);vx_range(j);vz_range(k);pd_range(l)];                
                    Theta = chol(M);
                    Theta_Bw = Theta*Bw_fcn(x);
                    sigma_ThBw(i,j,k,l) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));

                    L = chol(W);                
%                     f = f_part_fcn(x);

                    df_dx0 = df_dx_fcn(x);                
                    F = dW_dt_fcn(x) - df_dx0*W - (df_dx0*W)'-2*lambda*W;

                    L_inv = inv(L);                
                    delta_u_den = eig((L_inv)'*(B*B')*L_inv);
                    delta_u(i,j,k,l) = 0.5*max(eig((L_inv)'*F*L_inv))/...
                        sqrt(min(delta_u_den(delta_u_den>0)));

                    R_CCM = B_perp'*F*B_perp;
                    eig_CCM(i,j,k,l) = min(eig(R_CCM));
                end
            end
        end
        Bw0 = Bw_fcn(x);
        sigma_Bw(i) = max(sqrt(eig(Bw0'*Bw0)));
    end
    % compute the tube size 
    Wbar_inv = eye(n)/W_bar;
    controller.tube_gain = 1/lambda*sqrt(value(cond_num_W))*max(sigma_Bw); 
    alpha_w = max(sigma_ThBw(:));
    d_bar_gain = alpha_w/lambda; % just a gain: needs to consider the bound on w to determine the final bounds 
    %     disp('Tube gain'); disp(d_bar);
    controller.d_bar = d_bar_gain*w_lim;
    controller.tube_gain_u = d_bar_gain*max(delta_u(:));
    controller.u_bnd = controller.tube_gain_u*w_lim;
    Wbar_inv_eigs = eig(Wbar_inv(1:2,1:2));
    controller.tube_gain_states = controller.tube_gain;
    controller.tube_gain_xz = d_bar_gain/sqrt(min(Wbar_inv_eigs));
    controller.tube_gain_x = d_bar_gain/sqrt(Wbar_inv(1,1));
    controller.tube_gain_z = d_bar_gain/sqrt(Wbar_inv(2,2));
    controller.tube_gain_phi = d_bar_gain/sqrt(Wbar_inv(3,3));
    controller.tube_gain_phidot = d_bar_gain/sqrt(Wbar_inv(6,6));
    controller.tube_gain_vx = d_bar_gain/sqrt(Wbar_inv(4,4));
    controller.tube_gain_vz = d_bar_gain/sqrt(Wbar_inv(5,5));
    
    controller.d_bar_gain = d_bar_gain;
    fprintf(1,'CCM, lambda = %.2f, cond(W) = %.3f,  tube gain (xz) = %.3f\n',lambda,cond_num_W,controller.tube_gain_xz);
    fprintf(1,'Control: %.3f\n',controller.u_bnd);
    fprintf(1,'min and max eigenvalues of W: %.3e, %.3e\n', min(min(eig_W(:,:,1))),max(max(eig_W(:,:,2))));
    fprintf(1,'minimum eigenvalue of CCM matrix (should be positive): %.3e\n',min(eig_CCM(:)));
    
    if controller.ccm_mat_ineq_form == CtrlDesignOpts.ccm_mat_ineq_use_rho       
        % extract rho_fcn
        c_rho = value(c_rho);     
        x = x_store; % must ensure that v_W and s contain "x" instead of "x_store"
        rho_fcn = 0;
        for i=1:n_monos_rho
            rho_fcn = rho_fcn+ c_rho(i)*v_rho(i);
        end
        rho_fcn = clean(rho_fcn, 1e-8);
        s = sdisplay(rho_fcn);
        syms x [n 1]
        syms rho_fcn
        rho_fcn = eval(s{1});
        rho_fcn = matlabFunction(rho_fcn,'Vars',{x});
        controller.rho_fcn = rho_fcn;
    end
    if save_rsts == 1
        file_name = ['ccm_' num2str(lambda0) '_plim_' num2str(p_lim/pi,2) 'pi.mat'];
        save(file_name,'plant','controller','state_set');
    end
elseif controller.type == CtrlDesignOpts.rccm
     % --- check matrix inequality conditions 
    w_range = linspace(-w_lim, w_lim, ctrl_N);
    eig_CCM = zeros(ctrl_N, ctrl_N, ctrl_N, ctrl_N,ctrl_N);
    eig_W = zeros(ctrl_N,ctrl_N,3);
    sigma_Bw = zeros(ctrl_N,1);         % Note that Bw depends on only x(3)

    for i = 1:length(p_range)
        for j = 1:length(vx_range)
            x = [0;0;p_range(i);vx_range(j);0;0];                
            W = W_fcn(x);
            M = W\eye(n);
            Y = Y_fcn(x);
            eig_W(i,j,1) = min(eig(W));
            eig_W(i,j,2) = max(eig(W));

            F2 = [lambda*W  zeros(n,nw) (C*W+D*Y)';
                  zeros(nw,n)   (alpha_opt-mu_opt)*eye(nw)  zeros(nw,nz);
                  C*W+D*Y    zeros(nz,nw)    alpha_opt*eye(nz)];
            eig_W(i,j,3) = min(eig(F2));
            
            for k = 1:length(vz_range)
                for l = 1:length(pd_range)
                    x = [randn(2,1);p_range(i);vx_range(j);vz_range(k);pd_range(l)];                
                  
                    for m = 1:length(w_range)
                        w = w_range(m);
%                         tmp = A_fcn(x,w)*W + B*Y;   
%                         F1 = [dW_dt_fcn(x,w)-(tmp+tmp')-lambda*W  -Bw_fcn(x); -Bw_fcn(x)' mu_opt*eye(nw)];
%%%%%%%%%%%%%%% for testing with (Cheby. poly) approximated functions %%%%%
                        tmp = A_approx_fcn(x,w)*W + B*Y;
                        F1 = [dW_dt_approx_fcn(x,w)-(tmp+tmp')-lambda*W  -Bw_approx_fcn(x); -Bw_approx_fcn(x)' mu_opt*eye(nw)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        eig_CCM(i,j,k,l,m) = min(eig(F1));
                    end
                end
            end
        end
        Bw0 = Bw_fcn(x);
        sigma_Bw(i) = max(sqrt(eig(Bw0'*Bw0)));
    end 
    plant1 = plant;
    % get the bound on x
    if size(plant.C,1)~=n  % 
        plant1.C = eye(n); plant1.D = zeros(n,nu);  % consider all the states
        [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_grid(plant1,controller,state_set);  
        controller.tube_gain_states = alpha_opt;
    end
        
    % get the bound on px, pz
    plant1.C = [eye(2) zeros(2,4)]; plant1.D = zeros(2,nu);    % consider only position states
    [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_grid(plant1,controller,state_set);  
    controller.tube_gain_xz = alpha_opt;
    controller.tube_xz = controller.tube_gain_xz*w_lim;
    
    % get the bound on vx, vz
    plant1.C = [0 0 0 1 0 0; 0 0 0 0 1 0]; 
    plant1.D = zeros(2,nu);    % consider only position states
    [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_grid(plant1,controller,state_set);  
    controller.tube_gain_vxvz = alpha_opt;
    controller.tube_vxvz = controller.tube_gain_vxvz*w_lim;
    
    % get the bound on phi, phi_dot 
    plant1.C = [0 0 1 0 0 0]; plant1.D = zeros(1,nu);   
    [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_grid(plant1,controller,state_set);  
    controller.tube_gain_phi = alpha_opt;
    controller.tube_phi = controller.tube_gain_phi*w_lim;
    
    plant1.C = [0 0 0 0 0 1]; plant1.D = zeros(1,nu);    
    [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_grid(plant1,controller,state_set);  
    controller.tube_gain_phidot = alpha_opt;
    controller.tube_phidot = controller.tube_gain_phidot*w_lim;
        
    % get the bound on u.
    plant1.C = zeros(nu,n); plant1.D = eye(nu);
    [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_grid(plant1,controller,state_set);
    controller.tube_gain_u = alpha_opt;
    controller.u_bnd = controller.tube_gain_u*w_lim;
    
    fprintf('RCCM, lambda = %.2f, tube gain (xz) = %.4f\n',lambda_opt,controller.tube_gain_xz);
    fprintf(1,'Control: %.3f\n', controller.u_bnd);
    fprintf(1,'min and max eigenvalues of W: %.3e, %.3e\n', min(min(eig_W(:,:,1))),max(max(eig_W(:,:,2))));
    fprintf(1,'minimum eigenvalue of lhs of PPG condition 1 (should be non-negative): %.3e\n',min(eig_CCM(:)));
    fprintf(1,'mininum eigenvalue of lhs of PPG condition 2 (should be non-negative): %.3e\n', min(min(eig_W(:,:,3))));
    if save_rsts == 1
        if  size(plant.C,1)>=6
            file_name = ['rccm_' num2str(lambda0,3) '_wmax_' num2str(w_lim,1) '_plim_' num2str(p_lim/pi,2) 'pi.mat'];
        else
            file_name = ['rccm_' num2str(lambda0,3) '_wmax_' num2str(w_lim,1) '_plim_' num2str(p_lim/pi,2) 'pi_pos.mat'];
        end
        save(file_name,'plant','controller','state_set');
    end
end