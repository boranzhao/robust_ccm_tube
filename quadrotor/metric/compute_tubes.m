lambda0 = lambda;
if isfield(plant,'df_dx')
    plant= rmfield(plant,{'df_dx','u','wdist','sinx','cosx','A','dBw_dx_fcn','A_fcn','df_dx_approx_fcn'});
end
if isfield(controller,'rho')
    controller= rmfield(controller,{'rho','c_rho','v_rho'});
end
if isfield(state_set,'box_lim')
    state_set = rmfield(state_set,{'box_lim','w_states'});
end
plant.state_set = state_set;

disp('Checking CCM conditions ...');

lambda = 0.998*lambda; % slightly reduce lambda to avoid infeasible problems due to numerical issues    
ctrl_N = 10;
T_range = linspace(T_lim(1), T_lim(2), ctrl_N);
phi_range = linspace(-phi_lim, phi_lim, ctrl_N);
theta_range = linspace(-theta_lim, theta_lim, ctrl_N);
controller.df_dx_fcn = df_dx_fcn;
controller.A_fcn = plant.df_dx_fcn;
if controller.type == CtrlDesignOpts.ccm  
    % --------------- following the approach in the following paper -------
    % S.  Singh, et al. Robust  feedback  motion  planning  via
    % contraction  theory. IJRR, 2019.
    delta_u = zeros(ctrl_N,ctrl_N,ctrl_N);
    eig_CCM = zeros(ctrl_N, ctrl_N, ctrl_N);
    eig_W = zeros(ctrl_N,ctrl_N,ctrl_N,2);
    sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N);
    sigma_Bw = norm(plant.Bw_fcn(zeros(n,1))); % Bw is constant for this example 
    for i = 1:length(T_range)
        for j = 1:length(phi_range)    
            for k = 1:length(theta_range)
                x0 = [randn(6,1);T_range(i);phi_range(j);theta_range(k)];    
                W = W_fcn(x0); % note that W depends only on T, phi and theta
                M = W\eye(n);
                eig_W(i,j,k,1) = min(eig(W));
                eig_W(i,j,k,2) = max(eig(W));
                Theta = chol(M);
                Theta_Bw = Theta*Bw_fcn(x0);
                sigma_ThBw(i,j,k) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));
   
                L = chol(W);               

                df_dx0 = plant.df_dx_fcn(x0);                
                F = dW_dt_fcn(x0,[0 0 0]') - df_dx0*W - (df_dx0*W)'-2*lambda*W; %%%%% there might be some issues here when the strong CCM condition does not hold

                L_inv = inv(L);                
%                 delta_u_den = eig((L_inv)'*(B*B')*L_inv);
%                 delta_u(i,j,k) = 0.5*max(eig((L_inv)'*F*L_inv))/...
%                     sqrt(min(delta_u_den(delta_u_den>0)));

                R_CCM = B_perp'*F*B_perp;
                eig_CCM(i,j,k) = min(eig(R_CCM));
            end
        end
    end
    % compute the tube size 
    Wbar_inv = eye(n)/W_bar;
    controller.tube_gain = 1/lambda*sqrt(value(cond_num_W))*sigma_Bw; 
    alpha_w = max(sigma_ThBw(:));
    d_bar_gain = alpha_w/lambda; % just a gain: needs to consider the bound on w to determine the final bounds 
    controller.d_bar = d_bar_gain*state_set.wdist_lim;
%     controller.tube_gain_u = d_bar_gain*max(delta_u(:));
%     controller.u_bnd = controller.tube_gain_u*state_set.wdist_lim;
    Wbar_inv_eigs = eig(Wbar_inv(1:3,1:3));
    controller.tube_gain_states = controller.tube_gain;
    controller.tube_gain_xyz = d_bar_gain/sqrt(min(Wbar_inv_eigs));
    controller.tube_gain_xy = d_bar_gain/sqrt(min(eig(Wbar_inv(1:2,1:2))));
    controller.tube_gain_z = d_bar_gain/sqrt(Wbar_inv(3,3));
    controller.tube_gain_vxyz =  d_bar_gain/sqrt(min(eig(Wbar_inv(4:6,4:6))));
    controller.tube_gain_vxy =  d_bar_gain/sqrt(min(eig(Wbar_inv(4:5,4:5))));
    controller.tube_gain_vz = d_bar_gain/sqrt(Wbar_inv(6,6));
    controller.tube_gain_T = d_bar_gain/sqrt(Wbar_inv(7,7));
    controller.tube_gain_phitheata = d_bar_gain/sqrt(min(eig(Wbar_inv(8:9,8:9))));
    controller.d_bar_gain = d_bar_gain;
    fprintf(1,'CCM, lambda = %.2f, cond(W) = %.3f,  tube gain (all states) = %.3f, tube gain (xyz) = %.3f, tube gain (xy) = %.3f\n',...
        lambda,cond_num_W,controller.tube_gain_states,controller.tube_gain_xyz,controller.tube_gain_xy);
%     fprintf(1,'Control: %.3f\n',controller.u_bnd);
    fprintf(1,'min and max eigenvalues of W: %.3e, %.3e\n',min(vec(eig_W(:,:,:,1))),max(vec(eig_W(:,:,:,2))));
    fprintf(1,'minimum eigenvalue of CCM matrix (should be positive): %.3e\n',min(eig_CCM(:)));
    
    if controller.ccm_mat_ineq_form == CtrlDesignOpts.ccm_mat_ineq_use_rho       
        % extract rho_fcn
        c_rho = value(c_rho);     
        x = x_store; % must ensure that W_poly and s contain "x" instead of "x_store"
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
elseif controller.type == CtrlDesignOpts.rccm
    Tdot_range = linspace(-Tdot_lim,Tdot_lim,3);
    phidot_range = linspace(-phidot_lim,phidot_lim,3);
    thetadot_range = linspace(-thetadot_lim, thetadot_lim,3);
     % --- check matrix inequality conditions 
    eig_CCM = zeros(ctrl_N, ctrl_N, ctrl_N,3,3,3);
    eig_W = zeros(ctrl_N,ctrl_N,ctrl_N,3);
    sigma_Bw = zeros(ctrl_N,1);         % Note that Bw depends on only x(3)
    for i = 1:length(T_range)
        for j = 1:length(phi_range)
            for k = 1:length(theta_range)
                x0 = [randn(6,1);T_range(i);phi_range(j);theta_range(k)];           
                W = W_fcn(x0);
                M = W\eye(n);
                Y = Y_fcn(x0);
                eig_W_tmp = eig(W);
                eig_W(i,j,k,1) = min(eig_W_tmp);
                eig_W(i,j,k,2) = max(eig_W_tmp);
                F2 = [lambda*W  zeros(n,nw) (C*W+D*Y)';
                      zeros(nw,n)   (alpha_opt-mu_opt)*eye(nw)  zeros(nw,nz);
                      C*W+D*Y    zeros(nz,nw)    alpha_opt*eye(nz)];
                eig_W(i,j,k,3) = min(eig(F2));
                
                tmp = df_dx_fcn(x0)*W + B*Y;
                for p = 1:length(Tdot_range)
                    for q = 1:length(phidot_range)             
                        for r = 1:length(thetadot_range)
                            u0 = [Tdot_range(p);phidot_range(q);thetadot_range(r)];
                            F1 = [dW_dt_fcn(x0,u0)-(tmp+tmp')-lambda*W  -Bw; -Bw' mu_opt*eye(nw)]; %%%%%%%%%%%%% there are something wrong with dW_dt_fcn
                            eig_CCM(i,j,k,p,q,r) = min(eig(F1));
                        end
                    end
                end
            end
        end
    end    
    fprintf(1,'min and max eigenvalues of W: %.3e, %.3e\n', min(vec(eig_W(:,:,:,1))),max(vec(eig_W(:,:,:,2))));
    fprintf(1,'minimum eigenvalue of lhs of PPG condition 1 (should be non-negative): %.3e\n',min(eig_CCM(:)));
    fprintf(1,'mininum eigenvalue of lhs of PPG condition 2 (should be non-negative): %.3e\n', min(vec(eig_W(:,:,:,3))));
    
    %% refine the tube sizes
    plant1 = plant;    
    if size(plant1.C,1)~=n  % 
        disp('get the bound on all states...');
        plant1.C = eye(n); plant1.D = zeros(n,nu);  % consider all the states
        [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_quad(plant1,controller,state_set);  
        controller.tube_gain_states = alpha_opt;
        fprintf('RCCM, lambda = %.2f, tube gain (all states) = %.4f\n',lambda_opt,controller.tube_gain_states);
    end
    
    disp('get the bound on px, py, pz...')
    plant1.C = [eye(3) zeros(3,n-3)]; plant1.D = zeros(3,nu);    % consider only position states
    
%     % ------------------------  spot: 10 times slower --------------------
%     tic;
%     [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_quad_spot(prog,plant1,controller,state_set);
%     toc;
%     % --------------------------------------------------------------------
    [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_quad(plant1,controller,state_set);
    
    controller.tube_gain_xyz = alpha_opt;
    fprintf('RCCM, lambda = %.2f, tube gain (xyz) = %.4f\n',lambda_opt,controller.tube_gain_xyz);
%     disp('get the bound on vx, vy, vz...');
%     plant1.C = zeros(3,9); plant1.C(1,4) = 1; plant1.C(2,5) = 1;  plant1.C(3,6) = 1; 
%     plant1.D = zeros(3,nu);   
%     [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_quad(prog,plant1,controller,state_set);
%     controller.tube_gain_vxyz = alpha_opt;
%     controller.tube_vxyz = controller.tube_gain_vxyz*state_set.wdist_lim;
%     
    disp('get the bound on T, phi, theta...');
    plant1.C = zeros(1,9); plant1.C(1,7) = 1;   
    plant1.D = zeros(1,nu);    %
    [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_quad(plant1,controller,state_set);
    controller.tube_gain_T = alpha_opt;
    
    plant1.C = zeros(2,9); plant1.C(1,8) = 1; plant1.C(2,9) = 1;    
    plant1.D = zeros(2,nu);    %
    [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_quad(plant1,controller,state_set);
    controller.tube_gain_rp = alpha_opt;
    fprintf(1,'tube gain (T): %.3f tube gain (rp): %.3f \n', controller.tube_gain_T,controller.tube_gain_rp);
        
    disp('get the bound on Tdot...');
    plant1.C = zeros(1,n); plant1.D = [1 0 0];
    [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_quad(plant1,controller,state_set);
    controller.tube_gain_Tdot = alpha_opt;
    
    disp('get the bound on phidot and thetadot...');
    plant1.C = zeros(2,n); plant1.D = [0 1 0;0 0 1];
    [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_quad(plant1,controller,state_set);
    controller.tube_gain_rpdot = alpha_opt;
    
    fprintf(1,'tube gain (Control): Tdot: %.3f rpdot: %.3f \n', controller.tube_gain_Tdot,controller.tube_gain_rpdot);
end
