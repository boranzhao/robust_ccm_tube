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

lambda = 0.998*lambda; % slightly reduce lambda to avoid constraint violations due to numerical issues    
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
    % --------------- following the paper below -------
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
    d_bar_gain = alpha_w/lambda; % just a gain: need to consider the bound on w to determine the final bounds 
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
                        % use (Cheby. poly) approximated functions
                        tmp = A_approx_fcn(x,w)*W + B*Y;
                        F1 = [dW_dt_approx_fcn(x,w)-(tmp+tmp')-lambda*W  -Bw_approx_fcn(x); -Bw_approx_fcn(x)' mu_opt*eye(nw)];

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
end