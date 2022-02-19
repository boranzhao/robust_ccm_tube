lambda0 = lambda;
plant.state_set = state_set;
disp('Checking RCCM conditions ...');
lambda = 0.998*lambda; % slightly reduce lambda to avoid infeasible problems due to numerical issues    
ctrl_N = 10;
T_range = linspace(T_lim(1), T_lim(2), ctrl_N);
phi_range = linspace(-phi_lim, phi_lim, ctrl_N);
theta_range = linspace(-theta_lim, theta_lim, ctrl_N);

A_fcn = @(x,w)  plant.df_dx_fcn(x) + plant.dBw_dx_fcn(x)*w; 
controller.df_dx_fcn = df_dx_fcn;
controller.A_fcn = A_fcn;

% --------------- following the approach in the following paper -------
% S.  Singh, et al. Robust  feedback  motion  planning  via
% contraction  theory. IJRR, 2019.
eig_CCM = zeros(ctrl_N, ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,ctrl_N,ctrl_N,2);
sigma_Bw = norm(plant.Bw_fcn(zeros(n,1))); % Bw is constant for this example 
for i = 1:length(T_range)
    for j = 1:length(phi_range)    
        for k = 1:length(theta_range)
            x0 = [randn(6,1);T_range(i);phi_range(j);theta_range(k)];    
            W = W_fcn(x0); % note that W depends only on T, phi and theta
            M = W\eye(n);
            eig_W(i,j,k,1) = min(eig(W));
            eig_W(i,j,k,2) = max(eig(W));
            df_dx0 = df_dx_fcn(x0);   
           
            F = dW_dt_fcn(x0) - df_dx0*W - (df_dx0*W)'-lambda*W; %%%%% there might be some issues here when the strong CCM condition does not hold
            R_CCM = B_perp'*F*B_perp-1/mu_opt*B_perp'*Bw*Bw'*B_perp;
            eig_CCM(i,j,k) = min(eig(R_CCM));
        end
    end
end
% compute the tube size 
fprintf(1,'min and max eigenvalues of W: %.3e, %.3e\n',min(vec(eig_W(:,:,:,1))),max(vec(eig_W(:,:,:,2))));
fprintf(1,'minimum eigenvalue of RCCM matrix (should be positive): %.3e\n',min(eig_CCM(:)));
lambda = lambda0;

