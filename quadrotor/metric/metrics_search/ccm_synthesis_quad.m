function [cond_num_W_opt,W_bar,w_lower,w_upper,W_coef,w_poly] = ccm_synthesis_quad(prog,plant,controller,W_list,W_eval,w_poly,w_poly_fcn,dw_poly_dt_fcn,lambda,state_set)
% The following codes are based on the codes available at https://github.com/stanfordASL/RobustMP
norm_scale = 1e-7;
W_scale = (1e-5)*diag([12;12;24;0.1;0.1;0.2;5;3.5;3.5]);
ccm_eps =  0.01;
A_fcn = plant.A_fcn;
n = plant.n;nu = plant.nu;

[prog, w_lower] = prog.newPos(1);
[prog, w_upper] = prog.newPos(1);
[prog, W_bar] = prog.newSym(n);


%W uniform bounds
prog = prog.withPos(w_lower-controller.w_lower_bound);
prog = prog.withPSD(w_upper*eye(n)-W_bar);

ctrl_N = 8;

T_range = linspace(state_set.T_lim(1), state_set.T_lim(2), 4);
phi_range = linspace(-state_set.phi_lim, state_set.phi_lim, ctrl_N);
theta_range = linspace(-state_set.theta_lim, state_set.theta_lim, ctrl_N);

% Note that LMI1 is affine w.r.t to u; therefore, evaluation at the
% vertices is enough
Tdot_range = [state_set.Tdot_lim, -state_set.Tdot_lim];
phidot_range = [state_set.phidot_lim, -state_set.phidot_lim];
thetadot_range = [state_set.thetadot_lim, -state_set.thetadot_lim];

for i = 1:length(T_range)
    for j = 1:length(phi_range)    
        for k = 1:length(theta_range)
            x0 = [randn(6,1);T_range(i);phi_range(j);theta_range(k)];    
            w_poly0 = w_poly_fcn(x0);
            W = W_list{1}*w_poly0(1);
            for s=2:length(w_poly0)
                W = W+ W_list{s}*w_poly0(s);
            end
            %W pos def  & W upper bound
            prog = prog.withPSD(W - w_lower*eye(n));     
            prog = prog.withPSD(W_bar - W);                   
             
            
            %CCM condition
            tmp = A_fcn(x0)*W;
            % PSD constraint: equation (11)
            if controller.ccm_mat_ineq_form == CtrlDesignOpts.ccm_mat_ineq_use_B_perp
                dw_poly_dt0 = dw_poly_dt_fcn([x0;0;0;0]);  
                dW_dt = W_eval(dw_poly_dt0); % note that the dependence u will disappear after multiplication by B_perp' and its transpose if the top left block does not depends on any states
                CCM_pos = plant.B_perp'*(dW_dt-(tmp+tmp')-2*lambda*W)*plant.B_perp;
                prog = prog.withPSD(CCM_pos - ccm_eps*eye(n-nu));                
            elseif controller.ccm_mat_ineq_form == CtrlDesignOpts.ccm_mat_ineq_use_rho 
                for p = 1:length(Tdot_range)
                    for q = 1:length(phidot_range)
                        for r = 1:length(thetadot_range)
                            u0 = [Tdot_range(p);phidot_range(q);thetadot_range(r)];
                            dw_poly_dt0 = dw_poly_dt_fcn([x0;u0]);  
                            dW_dt = W_eval(dw_poly_dt0);                        
                            CCM_pos = dW_dt - (tmp+tmp')+ controller.rho*(plant.B*plant.B')-2*lambda*W;              
                            prog = prog.withPSD(CCM_pos - ccm_eps*eye(n));
                        end
                    end
                end
            end  
            
            % using for loop
%             tic;
%             W = W_list{1}*w_poly0(1);
%             dW_dt =  W_list{1}*dw_poly_dt0(1); 
%             Y = Y_list{1}*y_poly0(1);
%             for l=2:length(w_poly0)
%                 W = W+ W_list{l}*w_poly0(l);
% %                 dW_dt = dW_dt + W_list{l}*dw_poly_dt0(l);
% %                 Y = Y + Y_list{l}*y_poly0(l);
%             end  
%             toc;
            
%             % ---------- using pre-generated functions ---------------
%             tic;
%             W2 = W_eval(w_poly0);
%             toc;
%             Y2 = Y_eval(y_poly0);
%             dW_dt2 = W_eval(dw_poly_dt0);
%             toc;
%             % ---------------------------------------------------------
        end
    end
end

%Norm constraint
free_vars = [prog.coneVar(2:end); prog.freeVar];
len = length(free_vars);
[prog, a] = prog.newPos(len);
prog = prog.withPos(-free_vars + a);
prog = prog.withPos(free_vars + a);

options = spot_sdp_mosek_options();
% options.solver_options.mosek.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
% options.solver_options.mosek.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)
options.verbose = 1;
        
disp('LMI formulation finished! Solving...');
SOS_soln = prog.minimize(norm_scale*sum(a) + w_upper*1e-3, @spot_mosek, options);

try
    solved = (strcmp(SOS_soln.info.solverInfo.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE'));% && ...
%            strcmp(SOS_soln.info.solverInfo.itr.solsta, 'OPTIMAL'));
catch
    solved = 0;
end

%% parse
if (solved == 0)
    disp('There are some issues...');
    w_lower = nan;
    w_upper = nan;
    W_bar = nan(n);
    cond_num_W_opt = nan;
    W_coef = nan;
    return;   
else
    disp('feasible, getting results...');
    w_lower = double(SOS_soln.eval(w_lower));
    w_upper = double(SOS_soln.eval(w_upper));
    cond_num_W_opt = w_upper/w_lower;
    W_bar = clean(double(SOS_soln.eval(W_bar)),1e-4);

    
    len_w_poly = length(w_poly);
    W_coef = zeros(n,n,len_w_poly);
    NNZ_list = zeros(len_w_poly,1);
    for i = 1:len_w_poly
        W_coef(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),1e-7);
        if sum(sum(abs(W_coef(:,:,i)))) > 0
            NNZ_list(i) = 1;
        end
    end
%     
    w_poly = w_poly(NNZ_list==1);
    W_coef = W_coef(:,:,NNZ_list==1);

    fprintf('%d non-zero monomials\n',length(w_poly));

%     %save in coefficient form (for copying over to C++)
%     p = W_poly_mat(find(NNZ_list),:);
%     save('Quad_Metric_Coeffs.mat','W_sol','p');
end
end

