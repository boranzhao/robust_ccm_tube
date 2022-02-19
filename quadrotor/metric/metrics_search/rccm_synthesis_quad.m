
function [alpha_opt,mu_opt,W_coef,Y_coef,w_poly,y_poly] = rccm_synthesis_quad(prog,plant,controller,W_list,W_eval,Y_list,Y_eval,w_poly,w_poly_fcn,dw_poly_dt_fcn,y_poly,y_poly_fcn,lambda,state_set)
% x is the variables that appear in W and Y
norm_scale = 1e-7;
ccm_eps = 0.01;
n = plant.n; nu = plant.nu;
nw = plant.nw;
nz = plant.nz;
A_fcn = plant.df_dx_fcn; %df_dx_approx_fcn; 
B = plant.B;
Bw = plant.Bw;
C = plant.C;
D = plant.D;
  
[prog, mu] = prog.newPos(1);
[prog, alpha] = prog.newPos(1);
ctrl_N = 7;

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
            y_poly0 = y_poly_fcn(x0);
            W = W_list{1}*w_poly0(1);

            Y = Y_list{1}*y_poly0(1);
            for s=2:length(w_poly0)
                W = W+ W_list{s}*w_poly0(s);
%                 dW_dt = dW_dt + W_list{l}*dw_poly_dt0(l);
                Y = Y + Y_list{s}*y_poly0(s);
            end 
            tmp = A_fcn(x0)*W + B*Y;
            
            for p = 1:length(Tdot_range)
                for q = 1:length(phidot_range)
                    for r = 1:length(thetadot_range)
                        u0 = [Tdot_range(p);phidot_range(q);thetadot_range(r)];
                        dw_poly_dt0 = dw_poly_dt_fcn([x0;u0]);  
                        dW_dt = W_eval(dw_poly_dt0);                        
                        LMI1 = [dW_dt-(tmp+tmp')-lambda*W  -Bw; -Bw' mu*eye(nw)];                
                        prog = prog.withPSD(LMI1-ccm_eps*eye(n+nw));
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

            
            LMI2 = [lambda*W  zeros(n,nw) (C*W+D*Y)';
                  zeros(nw,n)   (alpha-mu)*eye(nw)  zeros(nw,nz);
                  C*W+D*Y    zeros(nz,nw)    alpha*eye(nz)];
            prog = prog.withPSD(LMI2);
              
            %W pos def
            prog = prog.withPSD(W- controller.w_lower_bound*eye(n));    
            % LMI1 for Linf gain
      % % ------------------------------------------------------------------------
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
SOS_soln = prog.minimize(norm_scale*sum(a) + alpha, @spot_mosek, options);

try
    solved = (strcmp(SOS_soln.info.solverInfo.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE'));% && ...
%            strcmp(SOS_soln.info.solverInfo.itr.solsta, 'OPTIMAL'));
catch
    solved = 0;
end

%% parse
if (solved == 0)
    disp('There are some issues...');
    alpha_opt = nan;
    mu_opt = nan;
    W_coef = nan;
    Y_coef = nan;
    return;   
else
    disp('feasible, getting results...');
    len_w_poly = length(w_poly);
    len_y_poly = length(y_poly);
    W_coef = zeros(n,n,len_w_poly);
    Y_coef = zeros(nu,n,len_y_poly);
    NZ_W_list = zeros(len_w_poly,1);
    NZ_Y_list = zeros(len_y_poly,1);
    for i = 1:len_w_poly
        W_coef(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),1e-8);
        if sum(sum(abs(W_coef(:,:,i)))) > 0
            NZ_W_list(i) = 1;
        end
    end
    for i = 1:len_y_poly
        Y_coef(:,:,i) = clean(double(SOS_soln.eval(Y_list{i})),1e-8);
        if sum(sum(abs(Y_coef(:,:,i)))) > 0
            NZ_Y_list(i) = 1;
        end
    end
    
    alpha_opt = double(SOS_soln.eval(alpha));    
    mu_opt = double(SOS_soln.eval(mu));  
    
    w_poly = w_poly(NZ_W_list==1);
    y_poly = y_poly(NZ_Y_list==1);
    W_coef = W_coef(:,:,NZ_W_list==1);
    Y_coef = Y_coef(:,:,NZ_Y_list==1);
    
    fprintf('%d (%d) non-zero monomials for W (Y)\n',length(w_poly),length(y_poly));
end
end