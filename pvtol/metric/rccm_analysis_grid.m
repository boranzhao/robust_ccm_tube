
function [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_grid(plant,controller,state_set)
%% compute a bound on the RCCM-based feedback control effort given the plant and W and Y. 
% x is the variables that appear in W and Y
tolerance = 1e-6;
n = plant.n; nw = plant.nw;  nu = plant.nu;
B = plant.B;
Bw_fcn = plant.Bw_fcn;
% C = zeros(nu,n); D = eye(nu); %To set z = u, we need
C = plant.C; D = plant.D;

nz = size(C,1);
ctrl_N = 10;
p_range = linspace(-state_set.p_lim, state_set.p_lim, ctrl_N);
vx_range = linspace(-state_set.vx_lim, state_set.vx_lim, ctrl_N);

% those matrix inequalities are bi-linear w.r.t vz and pd
vz_range = linspace(-state_set.vz_lim, state_set.vz_lim, 5);
pd_range = linspace(-state_set.pd_lim, state_set.pd_lim, 5);

% those matrix inequalities are linear w.r.t w
w_range = [-state_set.w_lim, state_set.w_lim];

W_fcn = controller.W_fcn;
dW_dt_fcn = controller.dW_dt_fcn;
Y_fcn = controller.Y_fcn;
A_approx_fcn = controller.A_approx_fcn;
A_fcn = controller.A_fcn;
% use a pre-selected
% mu = sdpvar;
% paras = [mu];

paras = [];
alpha = sdpvar;
lambda = sdpvar;
mu = sdpvar;

constraints =[lambda>=1e-3; alpha>=1e-3; mu>=1e-3];
for i = 1:length(p_range)
    for j = 1:length(vx_range)
        x = [0;0;p_range(i);vx_range(j);0;0];                
        W = W_fcn(x);
        Y = Y_fcn(x);

        F2 = [lambda*W  zeros(n,nw) (C*W+D*Y)';
              zeros(nw,n)   (alpha-mu)*eye(nw)  zeros(nw,nz);
              C*W+D*Y    zeros(nz,nw)    alpha*eye(nz)];
        constraints =[constraints F2>=0];

        for k = 1:length(vz_range)
            for l = 1:length(pd_range)
                x = [randn(2,1);p_range(i);vx_range(j);vz_range(k);pd_range(l)];                

                for m = 1:length(w_range)   
                    w = w_range(m);
                    %%%%%%%%%%%%%%%%%%% for testing %%%%%%%%%%%%%%%%%%%
    %                 tmp = A_approx_fcn(x,w)*W + B*Y;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
                    tmp = A_fcn(x,w)*W + B*Y;
                    F1 = [dW_dt_fcn(x,w)-(tmp+tmp')-lambda*W  -Bw_fcn(x); -Bw_fcn(x)' mu*eye(nw)];
                    %%%%%%%%%%%%%%%%%%% for testing %%%%%%%%%%%%%%%%%%%
%                     F1 = [dW_dt_approx_fcn(x,w)-(tmp+tmp')-lambda*W  -Bw_approx_fcn(x); -Bw_approx_fcn(x)' mu_opt*eye(nw)];
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    constraints =[constraints F1>=0];
                end
            end
        end
    end
end 

% solve the problem
ops = sdpsettings('solver','mosek','verbose',0);
sol = optimize(constraints,alpha,ops);
disp(sol.info);
if (sol.problem == 0 || sol.problem == 4)     
    alpha_opt = value(alpha);  
    mu_opt = value(mu);
    lambda_opt = value(lambda);
else
    alpha_opt = inf;
    mu_opt = inf;
    lambda_opt = inf;
end
end