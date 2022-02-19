
function [alpha_opt,mu_opt,lambda_opt] = rccm_analysis_quad(plant,controller,state_set)
%% compute a bound on the RCCM-based feedback control effort given the plant and W and Y. 
% x is the variables that appear in W and Y
tolerance = 1e-6;
n = plant.n; nw = plant.nw;  
B = plant.B; Bw = plant.Bw;
% C = zeros(nu,n); D = eye(nu); %To set z = u, we need
C = plant.C; D = plant.D;

nz = size(C,1);
ctrl_N = 7;

T_range = linspace(state_set.T_lim(1), state_set.T_lim(2), 4);
phi_range = linspace(-state_set.phi_lim, state_set.phi_lim, ctrl_N);
theta_range = linspace(-state_set.theta_lim, state_set.theta_lim, ctrl_N);

% Note that LMI1 is affine w.r.t to u; therefore, evaluation at the
% vertices is enough
Tdot_range = [state_set.Tdot_lim, -state_set.Tdot_lim];
phidot_range = [state_set.phidot_lim, -state_set.phidot_lim];
thetadot_range = [state_set.thetadot_lim, -state_set.thetadot_lim];

W_fcn = controller.W_fcn;
dW_dt_fcn = controller.dW_dt_fcn;
Y_fcn = controller.Y_fcn;
A_fcn = controller.A_fcn;
% use a pre-selected
% mu = sdpvar;
% paras = [mu];

alpha = sdpvar;
lambda = sdpvar;
mu = sdpvar;

cons =[lambda>=1e-3; alpha>=1e-3; mu>=1e-3];
for i = 1:length(T_range)
    for j = 1:length(phi_range) 
        for k = 1:length(theta_range)
            x0 = [randn(6,1);T_range(i);phi_range(j);theta_range(k)];    
            W = W_fcn(x0);
            Y = Y_fcn(x0);            
            tmp = A_fcn(x0)*W + B*Y;  
            for p = 1:length(Tdot_range)
                for q = 1:length(phidot_range)
                    for r = 1:length(thetadot_range)
                        u0 = [Tdot_range(p);phidot_range(q);thetadot_range(r)];
                        dW_dt = dW_dt_fcn(x0,u0); 
                        LMI1 = [dW_dt-(tmp+tmp')-lambda*W  -Bw; -Bw' mu*eye(nw)];  
                        cons = [cons LMI1>=0];
                    end
                end
            end
            LMI2 = [lambda*W  zeros(n,nw) (C*W+D*Y)';
                  zeros(nw,n)   (alpha-mu)*eye(nw)  zeros(nw,nz);
                  C*W+D*Y    zeros(nz,nw)    alpha*eye(nz)];
            cons = [cons LMI2>=0];
        end
    end
end

% solve the problem
ops = sdpsettings('solver','mosek','verbose',0);
% disp('LMI formulation finished! Solving...');
sol = optimize(cons,alpha,ops);
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