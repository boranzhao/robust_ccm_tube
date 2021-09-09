function [alpha_opt,mu_opt,w_lower,max_residual] = rccm_ppg(plant,controller,W,dW_dt,Y,vars_WY,lambda,state_set)
% x is the variables that appear in W and Y
tolerance = 1e-6;
n = plant.n;
nw = plant.nw;
nz = plant.nz;
A = plant.A;
B = plant.B;
Bw = plant.Bw;
C = plant.C;
D = plant.D;

mu = sdpvar;alpha = sdpvar;   

y = sdpvar(n+nw+nz,1);   % for converting a matrix-valuded SOS to a scalar-valued SOS
y1 = y(1:n+nw);
yW = y(1:n);

tmp = A*W + B*Y;
M_pos1 = [dW_dt-(tmp+tmp')-lambda*W  -Bw; -Bw' mu*eye(nw)];
% tmp =  df_dx*W + B*Y;
yMy1 = y1'*M_pos1*y1;
M_pos2 = [lambda*W  zeros(n,nw) (C*W+D*Y)';
      zeros(nw,n)   (alpha-mu)*eye(nw)  zeros(nw,nz);
      C*W+D*Y    zeros(nz,nw)    alpha*eye(nz)];
yMy2 = y'*M_pos2*y;

W_lower = sdpvar; 
y_W_Wlower_y =  yW'*(W-W_lower*eye(n))*yW;
paras = [vars_WY;mu;W_lower];

F = [];
if state_set.consider_state_set == 1
    % considering the compact set using S procedure
    for i=1:length(state_set.box_lim) 
        if i<= state_set.num_consts_4_W_states  
            % enforce a lower bound for W            
            [~,c_lower,v_lower] = polynomial([state_set.W_states;yW],state_set.lagrange_deg_W);
            
            if ~(isfield(state_set,'lagrange_quadratic_only') && state_set.lagrange_quadratic_only == 0)
                % only take the terms quadratic in y
                index = [];

                for k=1:length(v_lower)
                    if sum(degree(v_lower(k),yW)) == 2
                        index = [index k];
                    end
                end
                c_lower = c_lower(index); v_lower = v_lower(index);
            end
            
            L_lower = c_lower'*v_lower;            
            y_W_Wlower_y =  y_W_Wlower_y-L_lower*state_set.box_lim(i);
            
            % Note that the second inequality only depends on W_states
            [~,c_L2,v_L2] = polynomial([state_set.W_states;y],state_set.lagrange_deg_W);
            
            if ~(isfield(state_set,'lagrange_quadratic_only') && state_set.lagrange_quadratic_only == 0)
                % only take the terms quadratic in y
                index = [];
                for k=1:length(v_L2)
                    if sum(degree(v_L2(k),y)) == 2
                        index = [index k];
                    end
                end
                c_L2 = c_L2(index); v_L2 = v_L2(index);
            end

            L2 = c_L2'*v_L2;

            yMy2 = yMy2 - L2*state_set.box_lim(i);
            paras = [paras;vec(c_lower);vec(c_L2)];
            F = [F sos(L_lower) sos(L2)];
        end

        % Lagrange multiplier
        if class(Bw(state_set.W_states_index)) == "sdpvar" 
            % W depends on some states whose derivatives are directly influenced by w
            variables = [state_set.W_states;state_set.other_lim_states;plant.w;y1];
        else
            variables = [state_set.W_states;state_set.other_lim_states;y1];            
        end
        
        [~,c_L1,v_L1] = polynomial(variables,state_set.lagrange_deg_ccm);
        if ~(isfield(state_set,'lagrange_quadratic_only') && state_set.lagrange_quadratic_only == 0)
            % only take the terms quadratic in y
            index = [];
            for k=1:length(v_L1)
                if sum(degree(v_L1(k),y1)) == 2
                    index = [index k];
                end
            end
            c_L1 = c_L1(index); v_L1 = v_L1(index);
        end

        L1 = c_L1'*v_L1;

        yMy1 = yMy1 - L1*state_set.box_lim(i);        
        paras = [paras;vec(c_L1)];
        F = [F sos(L1)];
    end
end

% SOS constraints
F = [F  sos(yMy1) sos(yMy2) sos(y_W_Wlower_y) mu>=tolerance alpha>=tolerance W_lower >= controller.W_lower_bound];
% solve the SOS problem
ops = sdpsettings('solver','mosek','verbose',0);
disp('Problem formulation finished! Start solving...');
[sol,~,~,res] = solvesos(F,alpha,ops,paras);
max_residual = max(res)
disp(sol.info)
if (sol.problem == 0 || sol.problem == 4) 
    alpha_opt = value(alpha);  
    mu_opt = value(mu);
else
    alpha_opt = inf;
    mu_opt = inf;
end
w_lower = value(W_lower);
end