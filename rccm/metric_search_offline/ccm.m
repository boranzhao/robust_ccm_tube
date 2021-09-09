function [cond_num_W_opt,W_upper,W_lower,W_bar,max_residual] = ccm(plant,controller,W,dW_dt,paras_W,lambda,state_set)
n = plant.n;nu = plant.nu;
y = sdpvar(n,1); W_bar = sdpvar(n,n,'symmetric');
% PSD constraint: equation (11)
if controller.ccm_mat_ineq_form == CtrlDesignOpts.ccm_mat_ineq_use_rho    
    M_pos = dW_dt - plant.df_dx*W - (plant.df_dx*W)'+controller.rho*(plant.B*plant.B')-2*lambda*W;
    y1 = y;
    paras = [paras_W;vec(controller.c_rho)];
elseif controller.ccm_mat_ineq_form == CtrlDesignOpts.ccm_mat_ineq_use_B_perp
     M_pos = plant.B_perp'*(dW_dt - plant.df_dx*W - (plant.df_dx*W)'-2*lambda*W)*plant.B_perp;
    y1 = y(1:n-nu);
    paras = paras_W;
end

yMy = y1'*M_pos*y1;
W_lower = sdpvar;
W_upper = sdpvar;
ops =sdpsettings('solver','mosek');
obj = W_upper;

y_Wupper_W_y =  y'*(W_bar-W)*y;
y_W_Wlower_y =  y'*(W-W_lower*eye(n))*y;

paras = [paras;W_lower;W_bar(:)];
F = [];  
if state_set.consider_state_set == 1
    % considering the compact set for states using S procedure        
    for i=1:length(state_set.box_lim)               
        if i<= state_set.num_consts_4_W_states            
            % Lagrange multiliers for enforcing constraints
            [~,c_Ll,v_Ll] = polynomial([state_set.W_states;y],state_set.lagrange_deg_W);
            if ~(isfield(state_set,'lagrange_quadratic_only') && state_set.lagrange_quadratic_only == 0)
                % only take the terms quadratic in y
                index = [];
                for k=1:length(v_Ll)
                    if sum(degree(v_Ll(k),y)) == 2
                        index = [index k];
                    end
                end
                c_Ll = c_Ll(index); v_Ll = v_Ll(index);
            end
            
            Ll = c_Ll'*v_Ll;

            [~,c_Lu,v_Lu] = polynomial([state_set.W_states;y],state_set.lagrange_deg_W);     
            % only take the terms quadratic in y
            if ~(isfield(state_set,'lagrange_quadratic_only') && state_set.lagrange_quadratic_only == 0)
                index = [];
                for k=1:length(v_Lu)
                    if sum(degree(v_Lu(k),y)) == 2
                        index = [index k];
                    end
                end
                c_Lu = c_Lu(index); v_Lu = v_Lu(index);
            end
            Lu = c_Lu'*v_Lu;
            
            y_Wupper_W_y =  y_Wupper_W_y-Lu*state_set.box_lim(i);
            y_W_Wlower_y =  y_W_Wlower_y-Ll*state_set.box_lim(i);
           
            paras = [paras;vec(c_Ll);vec(c_Lu)];
            F = [F sos(Ll) sos(Lu)];
        end
                  
        [~,c_Lm,v_Lm] = polynomial([state_set.W_states;state_set.other_lim_states;y1],state_set.lagrange_deg_ccm); % for M_pos    
        if ~(isfield(state_set,'lagrange_quadratic_only') && state_set.lagrange_quadratic_only == 0)
             % only take the terms quadratic in y
            index = [];
            for k=1:length(v_Lm)
                if sum(degree(v_Lm(k),y)) == 2
                    index = [index k];
                end
            end
            c_Lm = c_Lm(index); v_Lm = v_Lm(index);
        end
        Lm = c_Lm'*v_Lm;
        yMy = yMy - Lm*state_set.box_lim(i);  

        paras = [paras;vec(c_Lm)];
        F = [F sos(Lm)];    
    end
end
F = [F sos(yMy) sos(y_W_Wlower_y) sos(y_Wupper_W_y) W_lower >= controller.W_lower_bound  W_upper*eye(n)>=W_bar]; % 

[sol,v,Q,res] = solvesos(F,obj,ops,paras);
max_residual = max(res)
disp(sol.info);
W_upper = value(W_upper); W_lower = value(W_lower);
cond_num_W_opt = W_upper/W_lower;
% if max_res <1e-3
%     w_upper = value(w_upper); w_lower = value(w_lower);
%     cond_num_W_opt = w_upper/w_lower;
% else
%     cond_num_W_opt = inf;
% end
W_bar = value(W_bar);
