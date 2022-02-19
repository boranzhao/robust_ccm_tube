function [u,Erem]= ccm_law(t,x_nom,u_nom,x,plant,controller)
geodesic = controller.geodesic; 
n = plant.n; N = geodesic.N; D = geodesic.D; 
persistent t_pre beq_pre copt_pre Erem_pre
if isempty(t_pre) || (t == 0)
    t_pre = -3;
    beq_pre = zeros(2*n,1);
    copt_pre = zeros(n*(D+1),1);
    Erem_pre = Inf;
end

w_nom = controller.w_nom;
u = u_nom;

% get the initial value of c corresponding to a straight line
c0 = zeros(n*(D+1),1);
%     for i =1:n    
%         c0((i-1)*(D+1)+1,1) = xStar(i);
%         c0((i-1)*(D+1)+2,1) = x(i)- xStar(i);    
%     end
% vectorized format to improve computational efficiency
i =1:n;    
c0((i-1)*(D+1)+1,1) = x_nom;
c0((i-1)*(D+1)+2,1) = x - x_nom; 

if isempty(geodesic.nlprob) 
    % indicates that the geodesic is (assumed to be) a straight line    
    c0 = transpose(reshape(c0,D+1,n)); % the ith row corresponds to the ith element
    gamma = c0*geodesic.T;
    gamma_s = c0*geodesic.Tdot;
    
    Erem = 0;
    for k=1:N+1  
        tmp = gamma_s(:,k)'*(controller.W_fcn(gamma(:,k))\gamma_s(:,k))*geodesic.w_cheby(k);
        Erem = Erem+ tmp ; % noite that W_fcn needs to be selected for each specific example. 
    end 
else
    beq = [x_nom;x];    

    % tic;
    if norm(beq-beq_pre)<1e-8 && ~isinf(Erem_pre)
        copt = copt_pre;
        Erem = Erem_pre;
    else
        % --------------- ipopt ----------------------------------------------
        % geodesic.nlprob.options.rl = beq;
        % geodesic.nlprob.options.ru = beq;
        % geodesic.nlprob.x0 = c0;
        % [copt,Erem,exitflag,info] = opti_ipopt(geodesic.nlprob,c0);
        % --------------------------------------------------------------

        % ---------------- matlab -----------------------
        geodesic.nlprob.beq = beq;
        geodesic.nlprob.x0 = c0;
        [copt,Erem,exitflag,info] = fmincon(geodesic.nlprob);
        if exitflag<0
            disp('geodesic optimization problem failed!');
        end
        % ------------------------------------------------
        beq_pre = beq;
        copt_pre = copt;
        Erem_pre = Erem;
    end
    % toc;
    % ----------------- compute the control law -----------------------
    %     tic;
    %     gamma = zeros(n,N+1);
    %     gamma_s = zeros(n,N+1);  
    %     for i = 1:n   
    %        gamma(i,:) = copt((i-1)*(D+1)+1:i*(D+1),:)'*T;       % gamma(i) is 1*(N+1); the ith elment of gamma on all the (N+1) nodes
    %        gamma_s(i,:) = copt((i-1)*(D+1)+1:i*(D+1),:)'*T_dot;
    %     end  
    %     toc;
    % vectorized format (more computationally efficient)
    copt = transpose(reshape(copt,D+1,n)); % the ith row corresponds to the ith element
    gamma = copt*geodesic.T;
    gamma_s = copt*geodesic.Tdot;
    % ----------------------------------------------------------------

    % % -------- verify whether the curve found is really a geodesic ----------
    % % according to equation (11) in Leung  & Manchester
    % error = 0;
    % for k=1:N+1
%         error = error + (gamma_s(:,k)'*(controller.W_fcn(gamma(:,k))\gamma_s(:,k))-Erem)^2*geodesic.w_cheby(k);
    % end
    % error = sqrt(error)/Erem;
    % if error>=1e-5
    % %     disp('The curve optimized is probably not a geodesic!');
    %     fprintf(1,'t= %.2e, Error = %.3e, the curve optimized is probably not a geodesic!\n',t,error);
    %     if error> 1e-2
    %         pause;
    %     end
    % end
    % % -----------------------------------------------------------------------
end
% tic;
plant_fx = plant.f_fcn(x);
plant_fx_nom =  plant.f_fcn(x_nom);
% numerical integral
if controller.type == CtrlDesignOpts.ccm % ccm
    if isfield(controller,'use_generated_code') && controller.use_generated_code == 1 && ~isempty(geodesic.nlprob) 
    % can further accelerate the computation
        u = compute_u_ccm_mex(x,x_nom,u_nom,Erem,gamma,gamma_s,int8(controller.ccm_law_form),...
                int8(controller.ccm_mat_ineq_form),controller.lambda,plant.B,plant_fx,plant_fx_nom,geodesic.w_cheby);
    else
        if controller.ccm_law_form == CtrlDesignOpts.ccm_integration && controller.ccm_mat_ineq_form == ...
            CtrlDesignOpts.ccm_mat_ineq_use_rho
            for k=1:N+1 
                u = u-0.5*controller.rho_fcn(gamma(:,k))*geodesic.w_cheby(k)*(plant.B'*((controller.W_fcn(gamma(:,k))+1e-6*eye(n))\gamma_s(:,k)));
            end
        else 
            gamma_s1_Mx = gamma_s(:,end)'/controller.W_fcn(x);
            phi0 = gamma_s1_Mx*(plant_fx + plant.B*u_nom) - ...
                gamma_s(:,1)'/controller.W_fcn(x_nom)*(plant.f_fcn(x_nom) + plant.B*u_nom) + ...
                controller.lambda*Erem;
            if phi0 <=0
                u = u_nom;
            else
                phi1 = gamma_s1_Mx*plant.B;
                phi1_sq = phi1*phi1';                
                if phi1_sq< 1e-8
                    u = u_nom;
                else
                    u = u_nom- phi0*phi1'/ (phi1*phi1'+eps);
                end
            end
        end
    end
elseif controller.type == CtrlDesignOpts.rccm  % robust ccm minimizing Linf gain
    if isfield(controller,'use_generated_code')&& controller.use_generated_code == 1 && ~isempty(geodesic.nlprob) && (isfield(plant,'Bw') && strcmp(plant.Bw,'double'))
        u = compute_u_rccm_mex(x,x_nom,u_nom,w_nom,Erem,int8(N),gamma,gamma_s,int8(controller.ccm_law_form),controller.lambda,controller.mu,plant.B,plant.Bw,plant_fx,plant_fx_nom,geodesic.w_cheby);
    else
        if controller.ccm_law_form == CtrlDesignOpts.ccm_integration
            for k=1:N+1 
                u = u + geodesic.w_cheby(k)*(controller.Y_fcn(gamma(:,k))*(controller.W_fcn(gamma(:,k))\gamma_s(:,k)));
            end
        elseif controller.ccm_law_form == CtrlDesignOpts.ccm_min_norm            
            gamma_s1_Mx = gamma_s(:,end)'/controller.W_fcn(x);
            w_hat =  w_nom + 1/controller.mu*(gamma_s1_Mx*plant.Bw_fcn(x))';
            x_nom_dot = plant.f_fcn(x_nom) + plant.B*u_nom + plant.Bw_fcn(x)*w_nom;
            phi0 = gamma_s1_Mx*(plant_fx + plant.B*u_nom+plant.Bw_fcn(x)*w_hat) - ...
                gamma_s(:,1)'/controller.W_fcn(x_nom)*x_nom_dot - controller.mu/2*norm(w_hat-w_nom)+ ...
                0.5*controller.lambda*Erem;
            if phi0 <=0
                u = u_nom;
            else
                phi1 = gamma_s1_Mx*plant.B;
                u = u_nom - phi0*phi1'/ (phi1*phi1'+eps);
            end            
        end  
    end        
end    

if (t-t_pre>= 0.4) && mod(t,1)< 0.1
    fprintf('t = %.1f s\n',t);
    t_pre = t;
end
end