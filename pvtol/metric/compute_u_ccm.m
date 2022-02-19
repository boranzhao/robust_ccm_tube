function u = compute_u_ccm(x,x_nom,u_nom,Erem,gamma,gamma_s,controller_ccm_law_form,...
    controller_ccm_mat_ineq_form,controller_lambda,plant_B,plant_fx,plant_f_fx_nom,geodesic_w_cheby)

if controller_ccm_law_form == CtrlDesignOpts.ccm_integration && controller_ccm_mat_ineq_form == ...
    CtrlDesignOpts.ccm_mat_ineq_use_rho
    u = u_nom;
    error('Code for using rho has not been generated yet!');
%     for k=1:N+1 
%         u = u-0.5*controller_rho_fcn(gamma(:,k))*geodesic_w_cheby(k)*(plant_B'*((W_fcn1(gamma(:,k))+1e-6*eye(n))\gamma_s(:,k)));
%     end
    
else 
    gamma_s1_x_Mx = gamma_s(:,end)'/W_fcn1(x);
    phi0 = gamma_s1_x_Mx*(plant_fx + plant_B*u_nom) - ...
        gamma_s(:,1)'/W_fcn1(x_nom)*(plant_f_fx_nom + plant_B*u_nom) + ...
        controller_lambda*Erem;
    if phi0 <=0
        u = u_nom;
    else
        phi1 = gamma_s1_x_Mx*plant_B;
        u = u_nom- phi0*phi1'/ (phi1*phi1');
    end
end