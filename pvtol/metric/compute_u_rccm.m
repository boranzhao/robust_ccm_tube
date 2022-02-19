function u = compute_u_rccm(x,x_nom,u_nom,w_nom,Erem,N,gamma,gamma_s,controller_ccm_law_form,controller_lambda,controller_mu,plant_B,plant_Bw,plant_fx,plant_f_fx_nom,geodesic_w_cheby)
u = u_nom;
if controller_ccm_law_form == CtrlDesignOpts.ccm_integration
    for k=1:N+1 
        u = u + geodesic_w_cheby(k)*(Y_fcn1(gamma(:,k))*(W_fcn1(gamma(:,k))\gamma_s(:,k)));
    end
elseif controller_ccm_law_form == CtrlDesignOpts.ccm_min_norm
    
    gamma_s1_x_Mx = gamma_s(:,end)'/W_fcn1(x);
    w_hat =  w_nom + 1/controller_mu*(gamma_s1_x_Mx*plant_Bw)';
    x_nom_dot = plant_f_fx_nom + plant_B*u_nom + plant_Bw*w_nom;
    phi0 = gamma_s1_x_Mx*(plant_fx + plant_B*u_nom+plant_Bw*w_hat) - ...
        gamma_s(:,1)'/W_fcn1(x_nom)*x_nom_dot - controller_mu/2*norm(w_hat-w_nom)+ ...
        0.5*controller_lambda*Erem;
    if phi0 <=0
        u = u_nom;
    else
        phi1 = gamma_s1_x_Mx*plant_B;
        u = u_nom - phi0*phi1'/ (phi1*phi1');
    end            
end        
