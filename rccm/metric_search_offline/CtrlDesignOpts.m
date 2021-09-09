classdef CtrlDesignOpts
   properties (Constant)
      open_loop = -1;           % an open-loop controller
      lqr = 0;                  % an LQR controller
      ccm = 1;                  % a CCM controller
      rccm = 2;                 % a RCCM controller
      ccm_min_norm = 1;         % using a min-norm type control law
      ccm_integration = 2;      % using a differential state feedback control law
      ccm_mat_ineq_use_B_perp = 1;    % using the LMI involving B_perp
      ccm_mat_ineq_use_rho = 2;       % using the LMI involving rho (which yields a differential feedback control law). 
   end
end
