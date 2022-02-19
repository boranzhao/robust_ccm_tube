function geodesic = setup_geo_opt_quad(n,D,N,controller)
% optimization variables: chebyshev coefficients for geodesics
% [c_10, c_11,..., c_1D, ..., c_n0, c_n1,..., c_nD]

% --------------------obtain chebyshev pseudospectral numerics-------------
% --------------------(to be used for computing the integral)--------------
[s,w_cheby] = clencurt(N); % t is 1 by N+1, with values lying b/t 0 and 1.


% compute Cheby basis at all points
[T, Tdot] = compute_cheby(N,D,s); % Both T and T_dot are D+1 by N+1
% for equality constraints
Aeq = [kron(eye(n),T(:,1)'); kron(eye(n),ones(1,D+1))];
Aeq = sparse(Aeq);

% --------- formulate and solve the NLP problem using OPTI --------------
ndec = n*(D+1);
if controller.use_generated_code == 1
    costf = @(c) RiemannEnergy1_mex(c,n,D,N,T,Tdot,w_cheby);
    grad = @(c) energyGradient1_mex(c,n,D,N,T,Tdot,w_cheby);     
else    
    costf = @(c) RiemannEnergy(c,n,D,N,T,Tdot,w_cheby,controller.W_fcn);
    grad = @(c) energyGradient(c,n,D,N,T,Tdot,w_cheby,controller.W_fcn,controller.dW_dxi_fcn); 
end

geodesic.D = D; geodesic.N = N; geodesic.ndec = ndec;
geodesic.T = T; geodesic.Tdot = Tdot;
geodesic.Aeq = Aeq; 
geodesic.costf = costf;
geodesic.grad = grad;
geodesic.w_cheby = w_cheby;


beq = zeros(2*n,1);c0 = zeros(n*(D+1),1);
% add some bounds to mitigate numerical issues when computing the geodesic
lb = -10*ones(n,D+1);
ub = 10*ones(n,D+1);
lb(7,:) = -20*ones(1,D+1);             % total thrust     
ub(7,:) = 50*ones(1,D+1);   % total thrust
% lb(8:9,:) = -5*ones(2,D+1);         % phi and theta
% ub(8:9,:) = 5*ones(2,D+1);          % phi and theta
lb = lb';lb = lb(:);
ub = ub';ub= ub(:);

% --------------- re-generate the code is necessary after change of ------
% controller or geodesic optimization settings: remember to re-generate the
% m-file functions, e.g., dW_dphi, dW_dvx, etc., first.------------------  
if controller.use_generated_code 
    answer = questdlg('Are the generated codes for this particular scenario?','Question for using C-code in simulation','Yes','No','No');
    switch answer 
        case 'Yes'
        case 'No'
            error('You cannot continue without including the generated codes for this scenario!');
    end
end

% ----------------------for using matlab fmincon solver--------------------
opts_opti = optiset('solver','matlab','maxiter',500,'tolrfun',1e-4,'tolafun',1e-4,'display','off','derivCheck','off'); 
Opt = opti('fun',costf,'grad',grad,'eq',Aeq,beq,'bounds',lb,ub,'ndec',ndec,'x0',c0,'options',opts_opti);
geodesic.nlprob = convMatlab(Opt.prob,opts_opti); 
geodesic.nlprob.options = ...
optimoptions(@fmincon,'Display','off','HessianApproximation','lbfgs',...
'MaxIterations',opts_opti.maxiter,'SpecifyObjectiveGradient',true,'CheckGradients',false,...
'OptimalityTolerance',opts_opti.tolrfun,'FunctionTolerance',opts_opti.tolrfun,'FunValCheck','on','StepTolerance',1.0e-8);
geodesic.opts_opti = opts_opti;
% -------------------------------------------------------------------------

end