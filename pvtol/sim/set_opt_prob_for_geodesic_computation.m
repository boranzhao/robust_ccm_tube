function geodesic = set_opt_prob_for_geodesic_computation(n,D,N,controller)
% optimization variables: chebyshev coefficients for geodesics
% [c_10, c_11,..., c_1D, ..., c_n0, c_n1,..., c_nD]

% --------------------obtain chebyshev pseudospectral numerics-------------
% --------------------(to be used for computing the integral)--------------
[s,w_cheby] = clencurt(N); % t is 1 by N+1, with values lying b/t 0 and 1.
% evaluate the value of the CCM
W = zeros(n,n,N+1);

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
lb = -20*ones(n,D+1);
ub = 20*ones(n,D+1);
lb(3:4,:) = -5*ones(2,D+1);  % phi and vx
ub(3:4,:) = 5*ones(2,D+1);   % phi and vx
lb = lb';lb = lb(:);
ub = ub';ub= ub(:);


%     [copt1,Erem,exitflag,info] = solve(Opt,c0);

% ---------- for using ipopt solver ---------------------------------------
% opts_opti = optiset('solver','ipopt','maxiter', 500,'display','iter');  %,,'derivCheck','on'
% Opt = opti('fun',geodesic.costf,'grad',geodesic.grad,'eq',geodesic.Aeq,beq,...
%             'bounds',lb,ub,'ndec',geodesic.ndec,'x0',c0,'options',geodesic.opts_opti);
% geodesic.nlprob = convIpopt(Opt.prob,geodesic.opts_opti); 
% -------------------------------------------------------------------------

% ----------for using matlab fmincon solver--------------------------------
opts_opti = optiset('solver','matlab','maxiter',500,'tolrfun',1e-4,'tolafun',1e-4,'display','off','derivCheck','off'); 
Opt = opti('fun',costf,'grad',grad,'eq',Aeq,beq,'bounds',lb,ub,'ndec',ndec,'x0',c0,'options',opts_opti);
geodesic.nlprob = convMatlab(Opt.prob,opts_opti); 
geodesic.nlprob.options = ...
optimoptions(@fmincon,'Display','off','HessianApproximation','lbfgs',...
'MaxIterations',opts_opti.maxiter,'SpecifyObjectiveGradient',true,'CheckGradients',false,...
'OptimalityTolerance',opts_opti.tolrfun,'FunctionTolerance',opts_opti.tolrfun,'FunValCheck','on','StepTolerance',1.0e-8);
% -------------------------------------------------------------------------

geodesic.opts_opti = opts_opti;
