function u = u_ccm(t,x, ccm, plant)
lambda = ccm.lambda;
Ts = ccm.Ts;
n  = ccm.n;
D = ccm.D;
ndec = n*(D+1);

x_nom = ccm.x_nom_fcn(t);
u_nom = ccm.u_nom_fcn(t);

beq = [x_nom;x];    
% get the initial value of c corresponding to a straight line
c0 = zeros(n*(D+1),1);
%     for i =1:n    
%         c0((i-1)*(D+1)+1,1) = xStar(i);
%         c0((i-1)*(D+1)+2,1) = x(i)- xStar(i);    
%     end

% vectorized format to improve computational efficiency
i =1:n;    
c0((i-1)*(D+1)+1,1) = x_nom;
c0((i-1)*(D+1)+2,1) = x-x_nom; 

%     tic;  
Opt = opti('fun',ccm.costf,'grad',ccm.grad,'eq',ccm.Aeq,beq,'ndec',ndec,'x0',c0,'options',ccm.opts_opti);
%     Opt.prob.Aeq = Aeq; Opt.prob.beq = beq; Opt.prob.x0 = x0;
[copt,fval,~,~] = solve(Opt,c0);
%     toc;    
energy = fval;  

% vectorized format (more computationally efficient)
copt = transpose(reshape(copt,D+1,n)); % the ith row of copt corresponds to the ith element of x
%     gamma = copt*T;
%     gamma_s = copt*Ts;    
%     for k=1:N+1    
%         u = u-0.5*rho_fcn(gamma(:,k))*w(k)*(B'*(W_fcn(gamma(:,k))\gamma_s(:,k))); 
%     end    

gams_0 = copt*Ts(:,1);
gams_1 = copt*Ts(:,end);
W_x = ccm.W_fcn(x);
W_xnom = ccm.W_fcn(x_nom);       
rhs = -2*(gams_1'*(W_x\plant.dynamics(x,u_nom))- ...
        gams_0'*(W_xnom\plant.dynamics(x_nom,u_nom))+...
        lambda*energy);  
lhs = 2*(W_x\plant.B_fcn(x))'*gams_1;
nlhs = norm(lhs);
if rhs >=0 || (nlhs ==0)
    u = u_nom;
else
    tmp = (rhs/nlhs^2)*lhs;
    if norm(tmp)>0.1
        disp('Large control signal from CCM!')
    end
    u = u_nom-(rhs/nlhs^2)*lhs;        
end  


% % ------------ just for testing -----------------
% energy = 0;
% u = u_nom;
% % -----------------------------------------------

ue = [u;energy];
if t>0.1 && mod(t,0.5)<= 2e-2
    fprintf('t= %.2f s\n',t);
end
end