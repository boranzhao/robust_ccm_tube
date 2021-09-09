
function E = RiemannEnergy1(c,n,D,N,T,T_dot,w)
% "1" means for code generation 
% persistent c_pre E_pre;
% if isempty(c_pre) % adding this may make the gradient calculatioin inaccurate
%     c_pre = zeros(n,D+1);
%     E_pre = 0;
% end
%#codegen
%% Compute the Riemann Energy using the pseudospectral method
%     gamma = zeros(n,N+1);
%     gamma_s = zeros(n,N+1);
%     for i = 1:n   
%        gamma(i,:) = c((i-1)*(D+1)+1:i*(D+1),:)'*T;       % gamma(i) is 1*(N+1); the ith elment of gamma on all the (N+1) nodes
%        gamma_s(i,:) = c((i-1)*(D+1)+1:i*(D+1),:)'*T_dot;
%     end   
    % vectorized format to improve computational efficiency
c = transpose(reshape(c,D+1,n)); % the ith row corresponds to the ith element

% if norm(c-c_pre)> 1e-8
gamma = c*T;
gamma_s = c*T_dot;
E = 0;
for k=1:N+1  
    tmp = gamma_s(:,k)'*(W_fcn1(gamma(:,k))\gamma_s(:,k))*w(k);
%         if tmp<0
%             pause;
%         end
    E = E+ tmp ; % noite that W_fcn needs to be selected for each specific example. 
end 
%     c_pre = c;
%     E_pre = E;
% else
%     E = E_pre;
% end
