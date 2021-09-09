function g = energyGradient1(c,n,D,N,T,T_dot,w) 
%% Compute the gradient of the Riemann Energy under the pseudospectral method
% for code generation
persistent c_pre g_pre;
if isempty(c_pre) % adding this may make the gradient calculatioin inaccurate
    c_pre = zeros(n,D+1);
    g_pre = zeros(1,(D+1)*n);
end
%     gamma = zeros(n,N+1);
%     gamma_s = zeros(n,N+1);
%     for i = 1:n   
%        gamma(i,:) = c((i-1)*(D+1)+1:i*(D+1),:)'*T;       % gamma(i) is 1*(N+1); the ith elment of gamma on all the (N+1) nodes
%        gamma_s(i,:) = c((i-1)*(D+1)+1:i*(D+1),:)'*T_dot;
%     end 

c = transpose(reshape(c,D+1,n)); % the ith row corresponds to the ith element
if norm(c-c_pre)> 1e-5
gamma = c*T;
gamma_s = c*T_dot;
g = zeros(1,(D+1)*n);  
%     M_x_gamma_s = zeros(n,N+1);

% vectorized format
for k = 1:N+1   
%     if norm(gamma(:,k))> 10
%         disp('gamma norm is out of range');
%     end
    W_fcn_eval = W_fcn1(gamma(:,k));
%         M_x_gamma_sk = (W_fcn(gamma(:,k))\gamma_s(:,k));
    M_x_gamma_sk = W_fcn_eval\gamma_s(:,k);
%         coder.extrinsic('dW_fcn');
    for i = 1:n
        dW_dxi = dW_dxi_fcn1(i,gamma(:,k)); 
%             dW_dxi = dW_fcn{i}(gamma(:,k)); 
        g((i-1)*(D+1)+(1:D+1)) = g((i-1)*(D+1)+(1:D+1))+M_x_gamma_sk'*([zeros(i-1,D+1);T_dot(1:D+1,k)';zeros(n-i,D+1)]*2-dW_dxi*M_x_gamma_sk*T(1:D+1,k)')*w(k);
    end
end
    c_pre = c;
    g_pre = g;
else
    g = g_pre;
end
end