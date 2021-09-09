function [T, T_dot] = compute_cheby(N,D,t)
T = zeros(D+1,N+1); % up to order D
U = zeros(D+1,N+1); 
T_dot = zeros(D+1,N+1);

T(1,:) = 1.0;   %D = 0
T(2,:) = t;     %D = 1

U(1,:) = 1.0;   %D = 0
U(2,:) = 2*t;   %D = 1

T_dot(2,:) = 1.0; %D = 1
for n = 2:D                     %order, posn = n+1
    T(n+1,:) = 2.0*t.*T(n,:) - T(n-1,:);
    U(n+1,:) = 2.0*t.*U(n,:) - U(n-1,:);
    T_dot(n+1,:) = n*U(n,:);    % from the fact that dT(n)/dt = U(n-1)
end
end