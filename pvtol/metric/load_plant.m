plant.g = 9.81;     % (m/s^2) gravity
plant.l = 0.25;     % (m) half-width of quad rotor
plant.m = 0.486;    % (m) mass of the quad rotor
plant.J = 0.00383;  % (kgm^2), moment of inertia
n = 6;              % #states
nu =2;              % #inputs 
nw = 1;             % #disturbances


x = sdpvar(n,1); 
x_store = x;
% approximating sin/cos with Chebshev polynomials (for applying SOS
% programming)
sinx = @(x) 0.9101*(x./(pi/3)) - 0.04466*(4*(x./(pi/3)).^3 - 3*(x./(pi/3))); 
cosx = @(x) 0.7441 -0.2499*(2*(x./(pi/3)).^2 -1);                            

%% system dynamics
w = sdpvar(nw,1);                   % disturbance
sin_p = sinx(x(3)); cos_p = cosx(x(3));
f = [x(4)*cos_p - x(5)*sin_p;       %px
    x(4)*sin_p + x(5)*cos_p;        %pz
    x(6);                           %phi
    x(6)*x(5)-plant.g*sin_p;        %vx
    -x(6)*x(4)-plant.g*cos_p;       %vz
    0];                             %phi_dot

% f written as a function handle: can also work when x has multiple columns
f_fcn = @(x) [x(4,:).*cos(x(3,:)) - x(5,:).*sin(x(3,:));    %px
            x(4,:).*sin(x(3,:)) + x(5,:).*cos(x(3,:));      %pz
            x(6,:);                                         %phi
            x(6,:).*x(5,:)-plant.g*sin(x(3,:));             %vx
            -x(6,:).*x(4,:)-plant.g*cos(x(3,:));            %vz
            zeros(1,size(x,2))];                 
B = [zeros(4,2); 1/plant.m 1/plant.m; plant.l/plant.J -plant.l/plant.J]; 
B_perp = [eye(4); zeros(2,4)];
Bw = [zeros(1,3),cosx(x(3)),-sinx(x(3)),0]'; 
Bw_fcn = @(x)[zeros(1,3),cos(x(3)),-sin(x(3)),0]';
df_dx = jacobian(f,x);
dBw_dx = jacobian(Bw,x);
A = df_dx + dBw_dx*w;



f_phi_fcn = @(x) x(6);
f_vx_fcn = @(x) x(6)*x(5)-plant.g*sin(x(3));
Bw_phi_fcn = @(x) 0;
Bw_vx_fcn = @(x) cos(x(3));


%% for testing with(Chebyshev polynomials) approximated fcns 
Bw_approx_fcn = @(x)[zeros(1,3),cosx(x(3)),-sinx(x(3)),0]';
f_vx_approx_fcn = @(x) x(6)*x(5)-plant.g*sinx(x(3));
Bw_vx_approx_fcn = @(x) cosx(x(3));
% approximated f_approx, df_dx_fcn
x = x_store;
s = sdisplay(f);
s2 = sdisplay(df_dx);
s3 = sdisplay(dBw_dx);
syms x [n 1]
syms f_approx_fcn [n 1]
syms df_dx_approx_fcn [n n]
syms dBw_dx_approx_fcn [n n]
for i=1:n    
    f_approx_fcn(i,1) = eval(s{i});    
    for j=1:n
        df_dx_approx_fcn(i,j) = eval(s2{i,j});
        dBw_dx_approx_fcn(i,j) =eval(s3{i,j}); 
    end
end
f_approx_fcn = matlabFunction(f_approx_fcn,'Vars',{x});
df_dx_approx_fcn = matlabFunction(df_dx_approx_fcn,'Vars',{x});
dBw_dx_approx_fcn = matlabFunction(dBw_dx_approx_fcn,'Vars',{x});
x = x_store;
%-------------------------------------------------------------------------

plant.sinx = sinx;
plant.cosx = cosx;
plant.df_dx = df_dx;
plant.f_fcn = f_fcn;
plant.A = A;
plant.B = B;
plant.B_fcn = @(x) B;
plant.dynamics = @(x,u) f_fcn(x)+ B*u;
plant.w = w;

plant.B_perp = B_perp;
plant.Bw = Bw;
plant.Bw_fcn = Bw_fcn;
plant.n = n; plant.nu=nu; plant.nw = nw; 


if controller.opt_pos_dev == 0
    C = [eye(n); zeros(nu,n)]; D =  [zeros(n,nu); 1*eye(nu)];  % consider both states and inputs
elseif controller.opt_pos_dev == 1
    C = [eye(2) zeros(2,4);zeros(nu,n)]; 
    D = [zeros(2,nu); 1*eye(nu)];                               % consider only position states and inputs
end
nz = size(C,1); 
plant.C = C; plant.D = D; plant.nz = nz;
