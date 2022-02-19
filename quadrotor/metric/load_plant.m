n = 9;                          % state dimension (ignoring psi)
nu = 3;                         % input dimension
nw = 3;                         % disturbance dimension
plant.g = 9.8066;                       % gravitational constant
x = msspoly('x',n); x_store = x;
u = msspoly('u',nu); plant.u= u;
      
% approximating sin_x/cos_x with Chebshev polynomials
sinx = @(x) 0.9101*(x/(pi/3)) - 0.04466*(4*(x/(pi/3))^3 - 3*(x/(pi/3))); % 0.8799 (close to 0.9101*3/pi, -0.03915
cosx = @(x) 0.7441-0.2499*(2*(x/(pi/3))^2 -1);                            % 0.7652, -0.2299

% ------------------------- system dynamics -------------------------------
wdist = msspoly('w',nw); plant.wdist = wdist;           % disturbance

sin_r = sinx(x(8));             % roll -> phi
cos_r = cosx(x(8));
sin_p = sinx(x(9));             % pitch -> theta
cos_p = cosx(x(9));

% directly approximate the product terms 
% cosrcosp = 1 - 0.5*x(8)^2 - 0.5*x(9)^2;
% sinrsinp = sin_r*sin_p;
% cosrcosp = cos_r*cos_p;
% sinrcosp = sin_r*cos_p;
% cosrsinp = cos_r*sin_p;


% up to 4th order
sinrsinp = (0.99405)*x(9)*x(8)+(-0.1551)*x(9)*x(8)^3+(-0.1551)*x(9)^3*x(8);
cosrcosp = (0.98804)+(-0.45303)*x(8)^2+(-0.45303)*x(9)^2+(0.20772)*x(9)^2*x(8)^2;
sinrcosp = (0.99104)*x(8)+(-0.15462)*x(8)^3+(-0.45441)*x(9)^2*x(8);
cosrsinp = (0.99104)*x(9)+(-0.15462)*x(9)^3+(-0.45441)*x(9)*x(8)^2;

% up to 3rd order
% sinrsinp = (0.99405)*x(9)*x(8)+(-0.1551)*x(9)*x(8)^3;
% cosrcosp = (0.98804)+(-0.45303)*x(8)^2+(-0.45303)*x(9)^2;
% sinrcosp = (0.99104)*x(8)+(-0.15462)*x(8)^3+(-0.45441)*x(9)^2*x(8);
% cosrsinp = (0.99104)*x(9)+(-0.15462)*x(9)^3+(-0.45441)*x(9)*x(8)^2;


b_T = [sin_p; -sinrcosp; cosrcosp];

f = [[x(4);x(5);x(6)];           %px,py,pz
     [0;0;plant.g] - b_T*x(7);   %vx,vy,vz
     zeros(3,1)];                %T, phi, theta

%gradients       
db_T_q = [0, cos_p;
         -cosrcosp, sinrsinp;
         -sinrcosp,-cosrsinp];
     
% f written as a function handle: can also work when x has multiple columns
f_fcn = @(x)[x(4:6,:);                           %px, py, pz
            -x(7,:).*sin(x(9,:));                %vx
            x(7,:).*cos(x(9,:)).*sin(x(8,:));    %vy
            plant.g-x(7,:).*cos(x(9,:)).*cos(x(8,:));  %vz
            zeros(3,size(x,2))];                 %T, phi, theta     

b_T_fcn = @(x) [sin(x(9)); -cos(x(9))*sin(x(8)); cos(x(9))*cos(x(8))];

%gradients
db_T_q_fcn =@(x) [0, cos(x(9));
    -cos(x(8))*cos(x(9)), sin(x(8))*sin(x(9));
    -sin(x(8))*cos(x(9)),-cos(x(8))*sin(x(9))];
B = [zeros(6,3); eye(3)]; 
B_perp = [eye(6); zeros(3,6)];
Bw = [zeros(3);eye(3);zeros(3,3)]; 
Bw_fcn = @(x)[zeros(3);eye(3);zeros(3,3)];

df_dx = [zeros(3), eye(3),zeros(3,3);
         zeros(3,6),-b_T, -db_T_q(:,1)*x(7), -db_T_q(:,2)*x(7);
         zeros(3,9)];
df_dx_fcn = @(x) [zeros(3), eye(3),zeros(3,3);
         zeros(3,6),-b_T_fcn(x), -(db_T_q_fcn(x)*[1;0])*x(7), -(db_T_q_fcn(x)*[0;1])*x(7);
         zeros(3,9)];

A = df_dx ; % note that A does not depend on u or w
plant.n = n; plant.nu=nu; plant.nw = nw; 
plant.sinx = sinx;
plant.cosx = cosx;
plant.df_dx = df_dx;
plant.f_fcn = f_fcn;
plant.A = A; plant.B = B;
plant.B_fcn = @(x) B;
plant.dynamics = @(x,u) f_fcn(x)+ B*u;
plant.B_perp = B_perp;
plant.Bw = Bw;
plant.Bw_fcn = Bw_fcn;
plant.df_dx_fcn = df_dx_fcn;
plant.A_fcn = df_dx_fcn;
plant.dBw_dx_fcn = @(x) zeros(n,3);
% choose the outputs whose deviation from their nominal trajs we want to minimize
if controller.opt_pos_dev == 0
      C= [eye(n); zeros(nu,n)]; D =  [zeros(n,nu); 0.01*diag([2 5 5])];  % consider both states and inputs
%        C= [eye(n)]; D =  [zeros(n,nu)];  % consider both states and inputs
elseif controller.opt_pos_dev == 1
%     C = [diag([1 1 1 0.1 0.1 0.1 0.01 0.1 0.1]);zeros(nu,n)]; 
%     D = [zeros(n,nu); 0.01*diag([2 5 5])];                               % consider only position states and inputs
    C = [eye(3) zeros(3,n-3);zeros(nu,n)]; 
    D = [zeros(3,nu); 0.01*diag([2 5 5])];                               % consider only position states and inputs
end
nz = size(C,1); plant.nz = nz; 
plant.C = C; plant.D = D;

df_dx_approx_fcn_vec = mss2fnc(df_dx(:),x,randn(length(x),2));
df_dx_approx_fcn = @(x) reshape(df_dx_approx_fcn_vec(x),9,9);
plant.df_dx_approx_fcn = df_dx_approx_fcn;
