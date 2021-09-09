% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% generate c codes from the matlab function files, accelerating the computation by a lot!!
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% variables
n = plant.n; nu = plant.nu; nw= plant.nw;
D = geodesic.D;
N = geodesic.N;
%% for CCM related functions
% codegen W_eval -args {zeros(n,1,'double')}
% codegen W_fcn1 -args {zeros(n,1,'double')}% 
% codegen dW_dxi_fcn1 -args {int8(1),zeros(n,1,'double')}


%% for online computing the geodesic and control signal
% note that N and D is defined in the main file.
codegen RiemannEnergy1 -args {zeros(n*(D+1),1,'double'),0,0,0,zeros(D+1,N+1,'double'),zeros(D+1,N+1,'double'),zeros(1,N+1)}
% movefile('RiemannEnergy_mex.mexw64','RiemannEnergy_pvtol_mex.mexw64', 'f');
codegen energyGradient1 -args {zeros(n*(D+1),1,'double'),0,0,0,zeros(D+1,N+1,'double'),zeros(D+1,N+1,'double'),zeros(1,N+1)}
% movefile('energyGradient_mex.mexw64','energyGradient_pvtol_mex.mexw64','f');

% for computing the control signal
codegen compute_u_ccm -args {zeros(n,1,'double'),zeros(n,1,'double'),zeros(nu,1,'double'),double(0),zeros(n,N+1,'double'),zeros(n,N+1,'double'),int8(1),int8(1),double(1),zeros(n,nu,'double'),zeros(n,1,'double'),zeros(n,1,'double'),zeros(1,N+1,'double')}

codegen compute_u_rccm -args {zeros(n,1,'double'),zeros(n,1,'double'),zeros(nu,1,'double'),zeros(nw,1,'double'),double(0),int8(1),zeros(n,N+1,'double'),zeros(n,N+1,'double'),int8(1),double(1),double(1),zeros(n,nu,'double'),zeros(n,nw,'double'),zeros(n,1,'double'),zeros(n,1,'double'),zeros(1,N+1,'double')}
