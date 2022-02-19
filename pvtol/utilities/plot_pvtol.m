function plot_pvtol(x,z,phi,c)
%% PVTOL cartoon
len = 0.25;
pvtol_bnds = [-len, len;
                0, 0]; %at 0 degrees

R = [cos(phi), -sin(phi); sin(phi), cos(phi)];
bnds = R*pvtol_bnds;
hold on
plot(bnds(1,:)'+x,bnds(2,:)'+z,'o-','color',c,'markersize',7,'markerfacecolor',c);
end