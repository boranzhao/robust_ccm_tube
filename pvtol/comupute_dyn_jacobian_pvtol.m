function [] = comupute_dyn_jacobian_pvtol()
syms x [6 1]; 
f  = [x(4,:).*cos(x(3,:)) - x(5,:).*sin(x(3,:));    %px
    x(4,:).*sin(x(3,:)) + x(5,:).*cos(x(3,:));     %pz
    x(6,:);                        %phi
    x(6,:).*x(5,:)-plant.g*sin(x(3,:));         %vx
    -x(6,:).*x(4,:)-plant.g*cos(x(3,:));        %vz
    zeros(1,size(x,2))]; 

df_dx = jacobian(f,x);
matlabFunction(df_dx,'File','df_dx','Vars',{x});
end

        