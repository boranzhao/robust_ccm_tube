function R_omg = R_omg(q)
%#codegen

%Euler rates to body rates:
R_omg =    [cos(q(2))*cos(q(3)), sin(q(3)), 0;
           -cos(q(2))*sin(q(3)), cos(q(3)), 0;
            sin(q(2)), 0, 1];

end