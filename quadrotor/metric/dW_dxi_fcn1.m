function dW_dxi = dW_dxi_fcn1(i,x) 


dW_dxi = (i==7)*dW_dT(x)+(i==8)*dW_dphi(x)+(i==9)*dW_dtheta(x);
