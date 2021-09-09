function dW_dxi = dW_dxi_fcn1(i,x) 

dW_dxi = (i==3)*dW_dphi(x)+(i==4)*dW_dvx(x);
