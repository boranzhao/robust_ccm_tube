% Create monomial functions
w_poly_fcn = mss2fnc(w_poly,x,randn(length(x),2));
dw_poly_T_fcn = mss2fnc(dw_poly_dx(:,1),x,randn(length(x),2));
dw_poly_phi_fcn = mss2fnc(dw_poly_dx(:,2),x,randn(length(x),2));
dw_poly_theta_fcn = mss2fnc(dw_poly_dx(:,3),x,randn(length(x),2));

% write a few functions to m files
W_exec = 'W_eval = @(ml)';
w_eval_str = '';
for i = 1:length(w_poly)
    if i<length(w_poly)
        W_exec = strcat(W_exec,sprintf('W_coef(:,:,%d)*ml(%d) +',i,i));
        if i == 1
            w_eval_str = strcat(w_eval_str,mat2str(W_coef(:,:,i),5),sprintf('*ml(%d) +',i),'...\n');
        else
            w_eval_str = [w_eval_str,'      ',mat2str(W_coef(:,:,i),5),sprintf('*ml(%d) +',i),'...\n'];
        end
    else
        W_exec = strcat(W_exec,sprintf('W_coef(:,:,%d)*ml(%d);',i,i));
        w_eval_str = [w_eval_str,'      ',mat2str(W_coef(:,:,i)),sprintf('*ml(%d);',i)];
    end
end
w_poly_fcn_str = func2str(w_poly_fcn);
dw_poly_T_fcn_str = func2str(dw_poly_T_fcn);
dw_poly_phi_fcn_str = func2str(dw_poly_phi_fcn);
dw_poly_theta_fcn_str = func2str(dw_poly_theta_fcn);

w_poly_fcn_str = strcat('ml = ', w_poly_fcn_str(5:end),';\n');
dw_poly_T_fcn_str = strcat('ml = ', dw_poly_T_fcn_str(5:end),';\n');
dw_poly_phi_fcn_str = strcat('ml = ', dw_poly_phi_fcn_str(5:end),';\n');
dw_poly_theta_fcn_str = strcat('ml = ', dw_poly_theta_fcn_str(5:end),';\n');

w_fcn_str = strcat('function W = W_fcn1(x)\n', w_poly_fcn_str,'W = ',w_eval_str,'\nend');
dW_dT_str = strcat('function dW_dT = dW_dT(x)\n', dw_poly_T_fcn_str,'dW_dT = ',w_eval_str,'\nend');
dW_dphi_str = strcat('function dW_dphi = dW_dphi(x)\n', dw_poly_phi_fcn_str,'dW_dphi = ',w_eval_str,'\nend');
dW_dtheta_str = strcat('function dW_dtheta = dW_dtheta(x)\n', dw_poly_theta_fcn_str,'dW_dtheta = ',w_eval_str,'\nend');

fid = fopen('W_fcn1.m','w');
fprintf(fid,w_fcn_str);
fclose(fid);

fid = fopen('dW_dT.m','w');
fprintf(fid,dW_dT_str);
fclose(fid);

fid = fopen('dW_dphi.m','w');
fprintf(fid,dW_dphi_str);
fclose(fid);

fid = fopen('dW_dtheta.m','w');
fprintf(fid,dW_dtheta_str);
fclose(fid);

% create the function handle
eval(W_exec);
W_fcn = @(x) W_eval(w_poly_fcn(x));

dW_dT = @(x) W_eval(dw_poly_T(x));
dW_dphi = @(x) W_eval(dw_poly_phi(x));
dW_dtheta = @(x) W_eval(dw_poly_theta(x));

dW_dxi_fcn = @(i,x) (i==7)*dW_dT(x)+(i==8)*dW_dphi(x)+(i==9)*dW_dtheta(x);

if controller.type == CtrlDesignOpts.ccm
    dW_dt_fcn = @(x,u) W_eval(dw_poly_dt_fcn([x;u])); 
elseif controller.type == CtrlDesignOpts.rccm
    dW_dt_fcn = @(x,u) W_eval(dw_poly_dt_fcn([x;u]));   
end


% ------------------------------ extract Y_fcn-----------------------------
if controller.type == CtrlDesignOpts.rccm 
    % write Y_fcn to m files
    Y_exec = 'Y_eval = @(ml)';
    y_eval_str = '';
    for i = 1:length(w_poly)
        if i<length(w_poly)
            Y_exec = strcat(Y_exec,sprintf('Y_coef(:,:,%d)*ml(%d) +',i,i));
            if i == 1
                y_eval_str = strcat(y_eval_str,mat2str(Y_coef(:,:,i),5),sprintf('*ml(%d) +',i),'...\n');
            else
                y_eval_str = [y_eval_str,'      ',mat2str(Y_coef(:,:,i),5),sprintf('*ml(%d) +',i),'...\n'];
            end
        else
            Y_exec = strcat(Y_exec,sprintf('Y_coef(:,:,%d)*ml(%d);',i,i));
            y_eval_str = [y_eval_str,'      ',mat2str(Y_coef(:,:,i)),sprintf('*ml(%d);',i)];
        end
    end 
    Y_exec = 'Y_eval = @(ml)';
    y_eval_str = '';
    for i = 1:length(w_poly)
        if i<length(w_poly)
            Y_exec = strcat(Y_exec,sprintf('Y_coef(:,:,%d)*ml(%d) +',i,i));
            if i == 1
                y_eval_str = strcat(y_eval_str,mat2str(Y_coef(:,:,i),5),sprintf('*ml(%d) +',i),'...\n');
            else
                y_eval_str = [y_eval_str,'      ',mat2str(Y_coef(:,:,i),5),sprintf('*ml(%d) +',i),'...\n'];
            end
        else
            Y_exec = strcat(Y_exec,sprintf('Y_coef(:,:,%d)*ml(%d);',i,i));
            y_eval_str = [y_eval_str,'      ',mat2str(Y_coef(:,:,i)),sprintf('*ml(%d);',i)];
        end
    end
    y_fcn_str = strcat('function Y = Y_fcn1(x)\n', w_poly_fcn_str,'Y = ',y_eval_str,'\nend');
    fid = fopen('Y_fcn1.m','w');
    fprintf(fid,y_fcn_str);
    fclose(fid);
    
    % create the function handle 'Y_fcn'
    eval(Y_exec);
    Y_fcn = @(x) Y_eval(w_poly_fcn(x));    
end




