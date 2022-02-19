function opt = spot_sdp_mosek_options()
evalc('[r,res]=mosekopt(''param'');');
mosek = res.param;
% mosek.MSK_IPAR_NUM_THREADS = 1;
% mosek.MSK_IPAR_INTPNT_MULTI_THREAD = 0;
solver_options.mosek = mosek;
opt = struct('verbose',0,...
             'dualize',0,...
             'solver_options',solver_options);
end