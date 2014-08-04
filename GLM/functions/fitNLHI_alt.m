function fit1 = fitNLHI_alt(fit1,kern_output,spkbs)

thresh = 1e-4;  
min_iter = 2; max_iter = 5;

LPlast = Inf;
iter = 0;


%extract matrix of STCcfs across modules
STCcf_mat = get_STCcf_mat(fit1);

%internal filter outputs of each module
g_mat = kern_output*STCcf_mat;

while (((LPlast-fit1.LP) >= thresh) && (iter < max_iter) && (iter < min_iter))
    LPlast = fit1.LP;  iter = iter + 1;
    
%     fprintf('const: %.5f\n',fit1.const);
    fit1 = fitNL_jmm(fit1,g_mat,spkbs,1);
%     [nll, pnll] = getLLGLM_STCBF(fit1,kern_output,spkbs,'none');
%     fit1.LL = nll; fit1.LP = pnll;
    
%     fprintf( '  Iter %2d 1: LL = %f\n', iter, fit1.LP)
    fit1 = fitHI_jmm(fit1,g_mat,spkbs,1);
    [nll, pnll] = getLLGLM_STCBF(fit1,kern_output,spkbs,'none');
    fit1.LL = nll; fit1.LP = pnll;   
    
%     fprintf( '  Iter %2d 2: LL = %f\n', iter, fit1.LP)
end