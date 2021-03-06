function fit1 = fitNLw_alt(fit1,kern_output,spkbs,fun_tol)

if nargin < 4
    fun_tol = 1e-4;
end

min_iter = 2; max_iter = 5;

LPlast = Inf;
iter = 0;


%extract matrix of STCcfs across modules
STCcf_mat = get_STCcf_mat(fit1);

%internal filter outputs of each module
g_mat = kern_output*STCcf_mat;

while (((LPlast-fit1.LP) >= fun_tol) && (iter < max_iter) && (iter < min_iter))
    LPlast = fit1.LP;  iter = iter + 1;
    
    fit1 = fitNL_nopsc(fit1,g_mat,spkbs,0);
    
    fit1 = fitWeights_stcb(fit1,g_mat,spkbs,0);
    
    [nll, pnll] = getLLGLM_STCBF(fit1,kern_output,spkbs,'none');
    fit1.LL = nll; fit1.LP = pnll;       
end