function fit1 = fitNLw_alt_lexp(fit1,WX,spkbs,silent)

if nargin < 4
    silent = 1;
end

thresh = 1e-4;  
min_iter = 2; max_iter = 5;

LPlast = Inf;
iter = 0;


%extract matrix of STCcfs across modules
k_mat = get_k_mat(fit1);

%internal filter outputs of each module
g_mat = WX*k_mat;

fit1 = fitWeights_lexp(fit1,g_mat,spkbs,silent);

while (((LPlast-fit1.LP) >= thresh) && (iter < max_iter) && (iter < min_iter))
    LPlast = fit1.LP;  iter = iter + 1;
    
    fit1 = fitNL_lexp(fit1,g_mat,spkbs,silent);
    
    fit1 = fitWeights_lexp(fit1,g_mat,spkbs,silent);
    
    [nll, pnll] = getLLGLM_lexp(fit1,WX,spkbs,'none');
    fit1.LL = nll; fit1.LP = pnll;       
end