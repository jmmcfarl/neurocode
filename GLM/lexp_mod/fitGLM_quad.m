function [myglm] = fitGLM_quad(myglm,X,spkbs,disp_type)



nmods = length(myglm.mods);
lambdaW = myglm.lambdaW;

[stimlen,klen]   = size(X); %stimulus dimension
if strcmp(myglm.basis,'pix')
    pixlen = klen;
elseif strcmp(myglm.basis,'white')
    pixlen = size(myglm.pix_conv_mat,2);
end

fsdim = myglm.mods(1).fsdim;
flen = pixlen/fsdim;
if strcmp(myglm.image_type,'2d');
    sdim = sqrt(fsdim);
elseif strcmp(myglm.image_type,'1d')
    sdim = myglm.mods(1).SDIM;
end

%initialize likelihood vecs
myglm.LP_seq = [];
myglm.LP_ax = []; %iteration number
myglm.LL_seq = [];

% %process localization penalties.  If > 0, need to compute filter COMs
% locLambdas = arrayfun(@(x) x.locLambda,myglm.mods);
% if max(locLambdas) > 0
%     if strcmp(myglm.image_type,'1d')
%         [myglm,COMs] = get_filter_coms_1d(myglm);
%     elseif strcmp(myglm.image_type,'2d')
%         [myglm,COMs] = get_filter_coms_2d(myglm);
%     else
%         error('Invalid image type')
%     end
% end

fprintf('INITIAL MODEL\n');
[ll0, ll0p] = getLLGLM_lexp(myglm,X,spkbs,disp_type);
myglm.LP = ll0p; myglm.LL = ll0;
myglm.LP_seq = [myglm.LP_seq ll0p];
myglm.LP_ax = [1];
myglm.LL_seq = [myglm.LL_seq ll0];

%compute binned spike vector
Robs = zeros(1,stimlen);
ftable = tabulate(spkbs);
Robs(ftable(:,1)) = ftable(:,2);

%extract matrix of internal filter coefs
k_mat = get_k_mat(myglm);
%internal filter outputs of each module
g_mat = X*k_mat;

%fit model weights
% myglm = fitWeights_lexp(myglm,g_mat,spkbs,1);

[ll0, ll0p] = getLLGLM_lexp(myglm,X,spkbs,disp_type);
myglm.LP_seq = [myglm.LP_seq ll0p];
myglm.LL_seq = [myglm.LL_seq ll0];
myglm.LP_ax = [myglm.LP_ax myglm.LP_ax(end)+1];

%do an initial set of iterations using Conj Gradient
% [myglm eflag] = fitFULLRF_lexp(myglm,X,Robs,50,'pcg',0);
eflag = 0;

%until we reach convergence keep trying LBGS optimizations
n_opt_rounds = 1;
can_exit = 0;
% while can_exit == 0
    while eflag == 0
        n_opt_rounds = n_opt_rounds + 1;
        fprintf('Optimization round %d\n',n_opt_rounds);
        [myglm eflag] = fitFULLRF_lexp(myglm,X,Robs,50,'lbgs',0,1e-4,1e-6);
    end
%     %check to see if we really art at a minimum
%     [myglm eflag] = fitFULLRF_lexp(myglm,X,Robs,3,'lbgs',0,1e-5,1e-8); 
%     if eflag ~= 0
%         can_exit = 1;
%     end    
% end

%%




%if any subunits have all coefs to 0, move the 0 to the subunit weight
for imod = 1:nmods
    if all(myglm.mods(imod).k == 0)
        myglm.mods(imod).w = 0;
    end
end

[ll0, ll0p] = getLLGLM_lexp(myglm,X,spkbs,disp_type);
myglm.LP_seq(end) = ll0p;
myglm.LL_seq = [myglm.LL_seq ll0];

% %if using localization recompute COMs
% if max(locLambdas) > 0
%     if strcmp(myglm.image_type,'1d')
%         [myglm,COMs] = get_filter_coms_1d(myglm);
%     elseif strcmp(myglm.image_type,'2d')
%         
%     else
%         error('Invalid image type')
%     end
% end

myglm.LL = ll0; myglm.LP = ll0p;

flen = length(myglm.mods(1).h);
for imod =1:length(myglm.mods);
    if max(abs(myglm.mods(imod).k)) > 0
        cur_inputs = X*myglm.mods(imod).k;
        [foy,fox] = ksdensity(cur_inputs);
        myglm.mods(imod).fox = fox;
        myglm.mods(imod).foy = foy;
        
        [foys,foxs] = ksdensity(cur_inputs(spkbs));
        myglm.mods(imod).foxs = foxs;
        myglm.mods(imod).foys = foys;
    else
        myglm.mods(imod).fox = linspace(-3.1,3.1,10);
        myglm.mods(imod).foy = zeros(1,10);
        myglm.mods(imod).foxs = myglm.mods(imod).fox;
        myglm.mods(imod).foys = myglm.mods(imod).foy;
    end
end

printf('\n\n\Model Fitting Complete! \n');
getLLGLM_lexp(myglm,X,spkbs,disp_type);


