function [myglm] = fitGLM_lexp_rot(myglm,kern_output,spkbs,disp_type,its_per_round)

optTol = 1e-4;
progTol = 1e-6;

nmods = length(myglm.mods);

if nargin < 5
    its_per_round = 200;
    check_user_in = 0;
else
    check_user_in = 1;
end

lambdaW = myglm.lambdaW;

[stimlen,nbvs]   = size(kern_output); %stimulus dimension
% if strcmp(myglm.basis,'pix')
%     pixlen = klen;
% elseif strcmp(myglm.basis,'white')
%     pixlen = size(myglm.pix_conv_mat,2);
% end

% fsdim = myglm.mods(1).fsdim;
% flen = pixlen/fsdim;
% if strcmp(myglm.image_type,'2d');
%     sdim = sqrt(fsdim);
% elseif strcmp(myglm.image_type,'1d')
%     sdim = myglm.mods(1).SDIM;
% end

%initialize likelihood vecs
myglm.LP_seq = [];
myglm.LP_ax = []; %iteration number
myglm.LL_seq = [];

% myglm = normalizeRFs_STCB(myglm,kern_output);

fprintf('INITIAL MODEL\n');
[ll0, ll0p] = getLLGLM_lexp_rot(myglm,kern_output,spkbs,disp_type);
myglm.LP = ll0p; myglm.LL = ll0;
myglm.LP_seq = [myglm.LP_seq ll0p];
myglm.LP_ax = [1];
myglm.LL_seq = [myglm.LL_seq ll0];

%compute binned spike vector
Robs = zeros(1,stimlen);
ftable = tabulate(spkbs);
Robs(ftable(:,1)) = ftable(:,2);

%extract matrix of internal filter coefs
STCcf_mat = get_STCcf_mat(myglm);
%internal filter outputs of each module
g_mat = kern_output*STCcf_mat;

%fit model weights
% myglm = fitWeights_lexp(myglm,g_mat,spkbs,0);
% [myglm eflag] = fitbetas_lexp_rot(myglm,g_mat,Robs,200,'lbgs',0);edit 

% [ll0, ll0p] = getLLGLM_lexp_rot(myglm,kern_output,spkbs,disp_type);
% myglm.LP_seq = [myglm.LP_seq ll0p];
% myglm.LL_seq = [myglm.LL_seq ll0];
% myglm.LP_ax = [myglm.LP_ax myglm.LP_ax(end)+1];
% 
% %do an initial set of iterations using Conj Gradient
% [myglm eflag] = fitFULLRF_lexp_rot(myglm,kern_output,Robs,50,'pcg',0);

% myglm = normalizeRFs_STCB(myglm,kern_output);
% myglm = fitWeights_lexp(myglm,g_mat,spkbs,0);
% [myglm eflag] = fitbetas_lexp_rot(myglm,g_mat,Robs,200,'lbgs',0);

%until we reach convergence keep trying LBGS optimizations
eflag = 0;
n_opt_rounds = 0;
can_exit = 0;
while eflag == 0
    n_opt_rounds = n_opt_rounds + 1;
    fprintf('Optimization round %d\n',n_opt_rounds);
    [myglm eflag] = fitFULLRF_lexp_rot(myglm,kern_output,Robs,its_per_round,'lbgs',0,optTol,progTol);
    if check_user_in == 1
        R = input('Continue? (y/n)','s');
        if R == 'n'
            break
        end
    else
        break
    end
    eflag = 0;
end

%% estimate hessian at current point
% nmods = length(myglm.mods);
% cur_params = [];
% for m = 1:nmods
%     cur_stc_kern = myglm.mods(m).STCcf';
%     cur_params = [cur_params cur_stc_kern]; %add STCBcoefs to initial param vector
% end
% cur_params(end+1) = myglm.const; %add constant offset term to params
% 
% % hess = 
% % [f,df] = rot_lexp_LLinternal(cur_params', Robs, kern_output, myglm)
% 


%%


[ll0, ll0p] = getLLGLM_lexp_rot(myglm,kern_output,spkbs,disp_type);
myglm.LP_seq(end) = ll0p;
myglm.LL_seq = [myglm.LL_seq ll0];

myglm.LL = ll0; myglm.LP = ll0p;

flen = length(myglm.mods(1).h);
% for imod =1:length(myglm.mods);
%     if max(abs(myglm.mods(imod).k)) > 0
%         cur_inputs = kern_output*myglm.mods(imod).STCcf;
%         [foy,fox] = ksdensity(cur_inputs);
%         myglm.mods(imod).fox = fox;
%         myglm.mods(imod).foy = foy;
%         
%         [foys,foxs] = ksdensity(cur_inputs(spkbs));
%         myglm.mods(imod).foxs = foxs;
%         myglm.mods(imod).foys = foys;
%     else
%         myglm.mods(imod).fox = linspace(-3.1,3.1,10);
%         myglm.mods(imod).foy = zeros(1,10);
%         myglm.mods(imod).foxs = myglm.mods(imod).fox;
%         myglm.mods(imod).foys = myglm.mods(imod).foy;
%     end
% end

printf('\n\n\Model Fitting Complete! \n');
% getLLGLM_lexp_rot(myglm,kern_output,spkbs,disp_type);


