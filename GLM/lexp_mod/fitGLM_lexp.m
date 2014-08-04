function [myglm] = fitGLM_lexp(myglm,X,spkbs,disp_type,its_per_round,optTol,progTol,sa_params)

%INPUTS:
% myglm: model structure
% X: X-matrix
% spkbs: vector of spike bin indices
% disp_type: display type: 'tots','full','min','none'
% its_per_round: number of iterations to fit model before checking for user input to continue

if nargin < 6
    optTol = 1e-4;
end
if nargin < 7
    progTol = 1e-6;
end
if nargin < 8
    sa_params = [];
end

nmods = length(myglm.mods);

if nargin < 5 | isempty(its_per_round)
    its_per_round = 500;
    check_user_in = 0;
else
    check_user_in = 1;
end

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

%process localization penalties.  If > 0, need to compute filter COMs
locLambdas = arrayfun(@(x) x.locLambda,myglm.mods);
if max(locLambdas) > 0
    if strcmp(myglm.image_type,'1d')
        [myglm,COMs] = get_filter_peaks_1d(myglm);
%         [myglm,COMs] = get_filter_coms_1d(myglm);
    elseif strcmp(myglm.image_type,'2d')
        [myglm,COMs] = get_filter_coms_2d(myglm);
    else
        error('Invalid image type')
    end
end

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

% [ll0, ll0p] = getLLGLM_lexp(myglm,X,spkbs,disp_type);
% myglm.LP_seq = [myglm.LP_seq ll0p];
% myglm.LL_seq = [myglm.LL_seq ll0];
% myglm.LP_ax = [myglm.LP_ax myglm.LP_ax(end)+1];

eflag = 0;
n_opt_rounds = 0;
can_exit = 0;
while eflag == 0
    n_opt_rounds = n_opt_rounds + 1;
    fprintf('Optimization round %d\n',n_opt_rounds);
    [myglm eflag] = fitFULLRF_lexp(myglm,X,Robs,its_per_round,'lbgs',0,optTol,progTol,sa_params);
    %     [myglm eflag] = fitFULLRF_lexp(myglm,X,Robs,its_per_round,'qnewton',0,optTol,progTol);
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

disp('Model Fitting Complete');
getLLGLM_lexp(myglm,X,spkbs,disp_type);


