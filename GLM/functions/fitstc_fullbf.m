function [myglm] = fitstc_fullbf(myglm,WX,spkbs,disp_type,max_iter)
% USAGE: [myglm ll0] = fitNLHI(myglm,ustim,spkbs,lltol)
%   iteratively fit nonlinearity and history term

if nargin < 5
    max_iter = Inf;
end

addpath('~/Dan_Matlab/')

[stimlen,klen]   = size(WX); %stimulus dimension
if strcmp(myglm.basis,'pix')
    pixlen = klen;
elseif strcmp(myglm.basis,'white')
    pixlen = size(myglm.pix_conv_mat,2);
end
fsdim = myglm.mods(1).fsdim;
flen = pixlen/fsdim;
if strcmp(myglm.image_type,'2d');
    sdim = sqrt(fsdim);
%     if sdim ~= round(sdim)
%         error('non-square array!')
%     end
elseif strcmp(myglm.image_type,'1d')
    sdim = myglm.mods(1).SDIM;
end

myglm.LP_seq = [];
myglm.LP_ax = [];
myglm.LL_seq = [];

locLambdas = arrayfun(@(x) x.locLambda,myglm.mods);
if max(locLambdas) > 0
    if strcmp(myglm.image_type,'1d')
        [myglm,COMs] = get_filter_coms_1d(myglm);      
    elseif strcmp(myglm.image_type,'2d')
         [myglm,COMs] = get_filter_coms_2d(myglm);             
    else
        error('Invalid image type')
    end
end

% %should we normalize here?
% myglm = normalizeRFs_full(myglm,WX);

fprintf('INITIAL MODEL\n');
[ll0, ll0p] = getLLGLM_FULL2d(myglm,WX,spkbs,disp_type);
myglm.LP = ll0p; myglm.LL = ll0;
myglm.LP_seq = [myglm.LP_seq ll0p];
myglm.LP_ax = [1];
myglm.LL_seq = [myglm.LL_seq ll0];

%compute binned spike vector
Robs = zeros(1,stimlen);
ftable = tabulate(spkbs);
Robs(ftable(:,1)) = ftable(:,2);

nmods = length(myglm.mods);

%determine if any fo the NL functions are unconstrained
for i = 1:nmods
    if strcmp(myglm.mods(i).nltype,'uncon')
        error('Cant fit NLs here')
    end
end

lambdaW = myglm.lambdaW;

%extract matrix of internal filter coefs
k_mat = get_k_mat(myglm);
%internal filter outputs of each module
g_mat = WX*k_mat;

myglm = fitWeights_full(myglm,g_mat,spkbs,1);

[ll0, ll0p] = getLLGLM_FULL2d(myglm,WX,spkbs,disp_type);
myglm.LP_seq = [myglm.LP_seq ll0p];
myglm.LL_seq = [myglm.LL_seq ll0];
myglm.LP_ax = [myglm.LP_ax myglm.LP_ax(end)+1];
    
thresh = 1e-3;
LLlast = Inf;
iter = 0;
% while ((LLlast-myglm.LL) >= thresh) && (iter < max_iter)
    lastglm = myglm;
    LLlast = lastglm.LL;  iter = iter + 1;
    
%     myglm = fitWeights_full(myglm,g_mat,spkbs,0);
%     [ll0, ll0p] = getLLGLM_FULL2d(myglm,WX,spkbs,disp_type);
%     myglm.LP_seq = [myglm.LP_seq ll0p];
%     myglm.LL_seq = [myglm.LL_seq ll0];
%     myglm.LP_ax = [myglm.LP_ax myglm.LP_ax(end)+1];
    
            
    fprintf('Iteration %d: Optimizing RF\n',iter);
    myglm = fitFULLRF2d_stc(myglm,WX,Robs,0);

    for imod = 1:nmods
        if all(myglm.mods(imod).k == 0)
            myglm.mods(imod).w = 0;
        end
    end
    
%     %not sure if we should normalize here??
%      myglm = normalizeRFs_full(myglm,WX);

     [ll0, ll0p] = getLLGLM_FULL2d(myglm,WX,spkbs,disp_type);
    myglm.LP_seq(end) = ll0p;
    myglm.LL_seq = [myglm.LL_seq ll0];
    
    if max(locLambdas) > 0
        if strcmp(myglm.image_type,'1d')
            [myglm,COMs] = get_filter_coms_1d(myglm);
        elseif strcmp(myglm.image_type,'2d')
            
        else
            error('Invalid image type')
        end
    end
    
    myglm.LL = ll0; myglm.LP = ll0p;
    
%     if isinf(max_iter)
%         to_cont = input('Continue (1=yes; 0=no)?');
%         if to_cont == 0
%             break
%         end
%     end
% end

flen = length(myglm.mods(1).h);
for imod =1:length(myglm.mods);
    if max(abs(myglm.mods(imod).k)) > 0
        cur_inputs = WX*myglm.mods(imod).k;
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

% fprintf('\n\n\Model Fitting Complete! \n');
getLLGLM_FULL2d(myglm,WX,spkbs,disp_type);


