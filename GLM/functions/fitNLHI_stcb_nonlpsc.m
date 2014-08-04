function [myglm] = fitNLHI_stcb_nonlpsc(myglm,ustim,spkbs,disp_type,max_iter,min_iter)
% USAGE: [myglm ll0] = fitNLHI(myglm,ustim,spkbs,lltol)
%   iteratively fit nonlinearity and history term

addpath('~/Dan_Matlab/')

if nargin < 5
    max_iter = 5;
end
if nargin < 6
    min_iter = 1;
end
w_TolFun = 1e-5;

stcfilts = myglm.STCbasis;

[stimlen,SDIM]   = size(ustim); %stimulus dimension

%precompute the stimulus filtered by each STC kernel
[klen,Nstcbf] = size(stcfilts);
flen = klen/SDIM;
if size(ustim,2) > SDIM %if stimulus is already time-embedded
    kern_output = ustim*stcfilts;
    Robs = zeros(1,stimlen);
    ftable = tabulate(spkbs);
    Robs(ftable(:,1)) = ftable(:,2);
else %if not
    kern_output = zeros(stimlen-flen+1,Nstcbf); %initialize filtered output to be (NT-N_tau+1)xNmods.
    for ikern = 1:Nstcbf;
        %for the given NL module, store the output of the stimulus filtered by
        %the internal kernel
        kern_output(:,ikern)= kfilterInput(ustim,stcfilts(:,ikern));
    end
    
    %compute binned spike vector
    foff = flen + length(myglm.mods(1).h) - 2;
    Robs = zeros(1,stimlen - foff);
    ftable = tabulate(spkbs);
    Robs(ftable(:,1)) = ftable(:,2);
end


myglm.LP_seq = [];
myglm.LP_ax = [];
myglm.LL_seq = [];

fprintf('INITIAL MODEL\n');
[ll0, ll0p] = getLLGLM_STCBF_nonlpsc(myglm,kern_output,spkbs,disp_type);
myglm.LP = ll0p; myglm.LL = ll0;
myglm.LP_seq = [myglm.LP_seq ll0p]; 
myglm.LP_ax = [1];
myglm.LL_seq = [myglm.LL_seq ll0];

%normalize kernel coefs to map onto consistent range of nonlinearity
% myglm = normalizeRFs_STCB(myglm,kern_output);

thresh = 1e-3;
LLlast = Inf;
iter = 0;
while ((((LLlast-myglm.LL) >= thresh) && (iter < max_iter))) || (iter < min_iter)
    lastglm = myglm;
    LLlast = lastglm.LL;  iter = iter + 1;
      
    %extract matrix of STCcfs across modules
    STCcf_mat = get_STCcf_mat(myglm);   
    %internal filter outputs of each module
    g_mat = kern_output*STCcf_mat;       

    fprintf('Iteration %d: Fitting Weights\n',iter);
    myglm = fitWeights_stcb_nonlpsc(myglm,g_mat,spkbs,1,w_TolFun);    
    [ll0 ll0p] = getLLGLM_STCBF_nonlpsc(myglm,kern_output,spkbs,disp_type);
    myglm.LP_seq = [myglm.LP_seq ll0p];
    myglm.LL_seq = [myglm.LL_seq ll0];
    myglm.LP_ax = [myglm.LP_ax myglm.LP_ax(end)+1];

    fprintf('Iteration %d: Optimizing RF\n',iter);
    myglm = fitSTCBF_nonlpsc(myglm,kern_output,Robs,1);
    [ll0, ll0p] = getLLGLM_STCBF_nonlpsc(myglm,kern_output,spkbs,disp_type);
    %     [ll0, ll0p] = getLLGLM_STCBF(myglm,kern_output,spkbs,disp_type);
    myglm.LP_seq = [myglm.LP_seq ll0p];
    myglm.LL_seq = [myglm.LL_seq ll0];
    myglm.LP_ax = [myglm.LP_ax myglm.LP_ax(end)+1];
   
    myglm.LL = ll0; myglm.LP = ll0p;
end

flen = length(myglm.mods(1).h);
for imod =1:length(myglm.mods);
    cur_inputs = kern_output*myglm.mods(imod).STCcf;
    cur_inputs(1:flen-1) = [];
    [foy,fox] = ksdensity(cur_inputs);
    myglm.mods(imod).fox = fox;
    myglm.mods(imod).foy = foy;
    
    [foys,foxs] = ksdensity(cur_inputs(spkbs));
    myglm.mods(imod).foxs = foxs;
    myglm.mods(imod).foys = foys;
end

% fprintf('\n\n\Model Fitting Complete! \n');
getLLGLM_STCBF_nonlpsc(myglm,kern_output,spkbs,disp_type);

%%
% drawnow;
%
% lldiff = ll0-llnow;
% if(lldiff > 0);
%     ll0   = llnow;
%     myglm = nextWfit;
% end;
% myglm.LL = llnow;
% hilen = length(myglm.mods(1).h);
% myglm = addFOdens(myglm,ustim,spkbs+hilen-1);

