function [myglm] = fitNLHI_stcb_nopsc(myglm,ustim,spkbs,disp_type,max_iter,min_iter)
% USAGE: [myglm ll0] = fitNLHI(myglm,ustim,spkbs,lltol)
%   iteratively fit nonlinearity and history term

addpath('~/Dan_Matlab/')

if nargin < 5
    max_iter = 5;
end

if nargin < 6
    min_iter = 0;
end

stcfilts = myglm.STCbasis;

[stimlen,SDIM]   = size(ustim); %stimulus dimension

%precompute the stimulus filtered by each STC kernel
[klen,Nstcbf] = size(stcfilts);
flen = klen/SDIM;
if size(ustim,2) > 1 %if stimulus is already time-embedded
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
[ll0, ll0p] = getLLGLM_STCBF(myglm,kern_output,spkbs,disp_type);
myglm.LP = ll0p; myglm.LL = ll0;
myglm.LP_seq = [myglm.LP_seq ll0p]; 
myglm.LP_ax = [1];
myglm.LL_seq = [myglm.LL_seq ll0];

%normalize kernel coefs to map onto consistent range of nonlinearity
myglm = normalizeRFs_STCB(myglm,kern_output);

thresh = 1e-3;
LPlast = Inf;
iter = 0;
while (((LPlast-myglm.LP) >= thresh) && (iter < max_iter)) || (iter < min_iter)
    lastglm = myglm;
    LPlast = lastglm.LP;  iter = iter + 1;
    
    fprintf('Iteration %d: Optimizing NL-HI\n',iter);
    myglm = fitNLw_alt(myglm,kern_output,spkbs);

    [ll0 ll0p] = getLLGLM_STCBF(myglm,kern_output,spkbs,disp_type);
    myglm.LP_seq = [myglm.LP_seq ll0p];
    myglm.LL_seq = [myglm.LL_seq ll0];
    myglm.LP_ax = [myglm.LP_ax myglm.LP_ax(end)+1];

    fprintf('Iteration %d: Optimizing RF\n',iter);
    myglm = fitSTCBF_nopsc(myglm,kern_output,Robs);
%     plotfo1d(myglm); 

    [ll0, ll0p] = getLLGLM_STCBF(myglm,kern_output,spkbs,disp_type);
    myglm.LL = ll0; myglm.LP = ll0p; 
    myglm.LP_seq = [myglm.LP_seq ll0p];
    myglm.LL_seq = [myglm.LL_seq ll0];
    myglm.LP_ax = [myglm.LP_ax myglm.LP_ax(end)+1];

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

fprintf('\n\nModel Fitting Complete! \n');
getLLGLM_STCBF(myglm,kern_output,spkbs,disp_type);

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

