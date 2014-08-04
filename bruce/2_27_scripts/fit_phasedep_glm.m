function [phase_filts,phase_offset,phase_kLin] = fit_phasedep_glm(init_mod,Xstim,XLin,Robs,phase_sig,n_phase_bins,lambdas,silent,kern_phasedep,L1)
% [phase_filts,phase_offset,phase_kLin] = fit_phasedep_glm(init_mod,Xstim,XLin,Robs,phase_sig,n_phase_bins,lambdas,silent,kern_phasedep,L1)
% INPUTS: 
%     init_mod: initial NIM fit
%     Xstim: stimulus matrix
%     XLin: Any additional linear predictors (no regularization)
%     Robs: binned spike counts
%     phase_sig: time series of phases (in radians)
%     n_phase_bins: number of phase bins
%     lambdas: regularization lambdas. 1x3 vector [time space phase]
%     silent: 1 to silence output
%     kern_phasedep: 1 to allow stimulus filter to be phase-dependent (default 1)
%     L1: L1 regularization (default 0)
% OUTPUTS:
%     phase_filts: n_phase_bins*K vector of filters for each phase bin
%     phase_offset: offset as a function of phase
%     phase_kLin: extra linear predictor

if nargin < 9
    kern_phasedep = 1;
end
if nargin < 10
    L1 = 0;
end
%%


    phase_bin_edges = linspace(-pi,pi,n_phase_bins + 1);
    phase_bin_edges(end) = [];
    PBx = phase_bin_edges;
    
    stim_params = init_mod.stim_params;
    
    n_bars = stim_params.stim_dims(2);
    nLags = stim_params.stim_dims(1);
    kLen = nLags*n_bars;
    filtLen = nLags*n_bars;
    
    %%
    %create L2 mat
    lambda = max(lambdas);
    mix_prop = lambdas/lambda;
    et = ones(nLags,1)*sqrt(mix_prop(1));
    ex = ones(n_bars,1)*sqrt(mix_prop(2));
    ep = ones(n_phase_bins,1)*sqrt(mix_prop(3));
    % et([1 end]) = 0;
    ex([1 end]) = 0;
    D1t = spdiags([et -2*et et], [-1 0 1], nLags, nLags)';
    D1x = spdiags([ex -2*ex ex], [-1 0 1], n_bars, n_bars)';
    D1p = spdiags([ep -2*ep ep], [-1 0 1], n_phase_bins, n_phase_bins)';
    %periodic bcs
    D1p(end,1) = sqrt(mix_prop(3)); D1p(1,end) = sqrt(mix_prop(3));
    
    It = speye(nLags); Ix = speye(n_bars); Ip = speye(n_phase_bins);
    L2_mat = kron(Ip, kron(Ix, D1t)) + kron(Ip, kron(D1x, It))...
        + kron(kron(D1p,Ix),It);
    
    %% CREATE TENT BASIS FUNCTION OUTPUTs
    nPBs = length(PBx);
    XPBx = zeros(length(phase_sig),nPBs);
    dp = median(diff(PBx));
    for j = 2:nPBs-1
        cur_set = find(phase_sig > PBx(j-1) & phase_sig < PBx(j)); %left side
        XPBx(cur_set,j) = XPBx(cur_set,j) + (1 - (PBx(j) - phase_sig(cur_set))/dp);
        cur_set = find(phase_sig >= PBx(j) & phase_sig < PBx(j+1)); %right side
        XPBx(cur_set,j) = XPBx(cur_set,j) + (1 - (phase_sig(cur_set)-PBx(j))/dp);
    end
    %first TB
    cur_set = find(phase_sig >= PBx(1) & phase_sig < PBx(2)); %right side
    XPBx(cur_set,1) = XPBx(cur_set,1) + (1 - (phase_sig(cur_set)-PBx(1))/dp);
    cur_set = find(phase_sig >= PBx(end)); %left side
    XPBx(cur_set,1) = XPBx(cur_set,1) + (1 - (pi-phase_sig(cur_set))/dp);
    %last TB
    cur_set = find(phase_sig > PBx(end-1) & phase_sig < PBx(end)); %left side
    XPBx(cur_set,end) = XPBx(cur_set,end) + (1 - (PBx(end)-phase_sig(cur_set))/dp);
    cur_set = find(phase_sig >= PBx(end)); %right side
    XPBx(cur_set,end) = XPBx(cur_set,end) + (1 - (phase_sig(cur_set)-PBx(end))/dp);
        
    %%
    clear K0
    K0(1:filtLen*nPBs) = repmat(init_mod.mods(1).filtK,n_phase_bins,1);
    K0((filtLen*nPBs+1):(filtLen*nPBs+nPBs)) = init_mod.spk_NL_params(1);
    K0((filtLen*nPBs+nPBs+1):(filtLen*nPBs+nPBs+stim_params.lin_dims)) = init_mod.kLin;
%     K0((filtLen*nPBs+1):(filtLen*nPBs+stim_params.lin_dims)) = init_mod.kLin;
%     K0(end+1) = init_mod.spk_NL_params(1);
%     XLin(:,end+1) = 1;
    K0 = K0';
    
    lambda_L1 = zeros(size(K0));
    lambda_L1(1:filtLen*nPBs) = L1;
    lambda_L1 = lambda_L1/sum(Robs); %since we are dealing with LL/spk
    
    %if using Mark Schmidt's optimization, some differences in option parameters
    optim_params.optTol = 1e-4;
    optim_params.progTol = 1e-6;
    optim_params.Method = 'lbfgs';
    if silent == 0
        optim_params.Display = 'iter';
    else
        optim_params.Display = 'off';
    end
    if max(lambda_L1) > 0
        [params] = L1General2_PSSas( @(K) phasedep_glmfit(K, Robs, Xstim,XPBx,XLin,L2_mat,lambda,stim_params,PBx), K0, lambda_L1,optim_params);
        phase_filts = params(1:filtLen*nPBs);
        phase_offset = params((filtLen*nPBs+1):(filtLen*nPBs+nPBs));
        phase_kLin = params((filtLen*nPBs+nPBs+1):end);
        
    else
    if kern_phasedep == 1
        [params] = minFunc( @(K) phasedep_glmfit(K, Robs, Xstim,XPBx,XLin,L2_mat,lambda,stim_params,PBx), K0, optim_params);
        phase_filts = params(1:filtLen*nPBs);
        phase_offset = params((filtLen*nPBs+1):(filtLen*nPBs+nPBs));
% phase_offset = [];
        phase_kLin = params((filtLen*nPBs+nPBs+1):end);
%         phase_kLin = params((filtLen*nPBs+1):end);
    else
        fixK = K0(1:filtLen*nPBs);
        K0 = K0(filtLen*nPBs+1:end);
        [params] = minFunc( @(K) phasedep_glmfit_offset(K, fixK, Robs, Xstim,XPBx,XLin,L2_mat,lambda,stim_params,PBx), K0, optim_params);
        phase_filts = fixK;
        phase_offset = params((1):(nPBs));
        phase_kLin = params((nPBs+1):end);
    end
    end
    
   