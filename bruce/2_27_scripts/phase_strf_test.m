for cc = 1:24;
    
%     phase_sig = squeeze(full_phasegrams(tr_inds,15,cc));
%     phase_sig = squeeze(full_alphaphase(tr_inds,cc));
    phase_sig = squeeze(full_deltaphase(tr_inds,cc));
%     phase_sig = squeeze(full_gammaphase(tr_inds,cc));
    n_phase_bins = 8;
    phase_bin_edges = linspace(-pi,pi,n_phase_bins + 1);
    phase_bin_edges(end) = [];
    PBx = phase_bin_edges;
    
    Robs = full_spkbinned(tr_inds,cc);
    
    mod_signs = [1]; %determines whether input is exc or sup (doesn't matter in the linear case)
    NL_types = {'lin'}; %define subunit as linear
    silent = 0; %display model optimization iterations
    
    n_bars = stim_params.stim_dims(2);
    kLen = nLags*n_bars;
    full_kLen = kLen*n_phase_bins;
    filtLen = nLags*n_bars;
    % K0 = rand(full_kLen,1)/full_kLen*0.1;
    % K0(end+1) = 0;
    
    %create L2 mat
    mix_prop = [1 0.2 5e-2];
    et = ones(nLags,1)*sqrt(mix_prop(1));
    ex = ones(n_bars,1)*sqrt(mix_prop(2));
    ep = ones(n_phase_bins,1)*sqrt(mix_prop(3));
    % et([1 end]) = 0;
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
    dp = unique(diff(PBx));
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
    
    stim_params.lin_dims = nPBs;
    
    %%
    clear K0
    K0(1:filtLen*nPBs) = repmat(fit0(cc).mods(1).filtK,n_phase_bins,1);
    K0((filtLen*nPBs+1):(filtLen*nPBs+stim_params.lin_dims)) = fit0(cc).spk_NL_params(1);
    % K0(end+1) = fit0(cc).spk_NL_params(1);
    K0 = K0';
    reg_lambda = 100;
    
    %if using Mark Schmidt's optimization, some differences in option parameters
    optim_params.optTol = 1e-4;
    optim_params.progTol = 1e-6;
    optim_params.Method = 'lbfgs';
    optim_params.Display = 'iter';
    [params] = minFunc( @(K) phasedep_glmfit(K, Robs, full_Xmat(tr_inds,:),XPBx,L2_mat,reg_lambda,stim_params,PBx), K0, optim_params);
    
    kfit(cc,:) = params(1:filtLen*nPBs);
    pmod(cc,:) = params((filtLen*nPBs+1):end);
end
%%
% kfit = params(1:filtLen*nPBs);
% pmod = params((filtLen*nPBs+1):end);

for cc = 1:24
karray = reshape(kfit(cc,:),[nLags n_bars nPBs]);
cm = max(abs(kfit(cc,:)));
cax = [-cm cm]*1;
for i = 1:8
    subplot(4,2,i)
    imagesc(squeeze(karray(:,:,i)));
    caxis(cax);
end
pause
clf
end