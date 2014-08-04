function [LL,prate] = fit_phasedep_glmeval(init_mod,fit_params,Xstim,XLin,Robs,phase_sig,n_phase_bins)

    phase_bin_edges = linspace(-pi,pi,n_phase_bins + 1);
    phase_bin_edges(end) = [];
    PBx = phase_bin_edges;
    
    stim_params = init_mod.stim_params;
    
    n_bars = stim_params.stim_dims(2);
    nLags = stim_params.stim_dims(1);
    kLen = nLags*n_bars;
    filtLen = nLags*n_bars;
        
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
   
    %% Compute predicted firing rate
    kmat = reshape(fit_params(1:filtLen*nPBs),filtLen,nPBs);
    KPLin = fit_params((filtLen*nPBs+1):(filtLen*nPBs + nPBs));
    KLin = fit_params((filtLen*nPBs + nPBs+1):(filtLen*nPBs + nPBs+stim_params.lin_dims));
    filt_outs = Xstim*kmat;
    
    G = XPBx*KPLin' + XLin*KLin';%initialize overall generating function G
    
    G = G + sum(filt_outs.*XPBx,2);
    max_gbeta = 50; %to prevent numerical overflow
    expg = exp(G);
    too_large = (G > max_gbeta);
    prate = log(1+expg); %alpha*log(1+exp(gbeta))
    prate(too_large) = G(too_large); %log(1+exp(x)) ~ x in limit of large x
    
    %enforce minimum predicted firing rate to avoid nan LLs
    min_pred_rate = 1e-50;
    if min(prate) < min_pred_rate
        prate(prate < min_pred_rate) = min_pred_rate; %minimum predicted rate
    end
    
    %% COMPUTE LL and LL gradient
    LL = sum(Robs.* log(prate) - prate); %up to an overall constant
    LL = LL/sum(Robs);