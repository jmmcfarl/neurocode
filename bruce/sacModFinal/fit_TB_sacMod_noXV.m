function [all_TB_sacmod] = fit_TB_sacMod_noXV(Robs,stimG,t_since_sac,fit_inds,poss_d2T,TB_params,Xsac,basemod_predrate)

xbuff = TB_params.xbuff; %buffer for time binning to minimize edge artifacts
backlag = TB_params.backlag;
forlag = TB_params.forlag;
Xtick = -(backlag+xbuff+1/2):(1):(forlag+xbuff+1/2);
n_Gbins = TB_params.n_Gbins;
n_Tbins = length(Xtick); %number of time bins (with buffer)
G_lambda = TB_params.G_lambdas; %d2T reg on G dimension
silent = 1;
n_slags = size(Xsac,2);

%initialize 2D TB data using t-since-sac and normalized G
TB_stim = [t_since_sac stimG];

%set G-bins based on prctiles
%         %equispaced binning
Ytick = my_prctile(TB_stim(fit_inds,2),linspace(0.5,99.5,n_Gbins)); %equipopulated binning

%in some cases there are an excess of G==0 values, which makes
%some adjacent bins identical. Need to cut these bin edges.
nd_bins = find(diff(Ytick) <= 0);
Ytick(nd_bins) = [];
cur_nGbins = length(Ytick);

%initialize TBs
TB = TentBasis2D(Xtick, Ytick);

%this is data within range of the TB centers
TB_used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
    TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
% TB_used_data = TB_used_data(ismember(TB_used_data,fit_inds));

null_rate = mean(Robs(TB_used_data));
nullLL = nansum(Robs(TB_used_data).*log2(null_rate) - null_rate)/sum(Robs(TB_used_data));

%process data with TBs
[TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(TB_used_data,:));
for ll = 1:length(poss_d2T)
    fprintf('Fitting TB gsac model lambda %d/%d\n',ll,length(poss_d2T));
    cur_sac_lambda = poss_d2T(ll);
    
    %this gives fixed reg strength in G-direction, but
    %cur_sac_lambda in T-direction
    L2_params = create_L2_params([],[1 n_Tbins*cur_nGbins],[n_Tbins cur_nGbins],2,3,[Inf Inf],[1 1/cur_sac_lambda*G_lambda]);
    TB_sacmod = regGLM_fit(TB_Xmat,Robs(TB_used_data),L2_params,cur_sac_lambda,[],[],silent);
    [LL, penLL, TB_pred_rate] = regGLM_eval(TB_sacmod,Robs(TB_used_data),TB_Xmat);
    
    LL = nansum(Robs(TB_used_data).*log2(TB_pred_rate) - TB_pred_rate)/sum(Robs(TB_used_data));
    TB_sacmod.LLimp = (LL - nullLL);
    TB_sacmod.nullLL = nullLL;
    TB_sacmod.LL = LL;
    TB_sacmod.ovInfo = mean(TB_pred_rate/mean(TB_pred_rate).*log2(TB_pred_rate/mean(TB_pred_rate)));
    
    %compute output of TB model
    TB_K = reshape(TB_sacmod.K,n_Tbins,cur_nGbins)';
    bin_areas = TB.GetBinAreas();
    TB_dist = TB_counts./bin_areas;
    normFac = trapz(trapz(TB_dist));
    TB_dist = TB_dist'/normFac;
    TB_rate = log(1 + exp(TB_K + TB_sacmod.theta));
    
    %INFO CALS
    cur_avg_rate = mean(Robs(TB_used_data)); %overall avg rate
    marg_gdist = trapz(TB_dist,2); %marginal distribution of G
    marg_sdist = trapz(TB_dist);
    marg_gsacrate = trapz(TB_dist.*TB_rate)./marg_sdist;
    marg_grate = trapz(TB_dist.*TB_rate,2)./marg_gdist;
    gsacdep_info = nan(1,n_Tbins);
    for tt = 1:n_Tbins
        gsacdep_info(tt) = trapz(TB_dist(:,tt).*TB_rate(:,tt)/marg_gsacrate(tt).*log2(TB_rate(:,tt)/marg_gsacrate(tt)))/trapz(TB_dist(:,tt));
    end
    
    TB_sacmod.anInfo = trapz(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate))/cur_avg_rate;
    
    TB_sacmod.an_sac_rate = marg_gsacrate(xbuff+1:end-xbuff);
    TB_sacmod.an_sac_info = gsacdep_info(xbuff+1:end-xbuff);
    TB_sacmod.marg_gdist = marg_gdist;
    TB_sacmod.marg_grate = marg_grate;
    TB_sacmod.lagX = Xtick(xbuff+1:end-xbuff);
    TB_sacmod.Gtick = Ytick;
    TB_sacmod.TB_rate = TB_rate(:,xbuff+1:end-xbuff);
    
    %calculate marginal g distribution using equispaced binning
    dYtick = median(diff(Ytick));
    equi_space_gX = linspace(Ytick(1)-dYtick,Ytick(end)+dYtick,51);
    equi_space_gdist = histc(TB_stim(TB_used_data,2),equi_space_gX);
    TB_sacmod.equi_gdist = equi_space_gdist(1:end-1)/sum(equi_space_gdist);
    TB_sacmod.equi_gX = 0.5*equi_space_gX(1:end-1) + 0.5*equi_space_gX(2:end);
    
    %% %compute response gain and offsets
    [sac_offset,sac_gain] = deal(nan(length(TB_sacmod.lagX),1));
    for ii = 1:length(TB_sacmod.lagX)
        temp = lscov([ones(cur_nGbins,1) marg_grate],TB_sacmod.TB_rate(:,ii),TB_sacmod.marg_gdist);
        sac_offset(ii) = temp(1); sac_gain(ii) = temp(2);
    end
    TB_sacmod.sac_gain = sac_gain;
    TB_sacmod.sac_offset = sac_offset;
    
    
    [sac_info,sac_LLimp] = deal(nan(n_slags,1));
    for ss = 1:n_slags
        temp = find(Xsac(TB_used_data,ss) == 1);
        cur_Robs = Robs(TB_used_data(temp));
        
        sac_info(ss) = nanmean(TB_pred_rate(temp).*log2(TB_pred_rate(temp)/mean(TB_pred_rate(temp))))/mean(TB_pred_rate(temp));
        
        %compute nullLL for data at this latency
        cur_nullLL = nansum(cur_Robs.*log2(mean(cur_Robs)) - mean(cur_Robs));
        cur_Nspks = sum(cur_Robs);
        
        sac_LLimp(ss) = (nansum(cur_Robs.*log2(TB_pred_rate(temp)) - TB_pred_rate(temp)) - cur_nullLL)/cur_Nspks;
    end
    
    TB_sacmod.sac_modinfo = sac_info;
    TB_sacmod.sac_LLimp = sac_LLimp;
    
    all_TB_sacmod{ll} = TB_sacmod;
    
end

