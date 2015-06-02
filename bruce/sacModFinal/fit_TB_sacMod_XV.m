function [all_TB_sacmod] = fit_TB_sacMod_XV(Robs,stimG,t_since_sac,TB_params,Xsac,basemod_predrate,TB_d2T_lambda,tr_inds,xv_inds)

%%
%this function is like the main TB-model estimation function, but it uses
%separate training and xval data, and evaluates the MSE of separate
%gain,offset, and gain+offset regression models

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

fit_inds = union(tr_inds,xv_inds);

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

cur_tr_inds = find(ismember(TB_used_data,tr_inds));
cur_xv_inds = find(ismember(TB_used_data,xv_inds));

%this gives fixed reg strength in G-direction, but
%cur_sac_lambda in T-direction
L2_params = create_L2_params([],[1 n_Tbins*cur_nGbins],[n_Tbins cur_nGbins],2,3,[Inf Inf],[1 1/TB_d2T_lambda*G_lambda]);
TB_sacmod = regGLM_fit(TB_Xmat(cur_tr_inds,:),Robs(TB_used_data(cur_tr_inds)),L2_params,TB_d2T_lambda,[],[],silent);
[~, ~, TB_pred_rate] = regGLM_eval(TB_sacmod,Robs(TB_used_data),TB_Xmat);

TB_sacmod.lagX = Xtick(xbuff+1:end-xbuff);

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
TB_sacmod.TB_rate = TB_rate(:,xbuff+1:end-xbuff);
TB_sacmod.marg_gdist = marg_gdist;
TB_sacmod.marg_grate = marg_grate;

%% %compute response gain and offsets
[saconly_pred,gainonly_pred,both_pred] = deal(nan(length(cur_xv_inds),1));
[sac_offset,sac_gain] = deal(nan(length(TB_sacmod.lagX),1));
for ii = 1:length(TB_sacmod.lagX)
    temp = lscov([ones(cur_nGbins,1) marg_grate],TB_sacmod.TB_rate(:,ii),TB_sacmod.marg_gdist);
    sac_offset(ii) = temp(1); sac_gain(ii) = temp(2);
end
TB_sacmod.sac_gain = sac_gain;
TB_sacmod.sac_offset = sac_offset;

%%
[nsac_offset,nsac_gain] = deal(nan(size(Xsac,2),1));
[offonly_pred,gainonly_pred,both_pred] = deal(nan(length(cur_xv_inds),1));
for ii = 1:size(Xsac,2)
    temp_tr = cur_tr_inds(Xsac(TB_used_data(cur_tr_inds),ii) == 1);
    temp_xv = find(Xsac(TB_used_data(cur_xv_inds),ii) == 1);
    base_Xmat = [ones(length(temp_tr),1) basemod_predrate(TB_used_data(temp_tr))];
    base_Xmat_xv = [ones(length(temp_xv),1) basemod_predrate(TB_used_data(cur_xv_inds(temp_xv)))];
    
    rr = regress(TB_pred_rate(temp_tr)-base_Xmat(:,2),base_Xmat(:,1));
    offonly_pred(temp_xv) = base_Xmat_xv(:,1)*rr+base_Xmat_xv(:,2);
    
    rr = regress(TB_pred_rate(temp_tr),base_Xmat(:,2));
    gainonly_pred(temp_xv) = base_Xmat_xv(:,2)*rr;
   
    rr = regress(TB_pred_rate(temp_tr),base_Xmat);
    both_pred(temp_xv) = base_Xmat_xv*rr;

end

offonly_R2 = 1-nanmean((TB_pred_rate(cur_xv_inds) - offonly_pred).^2)/nanvar(TB_pred_rate(cur_xv_inds));
gainonly_MSE = 1-nanmean((TB_pred_rate(cur_xv_inds) - gainonly_pred).^2)/nanvar(TB_pred_rate(cur_xv_inds));
both_MSE = 1-nanmean((TB_pred_rate(cur_xv_inds) - both_pred).^2)/nanvar(TB_pred_rate(cur_xv_inds));

