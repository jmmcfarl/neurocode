function [TB_rate_R2] = get_TB_gain_offset_R2(Robs,stimG,t_since_sac,TB_params,Xsac,basemod_predrate,TB_d2T_lambda,tr_inds,xv_inds)

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
TB_sacmod = regGLM_fit(TB_Xmat,Robs(TB_used_data),L2_params,TB_d2T_lambda,[],[],silent); %fit using all data
[~, ~, TB_pred_rate] = regGLM_eval(TB_sacmod,Robs(TB_used_data),TB_Xmat);

%% go through each saccade latency and fit a set of regression models (using training data only to fit regression model params, then evaluate R2 using xval data)
[offonly_pred,gainonly_pred,both_pred] = deal(nan(length(cur_xv_inds),1));
for ii = 1:size(Xsac,2)
    temp_tr = cur_tr_inds(Xsac(TB_used_data(cur_tr_inds),ii) == 1);
    temp_xv = find(Xsac(TB_used_data(cur_xv_inds),ii) == 1);
    base_Xmat = [ones(length(temp_tr),1) basemod_predrate(TB_used_data(temp_tr))];
    base_Xmat_xv = [ones(length(temp_xv),1) basemod_predrate(TB_used_data(cur_xv_inds(temp_xv)))];
    
    %r(g,tau) = r0(g) + c(tau)
    rr = regress(TB_pred_rate(temp_tr)-base_Xmat(:,2),base_Xmat(:,1));
    offonly_pred(temp_xv) = base_Xmat_xv(:,1)*rr+base_Xmat_xv(:,2);
    
    %r(g,tau) = a(tau)*r0(g)
    rr = regress(TB_pred_rate(temp_tr),base_Xmat(:,2));
    gainonly_pred(temp_xv) = base_Xmat_xv(:,2)*rr;
   
    %r(g,tau) = a(tau)*r0(g) + c(tau)
    rr = regress(TB_pred_rate(temp_tr),base_Xmat);
    both_pred(temp_xv) = base_Xmat_xv*rr;
end

sac_resid = TB_pred_rate(cur_xv_inds) - basemod_predrate(TB_used_data(cur_xv_inds)); %this is the total amount of saccade modulation variance (to be explained)

%compute xval R2
TB_rate_R2.offset_only = 1-nanmean((TB_pred_rate(cur_xv_inds) - offonly_pred).^2)/nanvar(sac_resid);
TB_rate_R2.gain_only = 1-nanmean((TB_pred_rate(cur_xv_inds) - gainonly_pred).^2)/nanvar(sac_resid);
TB_rate_R2.both = 1-nanmean((TB_pred_rate(cur_xv_inds) - both_pred).^2)/nanvar(sac_resid);

