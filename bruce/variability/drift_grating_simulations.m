clear all

global Expt_name monk_name bar_ori rec_type

Expt_name = 'G086';
bar_ori = 0; %bar orientation to use (only for UA recs)
monk_name = 'jbe';

mod_data_name = 'corrected_models2';
compact_data_name = 'packaged_data';
et_anal_name = 'full_eyetrack_Rinit';

et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
mod_data_dir = ['~/Analysis/bruce/' Expt_name '/models'];

%%
data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir')
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
end
if ~exist(data_dir,'dir');
    error('Couldnt find data directory');
end
Edata_file = strcat(data_dir,'/',monk_name,Expt_name,'Expts');
fprintf('Loading %s\n',Edata_file);
load(Edata_file);

%is this a laminar probe or utah array rec?
if strcmp(Expts{1}.Header.DataType,'GridData 96')
    rec_type = 'UA';
elseif strcmp(Expts{1}.Header.DataType,'Spike2')
    rec_type = 'LP';
end

compact_file = [data_dir sprintf('/%s_ori%d',compact_data_name,bar_ori)];
fprintf('Loading %s\n',compact_file);
load(compact_file);

%%
%if using coil info
if any(params.use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end
et_file = [et_dir sprintf('/%s_ori%d',et_anal_name,bar_ori)];
fprintf('Loading %s\n',et_file);
load(et_file,'it_fix_post_*','drift_post_*','et_params');

mod_file = [mod_data_dir sprintf('/%s_ori%d',mod_data_name,bar_ori)];
fprintf('Loading %s\n',mod_file);
load(mod_file);

%%
spatial_usfac = 2;
base_sp_dx = mode(expt_data.expt_dw);
if length(unique(expt_data.expt_dw)) > 1
    fprintf('Warning, multiple different dws detected, using %.3f\n',base_sp_dx);
end
sp_dx = base_sp_dx/spatial_usfac/params.scale_fac; %model dx in deg

%%
NT = length(used_inds);
fullNT = size(spike_data.binned_mua,1);
n_trials = length(time_data.trial_flip_ids);
n_blocks = length(expt_data.used_blocks);

all_t_axis = time_data.t_axis;
trial_start_inds = [1+time_data.trial_flip_inds];
trial_end_inds = [time_data.trial_flip_inds(2:end); fullNT];
all_trialvec = nan(fullNT,1);
for ii = 1:n_trials
    all_trialvec(trial_start_inds(ii):trial_end_inds(ii)) = time_data.trial_flip_ids(ii);
end

block_start_inds = [1+time_data.block_flip_inds];
block_end_inds = [time_data.block_flip_inds(2:end); fullNT];
all_blockvec = nan(fullNT,1);
for ii = 1:n_blocks
    all_blockvec(block_start_inds(ii):block_end_inds(ii)) = time_data.block_flip_ids(ii);
end
Xblock = zeros(fullNT,n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% CREATE SACCADE PREDICTOR MATS
corrected_eye_vals_interp = ET_data.interp_eye_pos;
sac_start_times = [ET_data.saccades(:).start_time];
sac_stop_times = [ET_data.saccades(:).stop_time];

interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));

saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);
used_is_blink = ET_data.is_blink(used_saccade_set);

sac_amps = [ET_data.saccades(:).amplitude];
is_micro = sac_amps(used_saccade_set) < 1;
big_sacs = find(~is_micro & ~used_is_blink');
micro_sacs = find(is_micro & ~used_is_blink');
sac_durs = [ET_data.saccades(:).duration];

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

fix_start_inds = sort([trial_start_inds; saccade_stop_inds]);
fix_stop_inds = sort([trial_end_inds; saccade_start_inds]);
fix_durs = fix_stop_inds-fix_start_inds;
fix_start_inds(fix_durs<=0) = [];
fix_stop_inds(fix_durs <= 0) = [];
fix_post_blink = ismember(fix_start_inds,saccade_stop_inds(used_is_blink));
n_fixs = length(fix_start_inds);

fix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
end

%%
sac_buff = round(0.12/params.dt);
in_sac_inds = zeros(NT,1);
for ii = 1:sac_buff+1
    cur_inds = saccade_stop_inds + (ii-1);
    cur_inds(cur_inds > NT) = [];
    in_sac_inds(cur_inds) = 1;
end
for ii = 1:length(saccade_start_inds)
    cur_inds = saccade_start_inds(ii):saccade_stop_inds(ii);
    in_sac_inds(cur_inds) = 1;
end

in_sac_inds = logical(in_sac_inds);

blink_inds = find(used_is_blink);
in_blink_inds = zeros(NT,1);
for ii = 1:length(blink_inds)
    cur_inds = saccade_start_inds(blink_inds(ii)):saccade_stop_inds(blink_inds(ii));
    in_blink_inds(cur_inds) = 1;
end
in_blink_inds = logical(in_blink_inds);

%%
sac_shift = round(0.05/params.dt);

cur_fix_post_mean = squeeze(it_fix_post_mean(end,:));
cur_fix_post_std = squeeze(it_fix_post_std(end,:));
cur_drift_post_mean = squeeze(drift_post_mean(end,:));
cur_drift_post_std = squeeze(drift_post_std(end,:));

[fin_tot_corr,fin_tot_std] = construct_eye_position(cur_fix_post_mean,cur_fix_post_std,...
    cur_drift_post_mean,cur_drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);
ET_pos_orth = fin_tot_corr*et_params.sp_dx;

%%
min_trial_dur = 4;
trial_durs = [trial_data(:).end_times] - [trial_data(:).start_times];
used_trials = find(trial_durs >= min_trial_dur & ismember(1:length(trial_durs),all_trialvec(used_inds)));
nf = (params.trial_dur - params.beg_buffer - params.end_buffer)/params.dt;
NusedTrials = length(used_trials);
tbt_EP_orth = nan(nf,NusedTrials);
tbt_in_sac = nan(nf,NusedTrials);
tbt_in_blink = nan(nf,NusedTrials);
for ii = 1:NusedTrials
    cur_inds = find(all_trialvec(used_inds) == used_trials(ii));
    tbt_EP_orth(:,ii) = ET_pos_orth(cur_inds);
    tbt_in_sac(:,ii) = in_sac_inds(cur_inds);
    tbt_in_blink(:,ii) = in_blink_inds(cur_inds);
end
tbt_in_sac = logical(tbt_in_sac);
tbt_in_blink = logical(tbt_in_blink);
%%
su_set = (params.n_probes+1):length(ModData);
bad_sus = [];
use_mods = [];
for ii = 1:length(su_set)
    if ~isempty(ModData(su_set(ii)).bestGQM)
        flen = ModData(su_set(ii)).bestGQM.stim_params.stim_dims(1);
        nPix = ModData(su_set(ii)).bestGQM.stim_params.stim_dims(2);

        use_mods = cat(1,use_mods,ModData(su_set(ii)).bestGQM);
        use_mods(end).mods([1 end]) = [];
    else
        bad_sus = cat(1,bad_sus,ii);
    end
end
su_set(bad_sus) = [];

%%
ex_cell = 103;
base_mod = ModData(ex_cell).bestGQM;
flen = base_mod.stim_params.stim_dims(1);
nPix = base_mod.stim_params.stim_dims(2);
% base_filt = reshape(base_mod.mods(2).filtK,[flen nPix]);
% 
% xax = (1:nPix)*sp_dx;
% tax  = (1:flen);
% tprofile = std(base_filt,[],2);
% 
gsf = ModData(ex_cell).tune_props.RF_gSF;
% gcent = mean(xax);
% gRF_SD = ModData(ex_cell).tune_props.RF_sigma;
% 
% [XX,TT] = meshgrid(xax,tax);
% gabor_RF = exp(-(XX-gcent).^2/(2*gRF_SD^2));
% gabor_RF = gabor_RF.*sin(2*pi*(XX*gsf));
% gabor_RF = bsxfun(@rdivide,gabor_RF,std(gabor_RF,[],2));
% gabor_RF = bsxfun(@times,gabor_RF,tprofile);
% 
% gabor_QP = exp(-(XX-gcent).^2/(2*gRF_SD^2));
% gabor_QP = gabor_QP.*sin(2*pi*(XX*gsf)+pi/2);
% gabor_QP = bsxfun(@rdivide,gabor_QP,std(gabor_QP,[],2));
% gabor_QP = bsxfun(@times,gabor_QP,tprofile);
% 
% base_mod.mods(2:end) = [];
% base_mod.mods(1).filtK = gabor_RF(:);
% base_mod.mods(1).NLtype = 'quad';
% base_mod.mods(2) = base_mod.mods(1);
% base_mod.mods(2).filtK = gabor_QP(:);
% 
% su_set = 1;
% use_mods = base_mod;
% 
% base_mod.mods(2) = [];
% base_mod.mods(1).NLtype = 'lin';
% use_mods(2) = base_mod;
% su_set = [1 2];
%%
stim_params = NMMcreate_stim_params([flen nPix]);

test_NT = 1e4;
test_dd = mode(expt_data.expt_dds);

test_stim = randi(2,test_NT,nPix);
test_stim(test_stim == 2) = -1;
is_zero = rand(test_NT,nPix) < test_dd/100;
test_stim(is_zero) = 0;

Xmat = create_time_embedding(test_stim,stim_params);

[test_gSD,test_meanrate] = deal(nan(length(su_set),1));
for ii = 1:length(su_set)
    [~,~,test_prate,~,gint] = NMMeval_model(use_mods(ii),[],Xmat);
    test_gSD(ii) = mean(std(gint));
    test_meanrate(ii) = mean(test_prate);
end
%%
poss_EP_scales = [0:0.1:2];
for ep = 1:length(poss_EP_scales)
    % EP_scale = 1;
    EP_scale = poss_EP_scales(ep)
    
    xax = (1:nPix)*sp_dx;
    tax  = (1:nf)*params.dt;
    
    gr_sf = ModData(ex_cell).tune_props.RF_gSF;
    gr_tf = 4;
    
    [XX,TT] = meshgrid(xax,tax);
    base_dg = sin(2*pi*(gr_sf*XX + gr_tf*TT));
    
    full_EP_orth = EP_scale*tbt_EP_orth(:);
    % full_EP_orth(tbt_in_blink(:)) = nan;
    % full_EP_orth(tbt_in_sac(:)) = nan;
    
    [XX,TT,RR] = meshgrid(xax,tax,1:NusedTrials);
    TT = reshape(permute(TT,[1 3 2]),[],length(xax));
    XX = reshape(permute(XX,[1 3 2]),[],length(xax));
    
    drift_grating = sin(2*pi*(gr_sf*(bsxfun(@plus,XX,full_EP_orth)) + gr_tf*TT));
    % drift_grating = reshape(drift_grating,[nf NusedTrials length(xax)]);
    
    %%
    throwout_win = round(0.15/params.dt);
    
    Xmat = create_time_embedding(drift_grating,stim_params);
    base_Xmat = create_time_embedding(base_dg,stim_params);
    
    [tot_resp_var,avg_acrosstrial_var,avg_acrosstrial_SD,tot_resp_mean,tot_resp_SD] = deal(nan(length(su_set),1));
    [all_base_psth,all_psth,all_psth_SD] = deal(nan(nf-throwout_win,length(su_set)));
    for ii = 1:length(su_set)
        [~,~,~,~,gint] = NMMeval_model(use_mods(ii),[],Xmat);
        contrast_scale = test_gSD(ii)/mean(nanstd(gint));
        
        [~,~,pred_rate] = NMMeval_model(use_mods(ii),[],Xmat*contrast_scale);
        pred_rate = reshape(pred_rate,nf,NusedTrials);
        pred_rate(tbt_in_blink) = nan;
        
        
        [~,~,base_pred_rate,~,gints] = NMMeval_model(use_mods(ii),[],base_Xmat*contrast_scale);
        epscale_data(ep).all_base_psth(:,ii) = base_pred_rate((throwout_win+1):end);
        
        ch_pred_rate = pred_rate((throwout_win+1):end,:);
        epscale_data(ep).all_psth(:,ii) = nanmean(ch_pred_rate,2);
        epscale_data(ep).all_psth_SD(:,ii) = nanstd(ch_pred_rate,[],2);
        
        epscale_data(ep).tot_resp_mean(ii) = nanmean(ch_pred_rate(:));
        epscale_data(ep).tot_resp_SD(ii) = nanstd(ch_pred_rate(:));
        epscale_data(ep).tot_resp_var(ii) = nanvar(ch_pred_rate(:));
        across_trial_var = nanvar(ch_pred_rate,[],2);
        epscale_data(ep).avg_acrosstrial_var(ii) = mean(nanvar(ch_pred_rate,[],2));
        epscale_data(ep).avg_acrosstrial_SD(ii) = mean(nanstd(ch_pred_rate,[],2));
        epscale_data(ep).psth_var(ii) = var(epscale_data(ep).all_psth(:,ii));
        epscale_data(ep).all_base_psth(ii) = var(epscale_data(ep).all_base_psth(:,ii));
        
        %%
        pred_rate(tbt_in_sac) = nan;
        
        [~,~,base_pred_rate,~,gints] = NMMeval_model(use_mods(ii),[],base_Xmat*contrast_scale);
        epscale_data_nosac(ep).all_base_psth(:,ii) = base_pred_rate((throwout_win+1):end);
        
        ch_pred_rate = pred_rate((throwout_win+1):end,:);
        epscale_data_nosac(ep).all_psth(:,ii) = nanmean(ch_pred_rate,2);
        epscale_data_nosac(ep).all_psth_SD(:,ii) = nanstd(ch_pred_rate,[],2);
        
        epscale_data_nosac(ep).tot_resp_mean(ii) = nanmean(ch_pred_rate(:));
        epscale_data_nosac(ep).tot_resp_SD(ii) = nanstd(ch_pred_rate(:));
        epscale_data_nosac(ep).tot_resp_var(ii) = nanvar(ch_pred_rate(:));
        across_trial_var = nanvar(ch_pred_rate,[],2);
        epscale_data_nosac(ep).avg_acrosstrial_var(ii) = mean(nanvar(ch_pred_rate,[],2));
        epscale_data_nosac(ep).avg_acrosstrial_SD(ii) = mean(nanstd(ch_pred_rate,[],2));
        epscale_data_nosac(ep).psth_var(ii) = var(epscale_data_nosac(ep).all_psth(:,ii));
        epscale_data_nosac(ep).all_base_psth(ii) = var(epscale_data_nosac(ep).all_base_psth(:,ii));
    end
    
end

%%
all_means = reshape([epscale_data(:).tot_resp_mean],length(su_set),[]);
all_SDs = reshape([epscale_data(:).avg_acrosstrial_SD],length(su_set),[]);
all_at_var = reshape([epscale_data(:).avg_acrosstrial_var],length(su_set),[]);
all_tot_var = reshape([epscale_data(:).tot_resp_var],length(su_set),[]);

all_means_ns = reshape([epscale_data_nosac(:).tot_resp_mean],length(su_set),[]);
all_SDs_ns = reshape([epscale_data_nosac(:).avg_acrosstrial_SD],length(su_set),[]);
all_at_var_ns = reshape([epscale_data_nosac(:).avg_acrosstrial_var],length(su_set),[]);
all_tot_var_ns = reshape([epscale_data_nosac(:).tot_resp_var],length(su_set),[]);

cmap = jet(length(su_set));
f1 = figure(); hold on
for ii = 1:length(su_set)
plot(poss_EP_scales*gsf,all_SDs(ii,:)./all_means(ii,:),'color',cmap(ii,:))
plot(poss_EP_scales*gsf,all_SDs_ns(ii,:)./all_means_ns(ii,:),'color',cmap(ii,:),'linestyle','--')
end

f2 = figure(); hold on
for ii = 1:length(su_set)
plot(poss_EP_scales*gsf,all_at_var(ii,:)./all_tot_var(ii,:),'color',cmap(ii,:))
plot(poss_EP_scales*gsf,all_at_var_ns(ii,:)./all_tot_var_ns(ii,:),'color',cmap(ii,:),'linestyle','--')
end

% plot(poss_EP_scales*gsf,all_SDs(2,:)./all_means(2,:),'r')
% plot(poss_EP_scales*gsf,all_SDs_ns(2,:)./all_means_ns(2,:),'r--')


