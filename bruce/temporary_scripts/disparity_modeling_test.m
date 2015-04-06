clear all
close all

% monName = 'jbe';
% Expt_name = 'M008';
monName = 'lem';
Expt_name = 'M312';

rec_type = 'LP';
data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
% data_dir = ['~/Data/bruce/' Expt_name];
% cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];


%%

min_trial_dur = 2;
trial_dur = 4;
use_nPix = 30;
dt = 0.03;
flen = 5;
Fr = 1;
n_probes = 24;
spatial_usfac = 1;

%%
cd(data_dir);
ExptFileName = strcat(monName,Expt_name,'Expts.mat');
fprintf('Loading existing file %s\n',ExptFileName);
load(ExptFileName);
cd stims/
load('./stim_data.mat');
load('./expt_data.mat');

%%
include_expts = {'rls.dpXceXFrRC'};
usable_expts = find(cellfun(@(x) ~isempty(x),Expts));
expt_names = cellfun(@(x) x.Header.expname,Expts(usable_expts),'uniformoutput',0);
        
cur_block_set = usable_expts(cellfun(@(x) any(strcmp(x,include_expts)),expt_names));
n_blocks = length(cur_block_set);

full_nPix = unique(expt_npix(cur_block_set));
if length(full_nPix) > 1
    warning('multiple npix detected');
end

%% PARSE SU CLUSTERS (ASSUME THEY STAY ON SAME PROBES)
load([data_dir '/CellList']);
temp_Clist = CellList(cur_block_set,:,1);
su_set = find(any(temp_Clist > 0));
su_numbers = unique(temp_Clist(temp_Clist > 0));
su_block_good = false(n_blocks,length(su_set));
[PP,BB] = meshgrid(1:n_probes,1:length(cur_block_set));
for ii = 1:length(su_set)
    bset = find(CellList(cur_block_set,:,1) == su_numbers(ii));
    su_block_good(BB(bset),ii) = true;
end
%% Spatial resolution
all_dws = cellfun(@(x) x.Stimvals.dw,Expts(cur_block_set));
base_sp_dx = mode(all_dws);
if length(unique(all_dws)) > 1
    fprintf('Warning, multiple different dws detected, using %.3f\n',base_sp_dx);
end
sp_dx = base_sp_dx/spatial_usfac; %model dx in deg
ov_RF_pos = Expts{cur_block_set(1)}.Stimvals.rf(1:2);

%%
buffer_pix = floor((full_nPix - use_nPix)/2);
[Xinds_up,~] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
cur_use_pix = (1/spatial_usfac:1/spatial_usfac:use_nPix) + buffer_pix;

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
all_Lstim_mat = [];
all_Rstim_mat = [];
all_frame_ce = [];
all_frame_dp = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_blk = [];
all_trial_Fr = [];
all_nframes = [];
all_nstims = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_trial_rptframes = [];

trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
cur_spkind_offset = 0;
for ee = 1:n_blocks;
    fprintf('Block %d of %d;  UNMATCHED EXPT TYPE\n',ee,n_blocks);
    cur_block = cur_block_set(ee);
    
    fname = sprintf('Expt%dClusterTimes.mat',cur_block_set(ee));
    load(fname);
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times' + cur_toffset);
%         all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},ones(size(Clusters{cc}.times')));
    end
    
    trial_start_times = [Expts{cur_block}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_block}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_block}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_block}.Trials(:).id];
    [un_ids,id_inds] = unique(trial_ids);
    rpt_trials = false;
    if length(un_ids) < length(trial_ids)
        rpt_trials = true;
        fprintf('Warning, repeat trial inds detected!\n');
        use_trials = [];
    else
        use_trials = find(trial_durs >= min_trial_dur);
    end
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    if dt == 0.03
        trial_Fr = [Expts{cur_block}.Trials(:).Fr];
        use_trials = use_trials(trial_Fr(use_trials) == 3);
    end
    
    trial_Se = [Expts{cur_block}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    trial_Fr = [Expts{cur_block}.Trials(:).Fr];
    trial_Fr = trial_Fr(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    all_trial_Fr = cat(1,all_trial_Fr,trial_Fr(use_trials)');
    all_trial_blk = cat(1,all_trial_blk,ones(length(use_trials),1)*ee);
        
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    
    n_trials = length(use_trials);
    cur_nrpt_frames = zeros(n_trials,1);
    [trial_nframes,trial_nstims] = deal(nan(n_trials,1));
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
        cur_stim_ce = Expts{cur_block}.Trials(use_trials(tt)).ce;
        cur_stim_dp = Expts{cur_block}.Trials(use_trials(tt)).dp;
        
        if trial_Fr(use_trials(tt)) == 3 && dt == 0.01
           cur_stim_times = repmat(cur_stim_times,[3 1]);
           cur_stim_times = bsxfun(@plus,cur_stim_times,[0; dt; 2*dt]);
           cur_stim_times = cur_stim_times(:); cur_stim_times = cur_stim_times(1:301)';
           cur_stim_ce = repmat(cur_stim_ce',[3 1]);
           cur_stim_ce = cur_stim_ce(:); cur_stim_ce = cur_stim_ce(1:301);
           cur_stim_dp = repmat(cur_stim_dp',[3 1]);
           cur_stim_dp = cur_stim_dp(:); cur_stim_dp = cur_stim_dp(1:301); 
        end
        
        if dt == 0.03
            left_stim_mats{use_trials(tt)} = left_stim_mats{use_trials(tt)}(1:3:size(left_stim_mats{use_trials(tt)},1),:);
            right_stim_mats{use_trials(tt)} = right_stim_mats{use_trials(tt)}(1:3:size(right_stim_mats{use_trials(tt)},1),:);
        end
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        trial_nframes(tt) = n_frames;
        trial_nstims(tt) = length(cur_stim_times);
        if isfield(Expts{cur_block}.Trials(use_trials(tt)),'rptframes')
            cur_nrpt_frames(tt) = length(Expts{cur_block}.Trials(use_trials(tt)).rptframes);
        end
        if n_frames > 0
            if length(cur_stim_times) == 1
                cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
            end
            cur_t_edges = [cur_stim_times'; cur_stim_times(end) + dt*Fr];
        end
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt & trial_nframes(tt) == trial_nstims(tt)
            use_frames = min(length(cur_stim_times),n_frames);
            cur_Lstim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            cur_Rstim_mat = double(right_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times+cur_toffset < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_nframes = [all_nframes; trial_nframes(:)];
            all_nstims = [all_nstims; trial_nstims(:)];
            all_stim_times = [all_stim_times; cur_stim_times' + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_Lstim_mat = [all_Lstim_mat; cur_Lstim_mat];
            all_Rstim_mat = [all_Rstim_mat; cur_Rstim_mat];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_frame_ce = [all_frame_ce; cur_stim_ce];
            all_frame_dp = [all_frame_dp; cur_stim_dp];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
    all_trial_rptframes = [all_trial_rptframes; cur_nrpt_frames];
    
    %need to keep track of block time offsets for LP recordings
    if strcmp(rec_type,'LP')
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
        cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
    end
end

%% BIN SPIKES FOR MU AND SU
clust_params.n_probes = n_probes;
all_binned_mua = get_quickbinned_mu(all_spk_times,all_clust_ids,all_t_axis,all_t_bin_edges,all_bin_edge_pts,clust_params);

all_binned_sua = nan(size(all_binned_mua,1),length(su_set));
for ii = 1:length(su_set)
    cur_inds = find(ismember(all_blockvec,find(su_block_good(:,ii))));
    all_binned_sua(cur_inds,ii) = all_binned_mua(cur_inds,su_set(ii));
end

%%
full_nPix_us = spatial_usfac*full_nPix;
stim_params_us = NMMcreate_stim_params([flen full_nPix_us*2],dt);

%%
% flip_frames = all_frame_ce == -1;
% all_Rstim_mat(flip_frames,:) = -all_Rstim_mat(flip_frames,:);
% % all_Lstim_mat(flip_frames,:) = -all_Lstim_mat(flip_frames,:);

%%
blank_frames = all_frame_ce == -1009;
all_Rstim_mat(blank_frames,:) = 0;
all_Lstim_mat(blank_frames,:) = 0;
%% CREATE STIMULUS MATRIX
stim_params = NIMcreate_stim_params([flen use_nPix*2],dt);
all_Xmat = create_time_embedding([all_Lstim_mat all_Rstim_mat],stim_params);

%%
all_Lstim_mat(all_Lstim_mat == -128) = nan;
all_Rstim_mat(all_Rstim_mat == -128) = nan;
all_Xmat(all_Xmat == -128) = nan;
bad_trials = [];
good_inds = [];
un_trials = unique(all_trialvec);
for tt = 1:length(un_trials)
    cur_set = find(all_trialvec == un_trials(tt));
    if any(isnan(reshape(all_Xmat(cur_set,:),[],1)))
        bad_trials = [bad_trials un_trials(tt)];
    else
        good_inds = [good_inds; cur_set(:)];
    end
end
%%
% DEFINE DATA USED FOR ANALYSIS
beg_buffer = 0.15;
end_buffer = 0.05;
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
NT = length(used_inds);

used_inds(ismember(all_trialvec(used_inds),bad_trials)) = [];

utrials = find(all_trial_Fr == dt/0.01);
% utrials = find(all_trial_Fr == 3);
used_inds(~ismember(all_trialvec(used_inds),utrials)) = [];

%%
silent = 0;
base_lambda_d2XT = 20;
base_lambda_L1 = 1;

if dt == 0.01 %use st smoothness
L2_params = create_L2_params([],[1; use_nPix*flen],[flen use_nPix],2,3,[0 0]);
L2_mat_left = generate_L2_mat(L2_params,2*use_nPix*flen);
L2_params = create_L2_params([],[(use_nPix*flen + 1); 2*use_nPix*flen],[flen use_nPix],2,3,[0 0]);
L2_mat_right = generate_L2_mat(L2_params,2*use_nPix*flen);
elseif dt == 0.03 %use spatial smoothness
L2_params = create_L2_params([],[1; use_nPix*flen],[flen use_nPix],2,2,[0 0]);
L2_mat_left = generate_L2_mat(L2_params,2*use_nPix*flen);
L2_params = create_L2_params([],[(use_nPix*flen + 1); 2*use_nPix*flen],[flen use_nPix],2,2,[0 0]);
L2_mat_right = generate_L2_mat(L2_params,2*use_nPix*flen);
end

L2_mat = L2_mat_left + L2_mat_right;

n_stim_filts = 2;
mod_signs = [1 1 1 -1 -1];
NL_types = [{'threshlin','threshlin','threshlin','threshlin','threshlin'}];
init_d2XT = [ones(n_stim_filts,1)];
init_L2 = [zeros(n_stim_filts,1);];
init_reg_params = NMMcreate_reg_params('lambda_custom',base_lambda_d2XT,'lambda_L1',base_lambda_L1);

optim_params.optTol = 1e-6;
optim_params.progTol = 1e-8;

mod_pred_rates = nan(length(used_inds),length(su_set));
% for ss = 1:length(su_set)
for ss = [6]
    fprintf('Fitting model for SU %d of %d\n',ss,length(su_set));
    Robs = all_binned_sua(used_inds,ss);
    cur_uinds = find(~isnan(Robs));
    
    if nansum(Robs) > 0
        gqm1 = NMMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params);
        gqm1 = NMMfit_filters(gqm1,Robs,all_Xmat(used_inds,:),[],cur_uinds,silent,optim_params,L2_mat);
        [LL, ~, pred_rate, G, gint] = NMMeval_model(gqm1,Robs,all_Xmat(used_inds,:),[],cur_uinds);
        
        gqm2 = NMMadjust_regularization(gqm1,1:length(mod_signs),'lambda_custom',base_lambda_d2XT./var(gint)');
        gqm2 = NMMadjust_regularization(gqm2,1:length(mod_signs),'lambda_L1',base_lambda_L1./std(gint)');
        gqm2 = NMMfit_filters(gqm2,Robs,all_Xmat(used_inds,:),[],cur_uinds,silent,[],L2_mat);
        
        all_gqm1(ss) = gqm1;
        all_gqm2(ss) = gqm2;
        
        [LL, penLL, pred_rate, G, gint, fgint, nullLL] = NMMeval_model(gqm2,Robs,all_Xmat(used_inds,:),[],cur_uinds);
        mod_pred_rates(cur_uinds,ss) = pred_rate;
    end
end

%%
% use_lag = 3;
% [XX,TT] = meshgrid(1:use_nPix*2,1:flen);
% 
% cur_use_Kinds = find(TT == use_lag);
% cur_stim_params = NMMcreate_stim_params([1 use_nPix*2]);
% 
% silent = 1;
% base_lambda_d2XT = 20;
% base_lambda_L1 = 0;
% 
% L2_params = create_L2_params([],[1; use_nPix],[1 use_nPix],2,2,[0 0]);
% L2_mat_left = generate_L2_mat(L2_params,2*use_nPix);
% L2_params = create_L2_params([],[(use_nPix + 1); 2*use_nPix],[1 use_nPix],2,2,[0 0]);
% L2_mat_right = generate_L2_mat(L2_params,2*use_nPix);
% 
% L2_mat = L2_mat_left + L2_mat_right;
% 
% n_stim_filts = 5;
% mod_signs = [1 1 1 -1 -1];
% NL_types = [{'lin' 'quad','quad','quad','quad'}];
% init_d2XT = [ones(n_stim_filts,1)];
% init_L2 = [zeros(n_stim_filts,1);];
% init_reg_params = NMMcreate_reg_params('lambda_custom',base_lambda_d2XT,'lambda_L1',base_lambda_L1);
% 
% pos_used_inds = used_inds(all_frame_ce(used_inds) == 0);
% neg_used_inds = used_inds(all_frame_ce(used_inds) == -1);
% 
% mod_pred_rates = nan(length(used_inds),length(su_set));
% % for ss = 1:length(su_set)
% for ss = [1 6]
%     fprintf('Fitting model for SU %d of %d\n',ss,length(su_set));
%     Robs = all_binned_sua(pos_used_inds,ss);
%     cur_uinds = find(~isnan(Robs));
%     
%     if nansum(Robs) > 0
%         gqm1 = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,init_reg_params);
%         gqm1 = NMMfit_filters(gqm1,Robs,all_Xmat(pos_used_inds,cur_use_Kinds),[],cur_uinds,silent,[],L2_mat);
%         [LL, ~, pred_rate, G, gint] = NMMeval_model(gqm1,Robs,all_Xmat(pos_used_inds,cur_use_Kinds),[],cur_uinds);
%         
%         gqm2 = NMMadjust_regularization(gqm1,1:length(mod_signs),'lambda_custom',base_lambda_d2XT./var(gint)');
%         gqm2 = NMMadjust_regularization(gqm2,1:length(mod_signs),'lambda_L1',base_lambda_L1./std(gint)');
%         gqm2 = NMMfit_filters(gqm2,Robs,all_Xmat(pos_used_inds,cur_use_Kinds),[],cur_uinds,silent,[],L2_mat);
%         
%         pos_ce_mod(ss) = gqm2;
%         
%         [LL, penLL, pred_rate, G, gint, fgint, nullLL] = NMMeval_model(gqm2,[],all_Xmat(used_inds,cur_use_Kinds));
%         mod_pred_rates(used_inds,ss) = pred_rate;
%     end
% %     
% %     fprintf('Fitting model for SU %d of %d\n',ss,length(su_set));
% %     Robs = all_binned_sua(neg_used_inds,ss);
% %     cur_uinds = find(~isnan(Robs));
% %     
% %     if nansum(Robs) > 0
% %         gqm1 = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,init_reg_params);
% %         gqm1 = NMMfit_filters(gqm1,Robs,all_Xmat(neg_used_inds,cur_use_Kinds),[],cur_uinds,silent,[],L2_mat);
% %         [LL, ~, pred_rate, G, gint] = NMMeval_model(gqm1,Robs,all_Xmat(neg_used_inds,cur_use_Kinds),[],cur_uinds);
% %         
% %         gqm2 = NMMadjust_regularization(gqm1,1:length(mod_signs),'lambda_custom',base_lambda_d2XT./var(gint)');
% %         gqm2 = NMMadjust_regularization(gqm2,1:length(mod_signs),'lambda_L1',base_lambda_L1./std(gint)');
% %         gqm2 = NMMfit_filters(gqm2,Robs,all_Xmat(neg_used_inds,cur_use_Kinds),[],cur_uinds,silent,[],L2_mat);
% %         
% %         neg_ce_mod(ss) = gqm2;
% %         
% % %         [LL, penLL, pred_rate, G, gint, fgint, nullLL] = NMMeval_model(gqm2,[],all_Xmat(used_inds,cur_use_Kinds));
% % %         mod_pred_rates(used_inds,ss) = pred_rate;
% %     end
%     
% end

%%
poss_latency = 1:10;
un_dp = unique(all_frame_dp);
dp_trig_avgs = nan(length(un_dp),length(poss_latency),length(su_set));
dp_mp = nan(length(un_dp),length(poss_latency),length(su_set));
use_dps = 3:length(un_dp);

for pp = 1:length(poss_latency)
    latency = poss_latency(pp);
    for ii = 1:length(un_dp)
        cur_set = find(all_frame_dp(used_inds) == un_dp(ii) & all_frame_ce(used_inds) == 1);
        cur_resp_set = cur_set + latency;
        bad = cur_resp_set > length(used_inds);
        cur_set(bad) = []; cur_resp_set(bad) = [];
        cur_resp_set(all_trialvec(used_inds(cur_resp_set)) ~= all_trialvec(used_inds(cur_set))) = [];
        dp_trig_avgs(ii,pp,:) = nanmean(all_binned_sua(used_inds(cur_resp_set),:));
        dp_mp(ii,pp,:) = nanmean(mod_pred_rates((cur_resp_set),:));
        
        cur_set = find(all_frame_dp(used_inds) == un_dp(ii) & all_frame_ce(used_inds) == -1);
        cur_resp_set = cur_set + latency;
        bad = cur_resp_set > length(used_inds);
        cur_set(bad) = []; cur_resp_set(bad) = [];
        cur_resp_set(all_trialvec(used_inds(cur_resp_set)) ~= all_trialvec(used_inds(cur_set))) = [];
        dp_trig_avgsR(ii,pp,:) = nanmean(all_binned_sua(used_inds(cur_resp_set),:));
        dp_mpR(ii,pp,:) = nanmean(mod_pred_rates((cur_resp_set),:));
    end
end


%%
un = 6;
ln = 1;
f1 = figure(); hold on
 plot(un_dp(3:end),dp_trig_avgs(3:end,ln,un)/dt,'b')
 plot(un_dp(3:end),dp_mp(3:end,ln,un)/dt,'r')
 plot(un_dp(3:end),dp_trig_avgsR(3:end,ln,un)/dt,'k')
 plot(un_dp(3:end),dp_mpR(3:end,ln,un)/dt,'g')
xlim([-1.1 1.1])
xlabel('Dp')
ylabel('Firing rate (Hz)')
legend('Forward corr ce1','Model-predicted ce1','Forward corr ce-1','Model-predicted ce-1')
%%
% nneg = 5; npos = 5;
% stc_thresh = -5e-3;
% 
% for ss = 1:length(su_set)
%     
%     fprintf('Computing STC for SU %d of %d\n',ss,length(su_set));
%     
%     Robs = all_binned_sua(used_inds,ss);
%     cur_used_inds = used_inds(~isnan(Robs));
%     if ~isempty(cur_used_inds)
%         Robs = all_binned_sua(cur_used_inds,ss);
%         avg_rate = mean(Robs);
%         
%         spikebins = convert_to_spikebins(Robs);
%         spike_cond_stim = all_Xmat(cur_used_inds(spikebins),:);
%         sta      = mean(spike_cond_stim) - mean(all_Xmat(cur_used_inds,:));
%         sta = sta/norm(sta);
%         proj_mat = sta'/(sta*sta')*sta;
%         stim_proj = all_Xmat(cur_used_inds,:) - all_Xmat(cur_used_inds,:)*proj_mat;
%         % stim_proj = stim_emb;
%         stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
%         [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
%         stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
%         
%         sua_data(ss).avg_rates = avg_rate;
%          sua_data(ss).var_rates = nanvar(Robs);
%        sua_data(ss).sta = sta;
%         sua_data(ss).stcs = stcs;
%         sua_data(ss).evals = diag(evals);
%         cur_evec_diff = diff(flipud(sua_data(ss).evals));
%         sua_data(ss).npos_stc = find(cur_evec_diff(1:250) < stc_thresh,1,'last');
%         sua_data(ss).nneg_stc = find(cur_evec_diff(end:-1:250) < stc_thresh,1,'last');
%         if isempty(sua_data(ss).npos_stc); sua_data(ss).npos_stc = 0; end;
%         if isempty(sua_data(ss).nneg_stc); sua_data(ss).nneg_stc = 0; end;
%         sua_data(ss).stc_use = true;
%         sua_data(ss).avg_rate = mean(Robs);
%         sua_data(ss).nspks = sum(Robs);
%         
%         %%
%     else
%         sua_data(ss).stc_use = false;
%     end
% end

%%
% close all
% for ss = 1:length(su_set)
% subplot(3,1,1)
% imagesc(reshape(sua_data(ss).sta,[flen use_nPix*2]))
% colormap(gray)
% subplot(3,1,2)
% imagesc(reshape(sua_data(ss).stcs(:,1),[flen use_nPix*2]))
% colormap(gray)
% subplot(3,1,3)
% imagesc(reshape(sua_data(ss).stcs(:,2),[flen use_nPix*2]))
% colormap(gray)
% pause
% clf
% end
%%
mu_set = 1:24;
mu_mod_pred_rates = nan(length(used_inds),length(mu_set));
for ss = 1:length(mu_set)
    fprintf('Fitting model for MU %d of %d\n',ss,n_probes);
    Robs = all_binned_mua(used_inds,ss);
    cur_uinds = find(~isnan(Robs));
    
    if nansum(Robs) > 0
        gqm1 = NMMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params);
        gqm1 = NMMfit_filters(gqm1,Robs,all_Xmat(used_inds,:),[],cur_uinds,silent,[],L2_mat);
        [LL, ~, pred_rate, G, gint] = NMMeval_model(gqm1,Robs,all_Xmat(used_inds,:),[],cur_uinds);
        
        gqm2 = NMMadjust_regularization(gqm1,1:length(mod_signs),'lambda_d2XT',base_lambda_d2XT./var(gint)');
        gqm2 = NMMadjust_regularization(gqm2,1:length(mod_signs),'lambda_L1',base_lambda_L1./std(gint)');
        gqm2 = NMMfit_filters(gqm2,Robs,all_Xmat(used_inds,:),[],cur_uinds,silent,[],L2_mat);
        
        all_mu_gqm1(ss) = gqm1;
        all_mu_gqm2(ss) = gqm2;
        
        [LL, penLL, pred_rate, G, gint, fgint, nullLL] = NMMeval_model(gqm2,Robs,all_Xmat(used_inds,:),[],cur_uinds);
        mu_mod_pred_rates(cur_uinds,ss) = pred_rate;
    end
end

%%
nneg = 5; npos = 5;
stc_thresh = -5e-3;
mu_set = 1:n_probes;
for ss = 1:length(mu_set)
    
    fprintf('Computing STC for SU %d of %d\n',ss,length(mu_set));
    
    Robs = all_binned_mua(used_inds,ss);
    cur_used_inds = used_inds(~isnan(Robs));
    if ~isempty(cur_used_inds)
        Robs = all_binned_mua(cur_used_inds,ss);
        avg_rate = mean(Robs);
        
        spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = all_Xmat(cur_used_inds(spikebins),:);
        sta      = mean(spike_cond_stim) - mean(all_Xmat(cur_used_inds,:));
        sta = sta/norm(sta);
        proj_mat = sta'/(sta*sta')*sta;
        stim_proj = all_Xmat(cur_used_inds,:) - all_Xmat(cur_used_inds,:)*proj_mat;
        % stim_proj = stim_emb;
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
        
        mua_data(ss).avg_rates = avg_rate;
         mua_data(ss).var_rates = nanvar(Robs);
       mua_data(ss).sta = sta;
        mua_data(ss).stcs = stcs;
        mua_data(ss).evals = diag(evals);
        cur_evec_diff = diff(flipud(mua_data(ss).evals));
        mua_data(ss).npos_stc = find(cur_evec_diff(1:250) < stc_thresh,1,'last');
        mua_data(ss).nneg_stc = find(cur_evec_diff(end:-1:250) < stc_thresh,1,'last');
        if isempty(mua_data(ss).npos_stc); mua_data(ss).npos_stc = 0; end;
        if isempty(mua_data(ss).nneg_stc); mua_data(ss).nneg_stc = 0; end;
        mua_data(ss).stc_use = true;
        mua_data(ss).avg_rate = mean(Robs);
        mua_data(ss).nspks = sum(Robs);
        
        %%
    else
        mua_data(ss).stc_use = false;
    end
end

%%
close all
for ss = 1:n_probes
subplot(3,1,1)
imagesc(reshape(mua_data(ss).sta,[flen use_nPix*2]))
colormap(gray)
subplot(3,1,2)
imagesc(reshape(mua_data(ss).stcs(:,1),[flen use_nPix*2]))
colormap(gray)
subplot(3,1,3)
imagesc(reshape(mua_data(ss).stcs(:,2),[flen use_nPix*2]))
colormap(gray)
pause
clf
end