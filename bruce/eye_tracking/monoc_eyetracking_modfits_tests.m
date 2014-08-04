clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 86;
Expt_name = sprintf('G0%d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

anal_name = 'monoc_eyecorr_hbar_mods';
use_right_eye = false;
%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 10;
use_nPix = 36;
stim_params = NIMcreate_stim_params([flen use_nPix],dt);
Fr = 1;
bar_ori = 90;
min_trial_dur = 0.75;

beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

n_use_blocks = 5;
%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat

%%
if strcmp(Expt_name,'G093')
    include_expts = {'rls.FaXwi','rls.FaXwiXimi'};
else
    include_expts = {'rls.Fa', 'rls.FaXimi'};
end
for ii = 1:length(Expts)
    expt_names{ii} = Expts{ii}.Header.expname;
    expt_dds(ii) = Expts{ii}.Stimvals.dd;
    expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
    expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
    expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
    expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
    included_type(ii) = any(strcmp(expt_names{ii},include_expts));
end
expt_has_ds(isnan(expt_has_ds)) = 0;
expt_has_ds(expt_has_ds == -1) = 0;
expt_binoc(isnan(expt_binoc)) = 0;
% cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori & expt_sac_dir == bar_ori);
cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);

if exist('./fin_aclust_data.mat','file')
    load('./fin_aclust_data.mat');
    [n_aclust_expts,n_aclust_probes] = size(autoclust);
else
    disp('No fin_aclust_data found.');
    autoclust = [];
    n_aclust_expts = 0; n_aclust_probes = 0;
end


if strcmp(Expt_name,'G087')
    cur_block_set(cur_block_set == 15) = []; %only 6 trials and causes problems
end

if strcmp(Expt_name,'G087')
    cur_block_set(cur_block_set == 15) = []; %only 6 trials and causes problems
end
if strcmp(Expt_name,'G093')
    cur_block_set(cur_block_set ==  28) = []; %only 6 trials and causes problems
end

sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

if length(cur_block_set) > n_use_blocks
    cur_block_set = cur_block_set(1:n_use_blocks);
end
n_blocks = length(cur_block_set);
%% load overall su data
cur_expt_id = find(su_data.expt_nums == Expt_num);
su_probes = find(su_data.is_su(cur_expt_id,:));
mua_probes = setdiff(1:96,su_probes); %probes with ONLY MU
aclust_probenums = [autoclust(cur_block_set(1),:).probenum];
autoclust = autoclust(:,ismember(aclust_probenums,su_probes));

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
% all_Xmat = [];
all_stim_mat = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_wi = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(96,1);
all_clust_ids = cell(length(su_probes),1);
for ee = 1:length(cur_block_set);
    if ismember(ee,grayback_gs_expts)
        fprintf('Expt %d Block %d of %d; grayback GS, ori:%d\n',Expt_num,ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,imback_gs_expts)
        fprintf('Expt %d Block %d of %d; imback GS, ori:%d\n',Expt_num,ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,sim_sac_expts)
        fprintf('Expt %d Block %d of %d; SimSac, ori:%d\n',Expt_num,ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    else
        fprintf('Expt %d Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_num,ee,length(cur_block_set));
    end
    cur_block = cur_block_set(ee);
    fname = sprintf('Expt%dClusterTimesDetails.mat',cur_block);
    load(fname);
    for cc = 1:96
        all_spk_times{cc} = cat(1,all_spk_times{cc},ClusterDetails{cc}.t');
        cur_ind = find(su_probes == cc);
        if ~isempty(cur_ind)
            all_clust_ids{cur_ind} = cat(1,all_clust_ids{cur_ind},autoclust(cur_block_set(ee),cur_ind).idx(:));
        end
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
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
    
    if strcmp(Expt_name,'G093')
        trial_wi = [Expts{cur_block}.Trials(:).wi];
        trial_wi = trial_wi(id_inds);
        all_trial_wi = cat(1,all_trial_wi,trial_wi(use_trials)');
    end
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    buffer_pix = floor((expt_npix(cur_block) - use_nPix)/2);
    cur_use_pix = (1:use_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        if n_frames > 0
            if length(cur_stim_times) == 1
                cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
                cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
            end
        end
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
%             bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
%             all_Xmat = [all_Xmat; bar_Xmat];
        all_stim_mat = [all_stim_mat; cur_stim_mat];
all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
end

%% BIN SPIKES FOR MU AND SU
mahal_thresh = su_data.mah_thresh;
all_binned_spikes = nan(length(all_t_axis),96);
su_used_blocks = false(length(cur_block_set),length(su_probes));
%for only-MU probes
for cc = 1:96
    if ~ismember(cc,su_probes) %if probe doesn't have an SU
        [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc},all_t_bin_edges);
        cur_spkhist(all_bin_edge_pts) = [];
        all_binned_spikes(:,cc) = cur_spkhist;
    end
end
%for SU probes
for ss = 1:length(su_probes)
    %for MUA
    cur_muahist = histc(all_spk_times{su_probes(ss)},all_t_bin_edges);
    cur_muahist(all_bin_edge_pts) = [];
    
    all_su_inds = find(all_clust_ids{ss} == 1);
    cur_suahist = histc(all_spk_times{su_probes(ss)}(all_su_inds),all_t_bin_edges);
    cur_suahist(all_bin_edge_pts) = [];
    
    spk_block_inds = round(interp1(all_t_axis,all_blockvec,all_spk_times{su_probes(ss)}));
    for ee = 1:length(cur_block_set);
        cur_block_inds = find(all_blockvec == ee);
        %if SU is isolated in this block
        if autoclust(cur_block_set(ee),ss).mahal_d > mahal_thresh || autoclust(cur_block_set(ee),ss).man_code == 4
            su_used_blocks(ee,ss) = true;
            all_binned_spikes(cur_block_inds,su_probes(ss)) = cur_suahist(cur_block_inds);
        else %otherwise use MUA
            all_binned_spikes(cur_block_inds,su_probes(ss)) = cur_muahist(cur_block_inds);
        end
    end
end
%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
if strcmp(Expt_name,'G093')
    un_wi_vals = unique(all_trial_wi);
    use_wi_trials = find(all_trial_wi == un_wi_vals(2));
    used_inds = used_inds(ismember(all_trialvec(used_inds),use_wi_trials));
end
NT = length(used_inds);

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
cd(data_dir)
emfile = ['jbe' Expt_name '.em.mat'];
load(emfile);

all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
all_eye_blockvec = [];
eye_smooth = 3;
for ee = 1:n_blocks;
    fprintf('Loading ET data for expt %d, block %d of %d\n',Expt_num,ee,n_blocks);
    cur_set = find(all_blockvec==ee);
    if ~isempty(cur_set)
        [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_t_axis(cur_set([1 end])));
        
        eye_dt = Expt.Header.CRrates(1);
        eye_fs = 1/eye_dt;
        lEyeXY = eye_vals_interp(:,1:2);
        rEyeXY = eye_vals_interp(:,3:4);
        
        %slight smoothing before computing speed
        sm_avg_eyepos = lEyeXY; eye_vel = lEyeXY; %initialization
        sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),eye_smooth);
        sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),eye_smooth);
        eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
        eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
        
        eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
        
        all_eye_vals = [all_eye_vals; lEyeXY rEyeXY];
        all_eye_speed = [all_eye_speed; eye_speed];
        all_eye_ts = [all_eye_ts; eye_ts_interp'];
        all_eye_blockvec = [all_eye_blockvec; ee*ones(size(eye_speed))];
    end
end

back_pts = 1 + find(diff(all_eye_ts) <= 0);
double_samples = [];
for i = 1:length(back_pts)
    next_forward = find(all_eye_ts > all_eye_ts(back_pts(i)-1),1,'first');
    double_samples = [double_samples back_pts(i):next_forward];
end
all_eye_ts(double_samples) = [];
all_eye_speed(double_samples) = [];
all_eye_vals(double_samples,:) = [];

interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_t_axis);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%% BIAS-CORRECTION ON EYE-TRACKING DATA AND REFINE USED DATA TO EXCLUDE DATA WITH EYES OUTSIDE FIXATION WINDOW

%correct for median offsets
corrected_interp_eyevals = interp_eye_vals;
corrected_interp_eyevals = bsxfun(@minus,corrected_interp_eyevals,nanmedian(corrected_interp_eyevals(used_inds,:)));
corrected_eyevals = all_eye_vals;
corrected_eyevals = bsxfun(@minus,corrected_eyevals,nanmedian(interp_eye_vals(used_inds,:)));

%now correct for dependence of orthoganol position on parallel position
par_thresh = 4;
orth_thresh = 1.5;
if bar_ori == 0
    if use_right_eye
    fit_inds = used_inds(abs(corrected_interp_eyevals(used_inds,4)) < orth_thresh & ...
        abs(corrected_interp_eyevals(used_inds,3)) < par_thresh);
    b = robustfit(corrected_interp_eyevals(fit_inds,3),corrected_interp_eyevals(fit_inds,4));
    pred_eyevals = corrected_interp_eyevals(used_inds,3)*b(2) + b(1);
    corrected_interp_eyevals(used_inds,4) = corrected_interp_eyevals(used_inds,4) - pred_eyevals;
    pred_eyevals = corrected_eyevals(:,3)*b(2) + b(1);
    corrected_eyevals(:,4) = corrected_eyevals(:,4) - pred_eyevals;
    end
    
    fit_inds = used_inds(abs(corrected_interp_eyevals(used_inds,2)) < orth_thresh & ...
        abs(corrected_interp_eyevals(used_inds,1)) < par_thresh);
    b = robustfit(corrected_interp_eyevals(fit_inds,1),corrected_interp_eyevals(fit_inds,2));
    pred_eyevals = corrected_interp_eyevals(used_inds,1)*b(2) + b(1);
    corrected_interp_eyevals(used_inds,2) = corrected_interp_eyevals(used_inds,2) - pred_eyevals;
    pred_eyevals = corrected_eyevals(:,1)*b(2) + b(1);
    corrected_eyevals(:,2) = corrected_eyevals(:,2) - pred_eyevals;    
else
    if use_right_eye
    fit_inds = used_inds(abs(corrected_interp_eyevals(used_inds,3)) < orth_thresh & ...
        abs(corrected_interp_eyevals(used_inds,4)) < par_thresh);
    b = robustfit(corrected_interp_eyevals(fit_inds,4),corrected_interp_eyevals(fit_inds,3));
    pred_eyevals = corrected_interp_eyevals(used_inds,4)*b(2) + b(1);
    corrected_interp_eyevals(used_inds,3) = corrected_interp_eyevals(used_inds,3) - pred_eyevals;
    pred_eyevals = corrected_eyevals(:,4)*b(2) + b(1);
    corrected_eyevals(:,3) = corrected_eyevals(:,3) - pred_eyevals;
    end
    
    fit_inds = used_inds(abs(corrected_interp_eyevals(used_inds,1)) < orth_thresh & ...
        abs(corrected_interp_eyevals(used_inds,2)) < par_thresh);
    b = robustfit(corrected_interp_eyevals(fit_inds,2),corrected_interp_eyevals(fit_inds,1));
    pred_eyevals = corrected_interp_eyevals(used_inds,2)*b(2) + b(1);
    corrected_interp_eyevals(used_inds,1) = corrected_interp_eyevals(used_inds,1) - pred_eyevals;
    pred_eyevals = corrected_eyevals(:,2)*b(2) + b(1);
    corrected_eyevals(:,1) = corrected_eyevals(:,1) - pred_eyevals;
end
if bar_ori == 0
   outside_window_left = find(abs(corrected_interp_eyevals(used_inds,1)) > par_thresh | ...
       abs(corrected_interp_eyevals(used_inds,2)) > orth_thresh);
   outside_window_right = find(abs(corrected_interp_eyevals(used_inds,3)) > par_thresh | ...
       abs(corrected_interp_eyevals(used_inds,4)) > orth_thresh);
else
   outside_window_left = find(abs(corrected_interp_eyevals(used_inds,2)) > par_thresh | ...
       abs(corrected_interp_eyevals(used_inds,1)) > orth_thresh);
   outside_window_right = find(abs(corrected_interp_eyevals(used_inds,4)) > par_thresh | ...
       abs(corrected_interp_eyevals(used_inds,3)) > orth_thresh);    
end
if use_right_eye
outside_window = union(outside_window_left,outside_window_right);
else
    outside_window = outside_window_left;
end

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
Xblock = zeros(length(all_stim_times),length(cur_block_set));
for i = 1:length(cur_block_set)
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end
block_mean_rates = nan(length(cur_block_set),96);
block_n_spikes = nan(length(cur_block_set),96);
for ee = 1:length(cur_block_set)
    cur_block_inds = used_inds(all_blockvec(used_inds) == ee);
    block_mean_rates(ee,:) = mean(all_binned_spikes(cur_block_inds,:));
    block_n_spikes(ee,:) = sum(all_binned_spikes(cur_block_inds,:));
end

%%
all_Xmat = create_time_embedding(all_stim_mat,stim_params);

sp_dx = 0.05;
max_shift = 20;
if bar_ori == 0
    measured_eyecors = round(corrected_interp_eyevals(:,2)/sp_dx);
elseif bar_ori == 90
    measured_eyecors = round(corrected_interp_eyevals(:,1)/sp_dx);
end
measured_eyecors(measured_eyecors < -max_shift) = -max_shift;
measured_eyecors(measured_eyecors > max_shift) = max_shift;
% measured_eyecors(measured_eyecors < -max_shift) = 0;
% measured_eyecors(measured_eyecors > max_shift) = 0;

all_stim_mat_cor = all_stim_mat;
for ii = 1:NT
    all_stim_mat_cor(used_inds(ii),:) = shift_matrix_Nd(all_stim_mat(used_inds(ii),:),...
        -measured_eyecors(used_inds(ii)),2);
end
all_Xmat_cor = create_time_embedding(all_stim_mat_cor,stim_params);
%% SUA STA/STC ANALYSIS
nneg = 5; npos = 5;
stc_thresh = -7.5e-4;

max_squared_filts = 2;

nmm_stim_params = NMMcreate_stim_params([flen use_nPix],dt);
nmm_stim_params(2) = NMMcreate_stim_params([length(cur_block_set),1],dt);

lambda_d2XT = 500;
lambda_L1 = 50;
block_L2 = 1;
sub_samp_fac = 10;
silent = 1;
optim_p.optTol = 1e-5; optim_p.progTol = 1e-9;

for ss = 1:length(su_probes)
    
    fprintf('Computing STC for SU %d of %d\n',ss,length(su_probes));
    
    cur_used_blocks = find(su_used_blocks(:,ss)); %blocks when SU
    cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
    cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
    cur_NT = length(cur_used_inds);
    
    if ~isempty(cur_used_inds)
        Robs = all_binned_spikes(cur_used_inds,su_probes(ss));
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
 
        
        
       spike_cond_stim = all_Xmat_cor(cur_used_inds(spikebins),:);
        sta_cor      = mean(spike_cond_stim) - mean(all_Xmat_cor(cur_used_inds,:));
        sta_cor = sta_cor/norm(sta_cor);
        proj_mat = sta_cor'/(sta_cor*sta_cor')*sta_cor;
        stim_proj = all_Xmat_cor(cur_used_inds,:) - all_Xmat_cor(cur_used_inds,:)*proj_mat;
        % stim_proj = stim_emb;
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs_cor  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs_cor  = stcs_cor(:,end:-1:1);
        
        sua_data(ss).block_avg_rates = block_mean_rates(cur_used_blocks,su_probes(ss));
        sua_data(ss).sta = sta;
        sua_data(ss).stcs = stcs;
        sua_data(ss).sta_cor = sta_cor;
        sua_data(ss).stcs_cor = stcs_cor;
        sua_data(ss).evecs = diag(evals);

        
        
        %%
%         n_squared_filts = min(sua_data(ss).npos_stc,max_squared_filts);
%         fprintf('Fitting NIM with %d exc squared\n',n_squared_filts);
%         mod_signs = [1 ones(1,n_squared_filts) 1];
%         NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
%         Xtargets = [ones(1,n_squared_filts+1) 2];
%         d2XT = [lambda_d2XT repmat(sqrt(lambda_d2XT),1,n_squared_filts) 0];
%         L2 = [zeros(1,n_squared_filts+1) block_L2];
%         reg_params = NMMcreate_reg_params('lambda_d2XT',d2XT,'lambda_L2',L2);        
%         init_filts = cell(length(mod_signs),1); init_filts{end} = zeros(size(Xblock,2),1);
%         fit0 = NMMinitialize_model(nmm_stim_params,mod_signs,NL_types,reg_params,Xtargets,init_filts); %initialize NIM
%         sub_sample = randperm(cur_NT);
%         sub_sample = sub_sample(1:round(cur_NT/sub_samp_fac));        
%         cur_Xmat{1} = all_Xmat(cur_used_inds(sub_sample),:); cur_Xmat{2} = Xblock(cur_used_inds(sub_sample),:);
%         fit0 = NMMfit_filters(fit0,Robs(sub_sample),cur_Xmat,[],[],silent,optim_p); %fit stimulus filters
%         
%         fit0 = NMMadjust_regularization(fit0,1:(1+n_squared_filts),'lambda_L1',[lambda_L1 sqrt(lambda_L1) sqrt(lambda_L1)]);
%         cur_Xmat{1} = all_Xmat(cur_used_inds,:); cur_Xmat{2} = Xblock(cur_used_inds,:);
%         fit0 = NMMfit_filters(fit0,Robs,cur_Xmat,[],[],silent); %fit stimulus filters
        
%         sua_data(ss).nimFit = fit0;
        
    else
        sua_data(ss).stc_use = false;
    end
end

%% MODEL FITS FOR MUA
for ss = 1:96
    
    fprintf('Computing STC for MU %d of %d\n',ss,96);
    su_probe_ind = find(su_probes == ss);
    if ~isempty(su_probe_ind)
        cur_used_blocks = find(su_used_blocks(:,su_probe_ind) == 0); %blocks when NO SU
        cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
    else
        cur_used_blocks = 1:length(cur_block_set); %blocks when NO SU
        cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
    end
    cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
    cur_NT = length(cur_used_inds);
    
    if ~isempty(cur_used_inds)
        Robs = all_binned_spikes(cur_used_inds,ss);
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
        
        mua_data(ss).block_avg_rates = block_mean_rates(cur_used_blocks,ss);
        mua_data(ss).sta = sta;
        mua_data(ss).stcs = stcs;
        mua_data(ss).evecs = diag(evals);
        cur_evec_diff = diff(flipud(mua_data(ss).evecs));
        mua_data(ss).npos_stc = find(cur_evec_diff(1:flen*use_nPix/2) < stc_thresh,1,'last');
        mua_data(ss).nneg_stc = find(cur_evec_diff(end:-1:end-flen*use_nPix/2) < stc_thresh,1,'last');
        if isempty(mua_data(ss).npos_stc); mua_data(ss).npos_stc = 0; end;
        if isempty(mua_data(ss).nneg_stc); mua_data(ss).nneg_stc = 0; end;
        mua_data(ss).stc_use = true;
        mua_data(ss).avg_rate = mean(Robs);
        mua_data(ss).nspks = sum(Robs);
        cur_n_grayback_blocks = sum(ismember(cur_used_blocks,grayback_gs_expts));
        cur_n_imback_blocks = sum(ismember(cur_used_blocks,imback_gs_expts));
        mua_data(ss).nbocks = [cur_n_grayback_blocks cur_n_imback_blocks];
        
        %%
        n_squared_filts = min(mua_data(ss).npos_stc,max_squared_filts);
        fprintf('Fitting NIM with %d exc squared\n',n_squared_filts);
        mod_signs = [1 ones(1,n_squared_filts) 1];
        NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
        Xtargets = [ones(1,n_squared_filts+1) 2];
        d2XT = [lambda_d2XT repmat(sqrt(lambda_d2XT),1,n_squared_filts) 0];
        L2 = [zeros(1,n_squared_filts+1) block_L2];
        reg_params = NMMcreate_reg_params('lambda_d2XT',d2XT,'lambda_L2',L2);        
        init_filts = cell(length(mod_signs),1); init_filts{end} = zeros(size(Xblock,2),1);
        fit0 = NMMinitialize_model(nmm_stim_params,mod_signs,NL_types,reg_params,Xtargets,init_filts); %initialize NIM
        sub_sample = randperm(cur_NT);
        sub_sample = sub_sample(1:round(cur_NT/sub_samp_fac));        
        cur_Xmat{1} = all_Xmat(cur_used_inds(sub_sample),:); cur_Xmat{2} = Xblock(cur_used_inds(sub_sample),:);
        fit0 = NMMfit_filters(fit0,Robs(sub_sample),cur_Xmat,[],[],silent,optim_p); %fit stimulus filters
        
        cur_Xmat{1} = all_Xmat(cur_used_inds,:); cur_Xmat{2} = Xblock(cur_used_inds,:);
        fit0 = NMMadjust_regularization(fit0,1:(1+n_squared_filts),'lambda_L1',[lambda_L1 sqrt(lambda_L1) sqrt(lambda_L1)]);
        fit0 = NMMfit_filters(fit0,Robs,cur_Xmat,[],[],silent); %fit stimulus filters

        mua_data(ss).nimFit = fit0;
        
    else
        mua_data(ss).stc_use = false;
    end
end


%%
cd(anal_dir)
modfit_used_blocks = cur_used_blocks;
save(anal_name,'sua_data','mua_data','modfit_used_blocks');

%% VIEW MODEL FITS COMPARED TO STA/STC
f1 = figure();
f2 = figure();
for ss = 1:length(su_probes)
    fprintf('SU %d of %d\n',ss,length(su_probes));
    
    figure(f1); clf
    cur_sta = sua_data(ss).sta; ca = max(abs(cur_sta));
    subplot(3,3,1)
    imagesc(reshape(cur_sta,flen,use_nPix)); caxis([-ca ca]);
    cur_stcs = sua_data(ss).stcs; ca = max(abs(cur_stcs(:)));
    for ii = 1:3
        subplot(3,3,3+ii)
        imagesc(reshape(cur_stcs(:,ii),flen,use_nPix)); caxis([-ca ca]);
    end
    for ii = 1:3
        subplot(3,3,6+ii)
        imagesc(reshape(cur_stcs(:,end-nneg+ii),flen,use_nPix)); caxis([-ca ca]);
    end
 
        figure(f2); clf
    cur_sta = sua_data(ss).sta_cor; ca = max(abs(cur_sta));
    subplot(3,3,1)
    imagesc(reshape(cur_sta,flen,use_nPix)); caxis([-ca ca]);
    cur_stcs = sua_data(ss).stcs_cor; ca = max(abs(cur_stcs(:)));
    for ii = 1:3
        subplot(3,3,3+ii)
        imagesc(reshape(cur_stcs(:,ii),flen,use_nPix)); caxis([-ca ca]);
    end
    for ii = 1:3
        subplot(3,3,6+ii)
        imagesc(reshape(cur_stcs(:,end-nneg+ii),flen,use_nPix)); caxis([-ca ca]);
    end

%     kmat = [sua_data(ss).nimFit.mods(1:end-1).filtK];
%     figure(f2); clf
%     subplot(2,2,1)
%     imagesc(reshape(kmat(:,1),flen,use_nPix));
%     ca = max(abs(kmat(:,1))); caxis([-ca ca]);
%     for ii = 1:(size(kmat,2)-1)
%         subplot(2,2,2+ii)
%         imagesc(reshape(kmat(:,ii+1),flen,use_nPix));
%         ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
%     end
    pause
    
end
%% FOR MUA
f1 = figure();
f2 = figure();
for ss = 1:96
    if mua_data(ss).stc_use
    fprintf('MU %d of %d\n',ss,96);
    
    figure(f1); clf
    cur_sta = mua_data(ss).sta; ca = max(abs(cur_sta));
    subplot(3,3,1)
    imagesc(reshape(cur_sta,flen,use_nPix)); caxis([-ca ca]);
    cur_stcs = mua_data(ss).stcs; ca = max(abs(cur_stcs(:)));
    for ii = 1:3
        subplot(3,3,3+ii)
        imagesc(reshape(cur_stcs(:,ii),flen,use_nPix)); caxis([-ca ca]);
    end
    for ii = 1:3
        subplot(3,3,6+ii)
        imagesc(reshape(cur_stcs(:,end-3+ii),flen,use_nPix)); caxis([-ca ca]);
    end
    
    kmat = [mua_data(ss).nimFit.mods(1:end-1).filtK];
    figure(f2); clf
    subplot(2,2,1)
    imagesc(reshape(kmat(:,1),flen,use_nPix));
    ca = max(abs(kmat(:,1))); caxis([-ca ca]);
    for ii = 1:(size(kmat,2)-1)
        subplot(2,2,2+ii)
        imagesc(reshape(kmat(:,ii+1),flen,use_nPix));
        ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
    end
    pause
    end
end

