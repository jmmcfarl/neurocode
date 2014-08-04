clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 88;
Expt_name = sprintf('G0%d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

anal_name = 'binoc_eyecorr_hbar';
mod_data_name = 'binoc_eyecorr_hbar_mods';
eye_data_name = 'binoc_eyecorr_hbar_eyeonly';
use_right_eye = false;
%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 10;
use_nPix = 24;
full_nPix = 36;
init_stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;
bar_ori = 0;
min_trial_dur = 0.75;

beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat

%%
if strcmp(Expt_name,'G093')
    include_expts = {'rls.FaXwi','rls.FaXwiXimi'};
else
    include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaRC'};
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
% cur_block_set = find(included_type & expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori & expt_sac_dir == bar_ori);
cur_block_set = find(included_type & expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);

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

n_blocks = length(cur_block_set);

sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

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
all_left_stim_mat = [];
all_right_stim_mat = [];
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
for ee = 1:n_blocks;
    if ismember(ee,grayback_gs_expts)
        fprintf('Expt %d Block %d of %d; grayback GS, ori:%d\n',Expt_num,ee,n_blocks,expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,imback_gs_expts)
        fprintf('Expt %d Block %d of %d; imback GS, ori:%d\n',Expt_num,ee,n_blocks,expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,sim_sac_expts)
        fprintf('Expt %d Block %d of %d; SimSac, ori:%d\n',Expt_num,ee,n_blocks,expt_bar_ori(cur_block_set(ee)));
    else
        fprintf('Expt %d Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_num,ee,n_blocks);
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
    buffer_pix = floor((expt_npix(cur_block) - full_nPix)/2);
    cur_use_pix = (1:full_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start/1e4;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        if n_frames > 0
            if length(cur_stim_times) == 1
                cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
            end
        end
                 cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
       cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
%             cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
%             bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
            cur_lstim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            cur_rstim_mat = double(right_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
%             left_Xmat = create_time_embedding(cur_lstim_mat,stim_params);
%             right_Xmat = create_time_embedding(cur_rstim_mat,stim_params);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
%             all_left_Xmat = [all_left_Xmat; left_Xmat];
%             all_right_Xmat = [all_right_Xmat; right_Xmat];
            all_left_stim_mat = [all_left_stim_mat; cur_lstim_mat];
            all_right_stim_mat = [all_right_stim_mat; cur_rstim_mat];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
end

%% select submatrix with central pixels
[Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix));

%% BIN SPIKES FOR MU AND SU
mahal_thresh = su_data.mah_thresh;
all_binned_spikes = nan(length(all_t_axis),96);
su_used_blocks = false(n_blocks,length(su_probes));
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
    for ee = 1:n_blocks;
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

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
Xblock = zeros(length(all_stim_times),n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

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
fprintf('Removing %d of %d samples, %.4f where eyes are outside the FW\n',length(outside_window),length(used_inds),length(outside_window)/length(used_inds));
used_inds(outside_window) = [];
NT = length(used_inds);

%% DETECT SACCADES
sac_thresh = 10;
peak_sig = [0; diff(sign(diff(all_eye_speed))); 0];
saccade_inds = find(peak_sig == -2 & all_eye_speed > sac_thresh);

peri_thresh = 3; %threshold eye speed for defining saccade boundary inds
thresh_cross_up = 1 + find(all_eye_speed(1:end-1) < peri_thresh & all_eye_speed(2:end) >= peri_thresh);
thresh_cross_down = 1 + find(all_eye_speed(1:end-1) >= peri_thresh & all_eye_speed(2:end) < peri_thresh);
sac_start_inds = nan(size(saccade_inds));
sac_stop_inds = nan(size(saccade_inds));
for ii = 1:length(saccade_inds)
    next_tc = find(thresh_cross_down > saccade_inds(ii),1,'first');
    if ~isempty(next_tc)
        sac_stop_inds(ii) = thresh_cross_down(next_tc);
    end
    prev_tc = find(thresh_cross_up < saccade_inds(ii),1,'last');
    if ~isempty(prev_tc)
        sac_start_inds(ii) = thresh_cross_up(prev_tc);
    end
    
end

%get rid of double-peaks
min_isi = 0.05; max_isi = Inf;
isis = [Inf; diff(sac_start_inds)]/eye_fs;
bad_isis = (isis < min_isi | isis > max_isi);
bad_sacs = find(isnan(sac_stop_inds) | isnan(sac_start_inds) | bad_isis);
saccade_inds(bad_sacs) = []; isis(bad_sacs) = []; sac_start_inds(bad_sacs) = []; sac_stop_inds(bad_sacs) = [];

saccade_times = all_eye_ts(saccade_inds);
sac_start_times = all_eye_ts(sac_start_inds);
sac_stop_times = all_eye_ts(sac_stop_inds);
sac_durs = sac_stop_times - sac_start_times;

% sac_dbuff = round(0.005/eye_dt);
% pre_inds = saccade_inds - sac_dbuff;
% pre_inds(pre_inds < 1) = 1;
% sac_pre_pos = all_eye_vals(pre_inds,:);
% post_inds = saccade_inds + sac_dbuff;
% post_inds(post_inds > length(all_eye_ts)) = length(all_eye_ts);
% sac_post_pos = all_eye_vals(post_inds,:);

sac_dbuff = round(0.05/eye_dt);
sac_pre_pos = nan(length(saccade_inds),4);
sac_post_pos = nan(length(saccade_inds),4);
for ii = 1:length(saccade_inds)
    pre_inds = (sac_start_inds(ii) - sac_dbuff):sac_start_inds(ii);
    pre_inds(pre_inds < 1) = 1;
    sac_pre_pos(ii,:) = median(all_eye_vals(pre_inds,:),1);
    post_inds = sac_stop_inds(ii):(sac_stop_inds(ii) + sac_dbuff);
    post_inds(post_inds > length(all_eye_ts)) = length(all_eye_ts);
    sac_post_pos(ii,:) = median(all_eye_vals(post_inds,:),1);
end

%use only left-eye signal here
sac_delta_pos = sac_post_pos(:,1:2) - sac_pre_pos(:,1:2);
sac_amps = sqrt(sum(sac_delta_pos.^2,2));
sac_dirs = atan2(sac_delta_pos(:,2),sac_delta_pos(:,1));

temp = ones(length(saccade_times),1);
saccades = struct('peak_time',mat2cell(saccade_times,temp),'start_time',mat2cell(sac_start_times,temp),...
    'stop_time',mat2cell(sac_stop_times,temp),'isi',mat2cell(isis,temp),...
    'duration',mat2cell(sac_durs,temp),'amplitude',mat2cell(sac_amps,temp),'direction',mat2cell(sac_dirs,temp),...
    'pre_Lx',mat2cell(sac_pre_pos(:,1),temp),'post_Lx',mat2cell(sac_post_pos(:,1),temp),...
    'pre_Ly',mat2cell(sac_pre_pos(:,2),temp),'post_Ly',mat2cell(sac_post_pos(:,2),temp),...
    'pre_Rx',mat2cell(sac_pre_pos(:,3),temp),'post_Rx',mat2cell(sac_post_pos(:,3),temp),...
    'pre_Ry',mat2cell(sac_pre_pos(:,4),temp),'post_Ry',mat2cell(sac_post_pos(:,4),temp));

sac_start_times = [saccades(:).start_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > dt);
sac_start_times(bad_sacs) = [];
saccades(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];

%% DEFINE POINTS TO ALLOW RAPID CHANGES (TRIAL STARTS AFTER SACCADES)
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];

saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps(used_saccade_set) < 1;
big_sacs = find(~is_micro);
micro_sacs = find(is_micro);

%push the effects of saccades forward in time
sac_shift = round(0.05/dt);

%define times when to resort to independent prior
use_prior = false(size(used_inds));
use_prior(trial_start_inds) = true;
for i = 1:length(saccade_start_inds)
    next_trial = trial_start_inds(find(trial_start_inds > saccade_start_inds(i),1,'first'));
    %mark a point (forward in time) as a saccade if it occurs within the
    %same trial as the saccade
    if next_trial > saccade_start_inds(i) + sac_shift
        use_prior(saccade_start_inds(i) + sac_shift) = 1;
    end
end

sac_backlag = round(0.05/dt);
sac_forlag = round(0.3/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));
Xsac = zeros(NT,n_sac_bins);
Xmsac = zeros(NT,n_sac_bins);
for ii = 1:n_sac_bins
    cur_sac_target = saccade_start_inds(big_sacs) + sac_bincents(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(big_sacs(uu))) = [];
    Xsac(cur_sac_target,ii) = 1;
    
    cur_sac_target = saccade_start_inds(micro_sacs) + sac_bincents(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(micro_sacs(uu))) = [];
    Xmsac(cur_sac_target,ii) = 1;
end

%% LOAD MODEL FITS
cd(anal_dir);
load(mod_data_name);

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
all_left_Xmat = create_time_embedding(all_left_stim_mat,init_stim_params);
all_right_Xmat = create_time_embedding(all_right_stim_mat,init_stim_params);
full_X{1} = [all_left_Xmat(used_inds,use_kInds) all_right_Xmat(used_inds,use_kInds)];
full_X{2} = Xblock(used_inds,:);
full_X{3} = Xsac;
full_X{4} = Xmsac;

full_stim_params = NMMcreate_stim_params([flen 2*use_nPix],dt);
full_stim_params(2) = NMMcreate_stim_params(n_blocks);
full_stim_params(3) = NMMcreate_stim_params(n_sac_bins,dt);
full_stim_params(4) = NMMcreate_stim_params(n_sac_bins,dt);

null_stim_params = full_stim_params(2:end);

sac_d2t = 100;
block_L2 = 1;
null_reg_params = NMMcreate_reg_params('lambda_d2T',[0; sac_d2t; sac_d2t],'lambda_L2',[block_L2; 0; 0]);

%create L2 mat for left eye
L2_params = create_L2_params([],[1 flen*use_nPix],[flen use_nPix],2,3,[Inf 0]);
L2_mat = generate_L2_mat(L2_params,2*flen*use_nPix);
%add L2 mat for sine term
L2_params = create_L2_params([],[flen*use_nPix+1 2*flen*use_nPix],[flen use_nPix],2,3,[Inf 0]);
L2_mat = L2_mat + generate_L2_mat(L2_params,2*flen*use_nPix);

tot_nUnits = length(su_probes) + 96;
all_mod_SU = zeros(tot_nUnits,1);

%% for MUs
for ss = 1:96
    fprintf('Computing base LLs for MU %d of %d\n',ss,96);
    if mua_data(ss).stc_use
        su_probe_ind = find(su_probes == ss);
        if ~isempty(su_probe_ind)
            cur_used_blocks = find(su_used_blocks(:,su_probe_ind) == 0); %blocks when NO SU
            cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        else
            cur_used_blocks = 1:n_blocks; %blocks when NO SU
            cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        end
        cur_used_inds = find(ismember(used_inds,cur_poss_inds));
        cur_NT = length(cur_used_inds);
        if ~isempty(cur_used_inds)
            Robs = all_binned_spikes(used_inds(cur_used_inds),ss);
            cur_X = get_Xcell_tInds(full_X,cur_used_inds);
            
            null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
            null_mod = NMMfit_filters(null_mod,Robs,cur_X(2:end));
            all_nullmod(ss) = null_mod;
            
            new_fit = mua_data(ss).nimFit;
            Xtargs = [new_fit.mods(:).Xtarget];
            stim_targs = find(Xtargs == 1);
            for ii = 1:length(new_fit.mods)
               new_fit.mods(ii).reg_params.boundary_conds = [Inf Inf Inf]; 
               new_fit.mods(ii).reg_params.mixing_prop = [1 1 1];
            end
            new_fit.stim_params = full_stim_params;
            new_fit = NMMadd_NLinput(new_fit,{'lin'},1,3,zeros(n_sac_bins,1)); %sac term
            new_fit = NMMadd_NLinput(new_fit,{'lin'},1,4,zeros(n_sac_bins,1)); %nsac term
            cur_nmods = length(new_fit.mods);
            new_fit = NMMadjust_regularization(new_fit,[cur_nmods-1 cur_nmods],'lambda_custom',[0 0]);
            new_fit = NMMadjust_regularization(new_fit,[cur_nmods-1 cur_nmods],'lambda_d2T',[sac_d2t sac_d2t]);
            Xtargs = [new_fit.mods(:).Xtarget];
            new_fit = NMMfit_filters(new_fit,Robs,cur_X,[],find(Xtargs ~= 1),1,[],L2_mat);
            all_mod_fits(ss) = new_fit;
            LL = NMMmodel_eval(new_fit,Robs,cur_X,[],L2_mat);
            all_mod_LLimp(ss) = (LL-null_mod.LL_seq(end))/log(2);
            all_mod_fits_withspkNL(ss) = NMMfit_logexp_spkNL(new_fit,Robs,cur_X);
            all_mod_nspks(ss) = sum(Robs); 
            all_mod_avgrate(ss) = mean(Robs);
        end
    else
        all_mod_LLimp(ss) = nan;
    end
end

%% for SUs
for ss = 1:length(su_probes)
    fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
    all_mod_SU(ss+96) = su_probes(ss);
    if sua_data(ss).stc_use
        cur_used_blocks = find(su_used_blocks(:,ss)); %blocks when SU
        cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        cur_used_inds = find(ismember(used_inds,cur_poss_inds));
        cur_NT = length(cur_used_inds);
        if ~isempty(cur_used_inds)
            Robs = all_binned_spikes(used_inds(cur_used_inds),su_probes(ss));
            cur_X = get_Xcell_tInds(full_X,cur_used_inds);
            
            null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
            null_mod = NMMfit_filters(null_mod,Robs,cur_X(2:end));
            all_nullmod(ss+96) = null_mod;
            
            %fit left-eye model
            new_fit = sua_data(ss).nimFit;
            new_fit.stim_params = full_stim_params;
            new_fit = NMMadd_NLinput(new_fit,{'lin'},1,3,zeros(n_sac_bins,1)); %sac term
            new_fit = NMMadd_NLinput(new_fit,{'lin'},1,4,zeros(n_sac_bins,1)); %nsac term
            cur_nmods = length(new_fit.mods);
            new_fit = NMMadjust_regularization(new_fit,[cur_nmods-1 cur_nmods],'lambda_custom',[0 0]);
            new_fit = NMMadjust_regularization(new_fit,[cur_nmods-1 cur_nmods],'lambda_L1',[0 0]);
            new_fit = NMMadjust_regularization(new_fit,[cur_nmods-1 cur_nmods],'lambda_d2T',[sac_d2t sac_d2t]);
            Xtargs = [new_fit.mods(:).Xtarget];
            for ii = 1:length(new_fit.mods)
               new_fit.mods(ii).reg_params.boundary_conds = [Inf Inf Inf]; 
               new_fit.mods(ii).reg_params.mixing_prop = [1 1 1];
            end
            n_sq_filt = sum(Xtargs==1)-1;
%             new_fit = NMMadjust_regularization(new_fit,find(Xtargs == 1),'lambda_custom',[lambda_custom ones(1,n_sq_filt)*sqrt(lambda_custom)]);
            new_fit = NMMfit_filters(new_fit,Robs,cur_X,[],find(Xtargs ~= 1),1,[],L2_mat);
%             new_fit = NMMfit_filters(new_fit,Robs,cur_X,[],[],1,[],L2_mat);
            all_mod_fits(ss+96) = new_fit;
            LL = NMMmodel_eval(new_fit,Robs,cur_X,[],L2_mat);
            all_mod_LLimp(ss+96) = (LL-null_mod.LL_seq(end))/log(2);
            all_mod_fits_withspkNL(ss+96) = NMMfit_logexp_spkNL(new_fit,Robs,cur_X);
            all_mod_nspks(ss+96) = sum(Robs);
            all_mod_avgrate(ss+96) = mean(Robs);
        end
    else
        all_mod_LLimp(ss+96) = nan;
    end
end

%% SELECT USABLE UNITS AND make Robs_mat
LL_imp_thresh = 5e-3;
min_avgrate = 2;
usable_units = find(all_mod_LLimp >= LL_imp_thresh & all_mod_avgrate/dt >= min_avgrate);
n_used_sus = sum(all_mod_SU(usable_units) ~= 0);
n_used_mus = sum(all_mod_SU(usable_units) == 0);
fprintf('Using %d SUs and %d MUs for analysis\n',n_used_sus,n_used_mus);
tr_set = usable_units;
n_tr_chs = length(tr_set);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = nan(length(used_inds),n_tr_chs);
for ss = 1:n_tr_chs
    if all_mod_SU(tr_set(ss)) == 0
        su_probe_ind = find(su_probes == tr_set(ss));
        if ~isempty(su_probe_ind)
            cur_used_blocks = find(su_used_blocks(:,su_probe_ind) == 0); %blocks when NO SU
        else
            cur_used_blocks = 1:n_blocks; %blocks when NO SU
        end
        cur_use = find(ismember(all_blockvec(used_inds),cur_used_blocks));
        Robs_mat(cur_use,ss) = all_binned_spikes(used_inds(cur_use),tr_set(ss));
    else
        su_ind = find(su_probes == all_mod_SU(tr_set(ss)));
        cur_used_blocks = find(su_used_blocks(:,su_ind)); %blocks when SU
        cur_use = find(ismember(all_blockvec(used_inds),cur_used_blocks));
        Robs_mat(cur_use,ss) = all_binned_spikes(used_inds(cur_use),all_mod_SU(tr_set(ss)));
    end
end

%% INITIALIZE TRANSITION PRIORS FOR HMM
sp_dx = 0.05;
max_shift = 12;
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);

%overall prior on shifts
eps_prior_sigma = 0.12; %0.125 start
leps_prior = -(shifts*sp_dx).^2/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior,n_shifts,1);

cdist = squareform(pdist(shifts'*sp_dx));
deps_sigma = 0.0175; %.015
lA = -cdist.^2/(2*deps_sigma^2);
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

%% PREPROCESS MODEL COMPONENTS
klen = 2*flen*use_nPix;
filt_bank = nan(n_tr_chs,klen,3);
lin_kerns = nan(n_tr_chs,n_blocks);
sac_kerns = nan(n_tr_chs,n_sac_bins);
msac_kerns = nan(n_tr_chs,n_sac_bins);
mod_spkNL_params = nan(n_tr_chs,3);
for ss = 1:n_tr_chs
    cur_Xtargs = [all_mod_fits(tr_set(ss)).mods(:).Xtarget];
    cur_k = [all_mod_fits(tr_set(ss)).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(ss,:,1:n_used_filts) = cur_k;
    mod_spkNL_params(ss,:) = all_mod_fits(tr_set(ss)).spk_NL_params;
    lin_kerns(ss,:) = all_mod_fits(tr_set(ss)).mods(cur_Xtargs == 2).filtK;
    sac_kerns(ss,:) = all_mod_fits(tr_set(ss)).mods(cur_Xtargs == 3).filtK;
    msac_kerns(ss,:) = all_mod_fits(tr_set(ss)).mods(cur_Xtargs == 4).filtK;
end
filt_bank = permute(filt_bank,[2 1 3]);
filt_bank = reshape(filt_bank,[flen 2*use_nPix n_tr_chs 3]);

left_filt_bank = filt_bank(:,1:use_nPix,:,:);
right_filt_bank = filt_bank(:,(1:use_nPix) + use_nPix,:,:);

%indicator predictions
block_out = Xblock(used_inds,:)*lin_kerns';
sac_out = Xsac*sac_kerns';
msac_out = Xmsac*msac_kerns';

%% ESTIMATE LL for each shift in each stimulus frame
cur_Xmat = [all_left_Xmat(used_inds,use_kInds) all_right_Xmat(used_inds,use_kInds)];

%precompute LL at all shifts for all units
LLs = nan(NT,length(tr_set),n_shifts);
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',xx,n_shifts);
    cur_Lfilt_shift = shift_matrix_Nd(left_filt_bank,shifts(xx),2);
    cur_Lfilt_shift = reshape(cur_Lfilt_shift,[flen*use_nPix n_tr_chs 3]);
    cur_Rfilt_shift = shift_matrix_Nd(right_filt_bank,shifts(xx),2);
    cur_Rfilt_shift = reshape(cur_Rfilt_shift,[flen*use_nPix n_tr_chs 3]);
    cur_filt_shift = cat(1,cur_Lfilt_shift,cur_Rfilt_shift);
    
    %outputs of stimulus models at current X-matrix shift
    gfuns = ones(NT,n_tr_chs);
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
    gfuns = gfuns + cur_Xmat*squeeze(cur_filt_shift(:,:,1));
    for ff = 2:3
        gfuns = gfuns + (cur_Xmat*squeeze(cur_filt_shift(:,:,ff))).^2;
    end
    
    %add contributions from extra lin kernels
    gfuns = gfuns + block_out + sac_out + msac_out;
    
    %incorporate beta
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,2)');
    
    too_large = gfuns > 50;
    pred_rate = log(1+exp(gfuns));
    pred_rate(too_large) = gfuns(too_large);
    
    %incorporate alpha
    pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(:,3)');
    
    pred_rate(pred_rate < 1e-20) = 1e-20;
    
    LLs(:,:,xx) = Robs_mat.*log(pred_rate) - pred_rate;
end

%% INFER EYE-POSITION USING ALL TR UNITS 
fprintf('Inferring eye-position using all cells\n');

frame_LLs = squeeze(nansum(LLs,2));
lalpha=zeros(NT,n_shifts);
lbeta = zeros(NT,n_shifts);
lscale=zeros(NT,1); %initialize rescaling parameters
%compute rescaled forward messages
lalpha(1,:) = leps_prior + frame_LLs(1,:);
lscale(1)=logsumexp(lalpha(1,:));
lalpha(1,:) = lalpha(1,:) - lscale(1);
for t=2:NT
    if use_prior(t)
        cur_lA = lA_tflip;
    else
        cur_lA = lA;
    end
    lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + frame_LLs(t,:);
    lscale(t) = logsumexp(lalpha(t,:));
    lalpha(t,:) = lalpha(t,:) - lscale(t);
end

%compute rescaled backward messages
lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
for t=NT-1:-1:1
    if use_prior(t+1)
        cur_lA = lA_tflip;
    else
        cur_lA = lA;
    end
    lf1 = lbeta(t+1,:) + frame_LLs(t+1,:);
    lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
end
lgamma= lalpha + lbeta;
lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));

gamma = exp(lgamma);
fpost_mean = sum(bsxfun(@times,gamma,shifts),2);
fpost_std = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - fpost_mean.^2);
[max_post,max_loc] = max(lgamma,[],2);
fshift_cor = shifts(max_loc);

%% FINALIZE AND SAVE INITIAL EYE-POSITION SEQUENCE
inferred_eyepos = fshift_cor*sp_dx;
inferred_eyepos_cor = inferred_eyepos;

%align jumps to saccade times
for i = length(saccade_start_inds):-1:1
    cur_inds = saccade_start_inds(i):(saccade_start_inds(i)+sac_shift);
    cur_inds(cur_inds > NT) = [];
    inferred_eyepos_cor(cur_inds) = inferred_eyepos(cur_inds(end));
    cur_back_inds = (saccade_start_inds(i)-(flen-sac_shift)):(saccade_start_inds(i)-1);
    cur_back_inds(cur_back_inds < 1) = [];
    inferred_eyepos_cor(cur_back_inds) = inferred_eyepos(saccade_start_inds(i));
end
inferred_eyepos_taxis = all_t_axis(used_inds);

cd(anal_dir)
save(eye_data_name,'inferred_eyepos_cor','inferred_eyepos_taxis');

%% RECONSTRUCT MAP STIMULUS
left_stim_mat_cor = all_left_stim_mat;
right_stim_mat_cor = all_right_stim_mat;
for ii = 1:NT
   left_stim_mat_cor(used_inds(ii),:) = shift_matrix_Nd(all_left_stim_mat(used_inds(ii),:),-fshift_cor(ii),2); 
   right_stim_mat_cor(used_inds(ii),:) = shift_matrix_Nd(all_right_stim_mat(used_inds(ii),:),-fshift_cor(ii),2); 
end
left_Xmat_cor = create_time_embedding(left_stim_mat_cor,init_stim_params);
right_Xmat_cor = create_time_embedding(right_stim_mat_cor,init_stim_params);
left_Xmat_cor = left_Xmat_cor(used_inds,use_kInds);
right_Xmat_cor = right_Xmat_cor(used_inds,use_kInds);

%% REFIT ALL CELLS THAT AREN'T IN POSS_XV_SET
silent = 1;
for ss = 1:length(tr_set)
    tr_cell = tr_set(ss);
    fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
    cur_unit_ind = find(tr_set == tr_cell);
    cur_uset = find(~isnan(Robs_mat(:,cur_unit_ind)));
    
    cur_X{1} = [left_Xmat_cor(cur_uset,:) right_Xmat_cor(cur_uset,:)];
    cur_X{2} = Xblock(used_inds(cur_uset),:);
    cur_X{3} = Xsac(cur_uset,:);
    cur_X{4} = Xmsac(cur_uset,:);
    
    ref_mod{1}(tr_cell) = all_mod_fits(tr_cell);
    ref_mod{1}(tr_cell) = NMMfit_filters(ref_mod{1}(tr_cell),Robs_mat(cur_uset,cur_unit_ind),...
        cur_X,[],[],silent,[],L2_mat); %fit stimulus filters
    
    ref_LL_imp(1,tr_cell) = (ref_mod{1}(tr_cell).LL_seq(end) - all_nullmod(tr_cell).LL_seq(end))/log(2);
    
    fprintf('Original: %.4f  New: %.4f\n',all_mod_LLimp(tr_cell),ref_LL_imp(1,tr_cell));

    %refit spk NL
    ref_mod_spkNL{1}(tr_cell) = NMMfit_logexp_spkNL(ref_mod{1}(tr_cell),Robs_mat(cur_uset,cur_unit_ind),cur_X);
end

%% REMOVE EYE-TRACKING BIASES
% corrected_interp_eyevals = interp_eye_vals;
% corrected_interp_eyevals = bsxfun(@minus,corrected_interp_eyevals,nanmedian(corrected_interp_eyevals(used_inds,:)));
% corrected_eyevals = all_eye_vals;
% corrected_eyevals = bsxfun(@minus,corrected_eyevals,nanmedian(interp_eye_vals(used_inds,:)));
% ver_thresh = 4;
% hor_thresh = 1;
% fit_inds = used_inds(abs(corrected_interp_eyevals(used_inds,3)) < hor_thresh & ...
%     abs(corrected_interp_eyevals(used_inds,4)) < ver_thresh);
% b = robustfit(corrected_interp_eyevals(fit_inds,4),corrected_interp_eyevals(fit_inds,3));
% pred_eyevals = corrected_interp_eyevals(used_inds,4)*b(2) + b(1);
% corrected_interp_eyevals(used_inds,3) = corrected_interp_eyevals(used_inds,3) - pred_eyevals;
% pred_eyevals = corrected_eyevals(:,4)*b(2) + b(1);
% corrected_eyevals(:,3) = corrected_eyevals(:,3) - pred_eyevals;
% 
% fit_inds = used_inds(abs(corrected_interp_eyevals(used_inds,1)) < hor_thresh & ...
%     abs(corrected_interp_eyevals(used_inds,2)) < ver_thresh);
% b = robustfit(corrected_interp_eyevals(fit_inds,2),corrected_interp_eyevals(fit_inds,1));
% pred_eyevals = corrected_interp_eyevals(used_inds,2)*b(2) + b(1);
% corrected_interp_eyevals(used_inds,1) = corrected_interp_eyevals(used_inds,1) - pred_eyevals;
% pred_eyevals = corrected_eyevals(:,2)*b(2) + b(1);
% corrected_eyevals(:,1) = corrected_eyevals(:,1) - pred_eyevals;
% 
eye_dsf = 4;
corrected_eyevals = downsample(corrected_eyevals,4);
all_eye_ts_ds = downsample(all_eye_ts,4);

% lcf = 0.02;
% [bb_f,aa_f] = butter(2,lcf/(eye_fs/2),'high');
% all_filt_eyevals = filtfilt(bb_f,aa_f,all_eye_vals);
% all_eye_ts_ds = downsample(all_eye_ts,4);
% all_filt_eyevals_ds = downsample(all_filt_eyevals,4);
% % all_filt_eyevals_ds = downsample(all_eye_vals,4);
%% COMPUTE DELTA POS CORRELATIONS
max_dur = 0.1;
sac_durs = [saccades(:).duration];
is_blink = sac_durs(used_saccade_set) > max_dur;

measured_pre_Ly = [saccades(:).pre_Ly];
measured_post_Ly = [saccades(:).post_Ly];
measured_pre_Lx = [saccades(:).pre_Lx];
measured_post_Lx = [saccades(:).post_Lx];
measured_pre_Ry = [saccades(:).pre_Ry];
measured_post_Ry = [saccades(:).post_Ry];
measured_pre_Rx = [saccades(:).pre_Rx];
measured_post_Rx = [saccades(:).post_Rx];
if bar_ori == 0
measured_delta_pos = [measured_post_Ly(used_saccade_set)' - measured_pre_Ly(used_saccade_set)'];
measured_delta_pos = [measured_delta_pos measured_post_Ry(used_saccade_set)' - measured_pre_Ry(used_saccade_set)'];
elseif bar_ori == 90
measured_delta_pos = [measured_post_Lx(used_saccade_set)' - measured_pre_Lx(used_saccade_set)'];
measured_delta_pos = [measured_delta_pos measured_post_Rx(used_saccade_set)' - measured_pre_Rx(used_saccade_set)'];
end
inferred_pos = fpost_mean;
% inferred_pos = mean(post_mean,2);
inferred_pre_pos = inferred_pos(saccade_start_inds-1);
inferred_post_pos = inferred_pos(saccade_start_inds + 6);
inferred_delta_pos = (inferred_post_pos - inferred_pre_pos)*sp_dx;

sac_delta_Lcorr =  corr(measured_delta_pos(~is_blink,1),inferred_delta_pos(~is_blink),'type','spearman');
sac_delta_Rcorr =  corr(measured_delta_pos(~is_blink,2),inferred_delta_pos(~is_blink),'type','spearman');
sac_delta_Acorr =  corr(mean(measured_delta_pos(~is_blink,:),2),inferred_delta_pos(~is_blink),'type','spearman');
sac_delta_LRcorr = corr(measured_delta_pos(~is_blink,1),measured_delta_pos(~is_blink,2),'type','spearman');


%% PLOT INFERRED MEASURED COMPARISON
inferred_pos = fpost_mean;
inferred_pre_pos = inferred_pos(saccade_start_inds-1);
inferred_post_pos = inferred_pos(saccade_start_inds + 6);
inferred_delta_pos = (inferred_post_pos - inferred_pre_pos)*sp_dx;
close all
hold on
plot(all_eye_ts_ds,corrected_eyevals(:,1),'r')
plot(all_eye_ts_ds,corrected_eyevals(:,3),'b')
% plot(all_t_axis(used_inds),inferred_eyepos_cor,'k.-')
% % stem(all_t_axis(used_inds(saccade_start_inds)),inferred_delta_pos,'ko','markersize',8,'linewidth',2);
% % stem(all_t_axis(used_inds(saccade_start_inds))-0.025,measured_delta_pos(:,1),'ro','markersize',8,'linewidth',2);
% % stem(all_t_axis(used_inds(saccade_start_inds))+0.025,measured_delta_pos(:,2),'o','markersize',8,'linewidth',2);
% 
% plot(all_t_axis(used_inds),left_fpost_mean*sp_dx,'k.-')
% plot(all_t_axis(used_inds),right_fpost_mean*sp_dx,'g.-')
plot(all_t_axis(used_inds),fpost_mean*sp_dx,'m.-')
ylim([-0.55 0.55]);
% 
%% NOW REFIT SEPARATE LEFT AND RIGHT EYE MODELS
half_stim_params = NMMcreate_stim_params([flen use_nPix],dt);
lambda_custom = 2000;
for ss = 1:length(tr_set)
    tr_cell = tr_set(ss);
    fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
    cur_unit_ind = find(tr_set == tr_cell);
    cur_uset = find(~isnan(Robs_mat(:,cur_unit_ind)));
    
    cur_X{1} = left_Xmat_cor(cur_uset,:);
    cur_X{2} = Xblock(used_inds(cur_uset),:);
    cur_X{3} = Xsac(cur_uset,:);
    cur_X{4} = Xmsac(cur_uset,:);
    
    Xtargs = [all_mod_fits(tr_cell).mods(:).Xtarget];
    left_ref_mod{1}(tr_cell) = all_mod_fits(tr_cell);
    stim_filts = find(Xtargs == 1);
    for ii =stim_filts
       left_ref_mod{1}(tr_cell).mods(ii).filtK = left_ref_mod{1}(tr_cell).mods(ii).filtK(1:flen*use_nPix);
    end
    left_ref_mod{1}(tr_cell) = NMMadjust_regularization(left_ref_mod{1}(tr_cell),find(Xtargs == 1),'lambda_custom',zeros(1,sum(Xtargs==1)));
    left_ref_mod{1}(tr_cell) = NMMadjust_regularization(left_ref_mod{1}(tr_cell),find(Xtargs == 1),'lambda_d2XT',lambda_custom*ones(1,sum(Xtargs==1)));
    left_ref_mod{1}(tr_cell) = NMMadjust_regularization(left_ref_mod{1}(tr_cell),find(Xtargs == 1),'lambda_L1',50*ones(1,sum(Xtargs==1)));
    left_ref_mod{1}(tr_cell).stim_params(1) = half_stim_params;
    
    left_ref_mod{1}(tr_cell) = NMMfit_filters(left_ref_mod{1}(tr_cell),Robs_mat(cur_uset,cur_unit_ind),...
        cur_X,[],[],silent); %fit stimulus filters
    left_ref_LL_imp(1,tr_cell) = (left_ref_mod{1}(tr_cell).LL_seq(end) - all_nullmod(tr_cell).LL_seq(end))/log(2);
    %refit spk NL
    left_ref_mod_spkNL{1}(tr_cell) = NMMfit_logexp_spkNL(left_ref_mod{1}(tr_cell),Robs_mat(cur_uset,cur_unit_ind),cur_X);
   
    %NOW DO RIGHT EYE
    cur_X{1} = right_Xmat_cor(cur_uset,:);
    right_ref_mod{1}(tr_cell) = all_mod_fits(tr_cell);
    for ii =stim_filts
       right_ref_mod{1}(tr_cell).mods(ii).filtK = right_ref_mod{1}(tr_cell).mods(ii).filtK((flen*use_nPix+1):end);
    end
    right_ref_mod{1}(tr_cell) = NMMadjust_regularization(right_ref_mod{1}(tr_cell),find(Xtargs == 1),'lambda_custom',zeros(1,sum(Xtargs==1)));
    right_ref_mod{1}(tr_cell) = NMMadjust_regularization(right_ref_mod{1}(tr_cell),find(Xtargs == 1),'lambda_d2XT',lambda_custom*ones(1,sum(Xtargs==1)));
    right_ref_mod{1}(tr_cell) = NMMadjust_regularization(right_ref_mod{1}(tr_cell),find(Xtargs == 1),'lambda_L1',50*ones(1,sum(Xtargs==1)));
    right_ref_mod{1}(tr_cell).stim_params(1) = half_stim_params;
    
    right_ref_mod{1}(tr_cell) = NMMfit_filters(right_ref_mod{1}(tr_cell),Robs_mat(cur_uset,cur_unit_ind),...
        cur_X,[],[],silent); %fit stimulus filters
    right_ref_LL_imp(1,tr_cell) = (right_ref_mod{1}(tr_cell).LL_seq(end) - all_nullmod(tr_cell).LL_seq(end))/log(2);
    %refit spk NL
    right_ref_mod_spkNL{1}(tr_cell) = NMMfit_logexp_spkNL(right_ref_mod{1}(tr_cell),Robs_mat(cur_uset,cur_unit_ind),cur_X);
    
end

%% FIRST DO LEFT EYE CORRECTIONS

%% SELECT USABLE UNITS AND make Robs_mat
left_usable_units = tr_set(left_ref_LL_imp(1,tr_set) >= LL_imp_thresh);
left_tr_set = left_usable_units;
n_tr_chs = length(left_tr_set);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = nan(length(used_inds),n_tr_chs);
for ss = 1:n_tr_chs
    if all_mod_SU(left_tr_set(ss)) == 0
        su_probe_ind = find(su_probes == left_tr_set(ss));
        if ~isempty(su_probe_ind)
            cur_used_blocks = find(su_used_blocks(:,su_probe_ind) == 0); %blocks when NO SU
        else
            cur_used_blocks = 1:n_blocks; %blocks when NO SU
        end
        cur_use = find(ismember(all_blockvec(used_inds),cur_used_blocks));
        Robs_mat(cur_use,ss) = all_binned_spikes(used_inds(cur_use),left_tr_set(ss));
    else
        su_ind = find(su_probes == all_mod_SU(left_tr_set(ss)));
        cur_used_blocks = find(su_used_blocks(:,su_ind)); %blocks when SU
        cur_use = find(ismember(all_blockvec(used_inds),cur_used_blocks));
        Robs_mat(cur_use,ss) = all_binned_spikes(used_inds(cur_use),all_mod_SU(left_tr_set(ss)));
    end
end


%% PREPROCESS MODEL COMPONENTS
klen = flen*use_nPix;
filt_bank = nan(n_tr_chs,klen,3);
lin_kerns = nan(n_tr_chs,n_blocks);
sac_kerns = nan(n_tr_chs,n_sac_bins);
msac_kerns = nan(n_tr_chs,n_sac_bins);
mod_spkNL_params = nan(n_tr_chs,3);
for ss = 1:n_tr_chs
    cur_Xtargs = [left_ref_mod{1}(left_tr_set(ss)).mods(:).Xtarget];
    cur_k = [left_ref_mod{1}(left_tr_set(ss)).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(ss,:,1:n_used_filts) = cur_k;
    mod_spkNL_params(ss,:) = left_ref_mod_spkNL{1}(left_tr_set(ss)).spk_NL_params;
    lin_kerns(ss,:) = left_ref_mod{1}(left_tr_set(ss)).mods(cur_Xtargs == 2).filtK;
    sac_kerns(ss,:) = left_ref_mod{1}(left_tr_set(ss)).mods(cur_Xtargs == 3).filtK;
    msac_kerns(ss,:) = left_ref_mod{1}(left_tr_set(ss)).mods(cur_Xtargs == 4).filtK;
end
filt_bank = permute(filt_bank,[2 1 3]);
filt_bank = reshape(filt_bank,[flen use_nPix n_tr_chs 3]);

%indicator predictions
block_out = Xblock(used_inds,:)*lin_kerns';
sac_out = Xsac*sac_kerns';
msac_out = Xmsac*msac_kerns';

%% ESTIMATE LL for each shift in each stimulus frame
cur_Xmat = all_left_Xmat(used_inds,use_kInds);

%precompute LL at all shifts for all units
LLs = nan(NT,length(left_tr_set),n_shifts);
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',xx,n_shifts);
    cur_filt_shift = shift_matrix_Nd(filt_bank,shifts(xx),2);
    cur_filt_shift = reshape(cur_filt_shift,[klen n_tr_chs 3]);
    
    %outputs of stimulus models at current X-matrix shift
    gfuns = ones(NT,n_tr_chs);
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
    gfuns = gfuns + cur_Xmat*squeeze(cur_filt_shift(:,:,1));
    for ff = 2:3
        gfuns = gfuns + (cur_Xmat*squeeze(cur_filt_shift(:,:,ff))).^2;
    end
    
    %add contributions from extra lin kernels
    gfuns = gfuns + block_out + sac_out + msac_out;
    
    %incorporate beta
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,2)');
    
    too_large = gfuns > 50;
    pred_rate = log(1+exp(gfuns));
    pred_rate(too_large) = gfuns(too_large);
    
    %incorporate alpha
    pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(:,3)');
    
    pred_rate(pred_rate < 1e-20) = 1e-20;
    
    LLs(:,:,xx) = Robs_mat.*log(pred_rate) - pred_rate;
end

%% INFER EYE-POSITION USING ALL TR UNITS 
fprintf('Inferring eye-position using all cells\n');
frame_LLs = squeeze(nansum(LLs,2));
lalpha=zeros(NT,n_shifts);
lbeta = zeros(NT,n_shifts);
lscale=zeros(NT,1); %initialize rescaling parameters
%compute rescaled forward messages
lalpha(1,:) = leps_prior + frame_LLs(1,:);
lscale(1)=logsumexp(lalpha(1,:));
lalpha(1,:) = lalpha(1,:) - lscale(1);
for t=2:NT
    if use_prior(t)
        cur_lA = lA_tflip;
    else
        cur_lA = lA;
    end
    lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + frame_LLs(t,:);
    lscale(t) = logsumexp(lalpha(t,:));
    lalpha(t,:) = lalpha(t,:) - lscale(t);
end

%compute rescaled backward messages
lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
for t=NT-1:-1:1
    if use_prior(t+1)
        cur_lA = lA_tflip;
    else
        cur_lA = lA;
    end
    lf1 = lbeta(t+1,:) + frame_LLs(t+1,:);
    lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
end
lgamma= lalpha + lbeta;
lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));

gamma = exp(lgamma);
[~,left_fpost_max] = max(gamma,[],2); 
left_fpost_mean = sum(bsxfun(@times,gamma,shifts),2);
left_fpost_std = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - left_fpost_mean.^2);

%% NOW DO RIGHT EYE CORRECTIONS

%% SELECT USABLE UNITS AND make Robs_mat
right_usable_units = tr_set(right_ref_LL_imp(1,tr_set) >= LL_imp_thresh);
right_tr_set = right_usable_units;
n_tr_chs = length(right_tr_set);


% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = nan(length(used_inds),n_tr_chs);
for ss = 1:n_tr_chs
    if all_mod_SU(right_tr_set(ss)) == 0
        su_probe_ind = find(su_probes == right_tr_set(ss));
        if ~isempty(su_probe_ind)
            cur_used_blocks = find(su_used_blocks(:,su_probe_ind) == 0); %blocks when NO SU
        else
            cur_used_blocks = 1:n_blocks; %blocks when NO SU
        end
        cur_use = find(ismember(all_blockvec(used_inds),cur_used_blocks));
        Robs_mat(cur_use,ss) = all_binned_spikes(used_inds(cur_use),right_tr_set(ss));
    else
        su_ind = find(su_probes == all_mod_SU(right_tr_set(ss)));
        cur_used_blocks = find(su_used_blocks(:,su_ind)); %blocks when SU
        cur_use = find(ismember(all_blockvec(used_inds),cur_used_blocks));
        Robs_mat(cur_use,ss) = all_binned_spikes(used_inds(cur_use),all_mod_SU(right_tr_set(ss)));
    end
end

%% PREPROCESS MODEL COMPONENTS
klen = flen*use_nPix;
filt_bank = nan(n_tr_chs,klen,3);
lin_kerns = nan(n_tr_chs,n_blocks);
sac_kerns = nan(n_tr_chs,n_sac_bins);
msac_kerns = nan(n_tr_chs,n_sac_bins);
mod_spkNL_params = nan(n_tr_chs,3);
for ss = 1:n_tr_chs
    cur_Xtargs = [right_ref_mod{1}(right_tr_set(ss)).mods(:).Xtarget];
    cur_k = [right_ref_mod{1}(right_tr_set(ss)).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(ss,:,1:n_used_filts) = cur_k;
    mod_spkNL_params(ss,:) = right_ref_mod_spkNL{1}(right_tr_set(ss)).spk_NL_params;
    lin_kerns(ss,:) = right_ref_mod{1}(right_tr_set(ss)).mods(cur_Xtargs == 2).filtK;
    sac_kerns(ss,:) = right_ref_mod{1}(right_tr_set(ss)).mods(cur_Xtargs == 3).filtK;
    msac_kerns(ss,:) = right_ref_mod{1}(right_tr_set(ss)).mods(cur_Xtargs == 4).filtK;
end
filt_bank = permute(filt_bank,[2 1 3]);
filt_bank = reshape(filt_bank,[flen use_nPix n_tr_chs 3]);

%indicator predictions
block_out = Xblock(used_inds,:)*lin_kerns';
sac_out = Xsac*sac_kerns';
msac_out = Xmsac*msac_kerns';

%% ESTIMATE LL for each shift in each stimulus frame
cur_Xmat = all_right_Xmat(used_inds,use_kInds);

%precompute LL at all shifts for all units
LLs = nan(NT,length(right_tr_set),n_shifts);
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',xx,n_shifts);
    cur_filt_shift = shift_matrix_Nd(filt_bank,shifts(xx),2);
    cur_filt_shift = reshape(cur_filt_shift,[klen n_tr_chs 3]);
    
    %outputs of stimulus models at current X-matrix shift
    gfuns = ones(NT,n_tr_chs);
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
    gfuns = gfuns + cur_Xmat*squeeze(cur_filt_shift(:,:,1));
    for ff = 2:3
        gfuns = gfuns + (cur_Xmat*squeeze(cur_filt_shift(:,:,ff))).^2;
    end
    
    %add contributions from extra lin kernels
    gfuns = gfuns + block_out + sac_out + msac_out;
    
    %incorporate beta
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,2)');
    
    too_large = gfuns > 50;
    pred_rate = log(1+exp(gfuns));
    pred_rate(too_large) = gfuns(too_large);
    
    %incorporate alpha
    pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(:,3)');
    
    pred_rate(pred_rate < 1e-20) = 1e-20;
    
    LLs(:,:,xx) = Robs_mat.*log(pred_rate) - pred_rate;
end

%% INFER EYE-POSITION USING ALL TR UNITS 
fprintf('Inferring eye-position using all cells\n');
frame_LLs = squeeze(nansum(LLs,2));
lalpha=zeros(NT,n_shifts);
lbeta = zeros(NT,n_shifts);
lscale=zeros(NT,1); %initialize rescaling parameters
%compute rescaled forward messages
lalpha(1,:) = leps_prior + frame_LLs(1,:);
lscale(1)=logsumexp(lalpha(1,:));
lalpha(1,:) = lalpha(1,:) - lscale(1);
for t=2:NT
    if use_prior(t)
        cur_lA = lA_tflip;
    else
        cur_lA = lA;
    end
    lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + frame_LLs(t,:);
    lscale(t) = logsumexp(lalpha(t,:));
    lalpha(t,:) = lalpha(t,:) - lscale(t);
end

%compute rescaled backward messages
lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
for t=NT-1:-1:1
    if use_prior(t+1)
        cur_lA = lA_tflip;
    else
        cur_lA = lA;
    end
    lf1 = lbeta(t+1,:) + frame_LLs(t+1,:);
    lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
end
lgamma= lalpha + lbeta;
lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));

gamma = exp(lgamma);
[~,right_fpost_max] = max(gamma,[],2);
right_fpost_mean = sum(bsxfun(@times,gamma,shifts),2);
right_fpost_std = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - right_fpost_mean.^2);


%% REFIT ALL CELLS THAT AREN'T IN POSS_XV_SET
silent = 1;
cur_fit_set = setdiff(tr_set,poss_xv_set);
for ss = 1:length(cur_fit_set)
    tr_cell = cur_fit_set(ss);
    fprintf('Refitting model for tr cell %d of %d\n',ss,length(cur_fit_set));
    cur_unit_ind = find(tr_set == tr_cell);
    cur_uset = find(~isnan(Robs_mat(:,cur_unit_ind)));
    
    cur_X{1} = X_sh(cur_uset,use_kInds);
    cur_X{2} = Xblock(used_inds(cur_uset),:);
    cur_X{3} = Xsac(cur_uset,:);
    cur_X{4} = Xmsac(cur_uset,:);
    
    ref_mod{1}(tr_cell) = all_mod_fits(tr_cell);
    ref_mod{1}(tr_cell) = NMMfit_filters(ref_mod{1}(tr_cell),Robs_mat(cur_uset,cur_unit_ind),...
        cur_X,[],[],silent); %fit stimulus filters
    
    ref_LL_imp(1,tr_cell) = (ref_mod{1}(tr_cell).LL_seq(end) - all_nullmod(tr_cell).LL_seq(end))/log(2);
    
    fprintf('Original: %.4f  New: %.4f\n',all_mod_LLimp(tr_cell),ref_LL_imp(1,tr_cell));

    %refit spk NL
    ref_mod_spkNL{1}(tr_cell) = NMMfit_logexp_spkNL(ref_mod{1}(tr_cell),Robs_mat(cur_uset,cur_unit_ind),cur_X);
    
end
