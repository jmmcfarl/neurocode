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

anal_name = 'binoc_eyecorr_vbar_v2';
mod_data_name = 'binoc_eyecorr_vbar_mods';
eye_data_name = 'binoc_eyecorr_vbar_eyeonly_v2';
use_right_eye = false;
%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 10;
use_nPix = 24;
full_nPix = 36;
init_stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;
bar_ori = 90;
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

%%
spatial_usfac = 2;
full_nPix_us = spatial_usfac*full_nPix;
all_left_stimmat_up = zeros(size(all_left_stim_mat,1),full_nPix_us);
all_right_stimmat_up = zeros(size(all_right_stim_mat,1),full_nPix_us);
for ii = 1:size(all_left_stim_mat,2)
    all_left_stimmat_up(:,2*(ii-1)+1) = all_left_stim_mat(:,ii);
    all_left_stimmat_up(:,2*(ii-1)+2) = all_left_stim_mat(:,ii);
    all_right_stimmat_up(:,2*(ii-1)+1) = all_right_stim_mat(:,ii);
    all_right_stimmat_up(:,2*(ii-1)+2) = all_right_stim_mat(:,ii);
end

stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);
% all_Xmat = create_time_embedding(all_stim_mat,stim_params);
all_left_Xmat_us = create_time_embedding(all_left_stimmat_up,stim_params_us);
all_right_Xmat_us = create_time_embedding(all_right_stimmat_up,stim_params_us);

%% select submatrix with central pixels
[Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix));

[Xinds_up,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-1/spatial_usfac & Xinds_up(:) <= cur_use_pix(end));
use_kInds_back = find(ismember(Xinds_up(use_kInds_up),cur_use_pix));

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
% fprintf('Removing %d of %d samples, %.4f where eyes are outside the FW\n',length(outside_window),length(used_inds),length(outside_window)/length(used_inds));
% used_inds(outside_window) = [];
% NT = length(used_inds);

%% DETECT TIMES WHEN TRIAL FAILURES OCCUR AND FORCE-END TRIALS
par_thresh = 4;
orth_thresh = 1.2;
back_samps = 1;

n_trials = length(all_trial_start_times);
out_of_trial = [];
for ii = 1:n_trials
    cur_trial_inds = used_inds(all_trialvec(used_inds) == ii);
    if ~isempty(cur_trial_inds)
        if bar_ori == 0
            cur_out = find(abs(corrected_interp_eyevals(cur_trial_inds,2)) >= orth_thresh | ...
                abs(corrected_interp_eyevals(cur_trial_inds,1)) >= par_thresh,1,'first');
        elseif bar_ori == 90
            cur_out = find(abs(corrected_interp_eyevals(cur_trial_inds,1)) >= orth_thresh | ...
                abs(corrected_interp_eyevals(cur_trial_inds,2)) >= par_thresh,1,'first');
        end
        if ~isempty(cur_out)
            if cur_out > back_samps
                cur_out = cur_out - back_samps;
            end
            cur_out = cur_out:length(cur_trial_inds);
        end
        out_of_trial = cat(1,out_of_trial,cur_trial_inds(cur_out(:)));
    end
end
fract_out = length(out_of_trial)/length(used_inds);
fprintf('Eliminating %.4f of data out of window\n',fract_out);
used_inds(ismember(used_inds,out_of_trial)) = [];
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

%% FIRST USE LEFT EYE ONLY
%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
full_X{1} = [all_left_Xmat_us(used_inds,use_kInds_up)];
full_X{2} = Xblock(used_inds,:);
full_X{3} = Xsac;
full_X{4} = Xmsac;

use_nPix_us = spatial_usfac*use_nPix;
full_stim_params = NMMcreate_stim_params([flen use_nPix_us],dt);
full_stim_params(2) = NMMcreate_stim_params(n_blocks);
full_stim_params(3) = NMMcreate_stim_params(n_sac_bins,dt);
full_stim_params(4) = NMMcreate_stim_params(n_sac_bins,dt);

null_stim_params = full_stim_params(2:end);

base_lambda_d2XT = 8;
base_lambda_L1 = 5;
silent = 0;

sac_d2t = 100;
block_L2 = 1;
null_reg_params = NMMcreate_reg_params('lambda_d2T',[0; sac_d2t; sac_d2t],'lambda_L2',[block_L2; 0; 0]);
sac_reg_params = NMMcreate_reg_params('lambda_d2T',sac_d2t);

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
            
            mod_signs = [new_fit.mods(:).sign];
            n_squared_filts = sum(Xtargs==1)-1;
            NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
            
            base_filts = [new_fit.mods(find(Xtargs==1)).filtK];
            base_filts = reshape(base_filts,[flen 2*use_nPix (n_squared_filts+1)]);
            base_filts = base_filts(:,1:use_nPix,:); %use only left-eye
            init_stim_filts = zeros(flen*use_nPix,n_squared_filts+1);
            init_stim_filts = reshape(init_stim_filts,[flen use_nPix (n_squared_filts+1)]);
            for ii = 1:use_nPix
                init_stim_filts(:,2*(ii-1)+1,:) = 0.5*base_filts(:,ii,:);
                init_stim_filts(:,2*(ii-1)+2,:) = 0.5*base_filts(:,ii,:);
            end
            init_stim_filts = reshape(init_stim_filts,use_nPix_us*flen,n_squared_filts+1);
            
            init_filts = cell(length(mod_signs),1);
            for ii = 1:n_squared_filts+1
                init_filts{ii} = init_stim_filts(:,ii);
            end
            init_filts{end} = new_fit.mods(find(Xtargs==2)).filtK;
            new_fit_up = NMMinitialize_model(full_stim_params,mod_signs,NL_types,[],Xtargs,init_filts); %initialize NIM
            
            new_fit_up.spk_NL_params(1) = new_fit.spk_NL_params(1);
            new_fit_up = NMMadd_NLinput(new_fit_up,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
            new_fit_up = NMMadd_NLinput(new_fit_up,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
            Xtargs = [new_fit_up.mods(:).Xtarget];
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(new_fit_up,Robs,cur_X);
            new_fit_up = NMMadjust_regularization(new_fit_up,find(Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
            new_fit_up = NMMadjust_regularization(new_fit_up,find(Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
            
            new_fit_up = NMMfit_filters(new_fit_up,Robs,cur_X,[],[],silent);
            
            all_left_mod_fits(ss) = new_fit_up;
            LL = NMMmodel_eval(new_fit_up,Robs,cur_X);
            all_left_mod_LLimp(ss) = (LL-null_mod.LL_seq(end))/log(2);
            all_left_mod_fits_withspkNL(ss) = NMMfit_logexp_spkNL(new_fit_up,Robs,cur_X);
            all_mod_nspks(ss) = sum(Robs);
            all_mod_avgrate(ss) = mean(Robs);
        end
    else
        all_left_mod_LLimp(ss) = nan;
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
            
            new_fit = sua_data(ss).nimFit;
            Xtargs = [new_fit.mods(:).Xtarget];
            stim_targs = find(Xtargs == 1);
            for ii = 1:length(new_fit.mods)
                new_fit.mods(ii).reg_params.boundary_conds = [Inf Inf Inf];
                new_fit.mods(ii).reg_params.mixing_prop = [1 1 1];
            end
            
            mod_signs = [new_fit.mods(:).sign];
            n_squared_filts = sum(Xtargs==1)-1;
            NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
            
            base_filts = [new_fit.mods(find(Xtargs==1)).filtK];
            base_filts = reshape(base_filts,[flen 2*use_nPix (n_squared_filts+1)]);
            base_filts = base_filts(:,1:use_nPix,:); %use only left-eye
            init_stim_filts = zeros(flen*use_nPix,n_squared_filts+1);
            init_stim_filts = reshape(init_stim_filts,[flen use_nPix (n_squared_filts+1)]);
            for ii = 1:use_nPix
                init_stim_filts(:,2*(ii-1)+1,:) = 0.5*base_filts(:,ii,:);
                init_stim_filts(:,2*(ii-1)+2,:) = 0.5*base_filts(:,ii,:);
            end
            init_stim_filts = reshape(init_stim_filts,use_nPix_us*flen,n_squared_filts+1);
            
            init_filts = cell(length(mod_signs),1);
            for ii = 1:n_squared_filts+1
                init_filts{ii} = init_stim_filts(:,ii);
            end
            init_filts{end} = new_fit.mods(find(Xtargs==2)).filtK;
            new_fit_up = NMMinitialize_model(full_stim_params,mod_signs,NL_types,[],Xtargs,init_filts); %initialize NIM
            
            new_fit_up.spk_NL_params(1) = new_fit.spk_NL_params(1);
            new_fit_up = NMMadd_NLinput(new_fit_up,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
            new_fit_up = NMMadd_NLinput(new_fit_up,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
            Xtargs = [new_fit_up.mods(:).Xtarget];
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(new_fit_up,Robs,cur_X);
            new_fit_up = NMMadjust_regularization(new_fit_up,find(Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
            new_fit_up = NMMadjust_regularization(new_fit_up,find(Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
            
            new_fit_up = NMMfit_filters(new_fit_up,Robs,cur_X,[],[],silent);
            
            all_left_mod_fits(ss+96) = new_fit_up;
            LL = NMMmodel_eval(new_fit_up,Robs,cur_X);
            all_left_mod_LLimp(ss+96) = (LL-null_mod.LL_seq(end))/log(2);
            all_left_mod_fits_withspkNL(ss+96) = NMMfit_logexp_spkNL(new_fit_up,Robs,cur_X);
            all_mod_nspks(ss+96) = sum(Robs);
            all_mod_avgrate(ss+96) = mean(Robs);
        end
    else
        all_left_mod_LLimp(ss+96) = nan;
    end
end

%% SELECT USABLE UNITS AND make Robs_mat
LL_imp_thresh = 4e-3;
min_avgrate = 2;
left_usable_units = find(all_left_mod_LLimp >= LL_imp_thresh & all_mod_avgrate/dt >= min_avgrate);
n_used_sus = sum(all_mod_SU(left_usable_units) ~= 0);
n_used_mus = sum(all_mod_SU(left_usable_units) == 0);
fprintf('Using %d SUs and %d MUs for analysis\n',n_used_sus,n_used_mus);
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

%% INITIALIZE TRANSITION PRIORS FOR HMM
sp_dx = 0.05/spatial_usfac;
max_shift = 24;
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);

%overall prior on shifts
eps_prior_sigma = 0.15; %0.125 start
leps_prior = -(shifts*sp_dx).^2/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior,n_shifts,1);

cdist = squareform(pdist(shifts'*sp_dx));
deps_sigma = 0.0125; %.015
lA = -cdist.^2/(2*deps_sigma^2);
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

%% PREPROCESS MODEL COMPONENTS
klen = flen*use_nPix_us;
filt_bank = nan(n_tr_chs,klen,3);
lin_kerns = nan(n_tr_chs,n_blocks);
sac_kerns = nan(n_tr_chs,n_sac_bins);
msac_kerns = nan(n_tr_chs,n_sac_bins);
mod_spkNL_params = nan(n_tr_chs,3);
for ss = 1:n_tr_chs
    cur_Xtargs = [all_left_mod_fits(left_tr_set(ss)).mods(:).Xtarget];
    cur_k = [all_left_mod_fits(left_tr_set(ss)).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(ss,:,1:n_used_filts) = cur_k;
    mod_spkNL_params(ss,:) = all_left_mod_fits_withspkNL(left_tr_set(ss)).spk_NL_params;
    lin_kerns(ss,:) = all_left_mod_fits(left_tr_set(ss)).mods(cur_Xtargs == 2).filtK;
    sac_kerns(ss,:) = all_left_mod_fits(left_tr_set(ss)).mods(cur_Xtargs == 3).filtK;
    msac_kerns(ss,:) = all_left_mod_fits(left_tr_set(ss)).mods(cur_Xtargs == 4).filtK;
end
filt_bank = permute(filt_bank,[2 1 3]);
filt_bank = reshape(filt_bank,[flen use_nPix_us n_tr_chs 3]);

%indicator predictions
block_out = Xblock(used_inds,:)*lin_kerns';
sac_out = Xsac*sac_kerns';
msac_out = Xmsac*msac_kerns';

%% ESTIMATE LL for each shift in each stimulus frame
cur_Xmat = [all_left_Xmat_us(used_inds,use_kInds_up)];

%precompute LL at all shifts for all units
LLs = nan(NT,length(left_tr_set),n_shifts);
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',xx,n_shifts);
    cur_filt_shift = shift_matrix_Nd(filt_bank,shifts(xx),2);
    cur_filt_shift = reshape(cur_filt_shift,[flen*use_nPix_us n_tr_chs 3]);
    
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
left_fpost_mean = sum(bsxfun(@times,gamma,shifts),2);
left_fpost_std = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - left_fpost_mean.^2);
[max_post,max_loc] = max(lgamma,[],2);
left_fshift_cor = shifts(max_loc);
%% VITERBI SEQ FOR LEFT EYE
delta=zeros(NT,n_shifts); %delta is the maximized log probability as a function of time
psi=zeros(NT,n_shifts); %Psi stores the most likely preceeding state
%initialization
delta(1,:) = leps_prior + frame_LLs(1,:);
%forward recursion
for t=2:NT
    if use_prior(t)
        cur_lA = lA_tflip;
    else
        cur_lA = lA;
    end
    temp = bsxfun(@plus,delta(t-1,:),cur_lA);
    [delta(t,:),psi(t,:)] = max(temp,[],2);
    delta(t,:) = delta(t,:) + frame_LLs(t,:);
end

% Backtracking Viterbi
[llik_best,state_seq(NT)] = max(delta(NT,:));
for t=(NT-1):-1:1
    state_seq(t) = psi(t+1,state_seq(t+1));
end
left_viterbi_shifts = shifts(state_seq);

%% NOW RIGHT-EYE ONLY
%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
full_X{1} = [all_right_Xmat_us(used_inds,use_kInds_up)];
full_X{2} = Xblock(used_inds,:);
full_X{3} = Xsac;
full_X{4} = Xmsac;

use_nPix_us = spatial_usfac*use_nPix;
full_stim_params = NMMcreate_stim_params([flen use_nPix_us],dt);
full_stim_params(2) = NMMcreate_stim_params(n_blocks);
full_stim_params(3) = NMMcreate_stim_params(n_sac_bins,dt);
full_stim_params(4) = NMMcreate_stim_params(n_sac_bins,dt);

null_stim_params = full_stim_params(2:end);

base_lambda_d2XT = 8;
base_lambda_L1 = 5;
silent = 0;

sac_d2t = 100;
block_L2 = 1;
null_reg_params = NMMcreate_reg_params('lambda_d2T',[0; sac_d2t; sac_d2t],'lambda_L2',[block_L2; 0; 0]);
sac_reg_params = NMMcreate_reg_params('lambda_d2T',sac_d2t);

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
            
            mod_signs = [new_fit.mods(:).sign];
            n_squared_filts = sum(Xtargs==1)-1;
            NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
            
            base_filts = [new_fit.mods(find(Xtargs==1)).filtK];
            base_filts = reshape(base_filts,[flen 2*use_nPix (n_squared_filts+1)]);
            base_filts = base_filts(:,(use_nPix+1):end,:); %use only right-eye
            init_stim_filts = zeros(flen*use_nPix,n_squared_filts+1);
            init_stim_filts = reshape(init_stim_filts,[flen use_nPix (n_squared_filts+1)]);
            for ii = 1:use_nPix
                init_stim_filts(:,2*(ii-1)+1,:) = 0.5*base_filts(:,ii,:);
                init_stim_filts(:,2*(ii-1)+2,:) = 0.5*base_filts(:,ii,:);
            end
            init_stim_filts = reshape(init_stim_filts,use_nPix_us*flen,n_squared_filts+1);
            
            init_filts = cell(length(mod_signs),1);
            for ii = 1:n_squared_filts+1
                init_filts{ii} = init_stim_filts(:,ii);
            end
            init_filts{end} = new_fit.mods(find(Xtargs==2)).filtK;
            new_fit_up = NMMinitialize_model(full_stim_params,mod_signs,NL_types,[],Xtargs,init_filts); %initialize NIM
            
            new_fit_up.spk_NL_params(1) = new_fit.spk_NL_params(1);
            new_fit_up = NMMadd_NLinput(new_fit_up,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
            new_fit_up = NMMadd_NLinput(new_fit_up,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
            Xtargs = [new_fit_up.mods(:).Xtarget];
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(new_fit_up,Robs,cur_X);
            new_fit_up = NMMadjust_regularization(new_fit_up,find(Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
            new_fit_up = NMMadjust_regularization(new_fit_up,find(Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
            
            new_fit_up = NMMfit_filters(new_fit_up,Robs,cur_X,[],[],silent);
            
            all_right_mod_fits(ss) = new_fit_up;
            LL = NMMmodel_eval(new_fit_up,Robs,cur_X);
            all_right_mod_LLimp(ss) = (LL-null_mod.LL_seq(end))/log(2);
            all_right_mod_fits_withspkNL(ss) = NMMfit_logexp_spkNL(new_fit_up,Robs,cur_X);
            all_mod_nspks(ss) = sum(Robs);
            all_mod_avgrate(ss) = mean(Robs);
        end
    else
        all_right_mod_LLimp(ss) = nan;
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
            
            new_fit = sua_data(ss).nimFit;
            Xtargs = [new_fit.mods(:).Xtarget];
            stim_targs = find(Xtargs == 1);
            for ii = 1:length(new_fit.mods)
                new_fit.mods(ii).reg_params.boundary_conds = [Inf Inf Inf];
                new_fit.mods(ii).reg_params.mixing_prop = [1 1 1];
            end
            
            mod_signs = [new_fit.mods(:).sign];
            n_squared_filts = sum(Xtargs==1)-1;
            NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
            
            base_filts = [new_fit.mods(find(Xtargs==1)).filtK];
            base_filts = reshape(base_filts,[flen 2*use_nPix (n_squared_filts+1)]);
            base_filts = base_filts(:,(use_nPix+1):end,:); %use only right-eye
            init_stim_filts = zeros(flen*use_nPix,n_squared_filts+1);
            init_stim_filts = reshape(init_stim_filts,[flen use_nPix (n_squared_filts+1)]);
            for ii = 1:use_nPix
                init_stim_filts(:,2*(ii-1)+1,:) = 0.5*base_filts(:,ii,:);
                init_stim_filts(:,2*(ii-1)+2,:) = 0.5*base_filts(:,ii,:);
            end
            init_stim_filts = reshape(init_stim_filts,use_nPix_us*flen,n_squared_filts+1);
            
            init_filts = cell(length(mod_signs),1);
            for ii = 1:n_squared_filts+1
                init_filts{ii} = init_stim_filts(:,ii);
            end
            init_filts{end} = new_fit.mods(find(Xtargs==2)).filtK;
            new_fit_up = NMMinitialize_model(full_stim_params,mod_signs,NL_types,[],Xtargs,init_filts); %initialize NIM
            
            new_fit_up.spk_NL_params(1) = new_fit.spk_NL_params(1);
            new_fit_up = NMMadd_NLinput(new_fit_up,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
            new_fit_up = NMMadd_NLinput(new_fit_up,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
            Xtargs = [new_fit_up.mods(:).Xtarget];
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(new_fit_up,Robs,cur_X);
            new_fit_up = NMMadjust_regularization(new_fit_up,find(Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
            new_fit_up = NMMadjust_regularization(new_fit_up,find(Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
            
            new_fit_up = NMMfit_filters(new_fit_up,Robs,cur_X,[],[],silent);
            
            all_right_mod_fits(ss+96) = new_fit_up;
            LL = NMMmodel_eval(new_fit_up,Robs,cur_X);
            all_right_mod_LLimp(ss+96) = (LL-null_mod.LL_seq(end))/log(2);
            all_right_mod_fits_withspkNL(ss+96) = NMMfit_logexp_spkNL(new_fit_up,Robs,cur_X);
            all_mod_nspks(ss+96) = sum(Robs);
            all_mod_avgrate(ss+96) = mean(Robs);
        end
    else
        all_right_mod_LLimp(ss+96) = nan;
    end
end

%% SELECT USABLE UNITS AND make Robs_mat
LL_imp_thresh = 5e-3;
min_avgrate = 2;
right_usable_units = find(all_right_mod_LLimp >= LL_imp_thresh & all_mod_avgrate/dt >= min_avgrate);
n_used_sus = sum(all_mod_SU(right_usable_units) ~= 0);
n_used_mus = sum(all_mod_SU(right_usable_units) == 0);
fprintf('Using %d SUs and %d MUs for analysis\n',n_used_sus,n_used_mus);
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

%% INITIALIZE TRANSITION PRIORS FOR HMM
% sp_dx = 0.05/spatial_usfac;
% max_shift = 24;
% dshift = 1;
% shifts = -max_shift:dshift:max_shift;
% n_shifts = length(shifts);
% zero_frame = find(shifts==0);
% 
% %overall prior on shifts
% eps_prior_sigma = 0.15; %0.125 start
% leps_prior = -(shifts*sp_dx).^2/(2*eps_prior_sigma^2);
% leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
% lA_tflip = repmat(leps_prior,n_shifts,1);
% 
% cdist = squareform(pdist(shifts'*sp_dx));
% deps_sigma = 0.0125; %.015
% lA = -cdist.^2/(2*deps_sigma^2);
% lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

%% PREPROCESS MODEL COMPONENTS
klen = flen*use_nPix_us;
filt_bank = nan(n_tr_chs,klen,3);
lin_kerns = nan(n_tr_chs,n_blocks);
sac_kerns = nan(n_tr_chs,n_sac_bins);
msac_kerns = nan(n_tr_chs,n_sac_bins);
mod_spkNL_params = nan(n_tr_chs,3);
for ss = 1:n_tr_chs
    cur_Xtargs = [all_right_mod_fits(right_tr_set(ss)).mods(:).Xtarget];
    cur_k = [all_right_mod_fits(right_tr_set(ss)).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(ss,:,1:n_used_filts) = cur_k;
    mod_spkNL_params(ss,:) = all_right_mod_fits_withspkNL(right_tr_set(ss)).spk_NL_params;
    lin_kerns(ss,:) = all_right_mod_fits(right_tr_set(ss)).mods(cur_Xtargs == 2).filtK;
    sac_kerns(ss,:) = all_right_mod_fits(right_tr_set(ss)).mods(cur_Xtargs == 3).filtK;
    msac_kerns(ss,:) = all_right_mod_fits(right_tr_set(ss)).mods(cur_Xtargs == 4).filtK;
end
filt_bank = permute(filt_bank,[2 1 3]);
filt_bank = reshape(filt_bank,[flen use_nPix_us n_tr_chs 3]);

%indicator predictions
block_out = Xblock(used_inds,:)*lin_kerns';
sac_out = Xsac*sac_kerns';
msac_out = Xmsac*msac_kerns';

%% ESTIMATE LL for each shift in each stimulus frame
cur_Xmat = [all_right_Xmat_us(used_inds,use_kInds_up)];

%precompute LL at all shifts for all units
LLs = nan(NT,length(right_tr_set),n_shifts);
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',xx,n_shifts);
    cur_filt_shift = shift_matrix_Nd(filt_bank,shifts(xx),2);
    cur_filt_shift = reshape(cur_filt_shift,[flen*use_nPix_us n_tr_chs 3]);
    
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
right_fpost_mean = sum(bsxfun(@times,gamma,shifts),2);
right_fpost_std = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - right_fpost_mean.^2);
[max_post,max_loc] = max(lgamma,[],2);
right_fshift_cor = shifts(max_loc);

%% VITERBI DECODING
delta=zeros(NT,n_shifts); %delta is the maximized log probability as a function of time
psi=zeros(NT,n_shifts); %Psi stores the most likely preceeding state
%initialization
delta(1,:) = leps_prior + frame_LLs(1,:);
%forward recursion
for t=2:NT
    if use_prior(t)
        cur_lA = lA_tflip;
    else
        cur_lA = lA;
    end
    temp = bsxfun(@plus,delta(t-1,:),cur_lA);
    [delta(t,:),psi(t,:)] = max(temp,[],2);
    delta(t,:) = delta(t,:) + frame_LLs(t,:);
end

% Backtracking Viterbi
[llik_best,state_seq(NT)] = max(delta(NT,:));
for t=(NT-1):-1:1
    state_seq(t) = psi(t+1,state_seq(t+1));
end
viterbi_shifts = shifts(state_seq);
