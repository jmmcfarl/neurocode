clear all
% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';

Expt_num = 94;

%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz
dt = 0.01;

min_trial_dur = 0.5;
beg_buffer = 0.3;

backlag = round(0.25/dt);
forwardlag = round(0.55/dt);

sua_sm_sig = (0.01/dt);
mua_sm_sig = (0.005/dt);

flen = 12;
use_nPix = 24;
full_nPix = 36;
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;
trial_dur = 4;
bar_ori = 0;

anal_name = 'hbar_stcprobmod_sacmod_analysis';

%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat


clear expt* included_type

Expt_name = sprintf('G0%d',Expt_num);
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
anal_dir = ['~/Analysis/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/G081/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
if ~strcmp(Expt_name,'G081')
    load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data
end

save_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/sac_mod'];
if ~exist(save_dir,'dir')
    system(['mkdir ' save_dir]);
end

if exist('./fin_aclust_data.mat','file')
    load('./fin_aclust_data.mat');
    [n_aclust_expts,n_aclust_probes] = size(autoclust);
else
    disp('No fin_aclust_data found.');
    autoclust = [];
    n_aclust_expts = 0; n_aclust_probes = 0;
end

%%
include_expts = {'rls.orXFaXFrRC'};
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
% cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1);
cur_block_set = find(included_type);

if Expt_num == 94
    cur_block_set(cur_block_set == 28) = [];
end
%% load overall su data
% cur_expt_id = find(su_data.expt_nums == Expt_num);
% su_probes = find(su_data.is_su(cur_expt_id,:));
% mua_probes = setdiff(1:96,su_probes); %probes with ONLY MU
% aclust_probenums = [autoclust(cur_block_set(1),:).probenum];
% autoclust = autoclust(:,ismember(aclust_probenums,su_probes));
su_probes = [];
mua_probes = 1:96;
%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
all_stim_ori = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_wi = [];
all_trial_Fr = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(96,1);
all_clust_ids = cell(length(su_probes),1);
for ee = 1:length(cur_block_set);
    fprintf('Expt %d Block %d of %d\n',Expt_num,ee,length(cur_block_set));
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
    trial_Fr = [Expts{cur_block}.Trials(:).Fr];
    
    [un_ids,id_inds] = unique(trial_ids);
    rpt_trials = false;
    if length(un_ids) < length(trial_ids)
        rpt_trials = true;
        fprintf('Warning, repeat trial inds detected!\n');
        use_trials = [];
    else
        use_trials = find(trial_durs >= min_trial_dur);
    end
    %         trial_durs = trial_durs(id_inds);
    %         trial_start_times = trial_start_times(id_inds);
    %         trial_end_times = trial_end_times(id_inds);
    
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
        
    all_trial_Fr = cat(1,all_trial_Fr,trial_Fr(use_trials)');
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    buffer_pix = floor((expt_npix(cur_block) - full_nPix)/2);
    cur_use_pix = (1:full_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start/1e4;
        cur_stim_ori = Expts{cur_block}.Trials(use_trials(tt)).or;
        if trial_Fr(use_trials(tt)) == 3
           cur_stim_ori = repmat(cur_stim_ori,1,3);
           cur_stim_ori = cur_stim_ori';
           cur_stim_ori = cur_stim_ori(:); 
           cur_stim_times = trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt));
           cur_stim_times = cur_stim_times(:);
           cur_stim_ori(length(cur_stim_times)+1:end) = [];
           cur_stim_times(length(cur_stim_ori)+1:end) = [];
        end
        
%         if length(cur_stim_times) == 1
%             cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
%             cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
%         end
        %             cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
            cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt];
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        
        if ~isempty(all_stim_times)
            if any(cur_stim_times < all_stim_times(end))
                fprintf('Warn trial %d\n',tt);
            end
        end
        all_stim_times = [all_stim_times; cur_stim_times];
        all_stim_ori = [all_stim_ori; cur_stim_ori];
        all_t_axis = [all_t_axis; cur_t_axis];
        all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
        all_tsince_start = [all_tsince_start; cur_tsince_start];
        all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
        all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
    end
    trial_cnt = trial_cnt + n_trials;
end
% bar_Xmat = create_time_embedding(all_stim_mat,stim_params);

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
used_inds = find(all_tsince_start >= beg_buffer);
if strcmp(Expt_name,'G093')
    un_wi_vals = unique(all_trial_wi);
    use_wi_trials = find(all_trial_wi == un_wi_vals(2) | all_trial_wi == un_wi_vals(3));
    used_inds = used_inds(ismember(all_trialvec(used_inds),use_wi_trials));
end

% used_inds = used_inds(all_trial_Fr(all_trialvec(used_inds)) == 1);
% used_inds = used_inds(all_trial_Fr(all_trialvec(used_inds)) == 3);

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)

emfile = ['jbe' Expt_name '.em.mat'];
load(emfile);

all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
all_eye_blockvec = [];
eye_smooth = 3;
for ee = 1:length(cur_block_set);
    fprintf('Loading ET data for expt %d, block %d of %d\n',Expt_num,ee,length(cur_block_set));
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


%% PICK OUT SACCADES FOR ANALYSIS
sac_start_times = [saccades(:).start_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > dt);
sac_start_times(bad_sacs) = [];
saccades(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];

sac_stop_times = [saccades(:).stop_time];
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));
interp_sac_stop_inds(isnan(interp_sac_stop_inds)) = 1;

sac_peak_times = [saccades(:).peak_time];
interp_sac_peak_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_peak_times));
interp_sac_peak_inds(isnan(interp_sac_peak_inds)) = 1;

max_dur = 0.1;
sac_durs = [saccades(:).duration];
is_blink = sac_durs > max_dur;

isis = [saccades(:).isi];
is_sacburst = false(length(saccades),1);
is_sacburst(isis(1:end-1) < 0.15 | isis(2:end) < 0.15) = true;

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps < 1;

sac_post_Lx = [saccades(:).post_Lx];
sac_post_Ly = [saccades(:).post_Ly];
sac_pre_Lx = [saccades(:).pre_Lx];
sac_pre_Ly = [saccades(:).pre_Ly];
sac_deltaX = sac_post_Lx - sac_pre_Lx;
sac_deltaY = sac_post_Ly - sac_pre_Ly;

sac_dirs = [0 90];
sim_sac_times = [0.7 1.4 2.1 2.8 3.5];
sac_thresh = 1.5;
sac_start_expts = all_blockvec(interp_sac_start_inds);
delta_sacpar = nan(size(saccades));
for bb = 1:length(sac_dirs)
    cur_bar_expts = find(expt_bar_ori(cur_block_set) == sac_dirs(bb));
    cur_sac_set = find(ismember(sac_start_expts,cur_bar_expts));
    
    cur_delta_sac_par = abs(sac_deltaX(cur_sac_set)*cos(sac_dirs(bb)*pi/180) + sac_deltaY(cur_sac_set)*sin(sac_dirs(bb)*pi/180));
    delta_sacpar(cur_sac_set) = cur_delta_sac_par;
end
% is_gsac = delta_sacpar' > sac_thresh;
is_gsac = sac_amps > sac_thresh;

%define which saccades to use
used_msacs = find(is_micro & ~is_blink);
used_gsacs = find(is_gsac & ~is_blink);
msac_start_inds = interp_sac_start_inds(used_msacs);
gsac_start_inds = interp_sac_start_inds(used_gsacs);

%% CREATE SACCADE PREDICTOR MATS
sac_backlag = 0.2;
sac_forlag = 0.5;
sac_bin_width = 1; %number of dt bins
sac_binspace = sac_bin_width*dt;
sac_bin_edges = -(sac_backlag-sac_binspace/2):sac_binspace:(sac_forlag+sac_binspace/2);
sac_bin_cents = 0.5*sac_bin_edges(1:end-1) + 0.5*sac_bin_edges(2:end);
n_sac_bins = length(sac_bin_cents);

%NOTE: Might want to enforce no inter-trial interaction criteria.
trial_gsac_mat = zeros(length(all_t_axis),n_sac_bins);
trial_msac_mat = zeros(length(all_t_axis),n_sac_bins);
for i = 1:n_sac_bins
    for bb = 1:sac_bin_width
        cur_inds = gsac_start_inds + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + bb;
        cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
        trial_gsac_mat(cur_inds,i) = 1;
        
        cur_inds = msac_start_inds + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + bb;
        cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
        trial_msac_mat(cur_inds,i) = 1;
    end
end


%% BLOCK PREDICTOR
Xblock = zeros(length(all_stim_times),length(cur_block_set));
for i = 1:length(cur_block_set)
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end
if Expt_num == 93
    Xwi = zeros(length(all_stim_times),2);
    cur_set = find(all_trial_wi == un_wi_vals(2));
    Xwi(cur_set,1) = 1;
    cur_set = find(all_trial_wi == un_wi_vals(3));
    Xwi(cur_set,2) = 1;
    
    Xblock = [Xblock Xwi];
end


%% CREATE STIM MATRIX
NT = length(used_inds);

un_oris = unique(all_stim_ori);
n_oris = length(un_oris);
all_stim_mat = zeros(length(all_stim_ori),n_oris);
for ii = 1:n_oris
    cur_set = find(all_stim_ori == un_oris(ii));
    all_stim_mat(cur_set,ii) = 1;
end
[Xinds,Tinds] = meshgrid(1:n_oris,1:flen);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix));

stim_params = NMMcreate_stim_params([flen n_oris],dt);
all_Xmat = create_time_embedding(all_stim_mat,stim_params);

%% FIT INITIAL MODELS
clear nmm_stim_params
nmm_stim_params(1) = stim_params;
% nmm_stim_params(2) = NMMcreate_stim_params(length(cur_block_set),dt);
% cur_X{1} = all_Xmat(used_inds,:);
% cur_X{2} = Xblock(used_inds,:);

L2_mat_params = create_L2_params([],[1 flen*n_oris],[flen n_oris],2,3,[Inf -1],[1 1]);
L2_mat = generate_L2_mat(L2_mat_params,flen*n_oris);

reg_params = NMMcreate_reg_params('lambda_custom',[1e2],'boundary_conds',[Inf Inf]);
silent = 1;
Fr1_inds = find(all_trial_Fr(all_trialvec(used_inds)) == 1);
Fr3_inds = find(all_trial_Fr(all_trialvec(used_inds)) == 3);
for ii = 1:96
    fprintf('Fitting model for unit %d of %d\n',ii,96);
    Robs = all_binned_spikes(used_inds,ii);
    
%     fit0(ii) = NMMinitialize_model(nmm_stim_params,[1 1],{'lin','lin'},reg_params,[1 2]);
%     fit0(ii) = NMMfit_filters(fit0(ii),Robs,cur_X,[],[],silent,[],L2_mat);    
    fit0(ii) = NMMinitialize_model(nmm_stim_params,[1],{'lin'},reg_params,[1]);
    fit0(ii) = NMMfit_filters(fit0(ii),Robs,all_Xmat(used_inds,:),[],[],silent,[],L2_mat);    
    fit0_Fr1(ii) = NMMinitialize_model(nmm_stim_params,[1],{'lin'},reg_params,[1]);
    fit0_Fr1(ii) = NMMfit_filters(fit0_Fr1(ii),Robs(Fr1_inds),all_Xmat(used_inds(Fr1_inds),:),[],[],silent,[],L2_mat);    
     
    fit0_Fr3(ii) = NMMinitialize_model(nmm_stim_params,[1],{'lin'},reg_params,[1]);
    fit0_Fr3(ii) = NMMfit_filters(fit0_Fr3(ii),Robs(Fr3_inds),all_Xmat(used_inds(Fr3_inds),:),[],[],silent,[],L2_mat);    

    avg_rate_Fr1(ii) = mean(Robs(Fr1_inds));
    avg_rate_Fr3(ii) = mean(Robs(Fr3_inds));
    
    [LL(ii), penLL, pred_rate, G, gint, fgint, nullLL(ii)] = NMMmodel_eval( fit0(ii), Robs, all_Xmat(used_inds,:));
    [LL_Fr1(ii), penLL, pred_rate, G, gint, fgint, nullLL_Fr1(ii)] = NMMmodel_eval( fit0_Fr1(ii), Robs(Fr1_inds), all_Xmat(used_inds(Fr1_inds),:));
    [LL_Fr3(ii), penLL, pred_rate, G, gint, fgint, nullLL_Fr3(ii)] = NMMmodel_eval( fit0_Fr3(ii), Robs(Fr3_inds), all_Xmat(used_inds(Fr3_inds),:));

end


%%
use_lag = 5;
use_tinds = find(Tinds(:) == use_lag);

stim_params = NMMcreate_stim_params([n_sac_bins n_oris],dt);
stim_params(2) = NMMcreate_stim_params([n_sac_bins 1],dt);
L2_mat_params = create_L2_params([],[1 n_sac_bins*n_oris],[n_sac_bins n_oris],2,3,[Inf -1],[1 .2]);
L2_mat = generate_L2_mat(L2_mat_params,flen*n_oris);

reg_params = NMMcreate_reg_params('lambda_custom',[30;0],'boundary_conds',[Inf -1;Inf Inf],'lambda_d2T',[0;200],'lambda_L2',[10;10]);

emb_stim = bsxfun(@times,reshape(all_Xmat(used_inds,use_tinds),...
    [length(used_inds) 1 length(use_tinds)]),trial_gsac_mat(used_inds,:));
emb_stim = reshape(emb_stim,size(emb_stim,1),n_sac_bins*n_oris);
cur_X{1} = emb_stim;
cur_X{2} = trial_gsac_mat(used_inds,:);

%%
for gg = 1:n_sac_bins
    cur_inds = find(trial_gsac_mat(used_inds,gg) == 1);
    sac_trg_mean_rates(:,gg) = mean(all_binned_spikes(used_inds(cur_inds),:));
end
%%
silent = 1;
Fr1_inds = find(all_trial_Fr(all_trialvec(used_inds)) == 1);
Fr3_inds = find(all_trial_Fr(all_trialvec(used_inds)) == 3);

for ii = 1:96;
    fprintf('Fitting model for unit %d of %d\n',ii,96);
    Robs = all_binned_spikes(used_inds,ii);
    
    %     fit1(ii) = NMMinitialize_model(stim_params,[1 1],{'lin','lin'},reg_params,[1 2]);
    %     fit1(ii) = NMMfit_filters(fit1(ii),Robs,cur_X,[],[],silent,[],L2_mat);
%     cur_X{1} = emb_stim(Fr1_inds,:);
%     cur_X{2} = trial_gsac_mat(used_inds(Fr1_inds),:);
%     fit1_Fr1(ii) = NMMinitialize_model(stim_params,[1 1],{'lin','lin'},reg_params,[1 2]);
%     fit1_Fr1(ii) = NMMfit_filters(fit1_Fr1(ii),Robs(Fr1_inds),cur_X,[],[],silent,[],L2_mat);
%     
%     cur_X{1} = emb_stim(Fr3_inds,:);
%     cur_X{2} = trial_gsac_mat(used_inds(Fr3_inds),:);
%     fit1_Fr3(ii) = NMMinitialize_model(stim_params,[1 1],{'lin','lin'},reg_params,[1 2]);
%     fit1_Fr3(ii) = NMMfit_filters(fit1_Fr3(ii),Robs(Fr3_inds),cur_X,[],[],silent,[],L2_mat);
end

%%
for ii = 1:96
    curfilt = reshape(fit0(ii).mods(1).filtK,flen,n_oris);
    [best_mag(ii),best_or(ii)] = max(curfilt(use_lag,:));
    [worst_mag(ii),worst_or(ii)] = min(curfilt(use_lag,:));
    
    curfilt = reshape(fit1(ii).mods(1).filtK,n_sac_bins,n_oris);
    best_tfilt(ii,:) = curfilt(:,best_or(ii));
    worst_tfilt(ii,:) = curfilt(:,worst_or(ii));
    
    sac_filt(ii,:) = fit1(ii).mods(2).filtK;
end

mod_tfilt = best_tfilt - worst_tfilt;