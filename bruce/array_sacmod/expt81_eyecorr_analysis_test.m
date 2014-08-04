clear all
% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';

Expt_num = 81;

%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz
dt = 0.01;

min_trial_dur = 0.5;
beg_buffer = 0.3;

backlag = round(0.25/dt);
forwardlag = round(0.55/dt);

sua_sm_sig = (0.01/dt);
mua_sm_sig = (0.005/dt);

flen = 10;
use_nPix = 17;
full_nPix = 17;
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;
trial_dur = 2;
bar_ori = 0;

anal_name = 'hbar_stcprobmod_sacmod_analysis';
eye_data_name = 'monoc_eyecorr_hbar_eyeonly';

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

%% LOAD DATA FOR EXPT 81
load ~/Data/bruce/G081/all_un_bar_pos
n_bar_pos = size(all_un_bar_pos,1);
all_un_bar_pos = all_un_bar_pos(:,[1 3]);

%%
for i = 1:length(Expts)
    if strcmp(Expts{i}.Header.expname,'grating.OpXseRC') | strcmp(Expts{i}.Header.expname,'grating.OpRC')
        is_bar_expt(i) = 1;
    else
        is_bar_expt(i) = 0;
    end
    
    if strcmp(Expts{i}.Stimvals.Bs,'image')
        expt_imback(i) = 1;
    else
        expt_imback(i) = 0;
    end
    
    expt_sim_sacs(i) = Expts{i}.Stimvals.ijump;
    expt_bar_ori(i) = Expts{i}.Stimvals.or;
end
expt_has_ds = (expt_sim_sacs == 0);
expt_bar_ori(expt_bar_ori == -45) = 135;
expt_binoc = zeros(size(expt_bar_ori));
expt_imback = expt_imback';
cur_block_set = find(is_bar_expt & expt_bar_ori == bar_ori);

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
    %         trial_durs = trial_durs(id_inds);
    %         trial_start_times = trial_start_times(id_inds);
    %         trial_end_times = trial_end_times(id_inds);
    
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
    
    if strcmp(Expt_name,'G093')
        trial_wi = [Expts{cur_block}.Trials(:).wi];
        trial_wi = trial_wi(id_inds);
        all_trial_wi = cat(1,all_trial_wi,trial_wi(use_trials)');
    end
    
    
    
    if expt_bar_ori(cur_block) == 0
        un_bar_pos = all_un_bar_pos(:,1);
    elseif expt_bar_ori(cur_block) == 90
        un_bar_pos = all_un_bar_pos(:,2);
    else
        error('Only supposed to have 0 and 90 here');
    end
    n_bar_pos = length(un_bar_pos);
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start/1e4;
        cur_Op = Expts{cur_block}.Trials(use_trials(tt)).Op;
        cur_phase = Expts{cur_block}.Trials(use_trials(tt)).ph;
        
        cur_t_edges = [cur_stim_times; Expts{cur_block}.Trials(use_trials(tt)).End(end)/1e4];
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        cur_bar_mat = zeros(length(cur_stim_times),n_bar_pos);
        for bb = 1:n_bar_pos
            cur_set = find(cur_Op==un_bar_pos(bb));
            pset = cur_phase(cur_set) == 0;
            nset = cur_phase(cur_set) == pi;
            cur_bar_mat(cur_set(pset),bb) = 1;
            cur_bar_mat(cur_set(nset),bb) = -1;
        end
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_t_axis = [all_t_axis; cur_t_axis];
        all_stim_mat = [all_stim_mat; cur_bar_mat];
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
sim_sac_times = [0.7 1.4];
sac_thresh = 1;
sac_start_expts = all_blockvec(interp_sac_start_inds);
delta_sacpar = nan(size(saccades));
for bb = 1:length(sac_dirs)
    cur_bar_expts = find(expt_bar_ori(cur_block_set) == sac_dirs(bb));
    cur_sac_set = find(ismember(sac_start_expts,cur_bar_expts));
    
    cur_delta_sac_par = abs(sac_deltaX(cur_sac_set)*cos(sac_dirs(bb)*pi/180) + sac_deltaY(cur_sac_set)*sin(sac_dirs(bb)*pi/180));
    delta_sacpar(cur_sac_set) = cur_delta_sac_par;
end
is_gsac = delta_sacpar' > sac_thresh;

all_sim_sacs = [];
sim_expt_inds = find(ismember(all_blockvec,sim_sac_expts));
sim_sacs = cell(length(sim_sac_times),1);
for ii = 1:length(sim_sac_times)
    sim_sacs{ii} = sim_expt_inds(all_tsince_start(sim_expt_inds(1:end-1)) < sim_sac_times(ii) & ...
        all_tsince_start(sim_expt_inds(2:end)) >= sim_sac_times(ii));
    all_sim_sacs = [all_sim_sacs; sim_sacs{ii}];
end

%define which saccades to use
used_msacs = find(is_micro & ~is_blink);
used_gsacs = find(is_gsac & ~is_blink);
msac_start_inds = interp_sac_start_inds(used_msacs);
gsac_start_inds = interp_sac_start_inds(used_gsacs);

%% CREATE SACCADE PREDICTOR MATS
sac_backlag = 0.2;
sac_forlag = 0.45;
sac_bin_width = 2; %number of dt bins
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

[Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix));

%% LOAD INFERRED EYE POSITIONS
NT = length(used_inds);
cd(anal_dir)
load(eye_data_name);

sac_shift = round(0.05/dt);

saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));

sp_dx = 0.05;
inferred_eyepos = inferred_eyepos_cor/sp_dx;
inferred_eyepos = round(interp1(inferred_eyepos_taxis,inferred_eyepos,all_t_axis(used_inds)));
inferred_eyepos_cor = inferred_eyepos;

%align jumps to saccade times
for i = length(saccade_start_inds):-1:1
    cur_inds = saccade_start_inds(i):(saccade_start_inds(i)+sac_shift);
    inferred_eyepos_cor(cur_inds) = inferred_eyepos(cur_inds(end));
    cur_back_inds = (saccade_start_inds(i)-(flen-sac_shift)):(saccade_start_inds(i)-1);
    cur_back_inds(cur_back_inds < 1) = [];
    inferred_eyepos_cor(cur_back_inds) = inferred_eyepos(saccade_start_inds(i));
end

interp_eyepos_cor = zeros(length(all_stim_times),1);
interp_eyepos_cor(used_inds) = inferred_eyepos_cor;
interp_eyepos_cor(isnan(interp_eyepos_cor)) = 0;

all_stim_mat_cor = all_stim_mat;
for ii = 1:NT
    all_stim_mat_cor(used_inds(ii),:) = shift_matrix_Nd(all_stim_mat(used_inds(ii),:),...
        -interp_eyepos_cor(used_inds(ii)),2);
end
all_Xmat_cor = create_time_embedding(all_stim_mat_cor,stim_params);
all_Xmat_cor = all_Xmat_cor(used_inds,use_kInds);
all_Xmat_uncor = create_time_embedding(all_stim_mat,stim_params);
all_Xmat_uncor = all_Xmat_uncor(used_inds,use_kInds);


%% SUA STA/STC ANALYSIS
addpath('~/James_scripts/TentBasis2D/')
nneg = 10; npos = 10;

nim_stim_params = NMMcreate_stim_params([flen use_nPix],dt);
nim_stim_params(2) = NMMcreate_stim_params([flen use_nPix],dt);
reg_params = NMMcreate_reg_params('lambda_d2XT',25);
for ss = 1:length(su_probes)
    %%
    fprintf('Computing STC for SU %d of %d\n',ss,length(su_probes));
    
    cur_used_blocks = find(su_used_blocks(:,ss)); %blocks when SU
    cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
    cur_used_inds = find(ismember(used_inds,cur_poss_inds));
    cur_NT = length(cur_used_inds);
        
    if ~isempty(cur_used_inds)
      
        Robs = all_binned_spikes(used_inds(cur_used_inds),ss);
        avg_rate = mean(Robs);
        
        spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = all_Xmat_uncor(cur_used_inds(spikebins),:);
        sta      = mean(spike_cond_stim) - mean(all_Xmat_uncor(cur_used_inds,:));
        sta_norm = norm(sta);
        sta = sta/norm(sta);
        spike_cond_stim = abs(all_Xmat_uncor(cur_used_inds(spikebins),:));
        abs_sta      = mean(spike_cond_stim) - mean(abs(all_Xmat_uncor(cur_used_inds,:)));
        abs_sta_norm = norm(abs_sta);
        abs_sta = abs_sta/norm(abs_sta);
        proj_mat = sta'/(sta*sta')*sta;
        stim_proj = all_Xmat_uncor(cur_used_inds,:) - all_Xmat_uncor(cur_used_inds,:)*proj_mat;
        % stim_proj = stim_emb;
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
        
        Xcell{1} = all_Xmat_cor(cur_used_inds,:);
        Xcell{2} = abs(all_Xmat_cor(cur_used_inds,:));
        fit0 = NMMinitialize_model(nim_stim_params,[1 1],{'lin','lin'},reg_params,[1 2]);
        fit0 = NMMfit_filters(fit0,Robs,Xcell);
        sua_data(ss).nimFit = fit0;
        
        sua_data(ss).sta_unc = sta;
        sua_data(ss).stc_unc = stcs;
        sua_data(ss).evals_unc = diag(evals);
        sua_data(ss).sta_norm_unc = sta_norm;
        sua_data(ss).abs_sta_unc = abs_sta;
        sua_data(ss).abs_sta_norm_unc = abs_sta_norm;
 
        spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = all_Xmat_cor(cur_used_inds(spikebins),:);
        sta      = mean(spike_cond_stim) - mean(all_Xmat_cor(cur_used_inds,:));
        sta_norm = norm(sta);
        sta = sta/norm(sta);
        spike_cond_stim = abs(all_Xmat_cor(cur_used_inds(spikebins),:));
        abs_sta      = mean(spike_cond_stim) - mean(abs(all_Xmat_cor(cur_used_inds,:)));
        abs_sta_norm = norm(abs_sta);
        abs_sta = abs_sta/norm(abs_sta);
        proj_mat = sta'/(sta*sta')*sta;
        stim_proj = all_Xmat_cor(cur_used_inds,:) - all_Xmat_cor(cur_used_inds,:)*proj_mat;
        % stim_proj = stim_emb;
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
        
        sua_data(ss).sta = sta;
        sua_data(ss).stc = stcs;
        sua_data(ss).evals = diag(evals);
        sua_data(ss).sta_norm = sta_norm;
        sua_data(ss).abs_sta = abs_sta;
        sua_data(ss).abs_sta_norm = abs_sta_norm;

    end
end

%%
ss = 7;
use_TInds = find(Tinds(:) == 6);
all_Xmat_cor = create_time_embedding(all_stim_mat_cor,stim_params);
all_Xmat_cor = all_Xmat_cor(used_inds,use_TInds);
% all_Xmat_cor = all_Xmat_cor(used_inds,use_kInds);
% all_Xmat_uncor = create_time_embedding(all_stim_mat,stim_params);
% all_Xmat_uncor = all_Xmat_uncor(used_inds,use_TInds);

cur_used_blocks = find(su_used_blocks(:,ss)); %blocks when SU
cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
cur_used_inds = find(ismember(used_inds,cur_poss_inds));
cur_NT = length(cur_used_inds);
Robs = all_binned_spikes(used_inds(cur_used_inds),ss);
avg_rate = mean(Robs);

% for gg = 1:length(sac_bin_cents)
%     gused_inds = find(trial_gsac_mat(used_inds(cur_used_inds),gg) == 1);
%     gsac_rate(gg) = mean(Robs(gused_inds));
%     
%     
%     spikebins = convert_to_spikebins(Robs(gused_inds));
%     spike_cond_stim = all_Xmat_cor(cur_used_inds(gused_inds(spikebins)),:);
%     sta      = mean(spike_cond_stim) - mean(all_Xmat_cor(cur_used_inds(gused_inds),:));
%     sta_norm = norm(sta);
%     sta = sta/norm(sta);
%     spike_cond_stim = abs(all_Xmat_cor(cur_used_inds(gused_inds(spikebins)),:));
%     abs_sta      = mean(spike_cond_stim) - mean(abs(all_Xmat_cor(cur_used_inds(gused_inds),:)));
%     abs_sta_norm = norm(abs_sta);
%     abs_sta = abs_sta/norm(abs_sta);
%    
% %     gsta(gg,:,:) = reshape(sta,flen,use_nPix);
% %     gsta_abs(gg,:,:) = reshape(abs_sta,flen,use_nPix);
%     gsta2(gg,:) = sta;
%     gsta2_abs(gg,:) = abs_sta;
% 
% 
%     gused_inds = find(trial_msac_mat(used_inds(cur_used_inds),gg) == 1);
%     msac_rate(gg) = mean(Robs(gused_inds));
%     
%     
%     spikebins = convert_to_spikebins(Robs(gused_inds));
%     spike_cond_stim = all_Xmat_cor(cur_used_inds(gused_inds(spikebins)),:);
%     sta      = mean(spike_cond_stim) - mean(all_Xmat_cor(cur_used_inds(gused_inds),:));
%     sta_norm = norm(sta);
%     sta = sta/norm(sta);
%     spike_cond_stim = abs(all_Xmat_cor(cur_used_inds(gused_inds(spikebins)),:));
%     abs_sta      = mean(spike_cond_stim) - mean(abs(all_Xmat_cor(cur_used_inds(gused_inds),:)));
%     abs_sta_norm = norm(abs_sta);
%     abs_sta = abs_sta/norm(abs_sta);
%    
% %     msta(gg,:,:) = reshape(sta,flen,use_nPix);
% %     msta_abs(gg,:,:) = reshape(abs_sta,flen,use_nPix);
%     msta2(gg,:) = sta;
%     msta2_abs(gg,:) = abs_sta;
% end
% 
any_msac = max(trial_msac_mat,[],2);
cur_used_inds = find(ismember(used_inds,cur_poss_inds));
cur_used_inds = cur_used_inds(any_msac(used_inds(cur_used_inds)) == 1);
temp = reshape(trial_msac_mat(used_inds(cur_used_inds),:),[length(cur_used_inds) 1 n_sac_bins]);
tempX = bsxfun(@times,all_Xmat_cor(cur_used_inds,:),temp);
tempX = reshape(tempX,length(cur_used_inds),use_nPix*n_sac_bins);
Robs = all_binned_spikes(used_inds(cur_used_inds),ss);

temp_stim_params = NMMcreate_stim_params([use_nPix n_sac_bins],dt);
temp_stim_params(2) = NMMcreate_stim_params([use_nPix n_sac_bins],dt);
% temp_stim_params(3) = NMMcreate_stim_params([n_sac_bins],dt);
temp_reg_params = NMMcreate_reg_params('lambda_d2XT',[20 20],'lambda_L2',[100 100]);
% temp_reg_params = NMMcreate_reg_params('lambda_d2T',[0 50]);
tempX_cell{1} = tempX; 
tempX_cell{2} = abs(tempX); 
% tempX_cell{3} = trial_msac_mat(used_inds(cur_used_inds),:);
% temp_nim = NMMinitialize_model(temp_stim_params,[1 1 1],{'lin','lin' 'lin'},temp_reg_params,[1 2 3]);
% init_filts{1} = zeros(size(tempX,2),1); 
% init_filts{1} = 0.001*randn(size(tempX,2),1); 
% init_filts{2} = [];
temp_nim = NMMinitialize_model(temp_stim_params,[1 1],{'lin','lin'},temp_reg_params,[1 2]);
% gtemp_nim = NMMadjust_regularization(gtemp_nim,[1],'spatial_boundaries','free');
% gtemp_nim = NMMadjust_regularization(gtemp_nim,[2],'lambda_d2XT',1);
optim_p.progTol = 1e-8;
optim_p.optTol = 1e-5;
optim_p.maxIter = 150;
temp_nim = NMMfit_filters(temp_nim,Robs,tempX_cell,[],[],0,optim_p);


any_gsac = max(trial_gsac_mat,[],2);
cur_used_inds = find(ismember(used_inds,cur_poss_inds));
cur_used_inds = cur_used_inds(any_gsac(used_inds(cur_used_inds)) == 1);
temp = reshape(trial_gsac_mat(used_inds(cur_used_inds),:),[length(cur_used_inds) 1 n_sac_bins]);
tempX = bsxfun(@times,all_Xmat_cor(cur_used_inds,:),temp);
tempX = reshape(tempX,length(cur_used_inds),use_nPix*n_sac_bins);

Robs = all_binned_spikes(used_inds(cur_used_inds),ss);
avg_rate = mean(Robs);


temp_stim_params = NMMcreate_stim_params([use_nPix n_sac_bins],dt);
temp_stim_params(2) = NMMcreate_stim_params([use_nPix n_sac_bins],dt);
temp_stim_params(3) = NMMcreate_stim_params([n_sac_bins],dt);
temp_reg_params = NMMcreate_reg_params('lambda_d2XT',[5 5 0],'lambda_d2T',[0 0 50]);
% temp_reg_params = NMMcreate_reg_params('lambda_d2T',[0 50]);
tempX_cell{1} = tempX; 
tempX_cell{2} = abs(tempX); 
tempX_cell{3} = trial_gsac_mat(used_inds(cur_used_inds),:);
% temp_nim = NMMinitialize_model(temp_stim_params,[1 1 1],{'lin','lin' 'lin'},temp_reg_params,[1 2 3]);
% init_filts{1} = zeros(size(tempX,2),1); 
% init_filts{1} = 0.001*randn(size(tempX,2),1); 
% init_filts{2} = [];
gtemp_nim = NMMinitialize_model(temp_stim_params,[1 1 1],{'lin','lin','lin'},temp_reg_params,[1 2 3]);
% gtemp_nim = NMMadjust_regularization(gtemp_nim,[1],'spatial_boundaries','free');
% gtemp_nim = NMMadjust_regularization(gtemp_nim,[2],'lambda_d2XT',1);
optim_p.progTol = 1e-8;
optim_p.optTol = 1e-5;
optim_p.maxIter = 150;
gtemp_nim = NMMfit_filters(gtemp_nim,Robs,tempX_cell,[],[],0,optim_p);

%%
close all
% f1 = figure('Name','unc');
f2 = figure('Name','cor');
f3 = figure('Name','NMM');
for ss = 1:length(su_probes)
    
    %%
    if ~isempty(sua_data(ss).sta)
        fprintf('SU %d of %d\n',ss,length(su_probes));
        
%         figure(f1); clf
%         cur_sta = sua_data(ss).sta_unc; ca = max(abs(cur_sta));
%         subplot(3,3,1)
%         imagesc(reshape(cur_sta,flen,use_nPix)); caxis([-ca ca]);
%         cur_sta = sua_data(ss).abs_sta_unc; ca = max(abs(cur_sta));
%         subplot(3,3,2)
%         imagesc(reshape(cur_sta,flen,use_nPix)); caxis([-ca ca]);
%         cur_stcs = sua_data(ss).stc_unc; ca = max(abs(cur_stcs(:)));
%         for ii = 1:3
%             subplot(3,3,3+ii)
%             imagesc(reshape(cur_stcs(:,ii),flen,use_nPix)); caxis([-ca ca]);
%         end
%         for ii = 1:3
%             subplot(3,3,6+ii)
%             imagesc(reshape(cur_stcs(:,end-nneg+ii),flen,use_nPix)); caxis([-ca ca]);
%         end
%         
        figure(f2); clf
        cur_sta = sua_data(ss).sta; ca = max(abs(cur_sta));
        subplot(3,3,1)
        imagesc(reshape(cur_sta,flen,use_nPix)); caxis([-ca ca]);
        cur_sta = sua_data(ss).abs_sta; ca = max(abs(cur_sta));
        subplot(3,3,2)
        imagesc(reshape(cur_sta,flen,use_nPix)); caxis([-ca ca]);
        cur_stcs = sua_data(ss).stc; ca = max(abs(cur_stcs(:)));
        for ii = 1:3
            subplot(3,3,3+ii)
            imagesc(reshape(cur_stcs(:,ii),flen,use_nPix)); caxis([-ca ca]);
        end
        for ii = 1:3
            subplot(3,3,6+ii)
            imagesc(reshape(cur_stcs(:,end-nneg+ii),flen,use_nPix)); caxis([-ca ca]);
        end
 
        figure(f3); clf
        subplot(1,2,1)
        cur_filt = sua_data(ss).nimFit.mods(1).filtK;ca = max(abs(cur_filt));
        imagesc(reshape(cur_filt,flen,use_nPix)); caxis([-ca ca]);
         subplot(1,2,2)
        cur_filt = sua_data(ss).nimFit.mods(2).filtK;ca = max(abs(cur_filt));
        imagesc(reshape(cur_filt,flen,use_nPix)); caxis([-ca ca]);
       pause
    end
end