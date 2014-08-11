clear all
close all

fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';


Expt_name = 'G099';
Expt_num = str2num(Expt_name(2:end));
data_dir = ['/media/NTlab_data2/Data/bruce/' Expt_name];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];

save_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
if ~exist(save_dir,'dir')
    system(['mkdir ' save_dir]);
end

if Expt_name(1) == 'G'
    n_probes = 96;
    expt_type = 'UA';
elseif Expt_name(1) == 'M'
    n_probes = 24;
    expt_type = 'LP';
    switch Expt_num
        case 266
            bar_ori = 80;
        case 270
            bar_ori = 60;
        case 275
            bar_ori = 135;
        case 277
            bar_ori = 70;
    end
    
else
    error('Unrecognized experiment name');
end

if Expt_num == 275 || Expt_num == 277
    rpt_seed = 1001; %M275 M277
else
    rpt_seed = 1e4; %m270 and 266
end

if strcmp(Expt_name,'M270')
    ignore_blocks = [5];
elseif strcmp(Expt_name,'G087')
    ignore_blocks = [15]; %only 6 trials and causes problems
elseif strcmp(Expt_name,'G093')
    ignore_blocks = [28]; %only 6 trials and causes problems
else
    ignore_blocks = [];
end

use_coils = [1 0];
good_coils = [1 0];

%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz
dt = 0.005;

min_trial_dur = 0.75;
beg_buffer = 0.2;
end_buffer = 0.05;

backlag = round(0.3/dt);
forwardlag = round(0.6/dt);

sua_sm_sig = (0.01/dt);
mua_sm_sig = (0.01/dt);
% mua_sm_sig = 0;

if strcmp(Expt_name,'G081') || ismember(Expt_num,[232 235 239])
    trial_dur = 2;
else
    trial_dur = 4;
end

EP_bounds = 1;%eye position boundary
micro_thresh = 1; %microsaccade amplitude threshold (deg)
max_sac_dur = 0.1; %maximum saccade duration (otherwise likely a blink)
sac_burst_isi = 0.15; %minimum inter-saccade interval before classifying sac as part of a 'burst'
gsac_thresh = 1;

%% LOAD EXPTS STRUCT
cd(data_dir)
if Expt_name(1) == 'G'
load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
else
load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
end

%% PARSE EXPTS STRUCT
is_sim_msac_expt = false;
include_expts = 'none.fc';

has_data = cellfun(@(x) length(x),Expts) > 0;
expt_names(has_data) = cellfun(@(x) x.Header.expname,Expts(has_data),'UniformOutput',false);
expt_backstim(has_data) = cellfun(@(x) x.Stimvals.Bs,Expts(has_data),'UniformOutput',false);
expt_imback = strcmp(expt_backstim,'image');
expt_dd(has_data) = cellfun(@(x) x.Stimvals.dd,Expts(has_data),'UniformOutput',true);
expt_bar_ori(has_data) = cellfun(@(x) x.Stimvals.or,Expts(has_data),'UniformOutput',true);
expt_Fr(has_data) = cellfun(@(x) x.Stimvals.Fr,Expts(has_data),'UniformOutput',true);
expt_sac_dir(has_data) = mod(cellfun(@(x) x.Stimvals.Fa,Expts(has_data),'UniformOutput',true),180);
expt_sac_amp(has_data) = cellfun(@(x) x.Stimvals.Fs,Expts(has_data),'UniformOutput',true);
expt_ce(has_data) = cellfun(@(x) x.Stimvals.ce,Expts(has_data),'UniformOutput',true);
expt_ijump(has_data) = cellfun(@(x) x.Stimvals.ijump,Expts(has_data),'UniformOutput',true);

included_type(has_data) = ismember(expt_names(has_data),include_expts);
cur_block_set = find(included_type);

cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];


%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_ttill_end = [];
all_blockvec = [];
all_trialvec = [];
all_trial_back = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_trial_blocknums = [];
all_bin_edge_pts = [];

all_spk_times = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);

trial_toffset = zeros(length(cur_block_set),1);
cur_spkind_offset = 0;
cur_toffset = 0;
for ee = 1:length(cur_block_set);
    cur_block = cur_block_set(ee);
    fprintf('Expt %s Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_name,ee,length(cur_block_set));
    
    fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    load(fname,'Clusters');
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
        all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
    end

    trial_start_times = [Expts{cur_block}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_block}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_block}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_block}.Trials(:).id];
    
    [un_ids,id_inds] = unique(trial_ids);
        use_trials = find(trial_durs >= min_trial_dur);
%     if length(un_ids) < length(trial_ids)
%         fprintf('Warning, repeat trial inds detected!\n');
%         use_trials = [];
%     else
%         use_trials = find(trial_durs >= min_trial_dur);
%     end
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    cur_t_edges = (trial_start_times(1):dt:trial_end_times(end))';
    cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
    all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
    all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
    all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
    all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
    all_trialvec = [all_trialvec; ones(size(cur_t_axis))*ee];
        
    if strcmp(expt_type,'LP')
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
        cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
    end
end

trial_start_inds = [1; find(diff(all_trialvec) > 0)+1];

%% BIN SPIKES FOR MU AND SU
clust_params.n_probes = n_probes;
clust_params.exclude_adjacent = false;

[all_binned_mua,all_binned_sua,Clust_data] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params);
SU_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;
SU_isodist = Clust_data.SU_isodists;
SU_Lratio = Clust_data.SU_Lratios;

%% SMOOTH AND NORMALIZE SPIKING DATA
all_mua_rate = nan(size(all_binned_mua));
mua_block_mean_rates = nan(length(cur_block_set),n_probes);
mua_block_n_spikes = nan(length(cur_block_set),n_probes);
for ee = 1:length(cur_block_set)
    cur_block_inds = find(all_blockvec==ee);
    if ~isempty(cur_block_inds)
        if mua_sm_sig > 0
            for cc = 1:n_probes
                all_mua_rate(cur_block_inds,cc) = jmm_smooth_1d_cor(all_binned_mua(cur_block_inds,cc),mua_sm_sig);
            end
        else
            all_mua_rate(cur_block_inds,:) = all_binned_mua(cur_block_inds,:);
        end
        mua_block_mean_rates(ee,:) = mean(all_mua_rate(cur_block_inds,:));
        mua_block_n_spikes(ee,:) = sum(all_binned_mua(cur_block_inds,:));
    end
end

all_mua_rate_norm = bsxfun(@rdivide,all_mua_rate,mean(all_mua_rate));

%%
niqf = 1/2/dt;
lcf = 0.1;
[b,a] = butter(2,[lcf/niqf],'high');

for ii = 1:n_probes
    all_mua_rate_norm(:,ii) = filtfilt(b,a,all_mua_rate_norm(:,ii));
end
% all_sua_rate_norm = filtfilt(b,a,all_sua_rate_norm);
%% DEFINE DATA USED FOR ANALYSIS
% used_inds = find(all_tsince_start >= beg_buffer & all_ttill_end >= end_buffer);
used_inds = (1:length(all_t_axis))';

%% PROCESS EYE TRACKING DATA
if isfield(Expts{cur_block_set(1)}.Header,'exptno')
    em_block_nums = cellfun(@(X) X.Header.exptno,Expts(cur_block_set),'uniformoutput',1); %block numbering for EM/LFP data sometimes isnt aligned with Expts struct
else
    em_block_nums = cur_block_set;
end

[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,em_block_nums,Expt_name,trial_toffset,good_coils);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data_v2(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori(cur_block_set),used_inds);

[saccades,et_params] = detect_saccades_v2(corrected_eye_vals,all_eye_vals,all_eye_speed,all_eye_ts,et_params);

is_blink = detect_blinks(all_eye_ts,all_eye_vals,saccades,et_params);

[saccades,is_blink] = merge_blinks(saccades,is_blink);

sac_start_times = [saccades(:).start_time];
sac_stop_times = [saccades(:).stop_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));
interp_sac_stop_inds(isnan(interp_sac_stop_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > dt);
sac_start_times(bad_sacs) = [];
sac_stop_times(bad_sacs) = [];
saccades(bad_sacs) = [];
is_blink(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];
interp_sac_stop_inds(bad_sacs) = [];

% out_of_range = (abs(corrected_eye_vals_interp(used_inds,1)) > 8 | abs(corrected_eye_vals_interp(used_inds,2) > 8));
% fract_out = sum(out_of_range)/length(used_inds);
% fprintf('Eliminating %.4f of data out of window\n',fract_out);

% used_inds(out_of_range) = [];
%% PROCESS SACCADE STATS

%interpolate saccade start times and get rid of saccades that aren't within
%the t-axis binning
sac_start_times = [saccades(:).start_time];
sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
sac_start_inds(isnan(sac_start_inds)) = 1;
sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),[saccades(:).stop_time]));
sac_error = abs(sac_start_times - all_t_axis(sac_start_inds)');
usac_set = find(ismember(sac_start_inds,used_inds) & ismember(sac_stop_inds,used_inds) & sac_error <= dt & sac_start_times >= all_t_axis(1) & [saccades(:).stop_time] <= all_t_axis(end));
saccades = saccades(usac_set);
used_is_blink = is_blink(usac_set);

sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),[saccades(:).start_time]));
sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),[saccades(:).stop_time]));
sac_peak_inds = round(interp1(all_t_axis,1:length(all_t_axis),[saccades(:).peak_time]));

%start and stop positions of each saccade
sac_prepos = reshape([saccades(:).pre_pos],[],length(saccades));
sac_postpos = reshape([saccades(:).post_pos],[],length(saccades));
sac_amp = [saccades(:).amplitude];

% %find saccades that start or end out of window
% out_bounds = abs(sac_prepos(2,:)) > EP_bounds | abs(sac_postpos(2,:)) > EP_bounds;

sacburst_set = find([saccades(:).isi] < sac_burst_isi | [saccades(:).next_isi] < sac_burst_isi);
micro_set = find([saccades(:).amplitude] < micro_thresh & ~used_is_blink');
msac_bursts = micro_set(ismember(micro_set,sacburst_set));
micro_set(ismember(micro_set,sacburst_set)) = []; %eliminate microsacs that are part of a 'burst'

%guided saccades are those whose parallel component is large enough and
%that aren't blinks
gsac_set = find(sac_amp > gsac_thresh & ~used_is_blink');


%%
dark_blocks = [3 4 6 7];
light_blocks = setdiff(cur_block_set,dark_blocks);

light_sacs = find(ismember(all_blockvec(sac_start_inds),light_blocks));
dark_sacs = find(ismember(all_blockvec(sac_start_inds),dark_blocks));

light_msacs = intersect(micro_set,light_sacs);
light_gsacs = intersect(gsac_set,light_sacs);
dark_msacs = intersect(micro_set,dark_sacs);
dark_gsacs = intersect(gsac_set,dark_sacs);

% n_amp_bins = 4;
% usac_set = find([saccades(:).duration] < max_sac_dur);
% sac_bin_edges = prctile([saccades(usac_set).amplitude],linspace(0,100,n_amp_bins+1));
% sac_bin_cents = 0.5*sac_bin_edges(1:end-1) + 0.5*sac_bin_edges(2:end);
% 
% pre_pos = reshape([saccades(:).pre_pos],2,[]);
% post_pos = reshape([saccades(:).post_pos],2,[]);
% pre_amp = sum(sqrt(pre_pos.^2));
% post_amp = sum(sqrt(post_pos.^2));
% in_bounds = find(pre_amp < 6 & post_amp < 6);
% 
% light_gsac_in = intersect(light_gsacs,in_bounds);
% light_msac_in = intersect(light_msacs,in_bounds);
% dark_gsac_in = intersect(dark_gsacs,in_bounds);


%% COMPUTE TRIG AVGS FOR MUA
nboot = [];

% %general averages
% [mua_data.msac_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(micro_set),backlag,forwardlag,nboot);
% [mua_data.gsac_avg,lags] = get_event_trig_avg(all_mua_rate_norm,sac_start_inds(gsac_set),backlag,forwardlag,nboot);

[light_msac_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(light_msacs),backlag,forwardlag,nboot);
[light_gsac_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(light_gsacs),backlag,forwardlag,nboot);
[dark_msac_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(dark_msacs),backlag,forwardlag,nboot);
[dark_gsac_avg,lags] = get_event_trig_avg_v3(all_mua_rate_norm,sac_start_inds(dark_gsacs),backlag,forwardlag,nboot);

%%
bad_units = [];
uset = setdiff(1:96,bad_units);
xl = [-0.3 0.6];

f1 = figure();
h1=shadedErrorBar(lags*dt,1+mean(dark_gsac_avg(:,uset),2),std(dark_gsac_avg(:,uset),[],2)/sqrt(length(uset)));
hold on
% h2=shadedErrorBar(lags*dt,mean(light_gsac_avg(:,uset)'),std(light_gsac_avg(:,uset)')/sqrt(length(uset)),{'color','r'});
ylim([0.9 1.1])
yl = ylim();
line(xl,[1 1],'color','k')
xlim(xl);
xlabel('Time (s)');
ylabel('Relative rate');

fig_width = 3.5; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'DarkSacMUA.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);
