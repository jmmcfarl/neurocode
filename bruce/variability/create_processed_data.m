
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements//');
addpath('~/James_scripts/bruce/processing/');
addpath('~/James_scripts/bruce/saccade_modulation/');

global Expt_name bar_ori monk_name rec_type

Expt_name = 'M270';
monk_name = 'lem';
bar_ori = 60; %bar orientation to use (only for UA or single-ori-LP recs)
rec_number = 1;

%if recording was split into two segments at some ori
if strcmp(Expt_name,'M011')
    rec_block_range =1:21; %M011
elseif strcmp(Expt_name,'M012')
    rec_block_range = 1:27;
else
    rec_block_range = nan;
end

Expt_num = str2num(Expt_name(2:end));

data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir')
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
end
if ~exist(data_dir,'dir');
    error('Couldnt find data directory');
end
Edata_file = strcat(data_dir,'/',monk_name,Expt_name,'Expts');
load(Edata_file);

%is this a laminar probe or utah array rec?
firste = find(cellfun(@(x) ~isempty(x),Expts),1);
if strcmp(Expts{firste}.Header.DataType,'GridData 96')
    rec_type = 'UA';
elseif strcmp(Expts{firste}.Header.DataType,'Spike2')
    rec_type = 'LP';
end

%% set params
params.micro_thresh = 1; %max amp of microsac (deg)
params.EP_bounds = 1;%eye position boundary (deg from central FP)
params.sac_burst_isi = 0.15; %inter-saccade interval for classifying microsaccade bursts
params.max_gsac_dur = 0.1; %maximum duration before we label a saccade a blink

params.min_trial_dur = 0.75; %in sec
params.stim_fs = 100; %in Hz
params.dt = 0.01; %time resolution (s)
params.Fr = 1; %number of repeats per bar pattern (frame rate of stim)

params.full_nPix=36; %default total number of bars per frame to store

%exclude data at beginning and end of each trial (in sec)
params.beg_buffer = 0.2;
params.end_buffer = 0.05;
params.trial_dur = 4;

%%
% [266-80 270-60 275-135 277-70 281-140 287-90 289-160 294-40 296-45 297-0/90 010-60 11-160 12-0 13-100 14-40 320-100]

%directory of bar orientations used for LP recs
if strcmp(rec_type,'LP')
    switch Expt_num
        case 266
            bar_ori = 80;
        case 270
            bar_ori = 60;
        case 275
            bar_ori = 135;
        case 277
            bar_ori = 70;
        case 281
            bar_ori = 140;
        case 287
            bar_ori = 90;
        case 289
            bar_ori = 160;
        case 294
            bar_ori = 40;
        case 296
            bar_ori = 45;
        case 297
            if ~ismember(bar_ori,[0 90])
                error('M297 is either 0 or 90 deg bars');
            end            
        case 309
            bar_ori = 120;
        case 5
            bar_ori = 50;
        case 9
            bar_ori = 0;
        case 10
            bar_ori = 60;
        case 11
            bar_ori = 160;
        case 12 
            bar_ori = 0;
        case 13
            bar_ori = 100;
        case 14
            bar_ori = 40;
        case 320
            bar_ori = 100;
    end
end

if strcmp(rec_type,'LP')
    params.n_probes = 24;
elseif strcmp(rec_type,'UA')
    params.n_probes = 96;
end
if strcmp(monk_name,'jbe')
    params.good_coils = [1 0];
        params.use_coils = [0 0];
elseif strcmp(monk_name,'lem');
    params.good_coils = [1 1];
    params.use_coils = [1 1];
end

cd(data_dir);
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
if ~exist(anal_dir,'dir')
    mkdir(anal_dir)
end

%set directory names
et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
mod_data_dir = ['~/Analysis/bruce/' Expt_name '/models'];

if rec_number > 1
    rec_block_range = setdiff(1:length(Expts),rec_block_range);
    cluster_dir = [cluster_dir sprintf('/rec%d',rec_number)];
end

%dont fit stim models using these blocks
ignore_blocks = [];
switch Expt_num
    case 270
        ignore_blocks = [5 19];
    case 289
        ignore_blocks = [25 26 27 28 29]; %25 and 26 have different dw, and 27-29 are clearly off somehow
    case 294
        ignore_blocks = [37 38 39]; %37-39 have slightly different dw used in these blocks
    case 86
        ignore_blocks = [16 17 28 30];
    case 87
        ignore_blocks = [15];
    case 93
        ignore_blocks = [28];
    case 11
        ignore_blocks = [29]; %different 
end

%problem with M270 where vd was wrong, need this correction factor to get
%correct units
if Expt_num==270
    params.scale_fac = 1.72;
else
    params.scale_fac = 1;
end

%some experiments have different conditions trial-interleaved
params.is_TBT_expt = false;
if Expt_num >= 275
    params.is_TBT_expt = true;
end

%if using a different # bar pixs in this expt
switch Expt_num
    case 270
        params.full_nPix=32;
    case  287
        params.full_nPix = 22;
    case 289
        params.full_nPix = 22;
    case 294
        params.full_nPix = 20;
        %     case 296
        %         full_nPix = 54;
end

%% SELECT BLOCKS FOR ANALYSIS
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs','rls.AllSac','rls.imiXFa','rls.FaXwi','rls.FaXwiXimi','rls.AllSacB','rls.froNoise'};
expt_names = cell(1,length(Expts));
expt_dds = nan(1,length(Expts));
expt_bar_ori = nan(1,length(Expts));
expt_sac_dir = nan(1,length(Expts));
expt_Fr = nan(1,length(Expts));
expt_sac_amp = nan(1,length(Expts));
expt_imback = nan(1,length(Expts));
included_type = false(1,length(Expts));
for ii = 1:length(Expts)
    if ~isempty(Expts{ii})
        expt_names{ii} = Expts{ii}.Header.expname;
        expt_dds(ii) = Expts{ii}.Stimvals.dd; %bar density
        expt_bar_ori(ii) = Expts{ii}.Stimvals.or; %bar ori
        expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180); %guided saccade dir
        expt_Fr(ii) = Expts{ii}.Stimvals.Fr; %bar frame rate
        expt_sac_amp(ii) = Expts{ii}.Stimvals.Fs; %saccade amplitude
        expt_imback(ii) = isfield(Expts{ii}.Trials,'imi'); %background type
        included_type(ii) = any(strcmp(expt_names{ii},include_expts));
    end
end
expt_has_ds(isnan(expt_has_ds)) = 0;
expt_has_ds(expt_has_ds == -1) = 0;
expt_binoc(isnan(expt_binoc)) = 0;

if strcmp(rec_type,'LP')
    expt_bar_ori(expt_bar_ori > 360) = bar_ori;
end

cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);
cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];
if length(unique(expt_dds(cur_block_set))) > 1
    fprintf('Warning, multiple dds detected!\n');
    main_dds = mode(expt_dds(cur_block_set));
    cur_block_set(expt_dds(cur_block_set) ~= main_dds) = [];
end

if ~isnan(rec_block_range)
cur_block_set(~ismember(cur_block_set,rec_block_range)) = []; %if specifying a usable block range
end

sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

n_blocks = length(cur_block_set);

all_nfs = cellfun(@(x) x.Stimvals.nf,Expts(cur_block_set));
if length(unique(all_nfs)) > 1
    fprintf('Warning, multiple different nfs detected: %.4f\n',all_nfs);
end

%set amp of guided saccades if there were any
if ~isempty(grayback_gs_expts) || ~isempty(imback_gs_expts) 
    params.gsac_amp = unique(expt_sac_amp(cur_block_set([grayback_gs_expts; imback_gs_expts])));
else
    params.gsac_amp = unique(expt_sac_amp(cur_block_set));
end
if length(params.gsac_amp) > 1
    fprintf('Multiple guided sac amps detected!\n');
end
%minimum (parallel) amplitude for a guided saccade to be included in
%analysis
params.gsac_thresh = mean(params.gsac_amp)/2;

%store expt data
expt_data.used_blocks = cur_block_set;
expt_data.expt_names = expt_names(cur_block_set);
expt_data.sim_sac_expts = sim_sac_expts;
expt_data.imback_gs_expts = imback_gs_expts;
expt_data.grayback_gs_expts = grayback_gs_expts;
expt_data.expt_dds = expt_dds(cur_block_set);
expt_data.expt_sac_dir = expt_sac_dir(cur_block_set);
expt_data.expt_Fr = expt_Fr(cur_block_set);
expt_data.expt_nf = all_nfs;
expt_data.expt_wi = cellfun(@(x) x.Stimvals.wi,Expts(cur_block_set));
expt_data.expt_dw = cellfun(@(x) x.Stimvals.dw,Expts(cur_block_set));

cur_expt_npix = unique(expt_npix(cur_block_set));
if length(cur_expt_npix) > 1
    warning('multiple Npix detected');
    cur_expt_npix = mode(cur_expt_npix);
end
if params.full_nPix > cur_expt_npix
    params.full_nPix = cur_expt_npix;
end

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

%this is a rediculous way of grabbing all this info...
all_stim_times = [];
all_stim_mat = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_wi = [];
all_trial_back = [];
all_trial_Ff = [];
all_trial_exvals = [];
all_trial_blk = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_trial_rptframes = [];
all_trial_nrptframes = [];
all_spk_times = cell(params.n_probes,1);
all_clust_ids = cell(params.n_probes,1);
all_spk_inds = cell(params.n_probes,1);
trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
cur_spkind_offset = 0;
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
    
    %get spiking data from this block
    fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    load(fname,'Clusters');
    for cc = 1:params.n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
        all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
    end
    
    %get trial start and end times
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
        use_trials = find(trial_durs >= params.min_trial_dur);
    end
    
    %add in time offset 
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    %get seed values for each used trial
    trial_Se = [Expts{cur_block}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    all_trial_blk = cat(1,all_trial_blk,ones(length(use_trials),1)*ee);
    
    %get stimulus area width
    if isfield(Expts{cur_block}.Trials,'wi')
        trial_wi = [Expts{cur_block}.Trials(:).wi];
        trial_wi = trial_wi(id_inds);
        all_trial_wi = cat(1,all_trial_wi,trial_wi(use_trials)');
    else
        all_trial_wi = cat(1,all_trial_wi,nan(length(use_trials),1));
    end
    
    %if trial-interleaved conditions, get the variable trial properties
    if params.is_TBT_expt
        if isfield(Expts{cur_block}.Trials,'Bs')
            trial_back = strcmp('image',{Expts{cur_block}.Trials(:).Bs});
        else
            trial_back = nan(1,length(use_trials));
        end
        trial_back = trial_back(id_inds); %background type
        
        %this stores all trial-varying experimental values
        if isfield(Expts{cur_block}.Trials,'exvals')
            exvals = reshape([Expts{cur_block}.Trials(:).exvals],length(Expts{cur_block}.Trials(1).exvals),[]);
            trial_exvals = exvals(:,id_inds)';
        else
            trial_exvals = nan(length(id_inds),3); 
        end
        
        if isfield(Expts{cur_block}.Trials,'Ff')
            trial_Ff = [Expts{cur_block}.Trials(:).Ff];
            trial_Ff = trial_Ff(id_inds);
        else
            trial_Ff = nan(1,length(id_inds));
        end
        all_trial_back = cat(1,all_trial_back,trial_back(use_trials)');
        all_trial_Ff = cat(1,all_trial_Ff,trial_Ff(use_trials)');
        all_trial_exvals = cat(1,all_trial_exvals,trial_exvals(use_trials,:));
    end
    
    %load in stimulus data
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    buffer_pix = floor((expt_npix(cur_block) - params.full_nPix)/2);
    if buffer_pix == -1 %in case there is slightly less pixels than desired by full_nPix, just buffer by a couple zeros
        for ii = 1:length(left_stim_mats)
            left_stim_mats{ii} = [zeros(size(left_stim_mats{ii},1),1) left_stim_mats{ii} zeros(size(left_stim_mats{ii},1),1)];
        end
        buffer_pix = 0;
    end
    cur_use_pix = (1:params.full_nPix) + buffer_pix; %use central pixels
    
    %cycle over trials
    n_trials = length(use_trials);
    cur_nrpt_frames = zeros(n_trials,1);
    cur_rpt_frames = cell(n_trials,1);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4; %stimulus frame onset times
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        if isfield(Expts{cur_block}.Trials(use_trials(tt)),'rptframes') %store any frames that got presented more than once
            cur_nrpt_frames(tt) = length(Expts{cur_block}.Trials(use_trials(tt)).rptframes);
            cur_rpt_frames{tt} = Expts{cur_block}.Trials(use_trials(tt)).rptframes;
        end
        if n_frames > 0
            if length(cur_stim_times) == 1 %if only the first stim time is stored
                cur_stim_times = (cur_stim_times:params.dt*params.Fr:(cur_stim_times + (n_frames-1)*params.dt*params.Fr))'; %create time axis
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
                cur_t_edges = [cur_stim_times; cur_stim_times(end) + params.dt*params.Fr];
            end
        end
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        %         cur_t_edges_up = cur_t_edges(1):up_dt:cur_t_edges(end);
        %         cur_t_axis_up = 0.5*cur_t_edges_up(1:end-1) + 0.5*cur_t_edges_up(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt)); %time since trial onset
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > params.min_trial_dur/params.dt %if using this trial
            use_frames = min(length(cur_stim_times),n_frames); %number of used frames
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix)); %stimulus matrix
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times+cur_toffset < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            
            %cat values
            all_stim_times = [all_stim_times; cur_stim_times + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_stim_mat = [all_stim_mat; cur_stim_mat];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
    all_trial_nrptframes = [all_trial_nrptframes; cur_nrpt_frames];
    all_trial_rptframes = cat(1,all_trial_rptframes,cur_rpt_frames);
    
    %need to keep track of block time offsets for LP recordings
    if strcmp(rec_type,'LP')
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
        cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
    end
end

if params.is_TBT_expt
    if any(~isnan(all_trial_exvals(:,3)))
        fprintf('Using exvals to define trial-by-trial conditions\n');
        all_trial_Ff(all_trial_exvals(:,3) == 1) = 0; %these are sim sac trials
        all_trial_Ff(all_trial_exvals(:,3) > 1) = 70;
        all_trial_back(all_trial_exvals(:,3) == 2) = 0; %these are gray-back trials
        all_trial_back(ismember(all_trial_exvals(:,3),[1 3])) = 1;
    end
end

if nanmax(all_trial_nrptframes) > 0
    warning('Some rpt frames detected');
end
%% BIN SPIKES FOR MU AND SU
clust_params.n_probes = params.n_probes;
if strcmp(rec_type,'LP')
    clust_params.exclude_adjacent = true;
else
    clust_params.exclude_adjacent = false;
end
[all_binned_mua,all_binned_sua,Clust_data,all_su_spk_times,~,all_mu_spk_times] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params,rec_block_range);
SU_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;

%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= params.beg_buffer & (params.trial_dur-all_tsince_start) >= params.end_buffer);
%for G093 use only data where stripe width is AT LEAST 2 deg
if strcmp(Expt_name,'G093')
    un_wi_vals = unique(all_trial_wi);
    use_wi_trials = find(all_trial_wi >= un_wi_vals(2));
    used_inds = used_inds(ismember(all_trialvec(used_inds),use_wi_trials));
end

%% PROCESS EYE TRACKING DATA
cd(data_dir)

if isfield(Expts{cur_block_set(1)}.Header,'exptno')
    em_block_nums = cellfun(@(X) X.Header.exptno,Expts(cur_block_set),'uniformoutput',1); %block numbering for EM/LFP data sometimes isnt aligned with Expts struct
else
    em_block_nums = cur_block_set;
end

%get raw ET data for these blocks
% [all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v3(all_t_axis,all_blockvec,em_block_nums,trial_toffset,params.good_coils);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori(cur_block_set),used_inds);

[saccades,et_params] = detect_saccades_v2(corrected_eye_vals,all_eye_vals,all_eye_speed,all_eye_ts,et_params); %detect saccades

is_blink = detect_blinks(all_eye_ts,all_eye_vals,saccades,et_params); %detect blinks

[saccades,is_blink] = merge_blinks(saccades,is_blink); %merge nearby blinks

%find saccade start/stop times aligned to expt time axis
sac_start_times = [saccades(:).start_time];
sac_stop_times = [saccades(:).stop_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));
interp_sac_stop_inds(isnan(interp_sac_stop_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > params.dt);
sac_start_times(bad_sacs) = [];
sac_stop_times(bad_sacs) = [];
saccades(bad_sacs) = [];
is_blink(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];
interp_sac_stop_inds(bad_sacs) = [];

%package ET data
ET_data.saccades = saccades;
ET_data.et_params = et_params;
ET_data.is_blink = is_blink;
ET_data.interp_eye_pos = corrected_eye_vals_interp;

%%

%in these two expts there was a problem with the repeat trials where the
%sequence would be shifted. We can still align these repeat trials by
%treating them as if there was some number of repeats of a 'zero-frame'
%before the trial started
if strcmp(Expt_name,'M012') || strcmp(Expt_name,'M013')
    rpt_trials = find(all_trial_Se < 1400);
    even_seeds = find(mod(all_trial_Se(rpt_trials),2) == 0);
    odd_seeds = find(mod(all_trial_Se(rpt_trials),2) == 1);
    frame_offsets = nan(length(rpt_trials),1);
    frame_offsets(even_seeds) = (all_trial_Se(rpt_trials(even_seeds)) - 1002)/2;
    frame_offsets(odd_seeds) = (all_trial_Se(rpt_trials(odd_seeds)) - 1001)/2;
    
    new_seeds = nan(length(rpt_trials),1);
    new_seeds(even_seeds) = 1002;
    new_seeds(odd_seeds) = 1001;
    all_trial_Se(rpt_trials) = new_seeds;
    
    for ii = 1:length(rpt_trials)
       all_trial_rptframes{rpt_trials(ii)} = cat(2,zeros(1,frame_offsets(ii)),all_trial_rptframes{rpt_trials(ii)});
       all_trial_nrptframes(rpt_trials(ii)) = length(all_trial_rptframes{rpt_trials(ii)});
    end
end


%% package trial-by-trial data
unique_seeds = unique(all_trial_Se);
seed_counts = hist(all_trial_Se,unique_seeds);
if ~params.is_TBT_expt
    [all_trial_Ff,all_trial_exvals,all_trial_back] = deal(nan(size(all_trial_start_times)));
end
trial_data = create_struct_from_mats('start_times',all_trial_start_times,'end_times',all_trial_end_times,'se',all_trial_Se,...
    'Ff',all_trial_Ff,'block_nums',all_trial_blocknums,'nrpt_frames',all_trial_nrptframes,'wi',all_trial_wi,'back_type',all_trial_back);
for ii = 1:length(trial_data); trial_data(ii).rpt_frames = all_trial_rptframes{ii}; end;

params.rpt_seeds = unique_seeds(seed_counts > 20);

%% package time data
time_data.t_axis = all_t_axis;
time_data.tsince_start = all_tsince_start;

time_data.trial_flip_inds = [1; find(all_trialvec(2:end) > all_trialvec(1:end-1))];
time_data.trial_flip_ids = all_trialvec(time_data.trial_flip_inds+1);
time_data.block_flip_inds = [1; find(all_blockvec(2:end) > all_blockvec(1:end-1))];
time_data.block_flip_ids = all_blockvec(time_data.block_flip_inds+1);

%% package stimulus
stimComp = compressTernNoise(all_stim_mat);

%% package binned spike data
spike_data.binned_mua = spikes_double2int8(all_binned_mua);
spike_data.binned_sua = spikes_double2int8(all_binned_sua);
spike_data.Clust_data = Clust_data;
spike_data.clust_params = clust_params;
spike_data.SU_spk_times = all_su_spk_times;

%% save packaged data
data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
if rec_number > 1
    data_name = strcat(data_name,sprintf('_r%d',rec_number));
end

save(data_name,'params','trial_data','expt_data','spike_data','stimComp','ET_data','time_data','used_inds');


%% quick check on SUA correlations to check for any missed double-units
n_sus = size(all_binned_sua,2);
norm_binned_spikes = nanzscore(all_binned_sua);
sua_corrmat = squeeze(nanmean(bsxfun(@times,reshape(norm_binned_spikes,[],1,n_sus),norm_binned_spikes)));

sua_corrmat(logical(eye(n_sus))) = nan;
f1 = figure();
imagescnan(sua_corrmat);

%% coarse xcorr test
% maxlags = 50;
% check_units = [1 3];
% uinds = find(~isnan(all_binned_sua(:,check_units(1))) & ~isnan(all_binned_sua(:,check_units(2))));
% [xc,xl] = xcov(all_binned_sua(uinds,check_units(1)),all_binned_sua(uinds,check_units(2)),maxlags,'coeff');
% f2 = figure();
% plot(xl*params.dt,xc);
    
%% fine timesclae xcorr check
% check_units = [1 3];
% rel_bins = (-0.05:0.001:0.05)';
% 
% joint_counts = zeros(size(rel_bins));
% spike_counts = nansum(all_binned_sua);
% [~,smaller_unit] = min(spike_counts(check_units));
% 
% if smaller_unit == 1
% ref_times = all_su_spk_times{check_units(1)};
% test_times = all_su_spk_times{check_units(2)};
% else
% ref_times = all_su_spk_times{check_units(2)};
% test_times = all_su_spk_times{check_units(1)};    
% end
% 
% nspks = length(ref_times);
% for ii = 1:nspks
%     joint_counts = joint_counts + histc(test_times - ref_times(ii),rel_bins);
% end
% 
% f1 = figure();
% plot(rel_bins(1:end-1),joint_counts(1:end-1));