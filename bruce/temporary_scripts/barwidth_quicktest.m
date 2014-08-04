clear all
% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_name = 'G086';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/G081/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/par_orth_sacs'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

recompute_prepdata = 1;
prep_data_name = [anal_dir '/sac_prepdata.mat'];
recompute_sacdata = 1;
sac_data_name = [anal_dir '/sac_sacdata.mat'];

%% PARSE TRIAL DATA STRUCTURES

stim_fs = 100; %in Hz
dt = 0.005;
beg_buffer = round(0.15/dt);
end_buffer = round(0.15/dt);

%% PARSE TRIAL DATA STRUCTURES
target_expt_name = {'rls.Fa' 'rls.seXFaRC'};
is_target_type = false(length(Expts),1);
expt_frame_dur = nan(length(Expts),1);
expt_sac_dir = nan(length(Expts),1);
expt_bar_ori = nan(length(Expts),1);
for i = 1:length(Expts)
    expt_names{i} = Expts{i}.Header.expname;
    if any(strcmp(expt_names{i},target_expt_name))
        is_target_type(i) = true;
    end
    expt_frame_dur(i) = Expts{i}.Stimvals.Fr;
    expt_bar_ori(i) = Expts{i}.Stimvals.or;
    expt_sac_dir(i) = mod(Expts{i}.Stimvals.Fa,180);
end

cur_expt_set = find(is_target_type & expt_bar_ori == 90 & expt_frame_dur == 1);
if strcmp(Expt_name,'G088')
    cur_expt_set(cur_expt_set == 22) = []; %only 4 trials and causes problems
end

%% COMPUTE TRIAL DATA
cd(data_dir);

 %%
 fprintf('Computing prep data\n');
    trial_cnt = 0;
    
    all_stim_times = [];
    all_t_axis = [];
    all_tsince_start = [];
    all_used_inds = false(0);
    all_phase_Xmat = [];
    all_binned_spks = [];
    all_exptvec = [];
    all_bar_or = [];
    all_frame_dur = [];
    all_sac_dir = [];
    all_trialvec = [];
    all_trial_start_times = [];
    all_trial_end_times = [];
    for ee = 1:length(cur_expt_set);
        fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
        cur_expt = cur_expt_set(ee);
        fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
        load(fname);
        
        trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
        trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
        trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
        trial_ids = [Expts{cur_expt}.Trials(:).id];
        
        [un_ids,id_inds] = unique(trial_ids);
        trial_durs = trial_durs(id_inds);
        trial_start_times = trial_start_times(id_inds);
        trial_end_times = trial_end_times(id_inds);
        
        use_trials = find(trial_durs >= 0.5);
        
        if isfield(Expts{cur_expt}.Trials(1),'Fr')
            trial_Fr = [Expts{cur_expt}.Trials(:).Fr];
            all_frame_dur = [all_frame_dur; trial_Fr(use_trials)'];
        elseif isfield(Expts{cur_expt}.Trials(1),'sl')
            trial_sl = [Expts{cur_expt}.Trials(:).sl];
            trial_sl(trial_sl == 0) = 1;
            all_frame_dur = [all_frame_dur; trial_sl(use_trials)'];            
        else
            all_frame_dur = [all_frame_dur; ones(length(use_trials),1)*expt_frame_dur(cur_expt)];
        end
        if isfield(Expts{cur_expt}.Trials(1),'or')
            trial_or = [Expts{cur_expt}.Trials(:).or];
            all_bar_or = [all_bar_or; trial_or(use_trials)'];
        else
            all_bar_or = [all_bar_or; ones(length(use_trials),1)*expt_bar_ori(cur_expt)];            
        end
        if isfield(Expts{cur_expt}.Trials(1),'Fa')
            trial_Fa = [Expts{cur_expt}.Trials(:).Fa];
            all_sac_dir = [all_sac_dir; trial_Fa(use_trials)'];
        else
            all_sac_dir = [all_sac_dir; ones(length(use_trials),1)*expt_sac_dir(cur_expt)];                        
        end
        
        all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
        all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
        n_trials = length(use_trials);
        for tt = 1:n_trials
            cur_stim_times = Expts{cur_expt}.Trials(use_trials(tt)).Start/1e4;
            if length(cur_stim_times) == 1
                cur_stim_times = trial_start_times(use_trials(tt)):1/stim_fs:trial_end_times(use_trials(tt));
            end
            cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
            cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
            cur_binned_spks = nan(length(cur_t_axis),96);
            for cc = 1:96
                cur_hist = histc(Clusters{cc}.times,cur_t_edges);
                cur_binned_spks(:,cc) = cur_hist(1:end-1);
                %             temp = convert_to_spikebins(cur_hist(1:end-1));
            end
            
            cur_used_inds = true(length(cur_t_axis),1);
            cur_used_inds(1:beg_buffer) = false;
            cur_used_inds((end-end_buffer+1):end) = false;
            
            cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
            
            all_stim_times = [all_stim_times; cur_stim_times(:)];
            all_t_axis = [all_t_axis; cur_t_axis];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_used_inds = [all_used_inds; cur_used_inds];
            all_binned_spks = [all_binned_spks; cur_binned_spks];
            all_exptvec = [all_exptvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
        end
        trial_cnt = trial_cnt + n_trials;
    end
    save(prep_data_name,'all_*');

%%
exptset1 = find(strcmp(expt_names(cur_expt_set),'rls.Fa'));
exptset2 = setdiff(1:length(cur_expt_set),exptset1);

for ee = 1:length(cur_expt_set)
    cur_expt_inds = find(all_exptvec == ee);
    expt_mean_rates(ee,:) = mean(all_binned_spks(cur_expt_inds,:))/dt;
end
ov_avg_rates = mean(expt_mean_rates);

