clear all
close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat
 
load ~/Data/bruce/grating_params_fin.mat
% grating_seq_id = [];
% grating_seq_ori = [];
% grating_seq_freq = [];
% grating_seq_phase = [];
for i = 1:3
%    grating_seq_id = [grating_seq_id ones(1,length(grating_params{i}))*i]; 
%    grating_seq_ori = [grating_seq_ori grating_params{i}.ori];
%    grating_seq_freq = [grating_seq_freq grating_params{i}.freq];
%    grating_seq_phase = [grating_seq_phase grating_params{i}.phase];
grating_seq_ori{i} = [grating_params{i}.ori];
grating_seq_freq{i} = [grating_params{i}.freq];
grating_seq_phase{i} = [grating_params{i}.phase];
end
%%
stim_fs = 1e4/117.5;
grating_blocks = [12 34 52] - 6;
all_stim_expt_ids = [];
all_stim_ids = [];
all_stim_start_times = [];
all_binned_spks = [];
for bb = 1:length(grating_blocks)
    fprintf('Analyzing block %d of %d\n',bb,length(grating_blocks));
    load(sprintf('Expt%dClusterTimes.mat',grating_blocks(bb)));
    
    n_trials(bb) = length(Expts{grating_blocks(bb)}.Trials);
    trial_completed{bb} = [Expts{grating_blocks(bb)}.Trials(:).Result];
    
%     cur_used = find(trial_completed{bb}==1);
    cur_used = 1:length(trial_completed{bb});
    cur_stim_ids = [];
    cur_stim_start_times = [];
    for i = 1:length(cur_used)
       cur_n_stims(i) = ceil(length(Expts{grating_blocks(bb)}.Trials(cur_used(i)).Seedseq)); 
       cur_trial_start_time = Expts{grating_blocks(bb)}.Trials(cur_used(i)).Start(1)/1e4;
       cur_seedseq = Expts{grating_blocks(bb)}.Trials(cur_used(i)).Seedseq(1:2:end);
       clen{bb}(i) = length(cur_seedseq);
%        cur_stim_ids = [cur_stim_ids; cur_seedseq(1:160)];
       cur_stim_ids = [cur_stim_ids; cur_seedseq];
       cur_t_edges = cur_trial_start_time + (0:2:2*length(cur_seedseq))/stim_fs;
%        if length(cur_seedseq) ~= length(cur_t_edges)-1
%            disp('unequal')
%            pause
%        end
       cur_stim_start_times = [cur_stim_start_times; cur_t_edges(1:end-1)'];
       
       binned_spks = nan(96,length(cur_t_edges)-1);
       for cc = 1:96
           temp = histc(Clusters{(cc)}.times,cur_t_edges);
           binned_spks(cc,:) = temp(1:end-1);
       end
       all_binned_spks = cat(1,all_binned_spks,binned_spks');
       
    end
    
    all_stim_ids = [all_stim_ids; cur_stim_ids];
    all_stim_start_times = [all_stim_start_times; cur_stim_start_times];
    all_stim_expt_ids = [all_stim_expt_ids; ones(length(cur_stim_ids),1)*grating_blocks(bb)];
end

%%
all_stim_ori = nan(size(all_stim_ids));
all_stim_freq = nan(size(all_stim_ids));
all_stim_phase = nan(size(all_stim_ids));

stim_offset = 4000;
cur_set = find(all_stim_ids >= stim_offset & all_stim_ids < stim_offset+1000); 
cur_stim_ids = all_stim_ids(cur_set)-stim_offset;
cur_set(cur_stim_ids > 160) = [];
cur_stim_ids(cur_stim_ids > 160) = [];
all_stim_ori(cur_set) = grating_seq_ori{1}(cur_stim_ids);
all_stim_freq(cur_set) = grating_seq_freq{1}(cur_stim_ids);
all_stim_phase(cur_set) = grating_seq_phase{1}(cur_stim_ids);

stim_offset = 5000;
cur_set = find(all_stim_ids >= stim_offset & all_stim_ids < stim_offset+1000); 
cur_stim_ids = all_stim_ids(cur_set)-stim_offset;
cur_set(cur_stim_ids > 160) = [];
cur_stim_ids(cur_stim_ids > 160) = [];
all_stim_ori(cur_set) = grating_seq_ori{2}(cur_stim_ids);
all_stim_freq(cur_set) = grating_seq_freq{2}(cur_stim_ids);
all_stim_phase(cur_set) = grating_seq_phase{2}(cur_stim_ids);

stim_offset = 6000;
cur_set = find(all_stim_ids >= stim_offset & all_stim_ids < stim_offset+1000); 
cur_stim_ids = all_stim_ids(cur_set)-stim_offset;
cur_set(cur_stim_ids > 160) = [];
cur_stim_ids(cur_stim_ids > 160) = [];
all_stim_ori(cur_set) = grating_seq_ori{3}(cur_stim_ids);
all_stim_freq(cur_set) = grating_seq_freq{3}(cur_stim_ids);
all_stim_phase(cur_set) = grating_seq_phase{3}(cur_stim_ids);

% all_stim_ori = mod(all_stim_ori+90,180);
all_stim_ori = mod(270-all_stim_ori,180);
%%
freq_bin_edges = linspace(48,195,8);
ori_bin_edges = linspace(0,180,12);
phase_bin_edges = linspace(0,360,2);
freq_bin_centers = 0.5*freq_bin_edges(1:end-1)+0.5*freq_bin_edges(2:end);
ori_bin_centers = 0.5*ori_bin_edges(1:end-1)+0.5*ori_bin_edges(2:end);
phase_bin_centers = 0.5*phase_bin_edges(1:end-1)+0.5*phase_bin_edges(2:end);

tot_n_stims = length(freq_bin_centers)*length(ori_bin_centers)*length(phase_bin_centers);
all_stim_code = nan(size(all_stim_ids));
cnt = 1;
code_freq = zeros(tot_n_stims,1);
code_ori = zeros(tot_n_stims,1);
code_phase = zeros(tot_n_stims,1);
code_cnt = zeros(tot_n_stims,1);
for i = 1:length(freq_bin_centers)
    for j = 1:length(ori_bin_centers)
        for k = 1:length(phase_bin_centers)
            cur_set = find(all_stim_ori >= ori_bin_edges(j) & all_stim_ori < ori_bin_edges(j+1) ...
                & all_stim_freq >= freq_bin_edges(i) & all_stim_freq < freq_bin_edges(i+1) ...
                & all_stim_phase >= phase_bin_edges(k) & all_stim_phase < phase_bin_edges(k+1));
           all_stim_code(cur_set) = cnt;
           code_freq(cnt) = freq_bin_centers(i);
           code_ori(cnt) = ori_bin_centers(j);
           code_phase(cnt) = phase_bin_centers(k);
           code_cnt(cnt) = length(cur_set);
           cnt = cnt + 1;
        end
    end
end

%%
dt = 2/stim_fs;
backlag = round(0);
forwardlag = round(0.15/dt);
clear stim_trg_avg n_used
stim_trg_avg = nan(tot_n_stims,96,backlag+forwardlag+1);

for i = 1:tot_n_stims
    cur_set = find(all_stim_code==i);
    cur_set(cur_set < backlag | cur_set > length(all_stim_ids)-forwardlag) = [];
    if length(cur_set) > 5
        for ll = 1:96
            [stim_trg_avg(i,ll,:),tlags] = get_event_trig_avg(all_binned_spks(:,ll),cur_set,backlag,forwardlag);
        end
    end
    n_used(i) = length(cur_set);
end

%%
ori_trg_avgs = nan(length(ori_bin_centers),96,backlag+forwardlag+1);
for i = 1:length(ori_bin_centers)
   cur_code = find(code_ori==ori_bin_centers(i));
   ori_trg_avgs(i,:,:) = nanmean(stim_trg_avg(cur_code,:,:),1);
end

%%
freq_trg_avgs = nan(length(freq_bin_centers),96,backlag+forwardlag+1);
for i = 1:length(freq_bin_centers)
   cur_code = find(code_freq==freq_bin_centers(i));
   freq_trg_avgs(i,:,:) = nanmean(stim_trg_avg(cur_code,:,:),1);
end

%%
use_lag = 3;
ori_tuning = squeeze(ori_trg_avgs(:,:,use_lag));
freq_tuning = squeeze(freq_trg_avgs(:,:,use_lag));
avg_rates = mean(all_binned_spks);
ori_tuning_n = bsxfun(@rdivide,ori_tuning,avg_rates);
freq_tuning_n = bsxfun(@rdivide,freq_tuning,avg_rates);

[~,pref_ori] = max(ori_tuning);
pref_ori = ori_bin_centers(pref_ori);
[~,pref_freq] = max(freq_tuning);
pref_freq = freq_bin_centers(pref_freq);

cd /home/james/Data/bruce/G075
save grating_unit_tuning *tuning* *centers avg_rates use_lag stim_fs dt pref_ori pref_freq 

%%
load ./ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

pref_ori_mat = nan(10,10);
pref_freq_mat = nan(10,10);
for i = 1:10
    for j = 1:10
        cur = find(X_pos==i & Y_pos==j);
        if ~isempty(cur)
            pref_ori_mat(i,j) = pref_ori(cur);
            pref_freq_mat(i,j) = pref_freq(cur);
        end
    end
end
