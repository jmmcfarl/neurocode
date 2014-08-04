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
for bb = 1:length(grating_blocks)
    n_trials(bb) = length(Expts{grating_blocks(bb)}.Trials);
    trial_completed{bb} = [Expts{grating_blocks(bb)}.Trials(:).Result];
    
    cur_used = find(trial_completed{bb}==1);
    cur_stim_ids = [];
    cur_stim_start_times = [];
    for i = 1:length(cur_used)
       cur_n_stims(i) = ceil(length(Expts{grating_blocks(bb)}.Trials(cur_used(i)).Seedseq)); 
       cur_trial_start_time = Expts{grating_blocks(bb)}.Trials(cur_used(i)).Start(1)/1e4;
       cur_seedseq = Expts{grating_blocks(bb)}.Trials(cur_used(i)).Seedseq(1:2:end);
       clen{bb}(i) = length(cur_seedseq);
       cur_stim_ids = [cur_stim_ids; cur_seedseq(1:160)];
       cur_stim_start_times = [cur_stim_start_times; cur_trial_start_time+(0:2:318)'/stim_fs];
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
all_stim_ori(cur_set) = grating_seq_ori{1}(cur_stim_ids);
all_stim_freq(cur_set) = grating_seq_freq{1}(cur_stim_ids);
all_stim_phase(cur_set) = grating_seq_phase{1}(cur_stim_ids);

stim_offset = 5000;
cur_set = find(all_stim_ids >= stim_offset & all_stim_ids < stim_offset+1000); 
cur_stim_ids = all_stim_ids(cur_set)-stim_offset;
all_stim_ori(cur_set) = grating_seq_ori{2}(cur_stim_ids);
all_stim_freq(cur_set) = grating_seq_freq{2}(cur_stim_ids);
all_stim_phase(cur_set) = grating_seq_phase{2}(cur_stim_ids);

stim_offset = 6000;
cur_set = find(all_stim_ids >= stim_offset & all_stim_ids < stim_offset+1000); 
cur_stim_ids = all_stim_ids(cur_set)-stim_offset;
all_stim_ori(cur_set) = grating_seq_ori{3}(cur_stim_ids);
all_stim_freq(cur_set) = grating_seq_freq{3}(cur_stim_ids);
all_stim_phase(cur_set) = grating_seq_phase{3}(cur_stim_ids);

freq_bin_edges = linspace(48,195,6);
ori_bin_edges = linspace(0,360,6);
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
Fs = 3e4;
dsf = 60;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1]/niqf,'high');
use_lfps = [1:16:96];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);

all_Vmat = [];
all_Vmatf = [];
all_t_ax = [];
all_stim_start_inds = [];
all_expt_inds = [];
for bb = 1:length(grating_blocks)
    fprintf('Analyzing block %d of %d\n',bb,length(grating_blocks));
    % bb = 1;
    
    filename = sprintf('Expt%dFullVmean.mat',grating_blocks(bb));
    load(filename);
    
    Vmat = [];
    Vmatf = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',grating_blocks(bb),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = V + sumv*FullV.sumscale;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dV = [];
        dVf = [];
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            dV = [dV decimate(V(cur_range),dsf)];
            dVf = [dVf filtfilt(filt_b,filt_a,decimate(V(cur_range),dsf))]; %do some high-pass filtering
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(ll,:) = dV;
        Vmatf(ll,:) = dVf;
   end
        
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,2)+1:end) = [];   
    
    cur_use_trials = find(all_stim_expt_ids==grating_blocks(bb));    
    use_stim_start_inds = round(interp1(t_ax,1:length(t_ax),all_stim_start_times(cur_use_trials)));
    all_stim_start_inds = [all_stim_start_inds; use_stim_start_inds + length(all_t_ax)];
    
    all_Vmat = [all_Vmat Vmat];
    all_Vmatf = [all_Vmatf Vmatf];
    all_t_ax = [all_t_ax t_ax];
    all_expt_inds = [all_expt_inds ones(1,length(t_ax))*bb];    
end

%%
backlag = round(Fsd*0.1);
forwardlag = round(Fsd*0.25);
clear stim_trg_avg
stim_trg_avg = nan(tot_n_stims,length(use_lfps),backlag+forwardlag+1);
stim_trg_sem = nan(tot_n_stims,length(use_lfps),backlag+forwardlag+1);

for i = 1:tot_n_stims
    cur_set = find(all_stim_code==i);
    cur_set(all_stim_start_inds(cur_set) > backlag & all_stim_start_inds(cur_set) < length(all_t_ax)-forwardlag);
    if length(cur_set) > 20
        for ll = 1:length(use_lfps)
            %         [stim_trg_avg(i,ll,:),tlags] = get_event_trig_avg(all_Vmatf(ll,:),all_stim_start_inds(cur_set),backlag,forwardlag);
            [stim_trg_mat,tlags] = get_event_trig_mat(all_Vmatf(ll,:),all_stim_start_inds(cur_set),backlag,forwardlag);
            stim_trg_avg(i,ll,:) = nanmean(stim_trg_mat);
            stim_trg_sem(i,ll,:) = nanstd(stim_trg_mat)/sqrt(length(cur_set));
        end
    end
    n_used(i) = length(cur_set);
end

%%
ov_avg = squeeze(mean(stim_trg_avg));
stim_trg_avgn = zeros(size(stim_trg_avg));
for i = 1:tot_n_stims
    stim_trg_avgn(i,:,:) = bsxfun(@minus,squeeze(stim_trg_avg(i,:,:)),ov_avg);    
end
%%
close all
use_elec = 1;
use_freq = 4;
% use_phases = phase_bin_centers([1 2 3 4]);
cmap = jet(length(ori_bin_centers));
for i = 1:length(ori_bin_centers)
%    cur_code = find(code_ori==ori_bin_centers(i) & code_freq==freq_bin_centers(use_freq) & code_phase == phase_bin_centers(use_phase)); 
   cur_code = find(code_ori==ori_bin_centers(i) & code_freq==freq_bin_centers(use_freq)); 
%    if ~isempty(cur_code)
%     plot(tlags/Fsd,squeeze(stim_trg_avgn(cur_code,use_elec,:)),'color',cmap(i,:));
    shadedErrorBar(tlags/Fsd,squeeze(stim_trg_avgn(cur_code,use_elec,:)),squeeze(stim_trg_sem(cur_code,use_elec,:)),{'color',cmap(i,:)});
%    end
    hold on
end
shg

%%
close all
use_elec = 1;
use_freq = 4;
figure
cmap = jet(length(freq_bin_centers));
for j = 1:4
    subplot(2,2,j)
for i = 1:length(freq_bin_centers)
   cur_code = find(code_freq==freq_bin_centers(i)); 
%    if ~isempty(cur_code)
%     plot(tlags/Fsd,squeeze(stim_trg_avg(cur_code,use_elec,:)),'color',cmap(i,:));
    shadedErrorBar(tlags/Fsd,squeeze(stim_trg_avgn(cur_code,j,:)),squeeze(stim_trg_sem(cur_code,j,:)),{'color',cmap(i,:)});
%    end
    hold on
end
xlim([-0.05 0.2])
xlabel('Time (s)')
ylabel('Amplitude (V)')
ylim([-5 5]*1e-5)
yl = ylim();
line([0 0],yl,'color','k')
title(sprintf('Channel %d',use_lfps(j)));
end
shg
%%
close all
use_elec = 1;
use_freq = 4;
cmap = jet(length(ori_bin_centers));
for j = 1:4
    subplot(2,2,j)
for i = 1:length(ori_bin_centers)
   cur_code = find(code_ori==ori_bin_centers(i)); 
%    if ~isempty(cur_code)
%     plot(tlags/Fsd,squeeze(stim_trg_avg(cur_code,use_elec,:)),'color',cmap(i,:));
    shadedErrorBar(tlags/Fsd,squeeze(stim_trg_avgn(cur_code,j,:)),squeeze(stim_trg_sem(cur_code,j,:)),{'color',cmap(i,:)});
%    end
    hold on
end
xlim([-0.05 0.15])
xlabel('Time (s)')
ylabel('Amplitude (V)')
ylim([-2 2]*1e-5)
yl = ylim();
line([0 0],yl,'color','k')
title(sprintf('Channel %d',use_lfps(j)));
end
shg