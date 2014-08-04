clear all
close all

stim_dt = 0.01;
usfac = 1;
tbspace = 1;
dt = stim_dt/usfac;

base_flen = 15;
use_nPix = 20;

min_trial_dur = 2;
trial_dur = 4;

n_probes = 24;

%%
monName = 'lem';
exp_name = 'M266';
block_nums = [3];
Expts = {};
for bb = 1:length(block_nums)
    dat_name = [pwd sprintf('/%s%s.%d.mat',monName,exp_name,block_nums(bb))];
    [a,cur_Expt] = APlaySpkFile(dat_name,'nospikes','noerrs');
    Expts = {Expts{:} cur_Expt{:}};
end
load('./stim_data.mat');
load('./expt_data.mat');

%%

fprintf('Computing prep data\n');
trial_cnt = 0;
cur_toffset = 0;

all_stim_times = [];
all_stim_mat = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);
trial_toffset = zeros(length(block_nums),1);
for ee = 1:length(Expts);
    fprintf('Block %d of %d\n',ee,length(Expts));
    cur_block_num = block_nums(ee);
    
    Fr = Expts{ee}.Stimvals.Fr;
        
    trial_start_times = [Expts{ee}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{ee}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{ee}.Trials(:).dur]/1e4;
    trial_ids = [Expts{ee}.Trials(:).id];
    
    [un_ids,id_inds] = unique(trial_ids);
    rpt_trials = false;
    if length(un_ids) < length(trial_ids)
        rpt_trials = true;
        fprintf('Warning, repeat trial inds detected!\n');
        use_trials = [];
    else
        use_trials = find(trial_durs >= min_trial_dur);
    end
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    
    trial_Se = [Expts{ee}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    
    fname = sprintf('Expt%d_stim.mat',cur_block_num);
    load(fname);
    buffer_pix = floor((expt_npix(ee) - use_nPix)/2);
    cur_use_pix = (1:use_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{ee}.Trials(use_trials(tt)).Start'/1e4;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        if length(cur_stim_times) == 1
            cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
            cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
        end
%         cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
%         if length(cur_t_edges) > trial_dur/dt + 1
%             cur_t_edges(round(trial_dur/dt+2):end) = [];
%         end
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            
            all_stim_mat = [all_stim_mat; cur_stim_mat];
            all_stim_times = [all_stim_times; cur_stim_times' + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
            
        end
    end
    trial_cnt = trial_cnt + n_trials;
    trial_toffset(ee) = all_t_bin_edges(end);
    cur_toffset = trial_toffset(ee);
end

%%
Fs = 1000;
dsf = 3;
Fsd = Fs/dsf;
niqf = Fs/2;
[bb,aa] = butter(4,[1 120]/niqf);

full_lfps = [];
full_lfp_taxis = [];
cur_toffset = 0;
for ee = 1:length(Expts);
    fprintf('Loading LFPs, Expt %d\n',block_nums(ee));
    fname = sprintf('lem%sA.%d.lfp.mat',exp_name,block_nums(ee));
    load(fname);
    
    Fs = 1/LFP.Header.CRsamplerate;
    
    n_trials(ee) = length(LFP.Trials);
    lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
    lfp_trial_ends = [LFP.Trials(:).End]/1e4;
    expt_lfp_t_axis = [];
    expt_lfps = [];
    for tt = 1:n_trials(ee)
        %         tt
        cur_npts = size(LFP.Trials(tt).LFP,1);
        cur_t_end(tt) = lfp_trial_starts(tt)+(cur_npts-1)/Fs;
        cur_t_axis = (lfp_trial_starts(tt):1/Fs:cur_t_end(tt)) + cur_toffset;
        
        if ~isempty(expt_lfp_t_axis)
            cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
        else
            cur_sp = 1;
        end
        cur_t_axis = cur_t_axis(cur_sp:end);
        cur_LFP = [LFP.Trials(tt).LFP];
        cur_LFP = cur_LFP(cur_sp:end,:);
        cur_LFP = filtfilt(bb,aa,cur_LFP);
        
        cur_LFP = downsample(cur_LFP,dsf);
        cur_t_axis = downsample(cur_t_axis,dsf);
        
        expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
        expt_lfps = [expt_lfps; cur_LFP];
    end
    
    cur_set = find(all_blockvec==ee);
    full_lfp_taxis = [full_lfp_taxis; expt_lfp_t_axis];
    full_lfps = [full_lfps; expt_lfps];
    cur_toffset = trial_toffset(ee);
end
full_lfps = full_lfps/std(full_lfps(:));

interp_lfps = interp1(full_lfp_taxis,full_lfps,all_t_axis);

%% CSD
%compute CSD
vars.Fs = Fsd;
vars.BrainBound = 1;
vars.ChanSep = 0.05;
vars.diam = 2; %0.5
full_CSDs = PettersenCSD(full_lfps','spline',vars)';
full_CSDs = full_CSDs/std(full_CSDs(:));

%%
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'));
params.Fs = Fsd;
params.tapers = [5 9];
params.fpass = [0.5 200];
params.segave = 1;
win = [3.5 3.5];

trial_starts = [1; 1+find(diff(all_trialvec) > 0)];
trial_ends = [find(diff(all_trialvec) > 0); length(all_trialvec)];
% lfp_trial_starts = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_t_axis(trial_starts)));
% lfp_trial_ends = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_t_axis(trial_ends)));
lfp_trial_starts = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_trial_start_times));
lfp_trial_ends = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_trial_end_times));

sMarkers = [lfp_trial_starts(:) lfp_trial_ends(:)];

clear P
for ii = 1:n_probes
[P(ii,:),f] = mtspectrumc_unequal_length_trials(full_lfps(:,ii),win,params,sMarkers);
end
lP = log10(P);
lPz = zscore(lP);

figure
subplot(2,1,1)
pcolor(f,1:24,lP);shading interp;
set(gca,'ydir','reverse')
% set(gca,'xscale','log');
xlabel('Frequency (Hz)','fontsize',12);
ylabel('Probe number','fontsize',12);
colorbar
title('Absolute power');
subplot(2,1,2)
pcolor(f,1:24,lPz);shading interp;
set(gca,'ydir','reverse','xscale','log');
caxis([-2 2]);
colorbar
xlabel('Frequency (Hz)','fontsize',12);
ylabel('Probe number','fontsize',12);
title('Relative power');
%%
forwardlag = round(Fsd*0.3);
backlag = round(0.1*Fsd);
[tstart_trig_lfps,tlags] = get_event_trig_avg(full_lfps,lfp_trial_starts,backlag,forwardlag);
[tstart_trig_csds,tlags] = get_event_trig_avg(full_CSDs,lfp_trial_starts,backlag,forwardlag);

figure
subplot(2,1,1)
pcolor(tlags/Fsd,1:24,tstart_trig_lfps'); shading flat
set(gca,'ydir','reverse');
ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
yl = ylim(); line([0 0],yl,'color','k');
xlabel('Time since trial onset (s)','fontsize',12);
ylabel('Probe number','fontsize',12);
title('Trial-trig avg LFP');
subplot(2,1,2)
pcolor(tlags/Fsd,1:24,tstart_trig_csds');shading flat
set(gca,'ydir','reverse');
ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
yl = ylim(); line([0 0],yl,'color','k');
xlabel('Time since trial onset (s)','fontsize',12);
ylabel('Probe number','fontsize',12);
title('Trial-trig avg CSD');


%%
forwardlag = round(Fsd*0.25);
backlag = round(0.05*Fsd);
n_lags = forwardlag + backlag + 1;
n_poss_bar_pos = size(all_stim_mat,2);
stim_trig_lfps = nan(n_lags,n_probes,n_poss_bar_pos);
stim_trig_csds = nan(n_lags,n_probes,n_poss_bar_pos);
for ii = 1:n_poss_bar_pos
    cur_stim_times = all_t_axis(all_stim_mat(:,ii) ~= 0);
    cur_stim_inds = round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),cur_stim_times));
    
    [stim_trig_lfps(:,:,ii),lags] = get_event_trig_avg(full_lfps,cur_stim_inds,backlag,forwardlag);
    [stim_trig_csds(:,:,ii),lags] = get_event_trig_avg(full_CSDs,cur_stim_inds,backlag,forwardlag);
end
all_stim_inds =  round(interp1(full_lfp_taxis,1:length(full_lfp_taxis),all_t_axis));
[allstim_trig_lfps,lags] = get_event_trig_avg(full_lfps,all_stim_inds,backlag,forwardlag);
[allstim_trig_csds,lags] = get_event_trig_avg(full_CSDs,all_stim_inds,backlag,forwardlag);

stim_trig_lfps = bsxfun(@minus,stim_trig_lfps,allstim_trig_lfps);
stim_trig_csds = bsxfun(@minus,stim_trig_csds,allstim_trig_csds);

lfp_bar_tuning = max(squeeze(std(stim_trig_lfps)));
csd_bar_tuning = max(squeeze(std(stim_trig_csds)));

figure
plot(1:n_poss_bar_pos,lfp_bar_tuning/max(lfp_bar_tuning),1:n_poss_bar_pos,csd_bar_tuning/max(csd_bar_tuning),'r');
xlabel('Bar position');
ylabel('Amp (au)');
legend('LFP','CSD');
%%
bar_pos = 10;
figure
subplot(2,1,1)
pcolor(lags/Fsd,1:n_probes,squeeze(stim_trig_lfps(:,:,bar_pos))');shading flat
set(gca,'ydir','reverse');
ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
yl = ylim(); line([0 0],yl,'color','k');
title('LFPs');
xlabel('Time since stim onset (s)','fontsize',12);
ylabel('Probe number','fontsize',12);
subplot(2,1,2)
pcolor(lags/Fsd,1:n_probes,squeeze(stim_trig_csds(:,:,bar_pos))');shading flat
set(gca,'ydir','reverse');
ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
yl = ylim(); line([0 0],yl,'color','k');
title('CSDs');
xlabel('Time since stim onset (s)','fontsize',12);
ylabel('Probe number','fontsize',12);

%%
cur_ch = 14;
figure
subplot(2,1,1)
pcolor(lags/Fsd,1:n_poss_bar_pos,squeeze(stim_trig_lfps(:,cur_ch,:))');shading flat
set(gca,'ydir','reverse');
ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
yl = ylim(); line([0 0],yl,'color','k');
title('LFPs');
xlabel('Time since stim onset (s)','fontsize',12);
ylabel('Bar position','fontsize',12);

subplot(2,1,2)
pcolor(lags/Fsd,1:n_poss_bar_pos,squeeze(stim_trig_csds(:,cur_ch,:))');shading flat
set(gca,'ydir','reverse');
ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
yl = ylim(); line([0 0],yl,'color','k');
title('CSDs');
xlabel('Time since stim onset (s)','fontsize',12);
ylabel('Bar position','fontsize',12);
