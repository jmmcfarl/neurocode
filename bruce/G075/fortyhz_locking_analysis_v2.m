clear all
% close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

%%
% loc_glob_blocks = [8 15 20 27 33 39 49] - 6;
forty_blocks = [7 13 17 19 23 26 29 31 36 38 48] - 6;
repeat_inds = [1:9]*1e3;
all_expt_id = [];
for bb = 1:length(forty_blocks)
    n_trials(bb) = length(Expts{forty_blocks(bb)}.Trials);
    trial_start_times{bb} = [Expts{forty_blocks(bb)}.Trials(:).Start]/1e4;
    trial_stop_times{bb} = [Expts{forty_blocks(bb)}.Trials(:).End]/1e4;
    trial_seof{bb} = [Expts{forty_blocks(bb)}.Trials(:).seof];
    trial_completed{bb} = [Expts{forty_blocks(bb)}.Trials(:).Result];
    trial_durs{bb} = trial_stop_times{bb} - trial_start_times{bb};
    all_expt_id = [all_expt_id forty_blocks(bb)*ones(1,n_trials(bb))];
end
all_trial_start_times = cell2mat(trial_start_times);
all_trial_stop_times = cell2mat(trial_stop_times);
all_trial_seof = cell2mat(trial_seof);
all_trial_completed = cell2mat(trial_completed);
all_trial_durs = cell2mat(trial_durs);

%%
Fs = 3e4;
% dsf = 120;Fsd = Fs/dsf;
dsf = 60;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1]/niqf,'high');
[filt_b2,filt_a2] = butter(2,[1 30]/niqf);
use_lfps = [1:8:96];
use_units = 1:96;
% use_lfps = [1 63];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);

spk_sm_win = round(Fsd*0.01);

scales = logspace(log10(8),log10(75),20);
scales = [scales 85 100 115 130 150 200 250];
scales = scales*60/dsf;
% scales = scales*2;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

for ss = 1:length(repeat_inds)
    all_block_stim_lfps{ss} = [];
    all_block_stim_lfpsf{ss} = [];
    all_block_stim_phase{ss} = [];
    all_block_stim_amp{ss} = [];
    all_block_stim_spks{ss} = [];
    all_block_stim_smspks{ss} = [];
end
for bb = 1:length(forty_blocks)
    fprintf('Analyzing block %d of %d\n',bb,length(forty_blocks));
    load(sprintf('Expt%dClusterTimes.mat',forty_blocks(bb)));
    
    filename = sprintf('Expt%dFullVmean.mat',forty_blocks(bb));
    load(filename);
    
    Vmat = [];
    Vmatf = [];
    phasegrams = [];
    ampgrams = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',forty_blocks(bb),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = V + FullV.sumscale*sumv;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dV = [];
        dVf = [];
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            dV = [dV filtfilt(filt_b,filt_a,decimate(V(cur_range),dsf))]; %do some high-pass filtering
            dVf = [dVf filtfilt(filt_b2,filt_a2,decimate(V(cur_range),dsf))]; %do some high-pass filtering
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(ll,:) = dV;
        Vmatf(ll,:) = dVf;
        temp = cwt(dV,scales,'cmor1-1');
        phasegrams(ll,:,:) = angle(temp)';
        ampgrams(ll,:,:) = abs(temp)';
    end
    
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,2)+1:end) = [];
    
    binned_spks = nan(length(use_units),length(t_ax));
    sm_binned_spks = nan(length(use_units),length(t_ax));
    for cc = 1:length(use_units)
        temp = histc(Clusters{use_units(cc)}.times,t_ax);
        binned_spks(cc,:) = temp;
        sm_binned_spks(cc,:) = jmm_smooth_1d_cor(temp,spk_sm_win);
    end
    
    cur_use_trials = find(all_expt_id==forty_blocks(bb) & all_trial_completed==1);
    use_trial_start_inds = round(interp1(t_ax,1:length(t_ax),all_trial_start_times(cur_use_trials)));
    use_trial_stop_inds = round(interp1(t_ax,1:length(t_ax),all_trial_stop_times(cur_use_trials)));
    
    for ss = 1:length(repeat_inds)
        cur_trial_set = find(all_trial_seof(cur_use_trials)==repeat_inds(ss));
        cur(ss) = length(cur_trial_set);
        if length(cur_trial_set) > 0
            block_stim_lfps = zeros(length(cur_trial_set),length(use_lfps),length(lags));
            block_stim_lfpsf = zeros(length(cur_trial_set),length(use_lfps),length(lags));
            block_stim_phase = zeros(length(cur_trial_set),length(use_lfps),length(lags),length(wfreqs));
            block_stim_amp = zeros(length(cur_trial_set),length(use_lfps),length(lags),length(wfreqs));
            block_stim_spks = zeros(length(cur_trial_set),length(use_units),length(lags));
            block_stim_smspks = zeros(length(cur_trial_set),length(use_units),length(lags));
            for tt = 1:length(cur_trial_set)
                cur_inds = (use_trial_start_inds(cur_trial_set(tt))-backlag):(use_trial_start_inds(cur_trial_set(tt))+forwardlag);
                used_inds = find(cur_inds > 0 & cur_inds < length(t_ax));
                block_stim_lfps(tt,:,used_inds) = Vmat(:,cur_inds(used_inds));
                block_stim_lfpsf(tt,:,used_inds) = Vmatf(:,cur_inds(used_inds));
                block_stim_phase(tt,:,used_inds,:) = phasegrams(:,cur_inds(used_inds),:);
                block_stim_amp(tt,:,used_inds,:) = ampgrams(:,cur_inds(used_inds),:);
                block_stim_spks(tt,:,used_inds) = binned_spks(:,cur_inds(used_inds));
                block_stim_smspks(tt,:,used_inds) = sm_binned_spks(:,cur_inds(used_inds));
            end
            all_block_stim_lfps{ss} = cat(1,all_block_stim_lfps{ss},block_stim_lfps);
            all_block_stim_lfpsf{ss} = cat(1,all_block_stim_lfpsf{ss},block_stim_lfpsf);
            all_block_stim_phase{ss} = cat(1,all_block_stim_phase{ss},block_stim_phase);
            all_block_stim_amp{ss} = cat(1,all_block_stim_amp{ss},block_stim_amp);
            all_block_stim_spks{ss} = cat(1,all_block_stim_spks{ss},block_stim_spks);
            all_block_stim_smspks{ss} = cat(1,all_block_stim_smspks{ss},block_stim_smspks);
        end
        
    end
end
%%
use_lags = find(lags/Fsd > 0.5 & lags/Fsd < 3.2);
for ss =1:length(repeat_inds)
    phase_locking(ss,:,:,:) = squeeze(nansum(exp(1i*all_block_stim_phase{ss})));
    phase_locking(ss,:,:,:) = abs(phase_locking(ss,:,:,:))/sum(~isnan(all_block_stim_amp{ss}(:,1,1,1)));
    avg_filt_avglfps(ss,:,:,:) = squeeze(mean(all_block_stim_amp{ss}.*cos(all_block_stim_phase{ss})))/2;
    avg_filt_semlfps(ss,:,:,:) = squeeze(std(all_block_stim_amp{ss}.*cos(all_block_stim_phase{ss})))/2/sqrt(size(all_block_stim_phase{ss},1));
    avg_lfps(ss,:,:) = squeeze(nanmean(all_block_stim_lfps{ss}));
    sem_lfps(ss,:,:) = squeeze(nanstd(all_block_stim_lfps{ss}))/sqrt(sum(~isnan(all_block_stim_amp{ss}(:,1,1,1))));
    avg_lfpsf(ss,:,:) = squeeze(nanmean(all_block_stim_lfpsf{ss}));
    sem_lfpsf(ss,:,:) = squeeze(nanstd(all_block_stim_lfpsf{ss}))/sqrt(sum(~isnan(all_block_stim_lfpsf{ss}(:,1,1))));
    std_lfpsf(ss,:,:) = squeeze(nanstd(all_block_stim_lfpsf{ss}));
    avg_spks(ss,:,:) = squeeze(nanmean(all_block_stim_spks{ss}));
    std_spks(ss,:,:) = squeeze(nanstd(all_block_stim_spks{ss}));
    sem_spks(ss,:,:) = squeeze(nanstd(all_block_stim_spks{ss}))/sqrt(sum(~isnan(all_block_stim_lfpsf{ss}(:,1,1))));
    avg_smspks(ss,:,:) = squeeze(nanmean(all_block_stim_smspks{ss}));
    sem_smspks(ss,:,:) = squeeze(nanstd(all_block_stim_smspks{ss}))/sqrt(sum(~isnan(all_block_stim_lfpsf{ss}(:,1,1))));
    std_smspks(ss,:,:) = squeeze(nanstd(all_block_stim_smspks{ss}));
end
phase_locking_time = squeeze(mean(phase_locking,2));
avg_phaselocking = squeeze(mean(phase_locking(:,:,use_lags,:),3));elavg_phaselock = squeeze(mean(avg_phaselocking,2));

%%
save fortyhz_lfp_avgs_v3 lags Fsd phase_locking* avg_filt* avg_lfps* sem_lfps* std_lfps* avg_*spks sem_*spks std_*spks all_block_stim_spks all_block_stim_phase

%%
stim_fs = 1e4/117.5;
stim_times = (160)/stim_fs;
cur_unit = 13;
% stim1 = 4;
% stim2 = 5;
% stim3 = 6;
stim1 = 7;
stim2 = 8;
stim3 = 9;
figure
shadedErrorBar(lags/Fsd,squeeze(avg_smspks(stim1,cur_unit,:)),squeeze(sem_smspks(stim1,cur_unit,:)),{'color','b'});
hold on
shadedErrorBar(lags/Fsd,squeeze(avg_smspks(stim2,cur_unit,:)),squeeze(sem_smspks(stim2,cur_unit,:)),{'color','r'});
shadedErrorBar(lags/Fsd,squeeze(avg_smspks(stim3,cur_unit,:)),squeeze(sem_smspks(stim3,cur_unit,:)),{'color','k'});

