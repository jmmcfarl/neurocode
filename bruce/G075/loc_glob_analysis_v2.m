clear all
close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

%%
loc_glob_blocks = [8 15 20 27 33 39 49] - 6;
movie_ids = [1e3 2e3 4e3 5e3 6e3 7e3 8e3 9e3 10e3 11e3 12e3 13e3];
all_expt_id = [];
for bb = 1:length(loc_glob_blocks)
    n_trials(bb) = length(Expts{loc_glob_blocks(bb)}.Trials);
    trial_start_times{bb} = [Expts{loc_glob_blocks(bb)}.Trials(:).Start]/1e4;
    trial_stop_times{bb} = [Expts{loc_glob_blocks(bb)}.Trials(:).End]/1e4;
    trial_movie_ids{bb} = [Expts{loc_glob_blocks(bb)}.Trials(:).backMov];
    trial_completed{bb} = [Expts{loc_glob_blocks(bb)}.Trials(:).Result];
    trial_durs{bb} = trial_stop_times{bb} - trial_start_times{bb};
    all_expt_id = [all_expt_id loc_glob_blocks(bb)*ones(1,n_trials(bb))];
end
all_trial_start_times = cell2mat(trial_start_times);
all_trial_stop_times = cell2mat(trial_stop_times);
all_trial_movie_ids = cell2mat(trial_movie_ids);
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

all_block_stim_lfps = [];
all_block_stim_lfpsf = [];
all_block_stim_phase = [];
all_block_stim_amp = [];
all_block_stim_spks = [];
all_block_stim_smspks = [];
for bb = 1:length(loc_glob_blocks)
    % bb = 1;
    Vmat = [];
    Vmatf = [];
    phasegrams = [];
    ampgrams = [];
    load(sprintf('Expt%dClusterTimes.mat',loc_glob_blocks(bb)));
    
    filename = sprintf('Expt%dFullVmean.mat',loc_glob_blocks(bb));
    load(filename);
    
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',loc_glob_blocks(bb),use_lfps(ll));
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
    
    cur_use_trials = find(all_expt_id==loc_glob_blocks(bb) & all_trial_completed==1);
    use_trial_start_inds = round(interp1(t_ax,1:length(t_ax),all_trial_start_times(cur_use_trials)));
    use_trial_stop_inds = round(interp1(t_ax,1:length(t_ax),all_trial_stop_times(cur_use_trials)));
    
    block_stim_lfps = nan(length(movie_ids),5,length(use_lfps),length(lags));
    block_stim_lfpsf = nan(length(movie_ids),5,length(use_lfps),length(lags));
    block_stim_phase = nan(length(movie_ids),5,length(use_lfps),length(lags),length(wfreqs));
    block_stim_amp = nan(length(movie_ids),5,length(use_lfps),length(lags),length(wfreqs));
    block_stim_spks = nan(length(movie_ids),5,length(use_units),length(lags));
    block_stim_smspks = nan(length(movie_ids),5,length(use_units),length(lags));
    for ss = 1:length(movie_ids)
       cur_trial_set = find(all_trial_movie_ids(cur_use_trials)==movie_ids(ss));
       for tt = 1:length(cur_trial_set)
           cur_inds = (use_trial_start_inds(cur_trial_set(tt))-backlag):(use_trial_start_inds(cur_trial_set(tt))+forwardlag);
           used_inds = find(cur_inds > 0 & cur_inds < length(t_ax));
           block_stim_lfps(ss,tt,:,used_inds) = Vmat(:,cur_inds(used_inds));
           block_stim_lfpsf(ss,tt,:,used_inds) = Vmatf(:,cur_inds(used_inds));
           block_stim_phase(ss,tt,:,used_inds,:) = phasegrams(:,cur_inds(used_inds),:);
           block_stim_amp(ss,tt,:,used_inds,:) = ampgrams(:,cur_inds(used_inds),:);
           block_stim_spks(ss,tt,:,used_inds) = binned_spks(:,cur_inds(used_inds));
           block_stim_smspks(ss,tt,:,used_inds) = sm_binned_spks(:,cur_inds(used_inds));
       end
    end
        
    all_block_stim_lfps = cat(2,all_block_stim_lfps,block_stim_lfps);
    all_block_stim_lfpsf = cat(2,all_block_stim_lfpsf,block_stim_lfpsf);
    all_block_stim_amp = cat(2,all_block_stim_amp,block_stim_amp);
    all_block_stim_phase = cat(2,all_block_stim_phase,block_stim_phase);
    all_block_stim_spks = cat(2,all_block_stim_spks,block_stim_spks);
    all_block_stim_smspks = cat(2,all_block_stim_smspks,block_stim_smspks);
end

%%
use_lags = find(lags/Fsd > 0.5 & lags/Fsd < 3.2);
for ss =1:length(movie_ids)
n_trials(ss) = sum(~isnan(all_block_stim_amp(ss,:,1,1,1)),2);
    phase_locking(ss,:,:,:) = squeeze(nansum(exp(1i*all_block_stim_phase(ss,:,:,:,:))));
    phase_locking(ss,:,:,:) = abs(phase_locking(ss,:,:,:))/n_trials(ss);
    avg_lfps(ss,:,:) = squeeze(nanmean(all_block_stim_lfps(ss,:,:,:,:)));
    sem_lfps(ss,:,:) = squeeze(nanstd(all_block_stim_lfps(ss,:,:,:,:)))/sqrt(n_trials(ss));
    avg_lfpsf(ss,:,:) = squeeze(nanmean(all_block_stim_lfpsf(ss,:,:,:,:)));
    sem_lfpsf(ss,:,:) = squeeze(nanstd(all_block_stim_lfpsf(ss,:,:,:,:)))/sqrt(n_trials(ss));
    std_lfpsf(ss,:,:) = squeeze(nanstd(all_block_stim_lfpsf(ss,:,:,:,:)));
    avg_smspk(ss,:,:) = squeeze(nanmean(all_block_stim_smspks(ss,:,:,:,:)));
    sem_smspk(ss,:,:) = squeeze(nanstd(all_block_stim_smspks(ss,:,:,:,:)))/sqrt(n_trials(ss));
    std_smspk(ss,:,:) = squeeze(nanstd(all_block_stim_smspks(ss,:,:,:,:)));
    avg_spk(ss,:,:) = squeeze(nanmean(all_block_stim_spks(ss,:,:,:,:)));
    sem_spk(ss,:,:) = squeeze(nanstd(all_block_stim_spks(ss,:,:,:,:)))/sqrt(n_trials(ss));
    std_spk(ss,:,:) = squeeze(nanstd(all_block_stim_spks(ss,:,:,:,:)));
end
phase_locking_time = squeeze(mean(phase_locking,2));
avg_phaselocking = squeeze(mean(phase_locking(:,:,use_lags,:),3));
elavg_phaselock = squeeze(mean(avg_phaselocking,2));
%%
save -v7.3 loc_glob_lfp_avgs_v2 lags Fsd phase_locking* avg_lfps* sem_lfps* std_lfps* wfreqs avg_smspk sem_smspk std_smspk all_block_stim_spks all_block_stim_phase

%%
temp = permute(avg_smspk,[1 3 2]);
temp = reshape(temp,[12*length(lags) 96]);
ov_avg_rate = mean(temp);
avg_nsmspk = bsxfun(@rdivide,avg_smspk,ov_avg_rate);
lg_lfp = squeeze(mean(avg_lfpsf,2));
lg_spks = squeeze(mean(avg_nsmspk,2));
lg_al_spks = avg_nsmspk;

for cur_un = 1:96
    shadedErrorBar(lags/Fsd,squeeze(avg_smspk(10,cur_un,:)),squeeze(sem_smspk(10,cur_un,:)),{'color','b'});
    hold on
    shadedErrorBar(lags/Fsd,squeeze(avg_smspk(11,cur_un,:)),squeeze(sem_smspk(11,cur_un,:)),{'color','r'});
    shadedErrorBar(lags/Fsd,squeeze(avg_smspk(12,cur_un,:)),squeeze(sem_smspk(12,cur_un,:)),{'color','k'});
    pause
    clf
end

%%
load ./fortyhz_lfp_avgs_v3
%%
temp = permute(avg_smspks,[1 3 2]);
temp = reshape(temp,[9*length(lags) 96]);
ov_avg_rate = mean(temp);
avg_nsmspk = bsxfun(@rdivide,avg_smspk,ov_avg_rate);
frty_lfp = squeeze(mean(avg_lfpsf,2));
frty_spks = squeeze(mean(avg_nsmspk,2));
frty_all_spks = avg_nsmspk;

%%
% close all
frty_stid = 5;
lg_stid = 11;

stim_fs = 1e4/117.5;
used_times = find(mod(lags(1:1880)/Fsd,0.47) > 0.15);
% used_times = find(mod(lags(1:1880)/Fsd,0.47) > 0);
used_times(1:235) = [];

ov_avg_rate_hz = ov_avg_rate*Fsd;
all_scfac = squeeze(lg_al_spks(lg_stid,:,:)./(frty_all_spks(frty_stid,:,:)+0.1));
% all_scfac = squeeze(lg_al_spks(lg_stid,:,:) - (frty_all_spks(frty_stid,:,:)+0.1));
maxval = 4;
minval = 0.25;
all_scfac(all_scfac > maxval) = maxval;
all_scfac(all_scfac < minval) = minval;
all_scfac(ov_avg_rate_hz < 1) = nan;

% scfac = lg_spks(lg_stid,:)./frty_spks(frty_stid,:);
scfac = squeeze(nanmean(all_scfac));
p = polyfit(frty_lfp(frty_stid,1:1880),lg_lfp(lg_stid,1:1880),1);
pval = polyval(p,frty_lfp(frty_stid,:));
% resid = lg_lfp(lg_stid,:) - pval;
% resid = lg_lfp(lg_stid,:) - frty_lfp(frty_stid,:);
resid = lg_lfp(lg_stid,:);
[a,b] = corrcoef(resid(used_times),scfac(used_times))
% [at,bt] = corrcoef(lg_lfp(lg_stid,1:1880),scfac(1:1880))

% figure
% subplot(2,1,1)
% % plot(lags/Fsd,frty_lfp(frty_stid,:))
% hold on
% % plot(lags/Fsd,lg_lfp(lg_stid,:),'r')
% plot(lags/Fsd,resid,'k')
% stim_times = (0:40:320)/stim_fs;
% yl = ylim();
% for i = 1:length(stim_times)
%     line(stim_times([i i]),yl,'color','k')
% end
% xlim([0.47 3.75])
% subplot(2,1,2)
% % plot(lags/Fsd,frty_spks(frty_stid,:))
% hold on
% % plot(lags/Fsd,lg_spks(lg_stid,:),'r')
% plot(lags/Fsd,scfac,'k')
% % shadedErrorBar(lags/Fsd,mean(all_scfac),std(all_scfac)/sqrt(96),'k')

figure
shadedErrorBar(lags/Fsd,nanmean(all_scfac),nanstd(all_scfac)/sqrt(96),'k')
hold on
plot(lags/Fsd,zscore(resid)/6+1,'r','linewidth',2)
plot(lags/Fsd,zscore(frty_lfp(frty_stid,:))/6+1,'b','linewidth',1)
% plot(lags/Fsd,zscore(lg_lfp(lg_stid,:))/6+1,'g','linewidth',1)
stim_times = (0:40:320)/stim_fs;
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
xlim([0.47 3.75])
hold on



