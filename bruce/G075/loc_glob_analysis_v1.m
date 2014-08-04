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
use_lfps = [1:8:96];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);

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
for bb = 1:length(loc_glob_blocks)
    % bb = 1;
    Vmat = [];
    Vmatf = [];
    phasegrams = [];
    ampgrams = [];
    
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
            dV = [dV decimate(V(cur_range),dsf)];
            dVf = [dVf filtfilt(filt_b,filt_a,decimate(V(cur_range),dsf))]; %do some high-pass filtering
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
    
    cur_use_trials = find(all_expt_id==loc_glob_blocks(bb) & all_trial_completed==1);
    use_trial_start_inds = round(interp1(t_ax,1:length(t_ax),all_trial_start_times(cur_use_trials)));
    use_trial_stop_inds = round(interp1(t_ax,1:length(t_ax),all_trial_stop_times(cur_use_trials)));
    
    block_stim_lfps = nan(length(movie_ids),5,length(use_lfps),length(lags));
    block_stim_lfpsf = nan(length(movie_ids),5,length(use_lfps),length(lags));
    block_stim_phase = nan(length(movie_ids),5,length(use_lfps),length(lags),length(wfreqs));
    block_stim_amp = nan(length(movie_ids),5,length(use_lfps),length(lags),length(wfreqs));
    for ss = 1:length(movie_ids)
       cur_trial_set = find(all_trial_movie_ids(cur_use_trials)==movie_ids(ss));
       for tt = 1:length(cur_trial_set)
           cur_inds = (use_trial_start_inds(cur_trial_set(tt))-backlag):(use_trial_start_inds(cur_trial_set(tt))+forwardlag);
           used_inds = find(cur_inds > 0 & cur_inds < length(t_ax));
           block_stim_lfps(ss,tt,:,used_inds) = Vmat(:,cur_inds(used_inds));
           block_stim_lfpsf(ss,tt,:,used_inds) = Vmatf(:,cur_inds(used_inds));
           block_stim_phase(ss,tt,:,used_inds,:) = phasegrams(:,cur_inds(used_inds),:);
           block_stim_amp(ss,tt,:,used_inds,:) = ampgrams(:,cur_inds(used_inds),:);
       end
    end
        
    all_block_stim_lfps = cat(2,all_block_stim_lfps,block_stim_lfps);
    all_block_stim_lfpsf = cat(2,all_block_stim_lfpsf,block_stim_lfpsf);
    all_block_stim_amp = cat(2,all_block_stim_amp,block_stim_amp);
    all_block_stim_phase = cat(2,all_block_stim_phase,block_stim_phase);
end

%%
use_lags = find(lags/Fsd > 0.5 & lags/Fsd < 3.2);
for ss =1:length(movie_ids)
n_trials(ss) = sum(~isnan(all_block_stim_amp(ss,:,1,1,1)),2);
    phase_locking(ss,:,:,:) = squeeze(nansum(exp(1i*all_block_stim_phase(ss,:,:,:))));
    phase_locking(ss,:,:,:) = abs(phase_locking(ss,:,:,:))/n_trials(ss);
    avg_lfps(ss,:,:) = squeeze(nanmean(all_block_stim_lfps(ss,:,:,:)));
    sem_lfps(ss,:,:) = squeeze(nanstd(all_block_stim_lfps(ss,:,:,:)))/sqrt(n_trials(ss));
    avg_lfpsf(ss,:,:) = squeeze(nanmean(all_block_stim_lfpsf(ss,:,:,:)));
    sem_lfpsf(ss,:,:) = squeeze(nanstd(all_block_stim_lfpsf(ss,:,:,:)))/sqrt(n_trials(ss));
    std_lfpsf(ss,:,:) = squeeze(nanstd(all_block_stim_lfpsf(ss,:,:,:)));
end
phase_locking_time = squeeze(mean(phase_locking,2));
avg_phaselocking = squeeze(mean(phase_locking(:,:,use_lags,:),3));
elavg_phaselock = squeeze(mean(avg_phaselocking,2));
%%
save loc_glob_lfp_avgs lags Fsd phase_locking* avg_lfps* sem_lfps* std_lfps*
%%
% for ss =1:length(movie_ids)
%     avg_lfps(ss,:,:) = squeeze(nanmean(all_block_stim_lfps(ss,:,:,:),2));
%     sem_lfps(ss,:,:) = squeeze(nanstd(all_block_stim_lfps(ss,:,:,:),[],2))/sqrt(sum(~isnan(all_block_stim_lfps(ss,:,1,1))));
%     avg_lfpsf(ss,:,:) = squeeze(nanmean(all_block_stim_lfpsf(ss,:,:,:),2));
%     sem_lfpsf(ss,:,:) = squeeze(nanstd(all_block_stim_lfpsf(ss,:,:,:),[],2))/sqrt(sum(~isnan(all_block_stim_lfps(ss,:,1,1))));
% end
%%
stim_fs = 1e4/117.5;
stim_times = (40:40:320)/stim_fs;
cur_el = 2;
close all
for i = 1:12
    if i==10
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','b'});
    elseif i==8
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','r'});
    elseif i==9
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','k'});
    elseif i==10
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','g'});
    end
    hold on
end

yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end

%%
stim_fs = 1e4/117.5;
stim_times = (40:40:320)/stim_fs;
cur_el = 2;
close all
for i = 1:12
    if i==10
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','b'});
    elseif i==8
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','r'});
    elseif i==9
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','k'});
    elseif i==10
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','g'});
    end
    hold on
end

yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end

%%
stim_fs = 1e4/117.5;
stim_times = (40:40:160)/stim_fs;
cur_el = 1;
cur_stim = 12;
% close all
figure
shadedErrorBar(lags(1:472)/Fsd,squeeze(avg_lfpsf(cur_stim,cur_el,1:472)),squeeze(sem_lfpsf(cur_stim,cur_el,1:472)),{'color','b'});
hold on
shadedErrorBar(lags(1:472)/Fsd,squeeze(avg_lfpsf(cur_stim,cur_el,473:944)),squeeze(sem_lfpsf(cur_stim,cur_el,473:944)),{'color','r'});

yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end

%%
use_lags = find(lags/Fsd > 0.5 & lags/Fsd < 3.2);

phase_locking = squeeze(nansum(exp(1i*all_block_stim_phase),2));
phase_locking = abs(phase_locking)/sum(~isnan(all_block_stim_amp(1,:,1,1,1)),2);
phase_locking_time = squeeze(mean(phase_locking,2));
avg_phaselocking = squeeze(mean(phase_locking(:,:,use_lags,:),3));
elavg_phaselock = squeeze(mean(avg_phaselocking,2));

amp_locking = 1./squeeze(nanstd(all_block_stim_amp,[],2));
avg_amplocking = squeeze(nanmean(amp_locking(:,:,use_lags,:),3));
elavg_amplock = squeeze(nanmean(avg_amplocking,2));

amp_avg = squeeze(nanmean(all_block_stim_amp,2));
avg_ampavg = squeeze(nanmean(amp_avg(:,:,use_lags,:),3));
elavg_ampavg = squeeze(nanmean(avg_ampavg,2));

%%
% phase_ax = linspace(-pi,pi,10);
% cur_el = 1;
% % use_freqs = [27 22 16];
% use_freqs = [27];
% for ss = 1:length(movie_ids)
%     for c = 1:length(use_freqs)
%     figure(1)
%         subplot(1,length(use_freqs),c)
%         pcolor(lags/Fsd,1:size(all_block_stim_phase,2),squeeze(all_block_stim_phase(ss,:,cur_el,:,use_freqs(c))));shading flat
%      figure(2)
%         subplot(1,length(use_freqs),c)
%         plot(lags/Fsd,squeeze(phase_locking(ss,cur_el,:,use_freqs(c))));
%      figure(3)
%         subplot(1,length(use_freqs),c)
%         temp = zeros(length(phase_ax),length(lags));
%         for i = 1:length(lags)
%             temp(:,i) = hist(squeeze(all_block_stim_phase(ss,:,cur_el,i,use_freqs(c))),phase_ax);
%         end
%         pcolor(lags/Fsd,phase_ax,temp);shading flat
%     end
%     ss
%     
%     pause
%     figure(1);clf;figure(2);clf;figure(3);clf;
%     
% end

%%
el_num = 1;
for ss = 1:length(movie_ids)
    ss
    pcolor(lags/Fsd,wfreqs,squeeze(phase_locking(ss,el_num,:,:))');shading flat
set(gca,'yscale','log')
caxis([0 0.7])
    pause
clf
end

%%
for el = 1:length(use_lfps)
    el
    imagesc(wfreqs,1:length(movie_ids),squeeze(avg_phaselocking(:,el,:)));shading flat
    pause
    clf
end