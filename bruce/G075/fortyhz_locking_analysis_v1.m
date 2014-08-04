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
use_lfps = [1:8:96];
% use_lfps = [1 63];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);

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
end
for bb = 1:length(forty_blocks)
    fprintf('Analyzing block %d of %d\n',bb,length(forty_blocks));
    
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
            for tt = 1:length(cur_trial_set)
                cur_inds = (use_trial_start_inds(cur_trial_set(tt))-backlag):(use_trial_start_inds(cur_trial_set(tt))+forwardlag);
                used_inds = find(cur_inds > 0 & cur_inds < length(t_ax));
                block_stim_lfps(tt,:,used_inds) = Vmat(:,cur_inds(used_inds));
                block_stim_lfpsf(tt,:,used_inds) = Vmatf(:,cur_inds(used_inds));
                block_stim_phase(tt,:,used_inds,:) = phasegrams(:,cur_inds(used_inds),:);
                block_stim_amp(tt,:,used_inds,:) = ampgrams(:,cur_inds(used_inds),:);
            end
            all_block_stim_lfps{ss} = cat(1,all_block_stim_lfps{ss},block_stim_lfps);
            all_block_stim_lfpsf{ss} = cat(1,all_block_stim_lfpsf{ss},block_stim_lfpsf);
            all_block_stim_phase{ss} = cat(1,all_block_stim_phase{ss},block_stim_phase);
            all_block_stim_amp{ss} = cat(1,all_block_stim_amp{ss},block_stim_amp);
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
end
phase_locking_time = squeeze(mean(phase_locking,2));
avg_phaselocking = squeeze(mean(phase_locking(:,:,use_lags,:),3));elavg_phaselock = squeeze(mean(avg_phaselocking,2));

%%
save fortyhz_lfp_avgs lags Fsd phase_locking* avg_filt* avg_lfps* sem_lfps* std_lfps*
%%
% stim_fs = 1e4/117.5;
% stim_times = (40:40:291)/stim_fs;
% stim_inds = 237:237:7*237;
% resh_avg_lfpsf = reshape(avg_lfpsf(:,:,1:7*237),length(repeat_inds),length(use_lfps),237,7);
% resh_avg_lfpsf = shiftdim(resh_avg_lfpsf,1);
% resh_avgs = zeros(length(use_lfps),237);
% use_j = [1 2 3 7 8 9];
% for i = 1:7
%     for j = 1:9
%         if ismember(j,use_j)
%             resh_avgs = resh_avgs + squeeze(resh_avg_lfpsf(:,:,i,j));
%         end
%     end
% end
% resh_avgs = resh_avgs/(7*length(use_j));

%%
stim_fs = 1e4/117.5;
stim_times = (160:160:320)/stim_fs;
cur_el = 3;
% close all
figure
for i = 1:9
    if i==7
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','b'});
    elseif i==8
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','r'});
%     elseif i==9
%         shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','k'});
    end
    hold on
end

yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
%%
stim_fs = 1e4/117.5;
stim_times = (160:160:320)/stim_fs;
cur_el = 3;
close all
for i = 1:9
    if i==7
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','b'});
    elseif i==8
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','r'});
%     elseif i==9
%         shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','k'});
    end
    hold on
end

yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
%%
stim_fs = 1e4/117.5;
stim_times = (160:160:320)/stim_fs;
cur_el = 1;
cur_el2 = 2;
close all
stim1 = 1;
stim2 = 2;
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(stim1,cur_el,:)),squeeze(sem_lfpsf(stim1,cur_el,:)),{'color','b'});
hold on
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(stim2,cur_el,:)),squeeze(sem_lfpsf(stim2,cur_el,:)),{'color','r'});
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(stim1,cur_el2,:)),squeeze(sem_lfpsf(stim1,cur_el2,:)),{'color','k'});
%     elseif i==9
%         shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','k'});

yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
%% TWO-ELECTRODES SAME STIM
stim_fs = 1e4/117.5;
stim_times = (160:160:320)/stim_fs;
cur_el1 = 1;
cur_el2 = 8;
close all
use_stim = 5;
figure
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(use_stim,cur_el1,:)),squeeze(sem_lfpsf(use_stim,cur_el1,:)),{'color','b'});
hold on
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(use_stim,cur_el2,:)),squeeze(sem_lfpsf(use_stim,cur_el2,:)),{'color','r'});

yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
%% FIRST-HALF SECOND-HALF COMPARE
stim_fs = 1e4/117.5;
stim_times = (160:160:320)/stim_fs;
cur_el = 1;
close all
use_stim = 5;
figure
shadedErrorBar(lags(1:942)/Fsd,squeeze(avg_lfpsf(use_stim,cur_el,1:942)),squeeze(sem_lfpsf(use_stim,cur_el,1:942)),{'color','b'});
hold on
shadedErrorBar(lags(1:942)/Fsd,squeeze(avg_lfpsf(use_stim,cur_el,943:1884)),squeeze(sem_lfpsf(use_stim,cur_el,943:1884)),{'color','r'});


%%
stim_fs = 1e4/117.5;
stim_times = (160)/stim_fs;
cur_el = 2;
stim1 = 4;
stim2 = 5;
% close all
figure
plot(lags/Fsd,squeeze(all_block_stim_lfpsf{stim1}(:,cur_el,:)),'b')
hold on
plot(lags/Fsd,squeeze(all_block_stim_lfpsf{stim2}(:,cur_el,:)),'r')
%     elseif i==9
%         plot(lags/Fsd,squeeze(all_block_stim_lfpsf{i}(:,cur_el,:)),'k')

% yl = ylim();

% for i = 1:length(stim_times)
%     line(stim_times([i i]),yl,'color','k')
% end

%%
stim_fs = 1e4/117.5;
stim_times = (160)/stim_fs;
el1 = 2;
el2 = 5;
stim = 1;

close all
plot(lags/Fsd,squeeze(all_block_stim_lfpsf{stim}(:,el1,:)),'b')
hold on
plot(lags/Fsd,squeeze(all_block_stim_lfpsf{stim}(:,el2,:)),'r')

%%
stim_fs = 1e4/117.5;
stim_times = (160)/stim_fs;
stim = 1;
trial1 = 1;
trial2 = 3;
close all
plot(lags/Fsd,squeeze(all_block_stim_lfpsf{stim}(trial1,:,:)),'b')
hold on
plot(lags/Fsd,squeeze(all_block_stim_lfpsf{stim}(trial2,:,:)),'r')