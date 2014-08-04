clear all
close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

%%
sim_sac_blocks = [10 14 18 22 24 28 32 37 40 42 50] - 6;
% sim_sac_blocks = [18 28 40] - 6; %sim sac
% sim_sac_blocks = [10 22 32 42] - 6; %sim sac a
% sim_sac_blocks = [14 24 37 50] - 6; %sim sac b

repeat_inds = [1e3 2e3 3e3 2.01e5 2.02e5 2.03e5 4.01e5 4.02e5 4.03e5];
all_expt_id = [];
for bb = 1:length(sim_sac_blocks)
    n_trials(bb) = length(Expts{sim_sac_blocks(bb)}.Trials);
    trial_start_times{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).Start]/1e4;
    trial_stop_times{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).End]/1e4;
    trial_seof{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).seof];
    trial_completed{bb} = [Expts{sim_sac_blocks(bb)}.Trials(:).Result];
    trial_durs{bb} = trial_stop_times{bb} - trial_start_times{bb};
    all_expt_id = [all_expt_id sim_sac_blocks(bb)*ones(1,n_trials(bb))];
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
for bb = 1:length(sim_sac_blocks)
    fprintf('Analyzing block %d of %d\n',bb,length(sim_sac_blocks));
    
            filename = sprintf('Expt%dFullVmean.mat',sim_sac_blocks(bb));
        load(filename);

    Vmat = [];
    Vmatf = [];
    phasegrams = [];
    ampgrams = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',sim_sac_blocks(bb),use_lfps(ll));
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
    
    cur_use_trials = find(all_expt_id==sim_sac_blocks(bb) & all_trial_completed==1);
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
    
    temp = squeeze(nansum(exp(1i*all_block_stim_phase{ss}),2));
    temp = abs(temp)/length(use_lfps);
    elcoh_locking(ss,:,:) = squeeze(mean(temp));
    avg_lfps(ss,:,:) = squeeze(nanmean(all_block_stim_lfps{ss}));
    sem_lfps(ss,:,:) = squeeze(nanstd(all_block_stim_lfps{ss}))/sqrt(sum(~isnan(all_block_stim_amp{ss}(:,1,1,1))));
    avg_lfpsf(ss,:,:) = squeeze(nanmean(all_block_stim_lfpsf{ss}));
    sem_lfpsf(ss,:,:) = squeeze(nanstd(all_block_stim_lfpsf{ss}))/sqrt(sum(~isnan(all_block_stim_lfpsf{ss}(:,1,1))));
    std_lfpsf(ss,:,:) = squeeze(nanstd(all_block_stim_lfpsf{ss}));
end
phase_locking_time = squeeze(mean(phase_locking,2));
avg_phaselocking = squeeze(mean(phase_locking(:,:,use_lags,:),3));
elavg_phaselock = squeeze(mean(avg_phaselocking,2));
%%
save sim_sac_lfp_avgs lags Fsd phase_locking* avg_lfps* sem_lfps* std_lfps* all_block_stim_lfps*

%%
stim_fs = 1e4/117.5;
stim_times = (40:40:291)/stim_fs;
stim_inds = [1 round(interp1(lags/Fsd,1:length(lags),stim_times)) length(lags)];

new_lags = 1:min(diff(stim_inds));
stim_cnt = 1;
clear stim_trig_avgs
for ss = 1:length(repeat_inds)
    cur_lfps = all_block_stim_lfpsf{ss};
    for ll =1:length(stim_inds)-1
        cur_inds = stim_inds(ll):stim_inds(ll+1);
        cur_inds(length(new_lags)+1:end) = [];
       cur_set = cur_lfps(:,:,cur_inds);
       stim_trig_avgs(stim_cnt,:,:) = mean(cur_set);
       stim_type(stim_cnt) = ss;
       stim_cnt = stim_cnt + 1;
    end
end
u1 = find(stim_type > 1 &stim_type <= 3);
u2 = find(stim_type > 3 & stim_type<= 6);
u3 = find(stim_type > 7 & stim_type<= 9);
avg_u1 = squeeze(mean(stim_trig_avgs(u1,:,:)));
avg_u2 = squeeze(mean(stim_trig_avgs(u2,:,:)));
avg_u3 = squeeze(mean(stim_trig_avgs(u3,:,:)));

mean_sigs = zeros(size(avg_lfpsf));
for ss = 1:length(repeat_inds)
    if ss == 1 | ss == 2 | ss == 3
        for ll = 1:length(stim_inds)-1
           cur_inds = stim_inds(ll):stim_inds(ll+1);
        cur_inds(length(new_lags)+1:end) = [];
           mean_sigs(ss,:,cur_inds) = avg_u1;
        end
    elseif ss == 4 | ss == 5 | ss == 6
        for ll = 1:length(stim_inds)-1
           cur_inds = stim_inds(ll):stim_inds(ll+1);
        cur_inds(length(new_lags)+1:end) = [];
           mean_sigs(ss,:,cur_inds) = avg_u2;
        end
    else
        for ll = 1:length(stim_inds)-1
           cur_inds = stim_inds(ll):stim_inds(ll+1);
        cur_inds(length(new_lags)+1:end) = [];
           mean_sigs(ss,:,cur_inds) = avg_u3;
        end
    end
end

for ss =1:length(repeat_inds)
    cur_block_stim_ms = bsxfun(@minus,all_block_stim_lfpsf{ss},mean_sigs(ss,:,:));
    avg_lfpsf_ms(ss,:,:) = squeeze(nanmean(cur_block_stim_ms));
    sem_lfpsf_ms(ss,:,:) = squeeze(nanstd(cur_block_stim_ms))/sqrt(sum(~isnan(cur_block_stim_ms(:,1,1))));
    
    for kk = 1:size(cur_block_stim_ms,1) %cycle over trials
        fprintf('Trial %d of %d\n',kk,size(cur_block_stim_ms,1));
        for jj = 1:size(cur_block_stim_ms,2) %cycle over electrodes
            cur_sig = squeeze(cur_block_stim_ms(kk,jj,:));
            temp = cwt(cur_sig,scales,'cmor1-1');
            all_ms_phases(kk,jj,:,:) = angle(temp)';
        end
    end
    phase_locking_ms(ss,:,:,:) = squeeze(nansum(exp(1i*all_ms_phases)));
    phase_locking_ms(ss,:,:,:) = abs(phase_locking_ms(ss,:,:,:))/sum(~isnan(all_ms_phases(:,1,1,1)));
end
%%
for ss =1:length(repeat_inds)
    cur_block_stim_ms = bsxfun(@minus,all_block_stim_lfpsf{ss},avg_lfpsf(ss,:,:));    
    for kk = 1:size(cur_block_stim_ms,1) %cycle over trials
        fprintf('Trial %d of %d\n',kk,size(cur_block_stim_ms,1));
        for jj = 1:size(cur_block_stim_ms,2) %cycle over electrodes
            cur_sig = squeeze(cur_block_stim_ms(kk,jj,:));
            temp = cwt(cur_sig,scales,'cmor1-1');
            all_ms_phases(kk,jj,:,:) = angle(temp)';
        end
    end
end
temp = squeeze(nansum(exp(1i*all_ms_phases),2));
temp = abs(temp)/length(use_lfps);
elcoh_locking_ms(ss,:,:) = squeeze(mean(temp));

%%
stim_fs = 1e4/117.5;
stim_times = (40:40:291)/stim_fs;
cur_el = 5;
close all
% cmap = jet(9);
% for i = 1:9
%     plot(lags/Fsd,squeeze(avg_lfps(i,cur_el,:)),'color',cmap(i,:))
%     hold on
% end
% for i = 1:9
%     if ismember(i,1:3)
%         plot(lags/Fsd,squeeze(avg_lfps(i,cur_el,:)),'b')
%     elseif ismember(i,4:6)
%         plot(lags/Fsd,squeeze(avg_lfps(i,cur_el,:)),'r')
%     else
%         plot(lags/Fsd,squeeze(avg_lfps(i,cur_el,:)),'k')
%     end
%     hold on
% end
for i = 1:9
    if i==7
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','b'});
    elseif i==8
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','r'});
    elseif i==9
        shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(i,cur_el,:)),squeeze(sem_lfpsf(i,cur_el,:)),{'color','k'});
    end
    hold on
end

yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end

%% TWO-ELECTRODES SAME STIM
stim_fs = 1e4/117.5;
stim_times = (40:40:320)/stim_fs;
cur_el1 = 5;
cur_el2 = 8;
% close all
% use_stim = 4;
% use_stim2 = 5;
% use_stim3 = 6;
use_stim = 7;
use_stim2 =8;
use_stim3 = 9;
figure
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(use_stim,cur_el1,:)),squeeze(sem_lfpsf(use_stim,cur_el1,:)),{'color','b'});
hold on
% shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(use_stim,cur_el2,:)),squeeze(sem_lfpsf(use_stim,cur_el2,:)),{'color','r'});
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(use_stim2,cur_el1,:)),squeeze(sem_lfpsf(use_stim2,cur_el1,:)),{'color','k'});
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(use_stim3,cur_el1,:)),squeeze(sem_lfpsf(use_stim3,cur_el1,:)),{'color','r'});
plot(lags/Fsd,squeeze(mean_sigs(use_stim,cur_el1,:)),'k','linewidth',2)
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end

%% TWO-ELECTRODES SAME STIM (MEAN SUBTRACTED)
stim_fs = 1e4/117.5;
stim_times = (40:40:320)/stim_fs;
cur_el1 = 1;
cur_el2 = 8;
% close all
use_stim = 7;
use_stim2 = 8;
use_stim3 = 9;
figure
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf_ms(use_stim,cur_el1,:)),squeeze(sem_lfpsf_ms(use_stim,cur_el1,:)),{'color','b'});
hold on
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf_ms(use_stim2,cur_el1,:)),squeeze(sem_lfpsf_ms(use_stim2,cur_el1,:)),{'color','k'});
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf_ms(use_stim3,cur_el1,:)),squeeze(sem_lfpsf_ms(use_stim3,cur_el1,:)),{'color','r'});
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
%%
cur_el = 2;
stim1 = 9;
% stim2 = 5;
% close all
figure
hold on
plot(lags/Fsd,squeeze(all_block_stim_lfpsf{stim1}(1,:,:)),'b')
plot(lags/Fsd,squeeze(all_block_stim_lfpsf{stim1}(2,:,:)),'r')
% plot(lags/Fsd,squeeze(all_block_stim_lfpsf{stim1}(3,:,:)),'k')
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(stim1,1,:)),squeeze(sem_lfpsf(stim1,1,:)),{'color','k'});

yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end

figure
hold on
plot(lags/Fsd,squeeze(all_block_stim_lfpsf{stim1}(:,1,:)),'b')
plot(lags/Fsd,squeeze(all_block_stim_lfpsf{stim2}(:,1,:)),'r')
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(stim1,1,:)),squeeze(sem_lfpsf(stim1,1,:)),{'color','k'});
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(stim2,1,:)),squeeze(sem_lfpsf(stim2,1,:)),{'color','g'});

yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end

%% TWO-ELECTRODES SAME STIM
Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
siz = [1280 1280];
Fs = 1/Pix2Deg;
xax = linspace(-Nxp/2,Nxp/2,Nxp)/Fs; yax = linspace(-Nyp/2,Nyp/2,Nyp)/Fs;
[XAX,YAX] = meshgrid(xax,yax);
x0_avg = 0.35;
y0_avg = -0.4;
xw = 1;yw=1;

close all

stim_fs = 1e4/117.5;
stim_times = (0:40:320)/stim_fs;
cur_el1 = 5;
use_stim = 9;
figure
shadedErrorBar(lags/Fsd,squeeze(avg_lfpsf(use_stim,cur_el1,:)),squeeze(sem_lfpsf(use_stim,cur_el1,:)),{'color','b'});
hold on
plot(lags/Fsd,squeeze(mean_sigs(use_stim,cur_el1,:)),'k','linewidth',2);
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
cd ~/Data/bruce/Expt_1_8_13_imfolder
se_offset = repeat_inds(use_stim);
figure
for i = 1:8
    cur_fname = sprintf('IM1%d.png',se_offset+i);
    cur_im = imread(cur_fname);
    imagesc(xax,yax,flipud(cur_im));colormap(gray);set(gca,'ydir','normal');
    hold on
    rectangle('Position',[x0_avg y0_avg xw,yw],'edgecolor','r');
    xlim([-2 2]); ylim([-2 2])

    figure(1)
    xlim([stim_times(i) stim_times(i)+0.75])
    
    pause
     figure(2)
   clf
end
