% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
Expt_name = 'G091';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
addpath([dir_prefix '/James_scripts/bruce/temporary_scripts/'])

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
%% PARSE TRIAL DATA STRUCTURES

stim_fs = 100; %in Hz
dt = 0.005;
beg_buffer = round(0.15/dt);
end_buffer = round(0.15/dt);

%%
% cur_expt_set = [44:49];
cur_expt_set = [44 45 46 47 48 49];
all_trial_start_times = [];
all_trial_stop_times = [];
all_trial_durs = [];
all_trial_nph = [];
% all_trial_or = [];
all_trial_tf = [];
for ee = 1:length(cur_expt_set)
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
%     trial_ids = [Expts{cur_expt}.Trials(:).id]; 
%     trial_or = [Expts{cur_expt}.Trials(:).or];
    trial_nph = [Expts{cur_expt}.Trials(:).nph];
    trial_tf = [Expts{cur_expt}.Trials(:).tf];
    
    use_trials = find(trial_durs >= 0.5);
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_stop_times = cat(1,all_trial_stop_times,trial_end_times(use_trials)');
%     all_trial_or = cat(1,all_trial_or,trial_or(use_trials)');
    all_trial_tf = cat(1,all_trial_tf,trial_tf(use_trials)');
    all_trial_nph = cat(1,all_trial_nph,trial_nph(use_trials)');
end

%%
% use_lfps = [4:4:96];
use_lfps = [4:2:96];
dsf = 60;
Fs = 1/3.333307000208032e-05;
Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[0.5]/niqf,'high');

all_V = [];
all_t_lfp = [];

for ee = 1:length(cur_expt_set)
    cur_expt = cur_expt_set(ee);
    
    filename = sprintf('Expt%dFullVmean.mat',cur_expt);
    load(filename);
    
    Vmat = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',cur_expt,use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        
        V = V + FullV.sumscale*sumv;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dV = [];
        %splice together multiple blocks
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            curV = decimate(V(cur_range),dsf);
            curV = filtfilt(filt_b,filt_a,curV);
            dV = [dV curV];
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(:,ll) = dV;
    end
    
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,1)+1:end) = [];
    
    fprintf('LFP len: %d\n',range(t_ax));
    
    all_V = cat(1,all_V,Vmat);
    all_t_lfp = cat(1,all_t_lfp,t_ax');
    
end

%%
use_trial_start_times = all_trial_start_times;
use_trial_end_times = all_trial_stop_times;

all_drifting_grate_trials = find(all_trial_nph == 0);
dg_tstart_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),use_trial_start_times(all_drifting_grate_trials)));
dg_tstop_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),all_trial_stop_times(all_drifting_grate_trials)));

all_rp_grate_trials = find(all_trial_nph > 0);
rp_tstart_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),use_trial_start_times(all_rp_grate_trials)));
rp_tstop_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),all_trial_stop_times(all_rp_grate_trials)));

%%

%WAVELET SCALES this gives log-freq spacing ~(2-100) hz
nwfreqs = 35;
min_freq = 2; max_freq = 100;
min_scale = 1/max_freq*Fsd;
max_scale = 1/min_freq*Fsd;
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

forlag = round(Fsd*4);
backlag = round(0.2*Fsd);

un_tfs = unique(all_trial_tf);
un_tfs(1) = [];

for ll = 1:length(use_lfps);
    fprintf('LFP %d of %d\n',ll,length(use_lfps));
    temp = cwt(all_V(:,ll),scales,'cmor1-1');
    cur_ampgram = abs(temp)';
    avg_pow = mean(cur_ampgram);
    std_pow = std(cur_ampgram);
    
    un_wfreqs = linspace(wfreqs(1),wfreqs(end),100);
    avg_pow_int = interp1(wfreqs,avg_pow,un_wfreqs);
    
    powfun = @(a,x)(a(1)*x.^a(2));
    BETA = nlinfit(un_wfreqs,avg_pow_int,powfun,[max(avg_pow) -2]);
    powfit = un_wfreqs.^BETA(2)*BETA(1);
    interp_powfit = interp1(un_wfreqs,powfit,wfreqs);
    
    cur_ampgram_norm = bsxfun(@rdivide,cur_ampgram,interp_powfit);
    % [trig_avg,lags] = get_event_trig_avg(cur_ampgram_norm,dg_tstart_inds,backlag,forlag,[]);
    % [trig_avg2,lags] = get_event_trig_avg(cur_ampgram_norm,rp_tstart_inds,backlag,forlag,[]);
    
    %%
    % clear or_trig*
    for oo = 1:length(un_tfs)
        fprintf('Orientation %d of %d\n',oo,length(un_tfs));
        cur_ori_trials = find(all_trial_tf == un_tfs(oo));
        cur_dg_trials = all_drifting_grate_trials(ismember(all_drifting_grate_trials,cur_ori_trials));
        cur_rp_trials = all_rp_grate_trials(ismember(all_rp_grate_trials,cur_ori_trials));
        
        cur_dg_tstart_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),use_trial_start_times(cur_dg_trials)));
        cur_dg_tstop_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),all_trial_stop_times(cur_dg_trials)));
        cur_rp_tstart_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),use_trial_start_times(cur_rp_trials)));
        cur_rp_tstop_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),all_trial_stop_times(cur_rp_trials)));
        
        [or_trig_avg(ll,oo,:,:),lags] = get_event_trig_avg(cur_ampgram_norm,cur_dg_tstart_inds,backlag,forlag,[]);
        [or_trig_avg2(ll,oo,:,:),lags] = get_event_trig_avg(cur_ampgram_norm,cur_rp_tstart_inds,backlag,forlag,[]);
        
        [or_trig_lfp(ll,oo,:),lags] = get_event_trig_avg(all_V(:,ll),cur_dg_tstart_inds,backlag,forlag,[]);
        [or_trig_lfp2(ll,oo,:),lags] = get_event_trig_avg(all_V(:,ll),cur_rp_tstart_inds,backlag,forlag,[]);
    end
    
end



%%
temp = squeeze(mean(or_trig_avg));
temp2 = squeeze(mean(or_trig_avg2));
% for ll = 1:length(use_lfps)
% temp = squeeze((or_trig_avg(ll,:,:,:)));
% temp2 = squeeze((or_trig_avg2(ll,:,:,:)));
ca2 = max(temp(:))*0.75;
ca1 = min(temp(:));
for oo = 1:5
   subplot(2,5,oo)
   pcolor(lags/Fsd,wfreqs,squeeze(temp(oo,:,:))');shading flat
   caxis([ca1 ca2])
   xlim([-0.2 2])
   title(sprintf('%2f Hz',un_tfs(oo)));
    subplot(2,5,5+oo)
   pcolor(lags/Fsd,wfreqs,squeeze(temp2(oo,:,:))');shading flat
   caxis([ca1 ca2])
   xlim([-0.2 2])
    title(sprintf('%2f Hz',un_tfs(oo)));
  
end
% pause
% clf
% end

%%
sm_win = round(0.005*Fsd);
temp = squeeze(mean(or_trig_avg));
temp2 = squeeze(mean(or_trig_avg2));
for oo = 1:5
    oo
    for ww = 1:length(wfreqs)
       temp(oo,:,ww) = jmm_smooth_1d_cor(temp(oo,:,ww),sm_win);
       temp2(oo,:,ww) = jmm_smooth_1d_cor(temp2(oo,:,ww),sm_win);
    end
end
ca2 = max(temp(:))*0.95;
ca1 = min(temp(:));
for oo = 1:5
   subplot(2,5,oo)
   pcolor(lags/Fsd,wfreqs,squeeze(temp(oo,:,:))');shading flat
   caxis([ca1 ca2])
   xlim([-0.2 2])
   title(sprintf('%.2f Hz Drifting',un_tfs(oo)));
    subplot(2,5,5+oo)
   pcolor(lags/Fsd,wfreqs,squeeze(temp2(oo,:,:))');shading flat
   caxis([ca1 ca2])
   xlim([-0.2 2])
    title(sprintf('%.2f Hz Jumping',un_tfs(oo)));
  
end
% pause
% clf
% end

%%
forlag = round(Fsd*1);
backlag = round(0.2*Fsd);

un_tfs = unique(all_trial_tf);
un_tfs(1) = [];
un_tfs(end) = [];

for oo = 1:length(un_tfs)
        cur_reltimeset = 1/un_tfs(oo):1/un_tfs(oo):(4-1/un_tfs(oo));
        
        cur_start_ind_set{oo} = [];
        cur_ori_trials = find(all_trial_tf == un_tfs(oo));
        cur_rp_trials = all_rp_grate_trials(ismember(all_rp_grate_trials,cur_ori_trials));
        for tt = 1:length(cur_rp_trials)
            cur_start_times = use_trial_start_times(cur_rp_trials(tt)) + cur_reltimeset;
            bad = find(cur_start_times(2:end) > use_trial_end_times(cur_rp_trials(tt)));
            cur_start_times(bad) = [];
            cur_start_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),cur_start_times));
            cur_start_ind_set{oo} = cat(1,cur_start_ind_set{oo},cur_start_inds');
            
        end
end

for ll = 1:length(use_lfps);
    fprintf('LFP %d of %d\n',ll,length(use_lfps));
    temp = cwt(all_V(:,ll),scales,'cmor1-1');
    cur_ampgram = abs(temp)';
    avg_pow = mean(cur_ampgram);
    std_pow = std(cur_ampgram);
    
    un_wfreqs = linspace(wfreqs(1),wfreqs(end),100);
    avg_pow_int = interp1(wfreqs,avg_pow,un_wfreqs);
    
    powfun = @(a,x)(a(1)*x.^a(2));
    BETA = nlinfit(un_wfreqs,avg_pow_int,powfun,[max(avg_pow) -2]);
    powfit = un_wfreqs.^BETA(2)*BETA(1);
    interp_powfit = interp1(un_wfreqs,powfit,wfreqs);
    
    cur_ampgram_norm = bsxfun(@rdivide,cur_ampgram,interp_powfit);
%     cur_ampgram_norm = bsxfun(@minus,cur_ampgram,avg_pow);
%     cur_ampgram_norm = bsxfun(@rdivide,cur_ampgram_norm,std_pow);
    
    %%
    % clear or_trig*
    for oo = 1:length(un_tfs)
        fprintf('Orientation %d of %d\n',oo,length(un_tfs));        
        
        [or_trig_avg_jump(ll,oo,:,:),lags] = get_event_trig_avg(cur_ampgram_norm,cur_start_ind_set{oo},backlag,forlag,[]);
        
        [or_trig_lfp_jump(ll,oo,:),lags] = get_event_trig_avg(all_V(:,ll),cur_start_ind_set{oo},backlag,forlag,[]);
    end
    
end

%%
close all
% temp = squeeze(mean(or_trig_avg_jump));
for ll = 1:length(use_lfps)
temp = squeeze((or_trig_avg_jump(ll,:,:,:)));
ca2 = max(temp(:))*0.85;
ca1 = min(temp(:));
for oo = 1:4
   subplot(2,2,oo)
   pcolor(lags/Fsd,wfreqs,squeeze(temp(oo,:,:))');shading flat
   caxis([ca1 ca2])
   xlim([-0.1 0.5])  
   ylim([2 70])
           cur_reltimeset = 0:1/un_tfs(oo):(4-1/un_tfs(oo));
           yl = ylim();
           for ii = 1:length(cur_reltimeset)
               line(cur_reltimeset([ii ii]),yl,'color','k')
           end
           title(sprintf('TF: %.2f',un_tfs(oo)));
   xlabel('Time since jump (s)','fontsize',16)
   ylabel('Frequency (Hz)','fontsize',16)
end
pause
clf
end

%%

%%
params.Fs = Fsd;
params.tapers = [4 7];
movingwin = [3.5 3.5];
sMarkers_rp = [rp_tstart_inds(:) rp_tstop_inds(:)];
sMarkers_dg = [dg_tstart_inds(:) dg_tstop_inds(:)];
% params.err = [2 .05];
params.err = [0];


clear S_rp_or S_dg_or
un_tfs = unique(all_trial_tf);
un_tfs(1) = [];
for oo = 1:length(un_tfs)
    fprintf('Orientation %d of %d\n',oo,length(un_tfs));
    cur_ori_trials = find(all_trial_tf == un_tfs(oo));
    cur_dg_trials = all_drifting_grate_trials(ismember(all_drifting_grate_trials,cur_ori_trials));
    cur_rp_trials = all_rp_grate_trials(ismember(all_rp_grate_trials,cur_ori_trials));
    
    cur_dg_tstart_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),use_trial_start_times(cur_dg_trials)));
    cur_dg_tstop_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),all_trial_stop_times(cur_dg_trials)));
    cur_rp_tstart_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),use_trial_start_times(cur_rp_trials)));
    cur_rp_tstop_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),all_trial_stop_times(cur_rp_trials)));
    
    cur_sMarkers_rp = [cur_rp_tstart_inds(:) cur_rp_tstop_inds(:)];
    cur_sMarkers_dg = [cur_dg_tstart_inds(:) cur_dg_tstop_inds(:)];
    
    for ll = 1:length(use_lfps)
        
        % [S_dg, f,Serr_dg]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers );
        [S_dg_or(oo,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, cur_sMarkers_dg );
        [S_rp_or(oo,ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, cur_sMarkers_rp );
        % [S_rp, f,Serr_rp]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers );
        
    end
end

cur_ori_trials = find(all_trial_tf == 0);
cur_tstart_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),use_trial_start_times(cur_ori_trials)));
cur_tstop_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),all_trial_stop_times(cur_ori_trials)));
cur_sMarkers = [cur_tstart_inds(:) cur_tstop_inds(:)];
clear S_null
for ll = 1:length(use_lfps)
    [S_null(ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, cur_sMarkers);
end

%%
cmap = jet(length(un_tfs));
figure
subplot(1,2,1);hold on
% shadedErrorBar(f,mean(log(S_null),1),std(log(S_null),[],1)/sqrt(length(use_lfps)),{'color','k'});
plot(f,mean(log(S_null),1),'color','k','linewidth',2);

subplot(1,2,2);hold on
% shadedErrorBar(f,mean(log(S_null),1),std(log(S_null),[],1)/sqrt(length(use_lfps)),{'color','k'});
plot(f,mean(log(S_null),1),'color','k','linewidth',2);

for ii = 1:length(un_tfs)
    subplot(1,2,1)
%     shadedErrorBar(f,squeeze(mean(log(S_dg_or(ii,:,:)),2)),squeeze(std(log(S_dg_or(ii,:,:)),[],2))/sqrt(length(use_lfps)),{'color',cmap(ii,:)});
    plot(f,squeeze(mean(log(S_dg_or(ii,:,:)),2)),'color',cmap(ii,:),'linewidth',2);
    hold on
    
    subplot(1,2,2)
%     shadedErrorBar(f,squeeze(mean(log(S_rp_or(ii,:,:)),2)),squeeze(std(log(S_rp_or(ii,:,:)),[],2))/sqrt(length(use_lfps)),{'color',cmap(ii,:)});
    plot(f,squeeze(mean(log(S_rp_or(ii,:,:)),2)),'color',cmap(ii,:),'linewidth',2);
    hold on
    
end

subplot(1,2,1)
xlim([1 150]); ylim([-30 -22])
set(gca,'xscale','log');
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Log power','fontsize',16)
legend('0','1','2','4','8.3','16.7');
title('Drifting grating','fontsize',16)

subplot(1,2,2)
xlim([1 150]); ylim([-30 -22])
set(gca,'xscale','log');
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Log power','fontsize',16)
legend('0','1','2','4','8.3','16.7');
title('Random-phase grating','fontsize',16)

%%
log_diffs = log(S_dg_or) - log(S_rp_or);
cmap = jet(length(un_tfs));
figure
for ii = 1:length(un_tfs)
%     shadedErrorBar(f,squeeze(mean(log_diffs(ii,:,:),2)),squeeze(std(log_diffs(ii,:,:),[],2))/sqrt(length(use_lfps)),{'color',cmap(ii,:)});
     plot(f,squeeze(mean(log_diffs(ii,:,:),2)),'color',cmap(ii,:),'linewidth',2);
   hold on
end
legend('1','2','4','8.3','16.7');
xlim([1 150]); 
xl = xlim();
line(xl,[0 0],'color','k');
yl = ylim(); yl = max(abs(yl));  ylim([-yl yl]);
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Log power difference (drifting - random-phase)','fontsize',16)

%%
close all
cmap = jet(length(un_oris));
for ll = 1:length(use_lfps)
figure
for ii = 1:length(un_oris)
    subplot(2,1,1)
    plot(f,squeeze(log(S_dg_ori(ii,ll,:))),'color',cmap(ii,:));
    hold on
    
    subplot(2,1,2)
    plot(f,squeeze(log(S_rp_or(ii,ll,:))),'color',cmap(ii,:));
    hold on
    pause
end

pause
close all
end
%%
[Cmn_dg,~,~,Smm_dg,f] = coherencyc_unequal_length_trials( all_V, movingwin, params, sMarkers_dg );

[Cmn_rp,~,~,Smm_rp,f] = coherencyc_unequal_length_trials( all_V, movingwin, params, sMarkers_rp);

figure

shadedErrorBar(f,mean(Cmn_dg,2),std(Cmn_dg,[],2)/(length(use_lfps)));
hold on
shadedErrorBar(f,mean(Cmn_rp,2),std(Cmn_rp,[],2)/(length(use_lfps)),{'color','r'});




