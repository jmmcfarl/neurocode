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
all_trial_or = [];
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
    all_trial_durs = cat(1,all_trial_durs,trial_durs(use_trials)');
end

%%
% use_lfps = [4:4:96];
use_lfps = [4:4:96];
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
% use_trial_start_times = all_trial_start_times + 0.1;
use_trial_start_times = all_trial_start_times + 0.6;
use_trial_end_times = all_trial_end_times;

all_drifting_grate_trials = find(all_trial_nph == 0);
dg_tstart_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),use_trial_start_times(all_drifting_grate_trials)));
dg_tstop_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),use_trial_end_times(all_drifting_grate_trials)));

all_rp_grate_trials = find(all_trial_nph > 0);
rp_tstart_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),use_trial_start_times(all_rp_grate_trials)));
rp_tstop_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),use_trial_end_times(all_rp_grate_trials)));

params.Fs = Fsd;
params.tapers = [2 3];
% movingwin = [3.5 3.5];
movingwin = [0.5 0.5];
% params.err = [2 .05];
params.err = [0];

cur_ori_trials = find(all_trial_tf == 4);
cur_dg_trials = all_drifting_grate_trials(ismember(all_drifting_grate_trials,cur_ori_trials));
cur_rp_trials = all_rp_grate_trials(ismember(all_rp_grate_trials,cur_ori_trials));

cur_dg_tstart_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),use_trial_start_times(cur_dg_trials)));
cur_dg_tstop_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),all_trial_stop_times(cur_dg_trials)));
cur_rp_tstart_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),use_trial_start_times(cur_rp_trials)));
cur_rp_tstop_inds = round(interp1(all_t_lfp,1:length(all_t_lfp),all_trial_stop_times(cur_rp_trials)));

cur_sMarkers_rp = [cur_rp_tstart_inds(:) cur_rp_tstop_inds(:)];
cur_sMarkers_dg = [cur_dg_tstart_inds(:) cur_dg_tstop_inds(:)];

sdur_rp = cur_sMarkers_rp(:,2) - cur_sMarkers_rp(:,1);
sdur_dg = cur_sMarkers_dg(:,2) - cur_sMarkers_dg(:,1);
cur_sMarkers_rp(sdur_rp < round(movingwin(1)*Fsd),:) = [];
cur_sMarkers_dg(sdur_dg < round(movingwin(1)*Fsd),:) = [];

clear S_dg S_rp
for ll = 1:length(use_lfps)
    
    % [S_dg, f,Serr_dg]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers );
    [S_dg(ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, cur_sMarkers_dg );
    [S_rp(ll,:), f]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, cur_sMarkers_rp );
    % [S_rp, f,Serr_rp]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers );
    
end

figure

shadedErrorBar(f,mean(log(S_dg)),std(log(S_dg))/sqrt(length(use_lfps)));
hold on
shadedErrorBar(f,mean(log(S_rp)),std(log(S_rp))/sqrt(length(use_lfps)),{'color','r'});

%%
params.Fs = Fsd;
params.tapers = [3 5];
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
% subplot(1,2,1);hold on
% % shadedErrorBar(f,mean(log(S_null),1),std(log(S_null),[],1)/sqrt(length(use_lfps)),{'color','k'});
% plot(f,mean(log(S_null),1),'color','k','linewidth',2);

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
% set(gca,'xscale','log');
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Log power','fontsize',16)
% legend('0','1','2','4','8.3','16.7');
title('Drifting grating','fontsize',16)

subplot(1,2,2)
xlim([1 150]); ylim([-30 -22])
% set(gca,'xscale','log');
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Log power','fontsize',16)
% legend('0','1','2','4','8.3','16.7');
title('Random-phase grating','fontsize',16)

cd ~/Analysis/bruce/
save grating_powspectra_tf un_tfs S_null S_dg_or S_rp_or f 

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




