clear all
% close all

addpath('~/James_scripts/general_functions/');

ExptNum = 239;
cd(['~/Data/bruce/M' num2str(ExptNum)]);
load ./random_bar_eyedata_ftime.mat

load(['lemM' num2str(ExptNum) 'Expts.mat']);
% 
% if ExptNum == 235
%     bar_expts(bar_expts==51) = []; %this has different set of stim positions
% end
% if ExptNum == 239
%     bar_expts(bar_expts==40) = []; %this has different set of stim positions
% end

%%
axis_or = Expts{bar_expts(1)}.Stimvals.or*pi/180; %in radians
new_dt = .005;
dsfrac = 2;

%%
load ./CellList.mat
good_sus = find(all(CellList(bar_expts,:,1) > 0));
use_sus = 1:24;
%%
flen = 15;
full_spkbinned = [];
full_exptvec = [];
full_exptvec_new = [];
full_trialvec_new = [];
full_taxis = [];
full_taxis_new = [];
full_trel = [];
for ee = 1:length(bar_expts);
    fprintf('Processing expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('Expt%dClusterTimes.mat',bar_expts(ee));
    load(fname);
    
    trial_durs{ee} = [Expts{bar_expts(ee)}.Trials(:).dur];
    n_trials(ee) = length(Expts{bar_expts(ee)}.Trials);
    all_t_axis = [];
    all_t_axis_new = [];
    all_trel = [];
    all_used_inds_new = [];
    all_expt_vec = [];
    all_trial_vec = [];
    all_binned_spikes = [];
    all_used_inds = [];
    for tt = 1:n_trials(ee)
        
        cur_t_axis = [Expts{bar_expts(ee)}.Trials(tt).Start]/1e4;
        cur_t_edges = [cur_t_axis; Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4];
        cur_t_edges_new = cur_t_edges(1):new_dt:cur_t_edges(end);
        cur_t_axis_new = cur_t_edges_new(1:end-1);
        
        cur_binned_spks = nan(length(cur_t_axis_new),length(use_sus));
        for cc = 1:length(use_sus)
            cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges_new);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        cur_used_inds_new = ones(length(cur_t_axis_new),1);
        cur_used_inds_new(1:flen) = 0;
        
        all_used_inds_new = [all_used_inds_new; cur_used_inds_new(:)];
        all_t_axis = [all_t_axis; cur_t_axis(:)];
        all_t_axis_new = [all_t_axis_new; cur_t_axis_new(:)];
        all_trel = [all_trel; cur_t_axis_new(:) - cur_t_axis(1)];
        all_binned_spikes = [all_binned_spikes; cur_binned_spks];
        all_expt_vec = [all_expt_vec; ones(length(cur_t_axis_new),1)*bar_expts(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_t_axis_new),1)*tt];
    end
    
    % %
    % cur_blink_start_times = all_eye_ts{ee}(all_blink_startinds{ee});
    %     cur_blink_stop_times = all_eye_ts{ee}(all_blink_stopinds{ee});
    %     cur_blink_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_start_times));
    %     cur_blink_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_stop_times));
    %     cur_poss_blinks = find(~isnan(cur_blink_start_inds) & ~isnan(cur_blink_stop_inds));
    %
    %     all_blink_inds = zeros(size(all_t_axis));
    %     for i = 1:length(cur_poss_blinks)
    %         all_blink_inds(cur_blink_start_inds(cur_poss_blinks(i)):cur_blink_stop_inds(cur_poss_blinks(i))) = 1;
    %     end
    %
    %     cur_blink_start_inds = round(interp1(all_t_axis_new,1:length(all_t_axis_new),cur_blink_start_times));
    %     cur_blink_stop_inds = round(interp1(all_t_axis_new,1:length(all_t_axis_new),cur_blink_stop_times));
    %     cur_poss_blinks = find(~isnan(cur_blink_start_inds) & ~isnan(cur_blink_stop_inds));
    %     all_blink_inds_new = zeros(size(all_t_axis_new));
    %     for i = 1:length(cur_poss_blinks)
    %         all_blink_inds_new(cur_blink_start_inds(cur_poss_blinks(i)):cur_blink_stop_inds(cur_poss_blinks(i))) = 1;
    %     end
    %     all_used_inds_new(all_blink_inds_new == 1) = 0;
    
    full_spkbinned = [full_spkbinned; all_binned_spikes];
    full_exptvec = [full_exptvec; ones(length(all_used_inds),1)*ee];
    full_exptvec_new = [full_exptvec_new; ones(length(all_used_inds_new),1)*ee];
    full_taxis = [full_taxis; all_t_axis];
    full_taxis_new = [full_taxis_new; all_t_axis_new];
    full_trel = [full_trel; all_trel];
    full_trialvec_new = [full_trialvec_new; all_trial_vec];
end

%% WITHIN BLOCK SPK RATE NORMALIZATION
full_spkbinned_norm = full_spkbinned;
for bb = 1:length(bar_expts)
    cur_set = find(full_exptvec_new==bb);
    full_spkbinned_norm(cur_set,:) = bsxfun(@rdivide,full_spkbinned(cur_set,:),mean(full_spkbinned(cur_set,:)));
end



%%
all_sac_inds = [];
all_microsac_inds = [];
all_firstsac_inds = [];
all_secondsac_inds = [];
all_trial_start = [];
for ee = 1:length(bar_expts)
    expt_trial_starts = [Expts{bar_expts(ee)}.Trials(:).TrialStart]/1e4;
    expt_trial_end = [Expts{bar_expts(ee)}.Trials(:).TrueEnd]/1e4;
    
    sac_start_times = all_eye_ts{ee}(all_sac_startinds{ee});
    sac_stop_times = all_eye_ts{ee}(all_sac_stopinds{ee});
    sac_amps = all_sac_peakamps{ee};
    
    cur_sac_set = find(all_expt_num==ee);
    if length(cur_sac_set) ~= length(sac_start_times)
        error('saccade mismatch')
    end
    intrial_sac_set = find(~isnan(all_trial_num(cur_sac_set)));
    
    infirst_sac_set = find(ismember(cur_sac_set,all_isfirst));
    insecond_sac_set = find(ismember(cur_sac_set,all_issecond));
    big_sac_set = [infirst_sac_set; insecond_sac_set];
    micro_sac_set = find(sac_amps' < 60 & ~isnan(all_trial_num(cur_sac_set)));
    micro_sac_set(ismember(micro_sac_set,big_sac_set)) = [];
    
    sac_times = sac_start_times(intrial_sac_set);
    first_sac_times = sac_start_times(infirst_sac_set);
    second_sac_times = sac_start_times(insecond_sac_set);
    micro_sac_times = sac_start_times(micro_sac_set);
    
    cur_e_set = find(full_exptvec_new == ee);
    
    sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),sac_times));
    sac_inds(isnan(sac_inds)) = [];
    micro_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),micro_sac_times));
    micro_sac_inds(isnan(micro_sac_inds)) = [];
    first_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),first_sac_times));
    first_sac_inds(isnan(first_sac_inds)) = [];
    second_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),second_sac_times));
    second_sac_inds(isnan(second_sac_inds)) = [];
    
    trial_start_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),expt_trial_starts));
    trial_start_inds(isnan(trial_start_inds)) = [];
    
    all_sac_inds = [all_sac_inds; cur_e_set(sac_inds(:))];
    all_microsac_inds = [all_microsac_inds; cur_e_set(micro_sac_inds(:))];
    all_firstsac_inds = [all_firstsac_inds; cur_e_set(first_sac_inds(:))];
    all_secondsac_inds = [all_secondsac_inds; cur_e_set(second_sac_inds(:))];
    all_trial_start = [all_trial_start; cur_e_set(trial_start_inds(:))];
end
all_bigsac_inds = sort([all_firstsac_inds; all_secondsac_inds]);
fprintf('Identified:\n %d msacs\n %d fsacs\n %d ssacs\n',length(all_microsac_inds),length(all_firstsac_inds),length(all_secondsac_inds));

all_firststim_inds = 1+find(full_trel(1:end-1) < 0.7 & full_trel(2:end) >= 0.7);
all_secondstim_inds = 1+find(full_trel(1:end-1) < 1.4 & full_trel(2:end) >= 1.4);
%%
clear *_avg_rate

ovmean_rate = mean(full_spkbinned)/new_dt;

% sm_sig = ceil(0.015/new_dt);
sm_sig = ceil(0.015/new_dt);
sm_spkbinned = zeros(size(full_spkbinned));
for cc =1:24
    sm_spkbinned(:,cc) = jmm_smooth_1d_cor(full_spkbinned_norm(:,cc),sm_sig);
end

%%
% for cc = 1:24
%     [sac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(:,cc),all_bigsac_inds,round(0.4/new_dt),round(0.5/new_dt));
%     [fsac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(:,cc),all_firstsac_inds,round(0.4/new_dt),round(0.5/new_dt));
%     [ssac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(:,cc),all_secondsac_inds,round(0.4/new_dt),round(0.5/new_dt));
%     [msac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(:,cc),all_microsac_inds,round(0.4/new_dt),round(0.5/new_dt));
%     [trial_trig_avg_rate(cc,:),cur_tlags] = get_event_trig_avg(sm_spkbinned(:,cc),all_trial_start,round(0/new_dt),round(2/new_dt));
%     [fstim_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(:,cc),all_firststim_inds,round(0.4/new_dt),round(0.5/new_dt));
%     [sstim_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(:,cc),all_secondstim_inds,round(0.4/new_dt),round(0.5/new_dt));
% end
% 
% %%
% close all
% figure
% imagesc(cur_lags*new_dt,1:24,sac_trig_avg_rate);colorbar;caxis([0.5 1.5]);
% 
% figure
% imagesc(cur_lags*new_dt,1:24,msac_trig_avg_rate);colorbar;caxis([0.75 1.25]);
% 
% figure
% subplot(2,1,1)
% imagesc(cur_lags*new_dt,1:24,fsac_trig_avg_rate);colorbar;caxis([0.5 1.5]);
% subplot(2,1,2)
% imagesc(cur_lags*new_dt,1:24,ssac_trig_avg_rate);colorbar;caxis([0.5 1.5]);
% 
% figure
% imagesc(cur_tlags*new_dt,1:24,trial_trig_avg_rate);colorbar;caxis([0.5 1.5]);
% 
% figure
% subplot(2,1,1)
% imagesc(cur_lags*new_dt,1:24,fstim_trig_avg_rate);colorbar;caxis([0.5 1.5]);
% subplot(2,1,2)
% imagesc(cur_lags*new_dt,1:24,sstim_trig_avg_rate);colorbar;caxis([0.5 1.5]);
% 
% %% RANDOMLY SPLIT TRIALS
% [c,ia,ic] = unique([full_exptvec_new full_trialvec_new],'rows');
% n_trials = length(ia);
% 
% rp = randperm(n_trials);
% rand_trial_vec = rp(ic);
% [~,ind_shuff] = sort(rand_trial_vec);
% 
% xv_frac = 0.2;
% n_xv_trials = floor(n_trials*xv_frac);
% rp_tset = randperm(n_trials);
% for ii = 1:5
%     cur = (ii-1)*n_xv_trials + (1:n_xv_trials);
%     cur_tset = rp_tset(cur);
%     xv_inds{ii} = find(ismember(ic,cur_tset));
% end

%%
% clear *_avg_rate
% for ii = 1:5
%     cur_bigsac_inds = find(ismember(xv_inds{ii},all_bigsac_inds));
%     cur_firstsac_inds = find(ismember(xv_inds{ii},all_firstsac_inds));
%     cur_secondsac_inds = find(ismember(xv_inds{ii},all_secondsac_inds));
%     cur_microsac_inds = find(ismember(xv_inds{ii},all_microsac_inds));
%     for cc = 1:24
%         [sac_trig_avg_rate(ii,cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(xv_inds{ii},cc),cur_bigsac_inds,round(0.4/new_dt),round(0.5/new_dt));
%         [fsac_trig_avg_rate(ii,cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(xv_inds{ii},cc),cur_firstsac_inds,round(0.4/new_dt),round(0.5/new_dt));
%         [ssac_trig_avg_rate(ii,cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(xv_inds{ii},cc),cur_secondsac_inds,round(0.4/new_dt),round(0.5/new_dt));
%         [msac_trig_avg_rate(ii,cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(xv_inds{ii},cc),cur_microsac_inds,round(0.4/new_dt),round(0.5/new_dt));
%     end
% end
%
% sac_trig_avg = squeeze(mean(sac_trig_avg_rate));
% sac_trig_sem = squeeze(std(sac_trig_avg_rate))/sqrt(5);
% msac_trig_avg = squeeze(mean(msac_trig_avg_rate));
% msac_trig_sem = squeeze(std(msac_trig_avg_rate))/sqrt(5);
%
% %%
% close all
% for cc = 1:24
%     fprintf('Cell %d. Avg rate: %.2f\n',cc,mean(full_spkbinned(:,cc))/new_dt);
%     shadedErrorBar(cur_lags*new_dt,sac_trig_avg(cc,:),sac_trig_sem(cc,:),{'color','k'});
%     hold on
%     shadedErrorBar(cur_lags*new_dt,msac_trig_avg(cc,:),msac_trig_sem(cc,:),{'color','r'});
%
%     pause
%     clf
% end

%%
use_expts = ones(24,length(bar_expts));

%for M235
if ExptNum == 235
use_expts(1,[1 2 5 8]) = 0;
use_expts(6,9) = 0;
use_expts(7,6) = 0;
end

%for M239
if ExptNum == 239
use_expts(14,7) = 0;
use_expts(15,[6 7]) = 0;
use_expts(16,2) = 0;
end

%%
clear *_avg_rate
for ii = 1:length(bar_expts)
    cur_inds = find(full_exptvec_new==ii);
    cur_bigsac_inds = find(ismember(cur_inds,all_bigsac_inds));
    cur_firstsac_inds = find(ismember(cur_inds,all_firstsac_inds));
    cur_secondsac_inds = find(ismember(cur_inds,all_secondsac_inds));
    cur_microsac_inds = find(ismember(cur_inds,all_microsac_inds));
    for cc = 1:24
        [sac_trig_avg_rate(ii,cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(cur_inds,cc),cur_bigsac_inds,round(0.4/new_dt),round(0.5/new_dt));
        [fsac_trig_avg_rate(ii,cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(cur_inds,cc),cur_firstsac_inds,round(0.4/new_dt),round(0.5/new_dt));
        [ssac_trig_avg_rate(ii,cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(cur_inds,cc),cur_secondsac_inds,round(0.4/new_dt),round(0.5/new_dt));
        [msac_trig_avg_rate(ii,cc,:),cur_lags] = get_event_trig_avg(sm_spkbinned(cur_inds,cc),cur_microsac_inds,round(0.4/new_dt),round(0.5/new_dt));
    end
end
for cc = 1:24
    sac_trig_avg(cc,:) = squeeze(mean(sac_trig_avg_rate(use_expts(cc,:)==1,cc,:)));
    sac_trig_sem(cc,:) = squeeze(std(sac_trig_avg_rate(use_expts(cc,:)==1,cc,:)))/sqrt(sum(use_expts(cc,:)==1));
    msac_trig_avg(cc,:) = squeeze(mean(msac_trig_avg_rate(use_expts(cc,:)==1,cc,:)));
    msac_trig_sem(cc,:) = squeeze(std(msac_trig_avg_rate(use_expts(cc,:)==1,cc,:)))/sqrt(sum(use_expts(cc,:)==1));
    fsac_trig_avg(cc,:) = squeeze(mean(fsac_trig_avg_rate(use_expts(cc,:)==1,cc,:)));
    fsac_trig_sem(cc,:) = squeeze(std(fsac_trig_avg_rate(use_expts(cc,:)==1,cc,:)))/sqrt(sum(use_expts(cc,:)==1));
    ssac_trig_avg(cc,:) = squeeze(mean(ssac_trig_avg_rate(use_expts(cc,:)==1,cc,:)));
    ssac_trig_sem(cc,:) = squeeze(std(ssac_trig_avg_rate(use_expts(cc,:)==1,cc,:)))/sqrt(sum(use_expts(cc,:)==1));
end
%%
% cd(sprintf('~/Analysis/bruce/M%d/sac_trig_avgs',ExptNum));
% to_print = 1;
% 
% close all
% for cc = 1:24
%     fprintf('Cell %d. Avg rate: %.2f\n',cc,mean(full_spkbinned(:,cc))/new_dt);
%     shadedErrorBar(cur_lags*new_dt,sac_trig_avg(cc,:),sac_trig_sem(cc,:),{'color','k'});
%     hold on
%     shadedErrorBar(cur_lags*new_dt,msac_trig_avg(cc,:),msac_trig_sem(cc,:),{'color','r'});
%     plot(cur_lags*new_dt,fsac_trig_avg(cc,:),'b')
%     plot(cur_lags*new_dt,ssac_trig_avg(cc,:),'g')
%     
%     if to_print == 1
%     fname = sprintf('Unit%d_sm%d',cc,sm_sig);
%     print('-dpng',fname);
%       close all  
%     else
%         pause
%         clf
%     end
% end

%%
cd(sprintf('~/Analysis/bruce/M%d',ExptNum));
sname = 'sacmod_data.mat';
save(sname,'*sac_trig*','cur_lags','new_dt','ovmean_rate')
%%
close all
use_units = 1:24;
shadedErrorBar(cur_lags*new_dt,mean(sac_trig_avg(use_units,:)),std(sac_trig_avg(use_units,:))/sqrt(length(use_units)),{'color','r'});
hold on;
shadedErrorBar(cur_lags*new_dt,mean(msac_trig_avg(use_units,:)),std(msac_trig_avg(use_units,:))/sqrt(length(use_units)),{'color','b'});
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');

%%
close all
use_units = 1:24;
min_rate =5;

cd('~/Analysis/bruce/M232');
load sacmod_data
% cd('~/Data/bruce/M232');
% load ./CellList.mat
% load ./random_bar_eyedata_ftime.mat bar_expts
% good_sus = find(all(CellList(bar_expts,:,1) > 0));
% cuse_units = good_sus;
cuse_units = use_units(ovmean_rate(use_units) > min_rate);
% figure
shadedErrorBar(cur_lags*new_dt,mean(sac_trig_avg(cuse_units,:)),std(sac_trig_avg(cuse_units,:))/sqrt(length(cuse_units)),{'color','r'});
hold on;
shadedErrorBar(cur_lags*new_dt,mean(msac_trig_avg(cuse_units,:)),std(msac_trig_avg(cuse_units,:))/sqrt(length(cuse_units)),{'color','b'});
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
title('M232');

cd('~/Analysis/bruce/M235');
load sacmod_data
% cd('~/Data/bruce/M235');
% load ./CellList.mat
% load ./random_bar_eyedata_ftime.mat bar_expts
good_sus = find(all(CellList(bar_expts,:,1) > 0));
% cuse_units = good_sus;
cuse_units = use_units(ovmean_rate(use_units) > min_rate);
figure
shadedErrorBar(cur_lags*new_dt,mean(sac_trig_avg(cuse_units,:)),std(sac_trig_avg(cuse_units,:))/sqrt(length(cuse_units)),{'color','r'});
hold on;
shadedErrorBar(cur_lags*new_dt,mean(msac_trig_avg(cuse_units,:)),std(msac_trig_avg(cuse_units,:))/sqrt(length(cuse_units)),{'color','b'});
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
title('M235');


cd('~/Analysis/bruce/M239');
load sacmod_data
% cd('~/Data/bruce/M235');
% load ./CellList.mat
% load ./random_bar_eyedata_ftime.mat bar_expts
% good_sus = find(all(CellList(bar_expts,:,1) > 0));
% cuse_units = good_sus;
cuse_units = use_units(ovmean_rate(use_units) > min_rate);
figure
shadedErrorBar(cur_lags*new_dt,mean(sac_trig_avg(cuse_units,:)),std(sac_trig_avg(cuse_units,:))/sqrt(length(cuse_units)),{'color','r'});
hold on;
shadedErrorBar(cur_lags*new_dt,mean(msac_trig_avg(cuse_units,:)),std(msac_trig_avg(cuse_units,:))/sqrt(length(cuse_units)),{'color','b'});
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
title('M239');


use_units = 1:24;
min_rate =5;

cd('~/Analysis/bruce/M232');
load sacmod_data
% cd('~/Data/bruce/M232');
% load ./CellList.mat
% load ./random_bar_eyedata_ftime.mat bar_expts
% good_sus = find(all(CellList(bar_expts,:,1) > 0));
% cuse_units = good_sus;
cuse_units = use_units(ovmean_rate(use_units) > min_rate);
% figure
shadedErrorBar(cur_lags*new_dt,mean(sac_trig_avg(cuse_units,:)),std(sac_trig_avg(cuse_units,:))/sqrt(length(cuse_units)),{'color','r'});
hold on;
shadedErrorBar(cur_lags*new_dt,mean(msac_trig_avg(cuse_units,:)),std(msac_trig_avg(cuse_units,:))/sqrt(length(cuse_units)),{'color','b'});
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
title('M232');

%%
use_units = 1:14;
min_rate = 2;
cur_savg = [];
cur_mavg = [];

cd('~/Analysis/bruce/M232');
load sacmod_data
cuse_units = use_units(ovmean_rate(use_units) > min_rate);
cur_savg = cat(1,cur_savg,sac_trig_avg(cuse_units,:));
cur_mavg = cat(1,cur_mavg,msac_trig_avg(cuse_units,:));

cd('~/Analysis/bruce/M239');
load sacmod_data
cuse_units = use_units(ovmean_rate(use_units) > min_rate);
cur_savg = cat(1,cur_savg,sac_trig_avg(cuse_units,:));
cur_mavg = cat(1,cur_mavg,msac_trig_avg(cuse_units,:));


figure
shadedErrorBar(cur_lags*new_dt,mean(cur_savg),std(cur_savg)/sqrt(size(cur_savg,1)),{'color','r'});
hold on;
shadedErrorBar(cur_lags*new_dt,mean(cur_mavg),std(cur_mavg)/sqrt(size(cur_mavg,1)),{'color','b'});
xl = xlim();yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlim([-0.2 0.4])
