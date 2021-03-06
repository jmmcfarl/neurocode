clear all
close all

ExptNum = 232;
dat_dir = ['~/Data/bruce/M' num2str(ExptNum)];
cd(dat_dir);
anal_dir = ['~/Analysis/bruce/M' num2str(ExptNum)];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

load ./random_bar_eyedata_ftime.mat

load(['lemM' num2str(ExptNum) 'Expts.mat']);
%%
rf_cent = [4.73 -2.2];
axis_or = 50*pi/180; %in radians

Fs = 1000;
dsf = 3;
Fsd = Fs/dsf;

niqf = Fs/2;
% [bb,aa] = butter(2,[1]/niqf,'high');
[bb,aa] = butter(2,[1 120]/niqf);
wname = 'cmor1-1';
min_freq = 1.5;
max_freq = 100;
n_freqs = 35;
max_scale = centfrq(wname)/min_freq*Fsd;
min_scale = centfrq(wname)/max_freq*Fsd;
scales = logspace(log10(min_scale),log10(max_scale),n_freqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

% dsfrac = 3;
% new_dt = .01/dsfrac;

%%
load ./CellList.mat
good_sus = find(all(CellList(bar_expts,:,1) > 0));
use_sus = 1:24;
%%
% flen = 30;
full_Xmat = [];
full_spkbinned = [];
full_exptvec = [];
full_exptvec_new = [];
full_trialvec_new = [];
full_taxis = [];
full_taxis_new = [];
full_old_t_inds = [];
full_bar_pos = [];
full_t_since_tstart = [];
full_used_inds = [];
full_used_inds_new = [];
for ee = 1:length(bar_expts);
    fprintf('Processing expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('Expt%dClusterTimes.mat',bar_expts(ee));
    load(fname);
    
    trial_durs{ee} = [Expts{bar_expts(ee)}.Trials(:).dur];
    n_trials(ee) = length(Expts{bar_expts(ee)}.Trials);
    all_t_axis = [];
    all_used_inds = [];
    all_t_axis_new = [];
    all_old_t_inds = [];
    all_used_inds_new = [];
    all_expt_vec = [];
    all_trial_vec = [];
    all_bar_Op = [];
    all_bar_Xmat = [];
    all_binned_spikes = [];
    all_used_inds = [];
    
    all_t_since_tstart = [];
    for tt = 1:n_trials(ee)
        
        cur_t_axis = [Expts{bar_expts(ee)}.Trials(tt).Start]/1e4;        
        cur_t_edges = [cur_t_axis; Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4];
        cur_t_edges_new = cur_t_edges(1):1/Fsd:cur_t_edges(end);
        cur_t_axis_new = cur_t_edges_new(1:end-1);
        
        cur_binned_spks = nan(length(cur_t_axis_new),length(use_sus));
        for cc = 1:length(use_sus)
            cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges_new);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
                
%         cur_bar_mat = zeros(length(cur_t_axis_new),length(un_bar_pos));
%         for b = 1:length(un_bar_pos)
%             cur_set = find(cur_bar_Op_new==un_bar_pos(b));
%             cur_bar_mat(cur_set,b) = 1;
%         end
%         bar_Xmat = makeStimRows(cur_bar_mat,flen);
        
%         cur_used_inds = ones(length(cur_t_axis),1);
%         cur_used_inds(1:flen) = 0;
        cur_used_inds_new = ones(length(cur_t_axis_new),1);
%         cur_used_inds_new(1:flen) = 0;
        cur_used_inds_new(1:round(0.15*Fsd)) = 0;
        
        cur_t_since_tstart = cur_t_axis_new - cur_t_axis(1);
        
        cur_used_inds_new(end-round(0.1*Fsd):end) = 0; %reduce phase estimation edge artifacts
        
        all_used_inds_new = [all_used_inds_new; cur_used_inds_new(:)];
%         all_t_axis = [all_t_axis; cur_t_axis(:)];
        all_t_since_tstart = [all_t_since_tstart; cur_t_since_tstart(:)];
%         all_used_inds = [all_used_inds; cur_used_inds(:)];
        all_t_axis_new = [all_t_axis_new; cur_t_axis_new(:)];
        all_binned_spikes = [all_binned_spikes; cur_binned_spks];
        all_expt_vec = [all_expt_vec; ones(length(cur_t_axis_new),1)*bar_expts(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_t_axis_new),1)*tt];
%         all_bar_Op = [all_bar_Op; cur_bar_Op_new(:)];
%         all_bar_Xmat = [all_bar_Xmat; bar_Xmat];
    end
    
    cur_blink_start_times = all_eye_ts{ee}(all_blink_startinds{ee});
    cur_blink_stop_times = all_eye_ts{ee}(all_blink_stopinds{ee});
%     cur_blink_start_inds = round(interp1(all_t_axis_new,1:length(all_t_axis),cur_blink_start_times));
%     cur_blink_stop_inds = round(interp1(all_t_axis_new,1:length(all_t_axis),cur_blink_stop_times));
%     cur_poss_blinks = find(~isnan(cur_blink_start_inds) & ~isnan(cur_blink_stop_inds));
%     
%     all_blink_inds = zeros(size(all_t_axis_new));
%     for i = 1:length(cur_poss_blinks)
%         all_blink_inds(cur_blink_start_inds(cur_poss_blinks(i)):cur_blink_stop_inds(cur_poss_blinks(i))) = 1;
%     end
%     
    cur_blink_start_inds = round(interp1(all_t_axis_new,1:length(all_t_axis_new),cur_blink_start_times));
    cur_blink_stop_inds = round(interp1(all_t_axis_new,1:length(all_t_axis_new),cur_blink_stop_times));
    cur_poss_blinks = find(~isnan(cur_blink_start_inds) & ~isnan(cur_blink_stop_inds));
    all_blink_inds_new = zeros(size(all_t_axis_new));
    for i = 1:length(cur_poss_blinks)
        all_blink_inds_new(cur_blink_start_inds(cur_poss_blinks(i)):cur_blink_stop_inds(cur_poss_blinks(i))) = 1;
    end
    
%     all_used_inds(all_blink_inds == 1) = 0;
%     all_used_inds = logical(all_used_inds);
    all_used_inds_new(all_blink_inds_new == 1) = 0;
    all_used_inds_new = logical(all_used_inds_new);
    
%     full_Xmat = [full_Xmat; all_bar_Xmat];
    full_spkbinned = [full_spkbinned; all_binned_spikes];
    full_exptvec = [full_exptvec; ones(length(all_used_inds),1)*ee];
    full_exptvec_new = [full_exptvec_new; ones(length(all_used_inds_new),1)*ee];
    full_taxis = [full_taxis; all_t_axis];
    full_taxis_new = [full_taxis_new; all_t_axis_new];
    full_t_since_tstart = [full_t_since_tstart; all_t_since_tstart];
    full_trialvec_new = [full_trialvec_new; all_trial_vec];
    full_bar_pos = [full_bar_pos; all_bar_Op];
    
    full_used_inds_new = [full_used_inds_new; all_used_inds_new];
    full_used_inds = [full_used_inds; all_used_inds];
    
end

%%
full_lfps = [];
full_phasegrams = [];
full_ampgrams = [];
for ee = 1:length(bar_expts);
    fprintf('Expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('lemM%dA.%d.lfp.mat',ExptNum,bar_expts(ee));
    load(fname);
    
    Fs = 1/LFP.Header.CRsamplerate;
    
    n_trials(ee) = length(LFP.Trials);
    %     lfp_trial_starts = [LFP.Trials(:).Start]/1e4;
    lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
    lfp_trial_ends = [LFP.Trials(:).End]/1e4;
    expt_lfp_t_axis = [];
    expt_lfps = [];
    for tt = 1:n_trials(ee)
        
        cur_npts = size(LFP.Trials(tt).LFP,1);
        cur_t_end(tt) = lfp_trial_starts(tt)+(cur_npts-1)/Fs;
        cur_t_axis = lfp_trial_starts(tt):1/Fs:cur_t_end(tt);
        
        if ~isempty(expt_lfp_t_axis)
            cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
        else
            cur_sp = 1;
        end
        cur_t_axis = cur_t_axis(cur_sp:end);
        cur_LFP = [LFP.Trials(tt).LFP];
        cur_LFP = filtfilt(bb,aa,cur_LFP); %high-pass filter
        
        cur_LFP = cur_LFP(cur_sp:end,:);
        
        expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
        expt_lfps = [expt_lfps; cur_LFP];
    end
    
    cur_set = find(full_exptvec_new==ee);
    interp_lfps = interp1(expt_lfp_t_axis,expt_lfps,full_taxis_new(cur_set));
    
    %     %compute CSD
    %     vars.Fs = Fs;
    %     vars.BrainBound = 1;
    %     vars.ChanSep = 0.05;
    %      vars.diam = 2; %0.5
    %     CSD = PettersenCSD(expt_lfps','spline',vars)';
    %     expt_lfps = CSD;
    %
    
    cur_phasegram = nan(length(cur_set),n_freqs,24);
    cur_ampgram = nan(length(cur_set),n_freqs,24);
    for ll = 1:24
        temp = cwt(interp_lfps(:,ll),scales,'cmor1-1');
        cur_phasegram(:,:,ll) = angle(temp)';
        cur_ampgram(:,:,ll) = abs(temp)';
    end    
    
    full_lfps = [full_lfps; interp_lfps];
    full_phasegrams = cat(1,full_phasegrams, cur_phasegram);
    full_ampgrams = cat(1,full_ampgrams, cur_ampgram);
end

%%
full_ampgrams_norm = full_ampgrams;
full_ampgrams_norm(full_used_inds_new==0,:,:) = nan;
full_ampgrams_norm = bsxfun(@minus,full_ampgrams_norm,nanmean(full_ampgrams_norm));
full_ampgrams_norm = bsxfun(@rdivide,full_ampgrams_norm,nanstd(full_ampgrams_norm));

full_phasegrams_used = full_phasegrams;
full_phasegrams_used(full_used_inds_new==0,:,:) = nan;

full_lfps_norm = full_lfps;
full_lfps_norm(full_used_inds_new==0,:) = nan;
full_lfps_norm = bsxfun(@minus,full_lfps_norm,nanmean(full_lfps_norm));
full_lfps_norm = bsxfun(@rdivide,full_lfps_norm,nanstd(full_lfps_norm));

%%
all_sac_inds = [];
all_microsac_inds = [];
all_firstsac_inds = [];
all_secondsac_inds = [];
% full_eye_ts = [];
full_eye_speed = [];
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
    
        interp_eye_speed = interp1(all_eye_ts{ee},all_eye_speed{ee},full_taxis_new(cur_e_set));
    full_eye_speed = [full_eye_speed; interp_eye_speed];
        
    sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),sac_times));
    sac_inds(isnan(sac_inds)) = [];
    micro_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),micro_sac_times));
    micro_sac_inds(isnan(micro_sac_inds)) = [];
    first_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),first_sac_times));
    first_sac_inds(isnan(first_sac_inds)) = [];
    second_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),second_sac_times));
    second_sac_inds(isnan(second_sac_inds)) = [];
    
    all_sac_inds = [all_sac_inds; cur_e_set(sac_inds(:))];
    all_microsac_inds = [all_microsac_inds; cur_e_set(micro_sac_inds(:))];
    all_firstsac_inds = [all_firstsac_inds; cur_e_set(first_sac_inds(:))];
    all_secondsac_inds = [all_secondsac_inds; cur_e_set(second_sac_inds(:))];
end
all_bigsac_inds = sort([all_firstsac_inds; all_secondsac_inds]);


%%
backlag = round(Fsd*0.4);
forwardlag = round(Fsd*0.4);

[first_trig_specgram,lags] = get_event_trig_avg(full_ampgrams_norm,all_firstsac_inds,backlag,forwardlag);
[second_trig_specgram,lags] = get_event_trig_avg(full_ampgrams_norm,all_secondsac_inds,backlag,forwardlag);
[bsac_trig_specgram,lags] = get_event_trig_avg(full_ampgrams_norm,all_bigsac_inds,backlag,forwardlag);
[micro_trig_specgram,lags] = get_event_trig_avg(full_ampgrams_norm,all_microsac_inds,backlag,forwardlag);

[second_trig_phaselock,lags] = get_event_trig_avg(exp(full_phasegrams_used*(1i)),all_secondsac_inds,backlag,forwardlag);
second_trig_phaselock = abs(second_trig_phaselock);
[first_trig_phaselock,lags] = get_event_trig_avg(exp(full_phasegrams_used*(1i)),all_firstsac_inds,backlag,forwardlag);
first_trig_phaselock = abs(first_trig_phaselock);
[bsac_trig_phaselock,lags] = get_event_trig_avg(exp(full_phasegrams_used*(1i)),all_bigsac_inds,backlag,forwardlag);
bsac_trig_phaselock = abs(bsac_trig_phaselock);
[micro_trig_phaselock,lags] = get_event_trig_avg(exp(full_phasegrams_used*(1i)),all_microsac_inds,backlag,forwardlag);
micro_trig_phaselock = abs(micro_trig_phaselock);

[micro_trig_avgs,lags] = get_event_trig_avg(full_lfps_norm,all_microsac_inds,backlag,forwardlag);
[bsac_trig_avgs,lags] = get_event_trig_avg(full_lfps_norm,all_bigsac_inds,backlag,forwardlag);
[first_trig_avgs,lags] = get_event_trig_avg(full_lfps_norm,all_firstsac_inds,backlag,forwardlag);
[second_trig_avgs,lags] = get_event_trig_avg(full_lfps_norm,all_secondsac_inds,backlag,forwardlag);
    

%%
sname = [anal_dir '/sac_lfpmod.mat'];
save(sname,'lags','Fsd','micro*','bsac*','first*','second*');


%%
close all
xax = [-0.2 0.4];
for cc = 1:24
plot(lags/Fsd,squeeze(bsac_trig_avgs(:,cc)))
hold on
plot(lags/Fsd,squeeze(micro_trig_avgs(:,cc)),'r')
xlim(xax)
pause
clf

end

%%
close all
xax = [-0.2 0.4];
for cc = 1:24
    subplot(2,2,1)
    pcolor(lags/Fsd,wfreqs,squeeze(bsac_trig_specgram(:,:,cc))');shading flat;colorbar
    xlim(xax)
    ca = caxis(); ca = max(abs(ca)); caxis([-ca ca]);
%     caxis([-0.15 0.15])
    subplot(2,2,3)
    pcolor(lags/Fsd,wfreqs,squeeze(bsac_trig_phaselock(:,:,cc))');shading flat;colorbar
    ylim([2 40])
    xlim(xax)
%     set(gca,'yscale','log')
    subplot(2,2,2)
    pcolor(lags/Fsd,wfreqs,squeeze(micro_trig_specgram(:,:,cc))');shading flat;colorbar
    xlim(xax)
    ca = caxis(); ca = max(abs(ca)); caxis([-ca ca]);
%     caxis([-0.15 0.15])
    subplot(2,2,4)
    pcolor(lags/Fsd,wfreqs,squeeze(micro_trig_phaselock(:,:,cc))');shading flat;colorbar
    ylim([2 40])
    xlim(xax)
%     set(gca,'yscale','log')
pause
clf
end

% for cc = 1:24
%     subplot(2,2,1)
%     pcolor(lags/Fsd,wfreqs,squeeze(first_trig_specgram(:,:,cc))');shading flat;colorbar
%     xlim(xax)
%     caxis([-0.15 0.15])
%     subplot(2,2,3)
%     pcolor(lags/Fsd,wfreqs,squeeze(first_trig_phaselock(:,:,cc))');shading flat;colorbar
%     ylim([2 40])
%     xlim(xax)
%     subplot(2,2,2)
%     pcolor(lags/Fsd,wfreqs,squeeze(second_trig_specgram(:,:,cc))');shading flat;colorbar
%     xlim(xax)
%     caxis([-0.15 0.15])
%     subplot(2,2,4)
%     pcolor(lags/Fsd,wfreqs,squeeze(second_trig_phaselock(:,:,cc))');shading flat;colorbar
%     ylim([2 40])
%     xlim(xax)
% pause
% clf
% end

% save sac_lfp_modulation lags Fsd micro_* first_* bsac_* second_*