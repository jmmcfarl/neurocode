clear all
close all

cd ~/Data/bruce/2_27_12/M232
% load ./random_bar_eyedata.mat
load ./random_bar_eyedata_ftime.mat

load ./lemM232Expts.mat
%%
rf_cent = [4.73 -2.2];
axis_or = 50*pi/180; %in radians

Fs = 1000;
dsf = 2.5;
Fsd = Fs/dsf;

niqf = Fs/2;
[bb,aa] = butter(2,[1]/niqf,'high');
scales = logspace(log10(1.7),log10(85),40);
scales = scales*dsf;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);


new_dt = .0025;
dsfrac = 4;

%%
load ./CellList.mat
good_sus = find(all(CellList(bar_expts,:,1) > 0));
use_sus = 1:24;
%%
load ./un_bar_pos.mat
un_bar_pos(1:2) = [];

flen = 30;
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
        cur_bar_Op = [Expts{bar_expts(ee)}.Trials(tt).Op];
        
        cur_bar_Op_new = repmat(cur_bar_Op,1,dsfrac)';
        cur_bar_Op_new = cur_bar_Op_new(:);
        cur_bar_Op_new = cur_bar_Op_new(2:end-1);
        
        cur_t_edges = [cur_t_axis; Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4];
        cur_t_edges_new = cur_t_edges(1):new_dt:cur_t_edges(end);
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
        cur_used_inds_new(1:round(0.2/new_dt)) = 0;
        
        cur_t_since_tstart = cur_t_axis_new - cur_t_axis(1);
        
%         cur_used_inds_new(end-round(0.2/new_dt):end) = 0; %reduce phase estimation edge artifacts
        
        all_used_inds_new = [all_used_inds_new; cur_used_inds_new(:)];
        all_t_axis = [all_t_axis; cur_t_axis(:)];
        all_t_since_tstart = [all_t_since_tstart; cur_t_since_tstart(:)];
%         all_used_inds = [all_used_inds; cur_used_inds(:)];
        all_t_axis_new = [all_t_axis_new; cur_t_axis_new(:)];
        all_binned_spikes = [all_binned_spikes; cur_binned_spks];
        all_expt_vec = [all_expt_vec; ones(length(cur_t_axis_new),1)*bar_expts(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_t_axis_new),1)*tt];
        all_bar_Op = [all_bar_Op; cur_bar_Op_new(:)];
%         all_bar_Xmat = [all_bar_Xmat; bar_Xmat];
    end
    
    cur_blink_start_times = all_eye_ts{ee}(all_blink_startinds{ee});
    cur_blink_stop_times = all_eye_ts{ee}(all_blink_stopinds{ee});
    cur_blink_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_start_times));
    cur_blink_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_stop_times));
    cur_poss_blinks = find(~isnan(cur_blink_start_inds) & ~isnan(cur_blink_stop_inds));
    
    all_blink_inds = zeros(size(all_t_axis));
    for i = 1:length(cur_poss_blinks)
        all_blink_inds(cur_blink_start_inds(cur_poss_blinks(i)):cur_blink_stop_inds(cur_poss_blinks(i))) = 1;
    end
    
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
    fname = sprintf('lemM232A.%d.lfp.mat',bar_expts(ee));
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
    
%     %compute CSD
%     vars.Fs = Fs;
%     vars.BrainBound = 1;
%     vars.ChanSep = 0.05;
%      vars.diam = 2; %0.5
%     CSD = PettersenCSD(expt_lfps','spline',vars)';
%     expt_lfps = CSD;
% 

    seg_start_inds = [1; 1+find(diff(expt_lfp_t_axis) > 2/Fs)];
    seg_stop_inds = [find(diff(expt_lfp_t_axis) > 2/Fs); length(expt_lfp_t_axis)];
    res_t_axis = [];
    res_phasegram = [];
    res_ampgram = [];
    res_lfps = [];
    for jj = 1:length(seg_start_inds)
        cur_res_inds = seg_start_inds(jj):seg_stop_inds(jj);
        
        down_lfps = resample(expt_lfps(cur_res_inds,:),Fsd,Fs);
        
        cur_res_t = expt_lfp_t_axis(seg_start_inds(jj)):1/Fsd:expt_lfp_t_axis(seg_stop_inds(jj));
        extra_len = size(down_lfps,1) - length(cur_res_t);
        cur_res_t = [cur_res_t cur_res_t(end) + (1:extra_len)/Fsd];
        res_t_axis = [res_t_axis; cur_res_t(:)];
        
        cur_phasegram = nan(length(cur_res_t),length(wfreqs),24);
        cur_ampgram = nan(length(cur_res_t),length(wfreqs),24);
        for ll = 1:24
            temp = cwt(down_lfps(:,ll),scales,'cmor1-1');
            cur_phasegram(:,:,ll) = angle(temp)';
            cur_ampgram(:,:,ll) = abs(temp)';
        end
        res_phasegram = cat(1,res_phasegram,cur_phasegram);
        res_ampgram = cat(1,res_ampgram,cur_ampgram);
        res_lfps = [res_lfps; down_lfps];
    end
    
    cur_set = find(full_exptvec_new==ee);
    interp_lfps = interp1(res_t_axis,res_lfps,full_taxis_new(cur_set));
    unwr_phasegram = unwrap(res_phasegram);
    interp_phasegrams = interp1(res_t_axis,unwr_phasegram,full_taxis_new(cur_set));
    interp_phasegrams = mod(interp_phasegrams+pi,2*pi)-pi;
    interp_ampgrams = interp1(res_t_axis,res_ampgram,full_taxis_new(cur_set));
    
    full_lfps = [full_lfps; interp_lfps];
    full_phasegrams = cat(1,full_phasegrams, interp_phasegrams);
    full_ampgrams = cat(1,full_ampgrams, interp_ampgrams);
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
full_lfps_norm = bsxfun(@minus,full_lfps_norm,nanmean(full_ampgrams_norm));
full_lfps_norm = bsxfun(@rdivide,full_ampgrams_norm,nanstd(full_ampgrams_norm));

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
[micro_trig_specgram,lags] = get_event_trig_avg(full_ampgrams_norm,all_microsac_inds,backlag,forwardlag);

[second_trig_phaselock,lags] = get_event_trig_avg(exp(full_phasegrams_used*(1i)),all_secondsac_inds,backlag,forwardlag);
second_trig_phaselock = abs(second_trig_phaselock);
[first_trig_phaselock,lags] = get_event_trig_avg(exp(full_phasegrams_used*(1i)),all_firstsac_inds,backlag,forwardlag);
first_trig_phaselock = abs(first_trig_phaselock);
[micro_trig_phaselock,lags] = get_event_trig_avg(exp(full_phasegrams_used*(1i)),all_microsac_inds,backlag,forwardlag);
micro_trig_phaselock = abs(micro_trig_phaselock);

[micro_trig_avgs,lags] = get_event_trig_avg(full_ampgrams_norm,all_microsac_inds,backlag,forwardlag);

%%
lin_X = zeros(length(full_taxis_new),length(bar_expts)-1);
for i = 1:length(bar_expts)-1
    cur_set = find(full_exptvec_new==i);
    lin_X(cur_set,i) = 1;
end

[c,ia,ic] = unique([full_exptvec_new full_trialvec_new],'rows');
n_trials = length(ia);

rp = randperm(n_trials);
rand_trial_vec = rp(ic);
[~,ind_shuff] = sort(rand_trial_vec);

xv_frac = 0.;
n_xv_trials = round(n_trials*xv_frac);
xv_tset = randperm(n_trials);
xv_tset(n_xv_trials+1:end) = [];
tr_tset = find(~ismember(1:n_trials,xv_tset));
xv_inds = find(ismember(ic,xv_tset));
tr_inds = find(ismember(ic,tr_tset))';

xv_inds(full_used_inds_new(xv_inds)==0) = [];
tr_inds(full_used_inds_new(tr_inds) == 0) = [];

% %  Probably have to debug the trial-shuffling procedure after making
% changes (to include edge samples in initial data parsing)
% tr_inds_shuffle = ind_shuff(tr_inds);
% xv_inds_shuffle = ind_shuff(xv_inds);












%%
cur_dt = 0.01;
% flen_t = 0.5;
tent_centers = [0:cur_dt:0.7];
tent_centers = round(tent_centers/new_dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

shift = round(0.2/new_dt);
tbmat = [zeros(ntents,tblen-2*shift-1) tbmat];
tent_centers = tent_centers-shift;

trial_sac_inds = zeros(size(full_taxis_new));
trial_sac_inds(all_sac_inds) = 1;
trial_sac_mat = zeros(length(full_taxis_new),ntents);
for i = 1:ntents
    trial_sac_mat(:,i) = conv(trial_sac_inds,tbmat(i,:),'same');
end

trial_msac_inds = zeros(size(full_taxis_new));
trial_msac_inds(all_microsac_inds) = 1;
trial_msac_mat = zeros(length(full_taxis_new),ntents);
for i = 1:ntents
    trial_msac_mat(:,i) = conv(trial_msac_inds,tbmat(i,:),'same');
end

trial_bsac_inds = zeros(size(full_taxis_new));
trial_bsac_inds(all_bigsac_inds) = 1;
trial_bsac_mat = zeros(length(full_taxis_new),ntents);
for i = 1:ntents
    trial_bsac_mat(:,i) = conv(trial_bsac_inds,tbmat(i,:),'same');
end

trial_fsac_inds = zeros(size(full_taxis_new));
trial_fsac_inds(all_firstsac_inds) = 1;
trial_fsac_mat = zeros(length(full_taxis_new),ntents);
for i = 1:ntents
    trial_fsac_mat(:,i) = conv(trial_fsac_inds,tbmat(i,:),'same');
end

trial_ssac_inds = zeros(size(full_taxis_new));
trial_ssac_inds(all_secondsac_inds) = 1;
trial_ssac_mat = zeros(length(full_taxis_new),ntents);
for i = 1:ntents
    trial_ssac_mat(:,i) = conv(trial_ssac_inds,tbmat(i,:),'same');
end

%%
tr_sac_inds = find(trial_sac_inds(tr_inds)==1);
tr_bsac_inds = find(trial_bsac_inds(tr_inds)==1);
tr_msac_inds = find(trial_msac_inds(tr_inds)==1);
tr_fsac_inds = find(trial_fsac_inds(tr_inds)==1);
tr_ssac_inds = find(trial_ssac_inds(tr_inds)==1);
% xv_sac_inds = find(trial_sac_inds(xv_inds)==1);
% xv_msac_inds = find(trial_msac_inds(xv_inds)==1);
% xv_fsac_inds = find(trial_fsac_inds(xv_inds)==1);
% xv_bsac_inds = find(trial_bsac_inds(xv_inds)==1);
% xv_ssac_inds = find(trial_ssac_inds(xv_inds)==1);


%%

fNT = length(full_taxis_new);
% new_phase_set_po = [reshape(cos(full_phasegrams),length(full_taxis_new),length(wfreqs)*24) reshape(sin(full_phasegrams),length(full_taxis_new),length(wfreqs)*24)];
new_phase_set = [reshape(full_ampgrams,fNT,length(wfreqs)*24).*reshape(cos(full_phasegrams),fNT,length(wfreqs)*24) ...
    reshape(full_ampgrams,fNT,length(wfreqs)*24).*reshape(sin(full_phasegrams),fNT,length(wfreqs)*24)];
% phase_elec_set = ones(length(wfreqs),1)*(1:24);
% phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];
% phase_freq_set = ones(length(wfreqs),1)*(1:24);
% phase_freq_set = [phase_freq_set(:); phase_freq_set(:)];

[phase_elec_set,phase_freq_set] = meshgrid(1:24,wfreqs);
phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];
phase_freq_set = [phase_freq_set(:); phase_freq_set(:)];

% gamma_freqs = wfreqs(wfreqs >= 25 & wfreqs <= 50);
% alpha_freqs = wfreqs(wfreqs >= 5 & wfreqs <= 15);

use_elecs = 1:24;
use_set = find(ismember(phase_elec_set,use_elecs));
% Xmat = [new_phase_set(tr_inds,use_set)];
% Xmatxv = [new_phase_set(xv_inds,use_set)];
% % Xmatsh = [new_phase_set(tr_inds_shuffle,use_set)];
% % Xmatxvsh = [new_phase_set(xv_inds_shuffle,use_set)];


NL_type = 0; % 0 = logexp || 1 = exp

% reg_params.dl2_freq = 5000;
% reg_params.dl2_ch = 5000;
% reg_params.dl2_freq = 2500; %20000
% reg_params.dl2_ch =  5000; %20000
% reg_params.dl2_freq = 30; %20000
% reg_params.dl2_ch =  1000; %20000
% reg_params.dl2_time = 100; %500
% reg_params.dl_time = 20;
% reg_params.dl_freq = 30; %100
% reg_params.dl_ch =  1000; %100
% reg_params.d2_phase = 20;
% reg_params.d2_time = 1;

reg_params.dl2_freq = 5; %20000
reg_params.dl2_ch =  200; %20000
reg_params.dl2_time = 100; %500
reg_params.dl_time = 20;
reg_params.dl_freq = 5; %100
reg_params.dl_ch =  200; %100
reg_params.d2_phase = 20;
reg_params.d2_time = 1;

% reg_params2.dl2_freq = 100; %20000
% reg_params2.dl2_ch =  0; %20000
% reg_params2.dl2_time = 100; %500
% reg_params2.dl_time = 20;
% reg_params2.dl_freq = 10; %100
% reg_params2.dl_ch =  0; %100
% reg_params2.d2_phase = 50;
% reg_params2.d2_time = 1;

silent = 1;

%%
for cc = 1:length(use_sus)
fprintf('Cell %d of %d\n',cc,length(use_sus));
    
    Robs = full_spkbinned(tr_inds,cc);
    tr_spkbns = convert_to_spikebins(Robs);
    
    glm_kern = get_k_mat(glm_fit(cc));
    stim_out = full_Xmat*glm_kern;
    stim_out_interp = stim_out(tr_inds);
    
    %NULL MODEL
    cur_Xmat = [lin_X(tr_inds,:)];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_null,grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, 1, [], [],[], [], [], [], NL_type);

    %STIM ONLY MODEL
    cur_Xmat = [stim_out_interp lin_X(tr_inds,:)];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_so(cc),grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, 1, [], [],[], [], [], [], NL_type);

%     cur_Xmat = [trial_sac_mat(tr_inds,:) stim_out_interp lin_X(tr_inds,:)];
%     lamrange = [reg_params.dl_time 1 ntents 0];
%     lamrange2 = [reg_params.dl2_time 1 ntents 0];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_sac,grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, silent, lamrange, lamrange2,[], [], [], [], NL_type);
%     so_sac_kern(cc,:) = fitp_sac.k(1:ntents);
%     predrate_out = cur_Xmat*fitp_sac.k(1:end-1) + fitp_sac.k(end);
%     [st_avg_sacmod(cc,:),cur_lags] = get_event_trig_avg(exp(predrate_out),tr_sac_inds,round(0.2/new_dt),round(0.5/new_dt));
    
%     %SAC ONLY MODEL
%     cur_Xmat = [trial_bsac_mat(tr_inds,:) trial_msac_mat(tr_inds,:) lin_X(tr_inds,:)];
%     lamrange = [reg_params.dl_time 1 ntents 0; reg_params.dl_time ntents+1 2*ntents 0;];
%     lamrange2 = [reg_params.dl2_time 1 ntents 0; reg_params.dl2_time ntents+1 2*ntents 0];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_sac,grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, silent, lamrange, lamrange2,[], [], [], [], NL_type);
%     bsac_kern(cc,:) = fitp_sac.k(1:ntents);
%     msac_kern(cc,:) = fitp_sac.k(ntents+1:2*ntents);

%     %SAC ONLY MODEL
%     cur_Xmat = [trial_sac_mat(tr_inds,:) lin_X(tr_inds,:)];
%     lamrange = [reg_params.dl_time 1 ntents 0];
%     lamrange2 = [reg_params.dl2_time 1 ntents 0];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_sac,grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, silent, lamrange, lamrange2,[], [], [], [], NL_type);
%     sac_kern(cc,:) = fitp_sac.k(1:ntents);
% 
%     %STIM + SAC MODEL
%     cur_Xmat = [trial_sac_mat(tr_inds,:) stim_out_interp lin_X(tr_inds,:)];
%     lamrange = [reg_params.dl_time 1 ntents 0];
%     lamrange2 = [reg_params.dl2_time 1 ntents 0];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_st_sac,grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, silent, lamrange, lamrange2,[], [], [], [], NL_type);
%     st_bac_kern(cc,:) = fitp_st_sac.k(1:ntents);
   

    use_elecs = 1:24;
%     if mod(cc,2)==0
%         use_elecs = 2:2:24;
%     else
%         use_elecs = 1:2:24;
%     end
    use_freqs = wfreqs(wfreqs < 35);
    use_set = find(ismember(phase_elec_set,use_elecs) & ismember(phase_freq_set,use_freqs));

%     %PHASE ONLY MODEL
%     cur_Xmat = [new_phase_set(tr_inds,use_set) lin_X(tr_inds,:)];
%     stim_params = [length(use_freqs),length(use_elecs)];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_ap(cc)] = fit_GLM_sinphase_model(cur_Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
%     ap_sinphase_cfilt(cc,:) = fitp_ap(cc).k(1:length(use_elecs)*length(use_freqs));
%     ap_sinphase_sfilt(cc,:) = fitp_ap(cc).k((length(use_elecs)*length(use_freqs)+1):length(use_elecs)*length(use_freqs)*2);
%     ap_ampkern(cc,:) = sqrt(ap_sinphase_cfilt(cc,:).^2 + ap_sinphase_sfilt(cc,:).^2);
% %     phasemod_out(cc,:) = cur_Xmat*fitp_ap(cc).k(1:end-1) + fitp_ap(cc).k(end);
%     phasemod_out(cc,:) = log(1+exp(phasemod_out(cc,:)));
%     [st_avg_phasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out(cc,:),tr_sac_inds,round(0.4/new_dt),round(0.5/new_dt));
%     [st_avg_mphasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out(cc,:),tr_msac_inds,round(0.4/new_dt),round(0.5/new_dt));
%     [st_avg_bphasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out(cc,:),tr_bsac_inds,round(0.4/new_dt),round(0.5/new_dt));


    %PHASE + STIM MODEL
    cur_Xmat = [new_phase_set(tr_inds,use_set) stim_out_interp lin_X(tr_inds,:)];
    stim_params = [length(use_freqs),length(use_elecs)];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_ap_st(cc)] = fit_GLM_sinphase_model(cur_Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    ap_st_sinphase_cfilt(cc,:) = fitp_ap_st(cc).k(1:length(use_elecs)*length(use_freqs));
    ap_st_sinphase_sfilt(cc,:) = fitp_ap_st(cc).k((length(use_elecs)*length(use_freqs)+1):length(use_elecs)*length(use_freqs)*2);
    ap_st_ampkern(cc,:) = sqrt(ap_st_sinphase_cfilt(cc,:).^2 + ap_st_sinphase_sfilt(cc,:).^2);
    phasemod_out(cc,:) = cur_Xmat*fitp_ap_st(cc).k(1:end-1) + fitp_ap_st(cc).k(end);
    phasemod_out(cc,:) = log(1+exp(phasemod_out(cc,:)));
    [st_avg_phasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out(cc,:),tr_sac_inds,round(0.4/new_dt),round(0.5/new_dt));
    [st_avg_mphasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out(cc,:),tr_msac_inds,round(0.4/new_dt),round(0.5/new_dt));
    [st_avg_bphasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out(cc,:),tr_bsac_inds,round(0.4/new_dt),round(0.5/new_dt));
% % 
%     %PHASE + SAC MODEL
%     cur_Xmat = [new_phase_set(tr_inds,use_set) trial_sac_mat(tr_inds,:) lin_X(tr_inds,:)];
%     stim_params = [length(wfreqs),length(use_elecs) length(tent_centers)];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_ap_sac(cc)] = fit_GLM_sinphase_model(cur_Xmat, Robs, K0,silent, stim_params, reg_params,1,NL_type);
%     ap_sac_sinphase_cfilt(cc,:) = fitp_ap_sac(cc).k(1:length(use_elecs)*length(wfreqs));
%     ap_sac_sinphase_sfilt(cc,:) = fitp_ap_sac(cc).k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
%     ap_sac_ampkern(cc,:) = sqrt(ap_sac_sinphase_cfilt(cc,:).^2 + ap_sac_sinphase_sfilt(cc,:).^2);

%     %PHASE + SAC + STIM MODEL
%     cur_Xmat = [new_phase_set(tr_inds,use_set) trial_sac_mat(tr_inds,:) stim_out_interp lin_X(tr_inds,:)];
%     stim_params = [length(wfreqs),length(use_elecs) length(tent_centers)];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_ap_st_sac(cc)] = fit_GLM_sinphase_model(cur_Xmat, Robs, K0,silent, stim_params, reg_params,1,NL_type);
%     ap_st_sac_sinphase_cfilt(cc,:) = fitp_ap_st_sac(cc).k(1:length(use_elecs)*length(wfreqs));
%     ap_st_sac_sinphase_sfilt(cc,:) = fitp_ap_st_sac(cc).k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
%     ap_st_sac_ampkern(cc,:) = sqrt(ap_st_sac_sinphase_cfilt(cc,:).^2 + ap_st_sac_sinphase_sfilt(cc,:).^2);

%     %SINGLE_ELECTRODE MODEL
%     use_singel = cc;
%     use_set_se = find(ismember(phase_elec_set,use_singel));
%     cur_Xmat = [new_phase_set(tr_inds,use_set_se) trial_sac_mat(tr_inds,:) stim_out_interp lin_X(tr_inds,:)];
%     stim_params = [length(wfreqs),length(use_singel) length(tent_centers)];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_apst_se(cc)] = fit_GLM_sinphase_model(cur_Xmat, Robs, K0,silent, stim_params, reg_params2,1,NL_type);
%     apst_se_sinphase_cfilt(cc,:) = fitp_apst_se(cc).k(1:length(use_singel)*length(wfreqs));
%     apst_se_sinphase_sfilt(cc,:) = fitp_apst_se(cc).k((length(use_singel)*length(wfreqs)+1):length(use_singel)*length(wfreqs)*2);
%     apst_se_ampkern(cc,:) = sqrt(apst_se_sinphase_cfilt(cc,:).^2 + apst_se_sinphase_sfilt(cc,:).^2);
% 
%     %PURE PHASE MODEL
%     cur_Xmat = [new_phase_set_po(tr_inds,use_set) trial_sac_mat(tr_inds,:) stim_out_interp lin_X(tr_inds,:)];
%     stim_params = [length(wfreqs),length(use_elecs) length(tent_centers)];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_pst(cc)] = fit_GLM_sinphase_model(cur_Xmat, Robs, K0,silent, stim_params, reg_params,1,NL_type);
%     pst_sinphase_cfilt(cc,:) = fitp_pst(cc).k(1:length(use_elecs)*length(wfreqs));
%     pst_sinphase_sfilt(cc,:) = fitp_pst(cc).k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
%     pst_ampkern(cc,:) = sqrt(pst_sinphase_cfilt(cc,:).^2 + pst_sinphase_sfilt(cc,:).^2);

    
% %     phasemod_out = Xmat(:,use_set)*fitp_po(cc).k(1:2*length(use_elecs)*length(wfreqs));
% %     phasemod_out = (Xmat(:,use_set)*fitp_post(cc).k(1:2*length(use_elecs)*length(wfreqs)));
%     phasemod_out = cur_Xmat*fitp_post(cc).k(1:end-1) + fitp_post(cc).k(end);
%     phasemod_out = log(1+exp(phasemod_out));
%     [st_avg_phasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_sac_inds,round(0.4/new_dt),round(0.5/new_dt));
%     [st_avg_mphasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_msac_inds,round(0.4/new_dt),round(0.5/new_dt));
%     [st_avg_bphasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_bsac_inds,round(0.4/new_dt),round(0.5/new_dt));
% %     [st_avg_fphasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_fsac_inds,round(0.4/new_dt),round(0.5/new_dt));
% %     [st_avg_sphasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_ssac_inds,round(0.4/new_dt),round(0.5/new_dt));
% %     [st_avg_allmod_po(cc,:),cur_lags] = get_event_trig_avg(exp(predrate_out),tr_sac_inds,round(0.2/new_dt),round(0.5/new_dt));

    
%     xv_Robs = full_spkbinned(xv_inds,cc);
%     stim_out_interp = stim_out(xv_inds);
%     
%     %STIM ONLY
%     cur_Xmat = [stim_out_interp lin_X(xv_inds,:)];
%     xv_so_pred_rate = cur_Xmat*fitp_so(cc).k(1:end-1) + fitp_so(cc).k(end);
%     if NL_type == 1
%         xv_so_pred_rate = exp(xv_so_pred_rate);
%     else
%         xv_so_pred_rate = log(1+exp(xv_so_pred_rate));
%     end
%     xv_so_LL(cc) = -sum(xv_Robs.*log(xv_so_pred_rate)-xv_so_pred_rate)/sum(xv_Robs);
    
%     %SAC ONLY
%     cur_Xmat = [trial_sac_mat(xv_inds,:) lin_X(xv_inds,:)];
%     xv_so_pred_rate = cur_Xmat*fitp_sac.k(1:end-1) + fitp_sac.k(end);
%     if NL_type == 1
%         xv_so_pred_rate = exp(xv_so_pred_rate);
%     else
%         xv_so_pred_rate = log(1+exp(xv_so_pred_rate));
%     end
%     xv_sac_LL(cc) = -sum(xv_Robs.*log(xv_so_pred_rate)-xv_so_pred_rate)/sum(xv_Robs);

%     cur_Xmat = [trial_bsac_mat(xv_inds,:) trial_msac_mat(xv_inds,:) stim_out_interp lin_X(xv_inds,:)];
%     xv_so_pred_rate = cur_Xmat*fitp_mulsac.k(1:end-1) + fitp_mulsac.k(end);
%     if NL_type == 1
%         xv_so_pred_rate = exp(xv_so_pred_rate);
%     else
%         xv_so_pred_rate = log(1+exp(xv_so_pred_rate));
%     end
%     xv_mulsac_LL(cc) = -sum(xv_Robs.*log(xv_so_pred_rate)-xv_so_pred_rate)/sum(xv_Robs);

% %STIM + SAC
% cur_Xmat = [trial_sac_mat(xv_inds,:) stim_out_interp lin_X(xv_inds,:)];
% xv_so_pred_rate = cur_Xmat*fitp_st_sac.k(1:end-1) + fitp_st_sac.k(end);
% if NL_type == 1
%     xv_so_pred_rate = exp(xv_so_pred_rate);
% else
%     xv_so_pred_rate = log(1+exp(xv_so_pred_rate));
% end
% xv_st_sac_LL(cc) = -sum(xv_Robs.*log(xv_so_pred_rate)-xv_so_pred_rate)/sum(xv_Robs);
  

% %PHASE ONLY MODEL
% cur_Xmat = [new_phase_set(xv_inds,use_set) lin_X(xv_inds,:)];
% xv_pred_rate = cur_Xmat*fitp_ap(cc).k(1:end-1) + fitp_ap(cc).k(end);
% if NL_type == 1
%     xv_pred_rate = exp(xv_pred_rate);
% else
%     xv_pred_rate = log(1+exp(xv_pred_rate));
% end
% xv_ap_LL(cc) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);
% 
% % 
% 
% %PHASE + STIM MODEL
%     cur_Xmat = [new_phase_set(xv_inds,use_set) stim_out_interp lin_X(xv_inds,:)];
% xv_pred_rate = cur_Xmat*fitp_ap_st(cc).k(1:end-1) + fitp_ap_st(cc).k(end);
% if NL_type == 1
%     xv_pred_rate = exp(xv_pred_rate);
% else
%     xv_pred_rate = log(1+exp(xv_pred_rate));
% end
% xv_ap_st_LL(cc) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);
% 
% phasemod_out = cur_Xmat*fitp_ap_st(cc).k(1:end-1) + fitp_ap_st(cc).k(end);
% phasemod_out = log(1+exp(phasemod_out));
% [xvst_avg_phasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,xv_sac_inds,round(0.4/new_dt),round(0.5/new_dt));
% [xvst_avg_mphasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,xv_msac_inds,round(0.4/new_dt),round(0.5/new_dt));
% [xvst_avg_bphasemod_post(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,xv_bsac_inds,round(0.4/new_dt),round(0.5/new_dt));
% 

% %PHASE + SAC MODEL
%     cur_Xmat = [new_phase_set(xv_inds,use_set) trial_sac_mat(xv_inds,:) lin_X(xv_inds,:)];
% xv_pred_rate = cur_Xmat*fitp_ap_sac(cc).k(1:end-1) + fitp_ap_sac(cc).k(end);
% if NL_type == 1
%     xv_pred_rate = exp(xv_pred_rate);
% else
%     xv_pred_rate = log(1+exp(xv_pred_rate));
% end
% xv_ap_sac_LL(cc) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);

% %PHASE + SAC MODEL
%     cur_Xmat = [new_phase_set(xv_inds,use_set) trial_sac_mat(xv_inds,:) stim_out_interp lin_X(xv_inds,:)];
% xv_pred_rate = cur_Xmat*fitp_ap_st_sac(cc).k(1:end-1) + fitp_ap_st_sac(cc).k(end);
% if NL_type == 1
%     xv_pred_rate = exp(xv_pred_rate);
% else
%     xv_pred_rate = log(1+exp(xv_pred_rate));
% end
% xv_ap_st_sac_LL(cc) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);

% %SINGLE_ELECTRODE
%     cur_Xmat = [new_phase_set(xv_inds,use_set_se) trial_sac_mat(xv_inds,:) stim_out_interp lin_X(xv_inds,:)];
%     xv_fphase_pred_rate = cur_Xmat*fitp_apst_se(cc).k(1:end-1) + fitp_apst_se(cc).k(end);
%     if NL_type == 0
%         xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
%     else
%         xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
%     end
%     xv_apst_se_LL(cc) = -sum(xv_Robs.*log(xv_fphase_pred_rate)-xv_fphase_pred_rate)/sum(xv_Robs);

%     %PHASE_ONLY 
%     cur_Xmat = [new_phase_set_po(xv_inds,use_set) trial_sac_mat(xv_inds,:) stim_out_interp lin_X(xv_inds,:)];
%     xv_fphase_pred_rate = cur_Xmat*fitp_pst(cc).k(1:end-1) + fitp_pst(cc).k(end);
%     if NL_type == 0
%         xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
%     else
%         xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
%     end
%     xv_pst_LL(cc) = -sum(xv_Robs.*log(xv_fphase_pred_rate)-xv_fphase_pred_rate)/sum(xv_Robs);
    



%     cur_Xmat = [lin_X(xv_inds,:)];
%     xv_null_pred_rate = cur_Xmat*fitp_null.k(1:end-1) + fitp_null.k(end);
%     if NL_type == 0
%         xv_null_pred_rate = log(1+exp(xv_null_pred_rate));
%     else
%         xv_null_pred_rate = exp(xv_null_pred_rate);
%     end
%     xv_null_LL(cc) = -sum(xv_Robs.*log(xv_null_pred_rate)-xv_null_pred_rate)/sum(xv_Robs);
%     
%     fprintf('Stim only: %.4f\n',(xv_null_LL(cc) - xv_so_LL(cc))/log(2));
%     fprintf('Sac only: %.4f\n',(xv_null_LL(cc) - xv_sac_LL(cc))/log(2));
%     fprintf('Sac + Stim: %.4f\n',(xv_null_LL(cc) - xv_st_sac_LL(cc))/log(2));
%     fprintf('Phase only: %.4f\n',(xv_null_LL(cc) - xv_ap_LL(cc))/log(2));
%     fprintf('Phase + sac: %.4f\n',(xv_null_LL(cc) - xv_ap_sac_LL(cc))/log(2));
%     fprintf('Phase + stim: %.4f\n',(xv_null_LL(cc) - xv_ap_st_LL(cc))/log(2));
%     fprintf('Phase + stim + sac: %.4f\n',(xv_null_LL(cc) - xv_ap_st_sac_LL(cc))/log(2));
%     fprintf('Se phase: %.4f\n',(xv_null_LL(cc) - xv_apst_se_LL(cc))/log(2));
%     fprintf('pure phase: %.4f\n',(xv_null_LL(cc) - xv_pst_LL(cc))/log(2));
   %% 
end
% po_phasekern = 180/pi*(atan2(po_sinphase_cfilt,po_sinphase_sfilt)+pi);
% post_phasekern = 180/pi*(atan2(post_sinphase_cfilt,post_sinphase_sfilt)+pi);

%%
close all
stim_imp = (xv_null_LL - xv_so_LL)/log(2);
sac_imp = (xv_null_LL - xv_sac_LL)/log(2);
sac_stim_imp = (xv_null_LL - xv_st_sac_LL)/log(2);
phase_imp = (xv_null_LL - xv_ap_LL)/log(2);
phase_sac_imp = (xv_null_LL - xv_ap_sac_LL)/log(2);
phase_stim_imp = (xv_null_LL - xv_ap_st_LL)/log(2);
phase_stim_sac_imp = (xv_null_LL - xv_ap_st_sac_LL)/log(2);
sephase_imp = (xv_null_LL - xv_apst_se_LL)/log(2);
purephase_imp = (xv_null_LL - xv_pst_LL)/log(2);

stim_imp = stim_imp./phase_stim_sac_imp;
sac_imp = sac_imp./phase_stim_sac_imp;
sac_stim_imp = sac_stim_imp./phase_stim_sac_imp;
phase_imp = phase_imp./phase_stim_sac_imp;
phase_sac_imp = phase_sac_imp./phase_stim_sac_imp;
phase_stim_imp = phase_stim_imp./phase_stim_sac_imp;
sephase_imp = sephase_imp./phase_stim_sac_imp;
purephase_imp = purephase_imp./phase_stim_sac_imp;
phase_stim_sac_imp = phase_stim_sac_imp./phase_stim_sac_imp;



sus = good_sus;
mus = setdiff(1:24,sus);
mus(mus==6) = [];
cur_sus = sus;

% figure
% boxplot([stim_imp(cur_sus)' sac_imp(cur_sus)' sac_stim_imp(cur_sus)' phase_imp(cur_sus)' phase_sac_imp(cur_sus)' phase_stim_imp(cur_sus)' phase_stim_sac_imp(cur_sus)'],...
%     {'stim','sac','stim+sac','LFP','LFP+sac','LFP+stim','LFP+stim+sac'});
% ylim([0 1.3])
% ylabel('LL improvement (bits/spk)','fontsize',16)
% % ylabel('Relative LL improvement','fontsize',16)

figure
boxplot([sephase_imp(cur_sus)' purephase_imp(cur_sus)' phase_stim_sac_imp(cur_sus)'],...
    {'single-electrode','phase only','Full model'});
ylim([0 1.05])
ylabel('LL improvement (bits/spk)','fontsize',16)
% ylabel('Relative LL improvement','fontsize',16)

%%
% stim_imp = (xv_null_LL - xv_so_LL)/log(2);
% sac_imp = (xv_null_LL - xv_sac_LL)/log(2);
% mulsac_imp = (xv_null_LL - xv_mulsac_LL)/log(2);
% phase_imp = (xv_null_LL - xv_pst_LL)/log(2);
% ampphase_se_imp = (xv_null_LL - xv_apst_se_LL)/log(2);
% ampphase_imp = (xv_null_LL - xv_apst_LL)/log(2);
% % stimphase_imp = (xv_null_LL - xv_post_LL)/log(2);
% 
% % figure
% % boxplot([stim_imp(:) mulsac_imp(:) stimphase_imp(:)]);
% 
% sus = good_sus;
% mus = setdiff(1:24,sus);
% mus(mus==6) = [];
% figure
% subplot(1,2,1)
% boxplot([stim_imp(sus)' mulsac_imp(sus)' ampphase_se_imp(sus)' phase_imp(sus)' ampphase_imp(sus)'],{'stim','sac','single-ch amp-phase','all-ch phase','all-ch amp-phase'});
% ylim([0 1.3])
% ylabel('LL improvement (bits/spk)','fontsize',16)
% title('SUs')
% subplot(1,2,2)
% boxplot([stim_imp(mus)' mulsac_imp(mus)' ampphase_se_imp(mus)' phase_imp(mus)' ampphase_imp(mus)'],{'stim','sac','single-ch amp-phase','all-ch phase','all-ch amp-phase'});
% ylim([0 1.3])
% ylabel('LL improvement (bits/spk)','fontsize',16)
% title('SUs')

%%
stim_imp = (xv_null_LL - xv_so_LL)/log(2);
sac_imp = (xv_so_LL - xv_sac_LL)/log(2);
mulsac_imp = (xv_so_LL - xv_mulsac_LL)/log(2);
phase_imp = (xv_so_LL - xv_post_LL)/log(2);

%%
sm_sig = (0.01/new_dt);
travg_spk_rates = mean(full_spkbinned(tr_inds,:));
xvavg_spk_rates = mean(full_spkbinned(xv_inds,:));
for cc =1:24
    norm_spkbinned(:,cc) = jmm_smooth_1d_cor(full_spkbinned(:,cc),sm_sig);
end
norm_spkbinnedxv = bsxfun(@rdivide,norm_spkbinned,xvavg_spk_rates);
norm_spkbinnedtr = bsxfun(@rdivide,norm_spkbinned,travg_spk_rates);
% norm_spkbinnedxv = norm_spkbinned;
% norm_spkbinnedtr = norm_spkbinned;

tr_trial_start_inds = 1 + find(diff(ic(tr_inds)) > 0);

clear *trig_avg_rate
for cc = 1:24
    [sac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(norm_spkbinnedtr(tr_inds,cc),tr_sac_inds,round(0.4/new_dt),round(0.5/new_dt));
    [bsac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(norm_spkbinnedtr(tr_inds,cc),tr_bsac_inds,round(0.4/new_dt),round(0.5/new_dt));
    [msac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(norm_spkbinnedtr(tr_inds,cc),tr_msac_inds,round(0.4/new_dt),round(0.5/new_dt));
    [fsac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(norm_spkbinnedtr(tr_inds,cc),tr_fsac_inds,round(0.4/new_dt),round(0.5/new_dt));
    [ssac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(norm_spkbinnedtr(tr_inds,cc),tr_ssac_inds,round(0.4/new_dt),round(0.5/new_dt));
    [trial_trig_avg_rate(cc,:),cur_tlags] = get_event_trig_avg(norm_spkbinnedtr(tr_inds,cc),tr_trial_start_inds,round(0/new_dt),round(2/new_dt));

    [trial_trig_avg_lfp(cc,:),cur_tlags] = get_event_trig_avg(full_lfps(tr_inds,cc),tr_trial_start_inds,round(0/new_dt),round(2/new_dt));

    %     [xvsac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(norm_spkbinnedxv(xv_inds,cc),xv_sac_inds,round(0.4/new_dt),round(0.5/new_dt));
%     [xvbsac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(norm_spkbinnedxv(xv_inds,cc),xv_bsac_inds,round(0.4/new_dt),round(0.5/new_dt));
%     [xvmsac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(norm_spkbinnedxv(xv_inds,cc),xv_msac_inds,round(0.4/new_dt),round(0.5/new_dt));
%     [xvfsac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(norm_spkbinnedxv(xv_inds,cc),xv_fsac_inds,round(0.4/new_dt),round(0.5/new_dt));
%     [xvssac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(norm_spkbinnedxv(xv_inds,cc),xv_ssac_inds,round(0.4/new_dt),round(0.5/new_dt));
end

% bad_uns = [12 23 24];
% msac_trig_avg_rate(bad_uns,:) = nan;
% bsac_trig_avg_rate(bad_uns,:) = nan;

%%
cd ~/Data/bruce/2_27_12/M232
% save random_bar_phase_models_lfp_allmods.mat fit* xv_* reg_params* *filt wfreqs *ampkern reg_params st* cur_lags *pst* *_kern tent_centers new_dt n_bar_pos un_bar_pos flen
% save random_bar_phase_models_lfp_ds2.mat fit* xv_* reg_params* *filt wfreqs *ampkern reg_params st* cur_lags *pst* *_kern tent_centers new_dt n_bar_pos un_bar_pos flen
% save random_bar_phase_models_csd2.mat fit* xv_* reg_params* *filt wfreqs *ampkern reg_params st* cur_lags *_kern tent_centers new_dt un_bar_pos flen


%%
for cc = 1:24
        pcolor(wfreqs,1:13,[reshape(ap_ampkern(cc,:)',length(wfreqs),12)'; zeros(1,length(wfreqs))]);shading interp; set(gca,'xscale','log');
pause
clf
end
%%
phasekern = 180/pi*(atan2(ap_sinphase_cfilt,ap_sinphase_sfilt)+pi);

[WF,CH] = meshgrid(wfreqs,(1:24)*0.05);
WF = WF'; CH = CH';
linw = linspace(wfreqs(end),wfreqs(1),50);
linch = (1:.05:24)*0.05;
[lWF,lCH] = meshgrid(linw,linch);

close all
for cc = 1:24
    cc
%         fprintf('Stim only: %.4f\n',(xv_null_LL(cc) - xv_so_LL(cc))/log(2));
%     fprintf('Sac + Stim: %.4f\n',(xv_null_LL(cc) - xv_sac_LL(cc))/log(2));
% %     fprintf('Stim + phase: %.4f\n',(xv_null_LL(cc) - xv_fphase_LL(cc))/log(2));
%     fprintf('Stim + phase: %.4f\n',(xv_null_LL(cc) - xv_post_LL(cc))/log(2));
% 
    subplot(3,2,1);hold on
    plot(tent_centers*new_dt,exp(so_sac_kern(cc,:)),'r')
%     plot(cur_lags*new_dt,exp(st_avg_phasemod_post(cc,:)),'b')
%     plot(cur_lags*new_dt,sac_trig_avg_rate(cc,:),'k')
% %     plot(cur_lags*new_dt,st_avg_phasemod_po(cc,:),'k')
% %     plot(tent_centers*new_dt,phase_sac_kern(cc,:),'k')
axis tight
    xlim([-0.2 0.4])
    yl = ylim(); xl = xlim();
    line([0 0],yl,'color','k')
    line(xl,[0 0],'color','k')
    xlabel('Time relative to saccade (s)')
%     title(sprintf('Sac imp: %.4f',sac_imp(cc)));
    
        subplot(3,2,2);
    k_mat = get_k_mat(glm_fit(cc));
    imagesc(un_bar_pos,(0:-1:-flen+1)*new_dt,flipud(reshape(k_mat,flen,n_bar_pos)));
    cca = max(abs(k_mat(:)));
    caxis([-cca cca]*0.85);
    xlabel('Bar position deg)')
    ylabel('Time lag (s)')
%     title(sprintf('Stim imp: %.4f',stim_imp(cc)));
    subplot(3,2,3)
%     F = TriScatteredInterp(WF(:),CH(:),po_ampkern(cc,:)');
%     Vq = F(lWF,lCH);
%     imagesc(linw,linch,Vq); set(gca,'ydir','normal');
    pcolor(wfreqs,1:25,[reshape(post_ampkern(cc,:)',length(wfreqs),24)'; zeros(1,length(wfreqs))]);shading interp; set(gca,'xscale','log');
%     pcolor(wfreqs,1:13,[reshape(post_ampkern(cc,:)',length(wfreqs),12)'; zeros(1,length(wfreqs))]);shading interp; set(gca,'xscale','log');
    ca1 = max(post_ampkern(cc,:));
    caxis([0 ca1]*0.75);
%     title(sprintf('Amplitude kernel, phase imp: %.4f',phase_imp(cc)))
    xlabel('Frequency (Hz)')
    ylabel('Depth (mm)')
    
    ca = max([max(abs(post_sinphase_cfilt(cc,:))) max(abs(post_sinphase_sfilt(cc,:)))]);
    subplot(3,2,4)
%     F = TriScatteredInterp(WF(:),CH(:),phasekern(cc,:)');
%     Vq = F(lWF,lCH);
%     imagesc(linw,linch,Vq); set(gca,'ydir','normal');
    pcolor(wfreqs,1:25,[reshape(post_phasekern(cc,:)',length(wfreqs),24)'; zeros(1,length(wfreqs))]);shading interp; set(gca,'xscale','log');
%     pcolor(wfreqs,1:13,[reshape(post_phasekern(cc,:)',length(wfreqs),12)'; zeros(1,length(wfreqs))]);shading interp; set(gca,'xscale','log');
    caxis([0 360]);
    title('Phase kernel')
    xlabel('Frequency (Hz)')
    ylabel('Depth (mm)')
    subplot(3,2,5)
%     F = TriScatteredInterp(WF(:),CH(:),sinphase_cfilt(cc,:)');
%     Vq = F(lWF,lCH);
%     imagesc(linw,linch,Vq); set(gca,'ydir','normal');
    pcolor(wfreqs,1:25,[reshape(post_sinphase_cfilt(cc,:)',length(wfreqs),24)'; zeros(1,length(wfreqs))]);shading interp; set(gca,'xscale','log');
%     pcolor(wfreqs,1:13,[reshape(post_sinphase_cfilt(cc,:)',length(wfreqs),12)'; zeros(1,length(wfreqs))]);shading interp; set(gca,'xscale','log');
    caxis(0.8*[-ca ca])
    title('Sine kernel')
    xlabel('Frequency (Hz)')
    ylabel('Depth (mm)')
    subplot(3,2,6)
%     F = TriScatteredInterp(WF(:),CH(:),sinphase_cfilt(cc,:)');
%     Vq = F(lWF,lCH);
%     imagesc(linw,linch,Vq); set(gca,'ydir','normal');
    pcolor(wfreqs,1:25,[reshape(post_sinphase_sfilt(cc,:)',length(wfreqs),24)'; zeros(1,length(wfreqs))]);shading interp; set(gca,'xscale','log');
%     pcolor(wfreqs,1:13,[reshape(post_sinphase_sfilt(cc,:)',length(wfreqs),12)'; zeros(1,length(wfreqs))]);shading interp; set(gca,'xscale','log');
    caxis(0.8*[-ca ca])
    title('Cosine kernel')
    xlabel('Frequency (Hz)')
    ylabel('Depth (mm)')
    
%     fillPage(gcf,'Papersize',[6 8])
%     gname = sprintf('lprobe_phasemodels_Unit%d_lfp',cc);
%     print(gname,'-dpng');
%     close
    
    pause
    clf
end

%%
for cc = 1:24
    temp_kmat = reshape(ap_ampkern(cc,:)',length(wfreqs),12)'; zeros(1,length(wfreqs));
    avg_fvec(cc,:) = mean(temp_kmat);
    cur_ch = floor(cc/2)+1;
    cur_ch(cur_ch > 12) = 12;
    same_fvec(cc,:) = temp_kmat(cur_ch,:);
end


%%

disp('Computing big sac stats')
n_boot = 100;
[sac_trg_mat,cur_lags] = get_event_trig_mat(norm_spkbinnedtr(tr_inds,:),tr_bsac_inds,round(0.4/new_dt),round(0.5/new_dt));
sac_trg_mat = reshape(sac_trg_mat,length(tr_bsac_inds),24*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 24]);
avg_bsac_trig_rate = squeeze(mean(boot_avgs))';
std_bsac_trig_rate = squeeze(std(boot_avgs))';

phasemod_out_norm = bsxfun(@rdivide,phasemod_out,travg_spk_rates');
[sac_trg_mat,cur_lags] = get_event_trig_mat(phasemod_out_norm',tr_bsac_inds,round(0.4/new_dt),round(0.5/new_dt));
sac_trg_mat = reshape(sac_trg_mat,length(tr_bsac_inds),24*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 24]);
avg_bsac_mod_pred = squeeze(mean(boot_avgs))';
std_bsac_mod_pred = squeeze(std(boot_avgs))';

disp('Computing first sac stats')
[sac_trg_mat,cur_lags] = get_event_trig_mat(norm_spkbinnedtr(tr_inds,:),tr_fsac_inds,round(0.4/new_dt),round(0.5/new_dt));
sac_trg_mat = reshape(sac_trg_mat,length(tr_fsac_inds),24*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 24]);
avg_fsac_trig_rate = squeeze(mean(boot_avgs))';
std_fsac_trig_rate = squeeze(std(boot_avgs))';

phasemod_out_norm = bsxfun(@rdivide,phasemod_out,travg_spk_rates');
[sac_trg_mat,cur_lags] = get_event_trig_mat(phasemod_out_norm',tr_fsac_inds,round(0.4/new_dt),round(0.5/new_dt));
sac_trg_mat = reshape(sac_trg_mat,length(tr_fsac_inds),24*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 24]);
avg_fsac_mod_pred = squeeze(mean(boot_avgs))';
std_fsac_mod_pred = squeeze(std(boot_avgs))';

disp('Computing second sac stats')
[sac_trg_mat,cur_lags] = get_event_trig_mat(norm_spkbinnedtr(tr_inds,:),tr_ssac_inds,round(0.4/new_dt),round(0.5/new_dt));
sac_trg_mat = reshape(sac_trg_mat,length(tr_ssac_inds),24*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 24]);
avg_ssac_trig_rate = squeeze(mean(boot_avgs))';
std_ssac_trig_rate = squeeze(std(boot_avgs))';

phasemod_out_norm = bsxfun(@rdivide,phasemod_out,travg_spk_rates');
[sac_trg_mat,cur_lags] = get_event_trig_mat(phasemod_out_norm',tr_ssac_inds,round(0.4/new_dt),round(0.5/new_dt));
sac_trg_mat = reshape(sac_trg_mat,length(tr_ssac_inds),24*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 24]);
avg_ssac_mod_pred = squeeze(mean(boot_avgs))';
std_ssac_mod_pred = squeeze(std(boot_avgs))';

% phasemod_out_norm2 = bsxfun(@rdivide,phasemod_out2,travg_spk_rates');
% [sac_trg_mat,cur_lags] = get_event_trig_mat(phasemod_out_norm2',tr_bsac_inds,round(0.4/new_dt),round(0.5/new_dt));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_bsac_inds),24*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 24]);
% avg_bsac_mod_pred2 = squeeze(mean(boot_avgs))';
% std_bsac_mod_pred2 = squeeze(std(boot_avgs))';

disp('Computing micro sac stats')
n_boot = 100;
[sac_trg_mat,cur_lags] = get_event_trig_mat(norm_spkbinnedtr(tr_inds,:),tr_msac_inds,round(0.4/new_dt),round(0.5/new_dt));
sac_trg_mat = reshape(sac_trg_mat,length(tr_msac_inds),24*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 24]);
avg_msac_trig_rate = squeeze(mean(boot_avgs))';
std_msac_trig_rate = squeeze(std(boot_avgs))';

[sac_trg_mat,cur_lags] = get_event_trig_mat(phasemod_out_norm',tr_msac_inds,round(0.4/new_dt),round(0.5/new_dt));
sac_trg_mat = reshape(sac_trg_mat,length(tr_msac_inds),24*length(cur_lags));
boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 24]);
avg_msac_mod_pred = squeeze(mean(boot_avgs))';
std_msac_mod_pred = squeeze(std(boot_avgs))';

% [sac_trg_mat,cur_lags] = get_event_trig_mat(phasemod_out_norm2',tr_msac_inds,round(0.4/new_dt),round(0.5/new_dt));
% sac_trg_mat = reshape(sac_trg_mat,length(tr_msac_inds),24*length(cur_lags));
% boot_avgs = bootstrp(n_boot,@mean,sac_trg_mat);
% boot_avgs = reshape(boot_avgs,[n_boot length(cur_lags) 24]);
% avg_msac_mod_pred2 = squeeze(mean(boot_avgs))';
% std_msac_mod_pred2 = squeeze(std(boot_avgs))';

%%
cd /home/james/Desktop
for cc = 1:24
    cc
    subplot(2,2,1)
    shadedErrorBar(cur_lags*new_dt,avg_bsac_trig_rate(cc,:),std_bsac_trig_rate(cc,:));
    hold on
    shadedErrorBar(cur_lags*new_dt,avg_bsac_mod_pred(cc,:),std_bsac_mod_pred(cc,:),{'r'});
%     shadedErrorBar(cur_lags*new_dt,avg_bsac_mod_pred2(cc,:),std_bsac_mod_pred2(cc,:),{'b'});
    xlim([-0.25 0.35])
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k')
    line([0 0],yl,'color','k')
    set(gca,'fontsize',14)
    xlabel('Time since saccade onset (s)','fontsize',16)
    ylabel('Relative firing rate','fontsize',16)
    title('Guided Saccades','fontsize',16)
    box off
    
     subplot(2,2,2)
    shadedErrorBar(cur_lags*new_dt,avg_msac_trig_rate(cc,:),std_msac_trig_rate(cc,:));
    hold on
    shadedErrorBar(cur_lags*new_dt,avg_msac_mod_pred(cc,:),std_msac_mod_pred(cc,:),{'r'});
%     shadedErrorBar(cur_lags*new_dt,avg_msac_mod_pred2(cc,:),std_msac_mod_pred2(cc,:),{'b'});
    xlim([-0.25 0.35])
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k')
    line([0 0],yl,'color','k')
    xlabel('Time since saccade onset (s)','fontsize',16)
    ylabel('Relative firing rate','fontsize',16)
    title('Microsaccades','fontsize',16)
    box off

         subplot(2,2,3)
    shadedErrorBar(cur_lags*new_dt,avg_fsac_trig_rate(cc,:),std_msac_trig_rate(cc,:));
    hold on
    shadedErrorBar(cur_lags*new_dt,avg_fsac_mod_pred(cc,:),std_msac_mod_pred(cc,:),{'r'});
    xlim([-0.25 0.35])
    xl = xlim(); 
    yl = ylim();
    line(xl,[1 1],'color','k')
    line([0 0],yl,'color','k')
    xlabel('Time since saccade onset (s)','fontsize',16)
    ylabel('Relative firing rate','fontsize',16)
    title('Outgoing sac','fontsize',16)
    box off

         subplot(2,2,4)
    shadedErrorBar(cur_lags*new_dt,avg_ssac_trig_rate(cc,:),std_msac_trig_rate(cc,:));
    hold on
    shadedErrorBar(cur_lags*new_dt,avg_ssac_mod_pred(cc,:),std_msac_mod_pred(cc,:),{'r'});
    xlim([-0.25 0.35])
    xl = xlim(); 
%     yl = ylim();
    line(xl,[1 1],'color','k')
    line([0 0],yl,'color','k')
    xlabel('Time since saccade onset (s)','fontsize',16)
    ylabel('Relative firing rate','fontsize',16)
    title('Returning sac','fontsize',16)
    box off

%     fillPage(gcf,'Papersize',[10 10]);
%     pname = sprintf('Obs_vs_pred_sacmod_unit%d',cc);
%     print(pname,'-dpdf');
%     close
    
    
    pause
    clf
end




%%
close all
for cc = 1:24
    subplot(3,2,1);hold on
%     plot(tent_centers*new_dt,exp(so_sac_kern(cc,:)),'r')
    plot(cur_lags*new_dt,(st_avg_bphasemod_post(cc,:)/travg_spk_rates(cc)),'r','linewidth',1)
%     plot(cur_lags*new_dt,zscore(xvst_avg_phasemod_po(cc,:)),'b','linewidth',1)
    plot(cur_lags*new_dt,(bsac_trig_avg_rate(cc,:)),'k')
%     plot(cur_lags*new_dt,zscore(xvsac_trig_avg_rate(cc,:)),'color',[0.2 0.8 0.2])
    axis tight
    xlim([-0.4 0.5])
    yl = ylim(); xl = xlim();
    line([0 0],yl,'color','k')
    line(xl,[0 0],'color','k')
    xlabel('Time relative to saccade (s)')
    
    subplot(3,2,2);hold on
%     plot(tent_centers*new_dt,exp(so_msac_kern(cc,:)),'r')
%    plot(cur_lags*new_dt,zscore(xvst_avg_mphasemod_po(cc,:)),'b','linewidth',1)
    plot(cur_lags*new_dt,(st_avg_mphasemod_post(cc,:)/travg_spk_rates(cc)),'r','linewidth',1)
    plot(cur_lags*new_dt,(msac_trig_avg_rate(cc,:)),'k')
%     plot(cur_lags*new_dt,zscore(xvmsac_trig_avg_rate(cc,:)),'color',[0.2 0.8 0.2])
    axis tight
    xlim([-0.4 0.5])
    yl = ylim(); xl = xlim();
    line([0 0],yl,'color','k')
    line(xl,[0 0],'color','k')
    xlabel('Time relative to saccade (s)')

    
        subplot(3,2,3);hold on
%     plot(tent_centers*new_dt,exp(so_sac_kern(cc,:)),'r')
%     plot(cur_lags*new_dt,zscore(st_avg_phasemod_po(cc,:)),'r','linewidth',1)
    plot(cur_lags*new_dt,zscore(xvst_avg_bphasemod_post(cc,:)),'b','linewidth',1)
%     plot(cur_lags*new_dt,zscore(sac_trig_avg_rate(cc,:)),'k')
    plot(cur_lags*new_dt,zscore(xvbsac_trig_avg_rate(cc,:)),'color',[0.2 0.8 0.2])
    axis tight
    xlim([-0.4 0.5])
    yl = ylim(); xl = xlim();
    line([0 0],yl,'color','k')
    line(xl,[0 0],'color','k')
    xlabel('Time relative to saccade (s)')
    
    subplot(3,2,4);hold on
%     plot(tent_centers*new_dt,exp(so_msac_kern(cc,:)),'r')
   plot(cur_lags*new_dt,zscore(xvst_avg_mphasemod_post(cc,:)),'b','linewidth',1)
%     plot(cur_lags*new_dt,zscore(st_avg_mphasemod_po(cc,:)),'r','linewidth',1)
%     plot(cur_lags*new_dt,zscore(msac_trig_avg_rate(cc,:)),'k')
    plot(cur_lags*new_dt,zscore(xvmsac_trig_avg_rate(cc,:)),'color',[0.2 0.8 0.2])
    axis tight
    xlim([-0.4 0.5])
    yl = ylim(); xl = xlim();
    line([0 0],yl,'color','k')
    line(xl,[0 0],'color','k')
    xlabel('Time relative to saccade (s)')

        subplot(3,2,5);hold on
%     plot(tent_centers*new_dt,exp(so_sac_kern(cc,:)),'r')
%     plot(cur_lags*new_dt,zscore(st_avg_phasemod_po(cc,:)),'r','linewidth',1)
%     plot(cur_lags*new_dt,zscore(xvst_avg_phasemod_po(cc,:)),'b','linewidth',1)
    plot(cur_lags*new_dt,zscore(sac_trig_avg_rate(cc,:)),'k')
    plot(cur_lags*new_dt,zscore(xvsac_trig_avg_rate(cc,:)),'color',[0.2 0.8 0.2])
    axis tight
    xlim([-0.4 0.5])
    yl = ylim(); xl = xlim();
    line([0 0],yl,'color','k')
    line(xl,[0 0],'color','k')
    xlabel('Time relative to saccade (s)')
    
    subplot(3,2,6);hold on
%     plot(tent_centers*new_dt,exp(so_msac_kern(cc,:)),'r')
%    plot(cur_lags*new_dt,zscore(xvst_avg_mphasemod_po(cc,:)),'b','linewidth',1)
%     plot(cur_lags*new_dt,zscore(st_avg_mphasemod_po(cc,:)),'r','linewidth',1)
    plot(cur_lags*new_dt,zscore(msac_trig_avg_rate(cc,:)),'k')
    plot(cur_lags*new_dt,zscore(xvmsac_trig_avg_rate(cc,:)),'color',[0.2 0.8 0.2])
    axis tight
    xlim([-0.4 0.5])
    yl = ylim(); xl = xlim();
    line([0 0],yl,'color','k')
    line(xl,[0 0],'color','k')
    xlabel('Time relative to saccade (s)')

    %         
%      subplot(2,2,3);hold on
% %     plot(tent_centers*new_dt,exp(so_fsac_kern(cc,:)),'r')
%     plot(cur_lags*new_dt,(st_avg_fphasemod_po(cc,:)),'r','linewidth',2)
%     plot(cur_lags*new_dt,(xvst_avg_fphasemod_po(cc,:)),'b','linewidth',2)
%     plot(cur_lags*new_dt,fsac_trig_avg_rate(cc,:),'k')
%     plot(cur_lags*new_dt,xvfsac_trig_avg_rate(cc,:),'color',[0.2 0.8 0.2])
%     axis tight
%     xlim([-0.4 0.5])
%     yl = ylim(); xl = xlim();
%     line([0 0],yl,'color','k')
%     line(xl,[0 0],'color','k')
%     xlabel('Time relative to saccade (s)')
% 
%      subplot(2,2,4);hold on
% %     plot(tent_centers*new_dt,exp(so_ssac_kern(cc,:)),'r')
%     plot(cur_lags*new_dt,(st_avg_sphasemod_po(cc,:)),'r','linewidth',2)
%     plot(cur_lags*new_dt,(xvst_avg_sphasemod_po(cc,:)),'b','linewidth',2)
%     plot(cur_lags*new_dt,ssac_trig_avg_rate(cc,:),'k')
%     plot(cur_lags*new_dt,xvssac_trig_avg_rate(cc,:),'color',[0.2 0.8 0.2])
%     axis tight
%     xlim([-0.4 0.5])
%     yl = ylim(); xl = xlim();
%     line([0 0],yl,'color','k')
%     line(xl,[0 0],'color','k')
%     xlabel('Time relative to saccade (s)')

    cc
    pause
    clf
end

%%
close all
for cc = 1:24
    subplot(1,2,1);hold on
    plot(cur_lags*new_dt,(st_avg_bphasemod_post(cc,:))/new_dt,'r','linewidth',1)
    plot(cur_lags*new_dt,(bsac_trig_avg_rate(cc,:))/new_dt,'k')
    axis tight
    xlim([-0.3 0.4])
    yl = ylim(); xl = xlim();
    line([0 0],yl,'color','k')
%     line(xl,[0 0],'color','k')
    xlabel('Time relative to saccade (s)')
    ylabel('Firing rate (Hz)','fontsize',16)
%     legend('Phase model prediction','Observed')
title('Guided saccades')
        line(xl,[avg_rate(cc) avg_rate(cc)]/new_dt,'color','k')
    
    subplot(1,2,2);hold on
    plot(cur_lags*new_dt,(st_avg_mphasemod_post(cc,:))/new_dt,'r','linewidth',1)
    plot(cur_lags*new_dt,(msac_trig_avg_rate(cc,:))/new_dt,'k')
%     plot(cur_lags*new_dt,zscore(xvmsac_trig_avg_rate(cc,:)),'color',[0.2 0.8 0.2])
    axis tight
    xlim([-0.3 0.4])
    yl = ylim(); xl = xlim();
    line([0 0],yl,'color','k')
%     line(xl,[0 0],'color','k')
    xlabel('Time relative to saccade (s)')
    ylabel('Firing rate (Hz)','fontsize',16)
%     legend('Phase model prediction','Observed')
title('Microsaccades')
        line(xl,[avg_rate(cc) avg_rate(cc)]/new_dt,'color','k')

    cc
    pause
    clf
end

%%
avg_rates = mean(full_spkbinned(:,:));
smooth_sig = 2;
close all
for cc = 1:24
    cc
    subplot(3,1,1); hold on
    plot(tent_centers*new_dt,so_sac_kern(cc,:),'b')
%     plot(tent_centers*new_dt,so_msac_kern(cc,:),'r')
%     plot(tent_centers*new_dt,so_fsac_kern(cc,:),'k')
%     plot(tent_centers*new_dt,so_ssac_kern(cc,:),'c')
    grid on
    subplot(3,1,2); hold on
    plot(cur_lags*new_dt,jmm_smooth_1d_cor(sac_trig_avg_rate(cc,:),smooth_sig)/avg_rates(cc),'b')
%     plot(cur_lags*new_dt,jmm_smooth_1d_cor(msac_trig_avg_rate(cc,:),smooth_sig)/avg_rates(cc),'r')
%     plot(cur_lags*new_dt,jmm_smooth_1d_cor(fsac_trig_avg_rate(cc,:),smooth_sig)/avg_rates(cc),'k')
%     plot(cur_lags*new_dt,jmm_smooth_1d_cor(ssac_trig_avg_rate(cc,:),smooth_sig)/avg_rates(cc),'c')
    grid on
    subplot(3,1,3); hold on
%      plot(cur_lags*new_dt,st_avg_phasemod_po(cc,:),'b')
%      plot(cur_lags*new_dt,st_avg_mphasemod_po(cc,:),'r')
%      plot(cur_lags*new_dt,st_avg_fphasemod_po(cc,:),'k')
%      plot(cur_lags*new_dt,st_avg_sphasemod_po(cc,:),'c')
     plot(cur_lags*new_dt,st_avg_phasemod_po(cc,:),'b')
%      plot(cur_lags*new_dt,exp(st_avg_mphasemod_po(cc,:)),'r')
%      plot(cur_lags*new_dt,exp(st_avg_fphasemod_po(cc,:)),'k')
%      plot(cur_lags*new_dt,exp(st_avg_sphasemod_po(cc,:)),'c')
     grid on
  
    pause
    clf
end

%%
use_bar_pos = full_bar_pos; 
use_bar_mat = zeros(length(full_taxis_new),n_bar_pos);
for b = 1:n_bar_pos
    cur_set = find(use_bar_pos == un_bar_pos(b));
    use_bar_mat(cur_set,b) = 1;
end
use_bar_pos(use_bar_pos < -10) = nan;

use_t_axis = (1:length(full_taxis_new))*new_dt;
%%
examp_cell = 17;
close all

Robs = full_spkbinned(:,examp_cell);
tr_spkbns = convert_to_spikebins(Robs);

glm_kern = get_k_mat(glm_fit(examp_cell));
stim_out = full_Xmat*glm_kern;

use_elecs = 1:24;
%     if mod(cc,2)==0
%         use_elecs = 2:2:24;
%     else
%         use_elecs = 1:2:24;
%     end
use_set = find(ismember(phase_elec_set,use_elecs));
use_set1 = find(ismember(phase_elec_set,use_elecs) & ismember(phase_freq_set,gamma_freqs));
use_set2 = find(ismember(phase_elec_set,use_elecs) & ismember(phase_freq_set,alpha_freqs));
phasemod_out = new_phase_set(:,use_set)*fitp_post(examp_cell).k(1:2*length(use_elecs)*length(wfreqs));
% gammod_out = new_phase_set(:,use_set1)*fitp_postgam(examp_cell).k(1:2*length(use_elecs)*length(gamma_freqs));
% alpmod_out = new_phase_set(:,use_set2)*fitp_postalp(examp_cell).k(1:2*length(use_elecs)*length(alpha_freqs));
gammod_out = new_phase_set(:,use_set1)*fitp_post(examp_cell).k(use_set1);
alpmod_out = new_phase_set(:,use_set2)*fitp_post(examp_cell).k(use_set2);
so_out = stim_out*fitp_so(examp_cell).k(1);

%  po_out = new_phase_set(:,use_set)*fitp_po(examp_cell).k(1:2*length(use_elecs)*length(wfreqs));

cur_Xmat = [new_phase_set(:,use_set) stim_out lin_X];
xv_fphase_pred_rate = cur_Xmat*fitp_post(examp_cell).k(1:end-1) + fitp_post(examp_cell).k(end);
if NL_type == 0
    xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
else
    xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
end
cur_Xmat = [stim_out lin_X];
xv_so_pred_rate = cur_Xmat*fitp_so(examp_cell).k(1:end-1) + fitp_so(examp_cell).k(end);
if NL_type == 0
    xv_so_pred_rate = log(1+exp(xv_so_pred_rate));
else
    xv_so_pred_rate = exp(xv_so_pred_rate);
end
sc_fac = fitp_post(examp_cell).k(2*length(use_elecs)*length(wfreqs)+1);

[gam_b,gam_a] = butter(2,[25 50]/(Fsd/2));
[alp_b,alp_a] = butter(2,[5 15]/(Fsd/2));
lfp_gam = filtfilt(gam_b,gam_a,full_lfps(:,examp_cell));
lfp_alp = filtfilt(alp_b,alp_a,full_lfps(:,examp_cell));

tick_y = [-20 -5];
f1 = figure;
for i = 246:n_trials
    cur_set = find(ic == i);
    
    subplot(5,1,[1])
    hold on
    plot(use_t_axis(cur_set),full_lfps(cur_set,examp_cell),'k')
    hold on
    plot(use_t_axis(cur_set),lfp_gam(cur_set),'r')
    plot(use_t_axis(cur_set),lfp_alp(cur_set),'b')
    xlabel('Time (s)','fontsize',16)
    xlim(use_t_axis(cur_set([1 end])))
    ylabel('LFP amplitude (z)')
    subplot(5,1,[2])
    hold on
    plot(use_t_axis(cur_set),phasemod_out(cur_set),'k')
    plot(use_t_axis(cur_set),gammod_out(cur_set),'r')
    plot(use_t_axis(cur_set),alpmod_out(cur_set),'b')
%     plot(use_t_axis(cur_set),sc_fac*stim_out(cur_set)+4.5,'b')
%     plot(use_t_axis(cur_set),so_out(cur_set),'g')
    xlabel('Time (s)','fontsize',16)
    xlim(use_t_axis(cur_set([1 end])))
%     legend('Phase model output','Stimulus model output')
    
    subplot(5,1,[3])
    plot(use_t_axis(cur_set),xv_fphase_pred_rate(cur_set)/new_dt,'r')
    hold on
    plot(use_t_axis(cur_set),xv_so_pred_rate(cur_set)/new_dt,'b')
    legend('Phase model prediction','Stimulus model prediction')
    cur_spks = tr_spkbns(ismember(tr_spkbns,cur_set));
    for ss = 1:length(cur_spks)
        line(use_t_axis(cur_spks([ss ss])),tick_y,'color','k')
    end
    % cur_sac_out = trial_bsac_mat*so_bsac_kern(examp_cell,:)' + trial_msac_mat*so_msac_kern(examp_cell,:)';
    % plot(use_t_axis,cur_sac_out,'c')
    xlim(use_t_axis(cur_set([1 end])))
    yl = ylim();
    ylim([-25 yl(2)])
    xl  =xlim();
    line(xl,[0 0],'color','k')
    xlabel('Time (s)','fontsize',16)
    ylabel('Firing rate (Hz)','fontsize',16)
    
    
    subplot(5,1,4)
    imagesc(use_t_axis(cur_set),un_bar_pos,use_bar_mat(cur_set,:)');colormap(gray);
    xlim(use_t_axis(cur_set([1 end])))
    xlabel('Time (s)','fontsize',16)
    ylabel('Bar position (deg)','fontsize',16)

    subplot(5,1,5)
    plot(use_t_axis(cur_set),full_eye_speed(cur_set),'k')
    xlim(use_t_axis(cur_set([1 end])))
     xlabel('Time (s)','fontsize',16)
    ylabel('Eye speed (deg/s)','fontsize',16)
   
    i
    pause;
    clf
end