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
dsf = 3;
Fsd = Fs/dsf;
niqf = Fs/2;
[bb,aa] = butter(2,120/niqf,'low');
% scales = logspace(log10(5),log10(60),30);
% scales = [scales 70 80 90];
% wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
% nwfreqs = length(wfreqs);

scales = logspace(log10(4),log10(60),35);
scales = [scales 70 80 90];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

new_dt = .0025;
dsfrac = 4;
%% COMPUTE PERPINDICULAR EYE POSITION
for ee = 1:length(bar_expts)
    
    avg_eyepos = 0.5*all_eye_vals_cal{ee}(:,1:2) + 0.5*all_eye_vals_cal{ee}(:,3:4);
    
    fix_r = sqrt(sum(avg_eyepos.^2,2));
    fix_angle = atan2(avg_eyepos(:,2),avg_eyepos(:,1));
    angle_diff = fix_angle - axis_or;
    fix_perp{ee} = fix_r.*sin(angle_diff);
    
end

%%
load ./CellList.mat
good_sus = find(all(CellList(bar_expts,:,1) > 0));
% use_sus = find(all(CellList(bar_expts,:,1) > 0));
use_sus = 1:24;
% used_expts{1} = 1:6;
% used_expts{2} = 1:5;
% used_expts{3} = 1:6;
% used_expts{4} =
%%
load ./un_bar_pos.mat
un_bar_pos(1:2) = [];

flen = 20;
full_Xmat = [];
full_spkbinned = [];
full_exptvec = [];
full_exptvec_new = [];
full_trialvec = [];
full_taxis = [];
full_taxis_new = [];
full_old_t_inds = [];
full_bar_pos = [];
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
    all_binned_spikes = [];
    all_used_inds = [];
    for tt = 1:n_trials(ee)
        
        cur_t_axis = [Expts{bar_expts(ee)}.Trials(tt).Start]/1e4;
        cur_bar_Op = [Expts{bar_expts(ee)}.Trials(tt).Op];
        %         cur_bar_e0 = [Expts{bar_expts(ee)}.Trials(tt).e0];
        
        cur_t_edges = [cur_t_axis; Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4];
        cur_t_edges_new = cur_t_edges(1):new_dt:cur_t_edges(end);
        cur_t_axis_new = cur_t_edges_new(1:end-1);
        
        cur_binned_spks = nan(length(cur_t_axis_new),length(use_sus));
        for cc = 1:length(use_sus)
            cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges_new);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        
        cur_used_inds = ones(length(cur_t_axis),1);
        cur_used_inds(1:flen) = 0;
        cur_used_inds_new = ones(length(cur_t_axis_new),1);
        
        cur_used_inds_new(end-round(0.1/new_dt):end) = 0; %reduce phase estimation edge artifacts
        
        all_used_inds_new = [all_used_inds_new; cur_used_inds_new(:)];
        all_t_axis = [all_t_axis; cur_t_axis(:)];
        all_used_inds = [all_used_inds; cur_used_inds(:)];
        all_t_axis_new = [all_t_axis_new; cur_t_axis_new(:)];
        %         all_old_t_inds = [all_old_t_inds; old_t_inds(:)];
        all_binned_spikes = [all_binned_spikes; cur_binned_spks];
        all_expt_vec = [all_expt_vec; ones(length(cur_t_axis_new),1)*bar_expts(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_t_axis),1)*tt];
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
    
    all_used_inds(all_blink_inds == 1) = 0;
    all_used_inds = logical(all_used_inds);
    all_used_inds_new(all_blink_inds_new == 1) = 0;
    all_used_inds_new = logical(all_used_inds_new);
    
    full_spkbinned = [full_spkbinned; all_binned_spikes(all_used_inds_new,:)];
    full_exptvec = [full_exptvec; ones(sum(all_used_inds),1)*ee];
    full_exptvec_new = [full_exptvec_new; ones(sum(all_used_inds_new),1)*ee];
    full_taxis = [full_taxis; all_t_axis(all_used_inds)];
    full_taxis_new = [full_taxis_new; all_t_axis_new(all_used_inds_new)];
    %     full_old_t_inds = [full_old_t_inds; all_old_t_inds(all_used_inds_new)];
    full_trialvec = [full_trialvec; all_trial_vec(all_used_inds)];
end

%%
full_lfps = [];
full_phasegrams = [];
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
        cur_LFP = filtfilt(bb,aa,cur_LFP);
        
        cur_LFP = cur_LFP(cur_sp:end,:);
        
        
        expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
        expt_lfps = [expt_lfps; cur_LFP];
    end
    
    seg_start_inds = [1; 1+find(diff(expt_lfp_t_axis) > 2/Fs)];
    seg_stop_inds = [find(diff(expt_lfp_t_axis) > 2/Fs); length(expt_lfp_t_axis)];
    res_t_axis = [];
    res_phasegram = [];
    res_lfps = [];
    for jj = 1:length(seg_start_inds)
        cur_res_inds = seg_start_inds(jj):seg_stop_inds(jj);
        cur_res_t = expt_lfp_t_axis(seg_start_inds(jj)):dsf/Fs:expt_lfp_t_axis(seg_stop_inds(jj));
        res_t_axis = [res_t_axis; cur_res_t(:)];
        
        interp_lfps = interp1(expt_lfp_t_axis(cur_res_inds),expt_lfps(cur_res_inds,:),cur_res_t);
        
        cur_phasegram = nan(length(cur_res_t),length(wfreqs),24);
        for ll = 1:24
            temp = cwt(interp_lfps(:,ll),scales,'cmor1-1');
            cur_phasegram(:,:,ll) = angle(temp)';
        end
        res_phasegram = cat(1,res_phasegram,cur_phasegram);
        res_lfps = [res_lfps; interp_lfps];
    end
    
    cur_set = find(full_exptvec_new==ee);
    interp_lfps = interp1(res_t_axis,res_lfps,full_taxis_new(cur_set));
    unwr_phasegram = unwrap(res_phasegram);
    interp_phasegrams = interp1(res_t_axis,unwr_phasegram,full_taxis_new(cur_set));
    interp_phasegrams = mod(interp_phasegrams+pi,2*pi)-pi;
    %     interp_phasegrams = interp1(res_t_axis,res_phasegram,full_taxis_new(cur_set));
    
    full_lfps = [full_lfps; interp_lfps];
    full_phasegrams = cat(1,full_phasegrams, interp_phasegrams);
end

%% create mapping to old taxis (with stimulus stored)
% full_old_inds = [];
% ind_cnt = 0;
% for ee = 1:length(bar_expts)
%     cur_set = find(full_exptvec==ee);
%     cur_set_new = find(full_exptvec_new==ee);
%     cur_inds = round(interp1(full_taxis(cur_set),1:length(cur_set),full_taxis_new(cur_set_new)));
%     full_old_inds = [full_old_inds; cur_inds + ind_cnt];
%     ind_cnt = ind_cnt + length(cur_set);
% end

%%
for cc = 1:24
    cur_spk_set = convert_to_spikebins(full_spkbinned(:,cc));
    phase_set{cc} = full_phasegrams(cur_spk_set,:,:);
    %%
   phase_locking(cc,:,:) = squeeze(nansum(exp(1i*phase_set{cc})));
   phase_locking(cc,:,:) = abs(phase_locking(cc,:,:))/sum(~isnan(phase_set{cc}(:,1,1)));


end

%%
avg_phaselocking = squeeze(mean(phase_locking,3));

use_probes = (1:24)+1;
use_probes(use_probes > 24) = 23;
best_phaselocking = nan(size(avg_phaselocking));
for c = 1:24
    best_phaselocking(c,:) = squeeze(phase_locking(c,:,use_probes(c)));
end
%%
save random_bar_phaselocking_data wfreqs phase_locking best_phaselocking avg_phaselocking

%%
close all
for cc = 1:24
    pcolor(wfreqs,1:24,squeeze(phase_locking(cc,:,:))');shading flat
    cc
    pause
    clf
end