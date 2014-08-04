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
scales = logspace(log10(5),log10(60),30);
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
    all_bar_e0 = [];
    all_bar_Xmat = [];
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
        
%         cur_interp_eye_perp = interp1(all_eye_ts{ee},fix_perp{ee},cur_t_axis);
%         cur_bar_Op_cor = cur_bar_Op - cur_interp_eye_perp;
        
        cur_bar_mat = zeros(length(cur_t_axis),length(un_bar_pos));
        for bb = 1:length(un_bar_pos)
            cur_set = find(cur_bar_Op==un_bar_pos(bb));
            cur_bar_mat(cur_set,bb) = 1;
%             cur_bar_mat(cur_set,bb) = cur_bar_e0(cur_set);
        end
        bar_Xmat = makeStimRows(cur_bar_mat,flen);
        
        cur_used_inds = ones(length(cur_t_axis),1);
        cur_used_inds(1:flen) = 0;
        cur_used_inds_new = ones(length(cur_t_axis_new),1);
        cur_used_inds_new(1:(flen*dsfrac)) = 0;
        
        cur_used_inds_new(end-round(0.1/new_dt):end) = 0; %reduce phase estimation edge artifacts
        
        all_used_inds_new = [all_used_inds_new; cur_used_inds_new(:)];
        all_t_axis = [all_t_axis; cur_t_axis(:)];
        all_used_inds = [all_used_inds; cur_used_inds(:)];
        all_t_axis_new = [all_t_axis_new; cur_t_axis_new(:)];
%         all_old_t_inds = [all_old_t_inds; old_t_inds(:)];
        all_binned_spikes = [all_binned_spikes; cur_binned_spks];
        all_expt_vec = [all_expt_vec; ones(length(cur_t_axis_new),1)*bar_expts(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_t_axis),1)*tt];
        all_bar_Op = [all_bar_Op; cur_bar_Op(:)];
%         all_bar_e0 = [all_bar_e0; cur_bar_e0(:)];
        all_bar_Xmat = [all_bar_Xmat; bar_Xmat];
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
        
    full_Xmat = [full_Xmat; all_bar_Xmat(all_used_inds,:)];
    full_spkbinned = [full_spkbinned; all_binned_spikes(all_used_inds_new,:)];
    full_exptvec = [full_exptvec; ones(sum(all_used_inds),1)*ee];
    full_exptvec_new = [full_exptvec_new; ones(sum(all_used_inds_new),1)*ee];
    full_taxis = [full_taxis; all_t_axis(all_used_inds)];
    full_taxis_new = [full_taxis_new; all_t_axis_new(all_used_inds_new)];
%     full_old_t_inds = [full_old_t_inds; all_old_t_inds(all_used_inds_new)];
    full_trialvec = [full_trialvec; all_trial_vec(all_used_inds)];
    full_bar_pos = [full_bar_pos; all_bar_Op(all_used_inds)];
end

%%
all_fix_inds = [];
all_first_inds = [];
all_second_inds = [];
all_small_inds = [];
all_tstart_inds = [];
all_tstop_inds = [];
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
    micro_sac_set = find(sac_amps < 50);
    
%     fix_times = sac_stop_times(intrial_sac_set);
%     first_fix_times = sac_stop_times(infirst_sac_set);
%     second_fix_times = sac_stop_times(insecond_sac_set);
    fix_times = sac_start_times(intrial_sac_set);
    first_fix_times = sac_start_times(infirst_sac_set);
    second_fix_times = sac_start_times(insecond_sac_set);
    small_fix_times = sac_start_times(micro_sac_set);
    
    cur_e_set = find(full_exptvec_new == ee);
    
    fix_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),fix_times));
    first_fix_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),first_fix_times));
    second_fix_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),second_fix_times));
    small_fix_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),small_fix_times));

    trial_start_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),expt_trial_starts));
    trial_stop_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),expt_trial_end));
    
    
    fix_inds(isnan(fix_inds)) = [];
    first_fix_inds(isnan(first_fix_inds)) = [];
    second_fix_inds(isnan(second_fix_inds)) = [];
    small_fix_inds(isnan(small_fix_inds)) = [];
    rid = find(isnan(trial_start_inds) | isnan(trial_stop_inds));
    trial_start_inds(rid) = [];
    trial_stop_inds(rid) = [];
    
    all_fix_inds = [all_fix_inds; cur_e_set(fix_inds(:))];
    all_first_inds = [all_first_inds; cur_e_set(first_fix_inds(:))];
    all_second_inds = [all_second_inds; cur_e_set(second_fix_inds(:))];
    all_small_inds = [all_small_inds; cur_e_set(small_fix_inds(:))];
    all_tstart_inds = [all_tstart_inds; cur_e_set(trial_start_inds(:))];
    all_tstop_inds = [all_tstop_inds; cur_e_set(trial_stop_inds(:))];
end

used_sac_inds = sort([all_first_inds; all_second_inds]);
% used_sac_inds = sort(all_fix_inds);
% small_set = setdiff(all_fix_inds,used_sac_inds);

%%
cur_dt = 0.01;
% flen_t = 0.5;
tent_centers = [0:cur_dt:0.8];
tent_centers = round(tent_centers/new_dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

shift = round(0.3/new_dt);
tbmat = [zeros(ntents,tblen-2*shift-1) tbmat];
tent_centers = tent_centers-shift;

trial_inds = zeros(size(full_taxis_new));
trial_inds(used_sac_inds) = 1;
trial_Tmat = zeros(length(full_taxis_new),ntents);
for i = 1:ntents
    trial_Tmat(:,i) = conv(trial_inds,tbmat(i,:),'same');
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
full_old_inds = [];
ind_cnt = 0;
for ee = 1:length(bar_expts)
   cur_set = find(full_exptvec==ee);
   cur_set_new = find(full_exptvec_new==ee);
   cur_inds = round(interp1(full_taxis(cur_set),1:length(cur_set),full_taxis_new(cur_set_new)));
   full_old_inds = [full_old_inds; cur_inds + ind_cnt];
   ind_cnt = ind_cnt + length(cur_set);
end

stim_Xmat = full_Xmat;
%%

load ./random_bar_unit_models_ftime.mat
% load ./bar_phase_models_ftime_v2
load ./random_bar_phase_models_v2
new_phase_set = [reshape(cos(full_phasegrams),length(full_taxis_new),length(wfreqs)*24) reshape(sin(full_phasegrams),length(full_taxis_new),length(wfreqs)*24)];
% phase_elec_set = [repmat(1:24,1,length(wfreqs)) repmat(1:24,1,length(wfreqs))];
% phase_elec_set = phase_elec_set(:);
phase_elec_set = ones(length(wfreqs),1)*(1:24);
phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];
NT = length(full_old_inds);

%%
for cur_cell = 1:24;
cur_cell
Robs = full_spkbinned(:,cur_cell);
% if cc < 24
%     nearest_probe = cc+1;
% else
%     nearest_probe = cc-1;
% end
nearest_probe = cc;
cur_alpha_phase = squeeze(full_phasegrams(:,:,nearest_probe));
Pmat = nan(NT,length(wfreqs)*nbins);
for ww = 1:length(wfreqs)
    cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
    cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
    cur_tb(:,end) = [];
    cur_set = ((ww-1)*nbins+1):ww*nbins;
    Pmat(:,cur_set) = cur_tb;
end

% se_phasemod_out(cur_cell,:) = Pmat*stim_ind_phase_pfilt(cur_cell,:)';
se_phasemod_out(cur_cell,:) = Pmat*phase_pfilt(cur_cell,:)';
se_phasemod_out(cur_cell,:) = se_phasemod_out(cur_cell,:) - mean(se_phasemod_out(cur_cell,:));   

    use_elecs = 1:24;
%     use_elecs(cc) = [];
    use_set = find(ismember(phase_elec_set,use_elecs)); 
% full_phase_filt = [stim_ind_sinphase_cfilt(cur_cell,:) stim_ind_sinphase_sfilt(cur_cell,:)];
full_phase_filt = [sinphase_cfilt(cur_cell,:) sinphase_sfilt(cur_cell,:)];
full_phasemod_out(cur_cell,:) = new_phase_set(:,use_set)*full_phase_filt'; 

full_phasemod_out(cur_cell,:) = full_phasemod_out(cur_cell,:) - mean(full_phasemod_out(cur_cell,:));

% glm_kern = get_k_mat(glm_fit(cur_cell));
% stim_out = stim_Xmat*glm_kern;
% % stim_out_interp(cur_cell,:) = stim_out(full_old_inds(tr_inds_new));
% 
rate_sm = 1;
sm_rate(cur_cell,:) = jmm_smooth_1d_cor(Robs,rate_sm);
end
%%
for i = 1:length(all_tstart_inds)
    cur_set = all_tstart_inds(i):all_tstop_inds(i);
    cur_msacs = all_fix_inds(all_fix_inds >= cur_set(1) & all_fix_inds < cur_set(end));
    cur_sacs = used_sac_inds(used_sac_inds >= cur_set(1) & used_sac_inds < cur_set(end));
    subplot(2,1,1)
    imagesc(full_taxis_new(cur_set),1:24,full_phasemod_out(:,cur_set));
    for j = 1:length(cur_msacs)
        line(full_taxis_new(cur_msacs([j j])),[1 24],'color','r','linewidth',2)
    end
    for j = 1:length(cur_sacs)
        line(full_taxis_new(cur_sacs([j j])),[1 24],'color','w','linewidth',2)
    end
    caxis([-2 2])
    subplot(2,1,2)
    imagesc(full_taxis_new(cur_set),1:24,se_phasemod_out(:,cur_set));
    for j = 1:length(cur_msacs)
        line(full_taxis_new(cur_msacs([j j])),[1 24],'color','r','linewidth',2)
    end
    for j = 1:length(cur_sacs)
        line(full_taxis_new(cur_sacs([j j])),[1 24],'color','w','linewidth',2)
    end
    caxis([-2 2])
    pause
    clf
end


% expt = 1;
% cur_set = find(full_exptvec_new==expt);

% figure
% plot(full_taxis_new(cur_set),zscore(full_lfps(cur_set,nearest_probe)))
% hold on
% plot(full_taxis_new(cur_set),se_phasemod_out(cur_set),'r')
% plot(full_taxis_new(cur_set),full_phasemod_out(cur_set)-4,'k')
% plot(full_taxis_new(cur_set),stim_out_interp(cur_set),'g','linewidth',2)
% 
% plot(full_taxis_new(cur_set),2+zscore(sm_rate(cur_set)),'k')

%%
backlag = round(0.1/new_dt);
forwardlag = round(0.4/new_dt);

used_msacs = all_fix_inds;
used_msacs(used_msacs < backlag | used_msacs > NT-forwardlag) = [];
[msac_trig_avg_se,lags] = get_event_trig_avg(se_phasemod_out,used_msacs,backlag,forwardlag);
[msac_trig_avg_full,lags] = get_event_trig_avg(full_phasemod_out,used_msacs,backlag,forwardlag);
[msac_trig_avg_rate,lags] = get_event_trig_avg(sm_rate,used_msacs,backlag,forwardlag);

used_sacs = used_sac_inds;
used_sacs(used_sacs < backlag | used_sacs > NT-forwardlag) = [];
[sac_trig_avg_se,lags] = get_event_trig_avg(se_phasemod_out,used_sacs,backlag,forwardlag);
[sac_trig_avg_full,lags] = get_event_trig_avg(full_phasemod_out,used_sacs,backlag,forwardlag);
[sac_trig_avg_rate,lags] = get_event_trig_avg(sm_rate,used_sacs,backlag,forwardlag);