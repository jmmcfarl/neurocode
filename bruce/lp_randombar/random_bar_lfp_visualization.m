clear all
close all

cd ~/Data/bruce/2_27_12/M232
load ./random_bar_eyedata.mat

load ./lemM232Expts.mat
%%
rf_cent = [4.73 -2.2];
axis_or = 50*pi/180; %in radians

Fs = 1000;
dsf = 4;
Fsd = Fs/dsf;
niqf = Fs/2;
[bb,aa] = butter(2,15/niqf,'low');

new_dt = .005;
dsfrac = 4;

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
% for ee = 1:length(bar_expts);
for ee = 1
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
load ./random_bar_eyedata.mat
all_fix_inds = [];
all_first_inds = [];
all_second_inds = [];
all_small_inds = [];
% for ee = 1:length(bar_expts)
for ee = 1
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

    fix_inds(isnan(fix_inds)) = [];
    first_fix_inds(isnan(first_fix_inds)) = [];
    second_fix_inds(isnan(second_fix_inds)) = [];
    small_fix_inds(isnan(small_fix_inds)) = [];
    
    all_fix_inds = [all_fix_inds; cur_e_set(fix_inds(:))];
    all_first_inds = [all_first_inds; cur_e_set(first_fix_inds(:))];
    all_second_inds = [all_second_inds; cur_e_set(second_fix_inds(:))];
    all_small_inds = [all_small_inds; cur_e_set(small_fix_inds(:))];
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
% for ee = 1:length(bar_expts);
 for ee = 1
   fprintf('Expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('lemM232A.%d.lfp.mat',bar_expts(ee));
    load(fname);
    
    Fs = 1/LFP.Header.CRsamplerate;
    
    n_trials(ee) = length(LFP.Trials);
    lfp_trial_starts = [LFP.Trials(:).Start]/1e4;
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
        
        res_lfps = [res_lfps; interp_lfps];
    end
    
    cur_set = find(full_exptvec_new==ee);
    interp_lfps = interp1(res_t_axis,res_lfps,full_taxis_new(cur_set));
    
    full_lfps = [full_lfps; interp_lfps];
 end

 %%
 close all
 n_trials = length(expt_trial_end);
 plot(expt_lfp_t_axis,expt_lfps(:,2))
 hold on
 yl = ylim();
 for i = 1:n_trials
     line([expt_trial_starts(i) expt_trial_starts(i)],yl,'color','g')
      line([expt_trial_end(i) expt_trial_end(i)],yl,'color','r')
 end
 
 all_first_times = full_taxis_new(all_first_inds);
 all_second_times = full_taxis_new(all_second_inds);
 n_first = length(all_first_inds);
 %  for i = 1:n_first
 %      line([all_first_times(i) all_first_times(i)],yl,'color','k')
 %  end
 plot(all_first_times,zeros(size(all_first_times)),'ko','linewidth',6)
 plot(all_second_times,zeros(size(all_second_times)),'ro','linewidth',6)
 
  all_small_times = full_taxis_new(all_small_inds);
 plot(all_small_times,zeros(size(all_small_times)),'g*','linewidth',4)
