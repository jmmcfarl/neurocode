clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
addpath('~/James_scripts/bruce/G081/');
%%
stim_fs = 100; %in Hz
Fs = 3e4;
dsf = 80;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1 40]/niqf);
use_lfps = [1:96];

backlag = round(Fsd*0.3);
forwardlag = round(Fsd*0.5);
lags = -backlag:forwardlag;

fmaxlag = round(Fsd*2);
flags = (0:round(2*Fsd));

%%
for i = 1:length(Expts)
    if strcmp(Expts{i}.Header.expname,'grating.OpXseRC') | strcmp(Expts{i}.Header.expname,'grating.OpRC')
        is_bar_expt(i) = 1;
    else
        is_bar_expt(i) = 0;
    end
    
    if strcmp(Expts{i}.Stimvals.Bs,'image')
        expt_image_back(i) = 1;
    else
        expt_image_back(i) = 0;
    end
    
    expt_sim_sacs(i) = Expts{i}.Stimvals.ijump;
    expt_bar_ori(i) = Expts{i}.Stimvals.or;
    
end
expt_bar_ori(expt_bar_ori == -45) = 135;

%%
load ./jbeG081.em.mat

%%
bar_oris = [0 45 90 135];
for EE = 1:length(bar_oris)
    fprintf('Analyzing %d ori expts\n',bar_oris(EE));
    
            cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE) & expt_sim_sacs > 0);
%     cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 0 & expt_bar_ori == bar_oris(EE));
%             cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE) & expt_sim_sacs > 0);

cur_expt_set(cur_expt_set<= 8) = []; %this expt had continuous bar luminance
cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?
cur_expt_set(ismember(cur_expt_set,[46 48 49 51])) = []; %problem with image background
cur_expt_set(cur_expt_set > 60) = []; %no rect

%
    all_stim_times = [];
    all_rel_stime = [];
    all_rel_etime = [];
    all_Op = [];
    all_phase = [];
    all_trial_start_times = [];
    all_trial_end_times = [];
    all_trial_durs = [];
    all_trial_expt = [];
    for ee = 1:length(cur_expt_set);
        fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
        cur_expt = cur_expt_set(ee);
        
        n_trials = length(Expts{cur_expt}.Trials);
        trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
        trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
        trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
        
        for tt = 1:n_trials
            cur_stim_times = Expts{cur_expt}.Trials(tt).Start/1e4;
            cur_Op = Expts{cur_expt}.Trials(tt).Op;
            cur_phase = Expts{cur_expt}.Trials(tt).ph;
            
            all_stim_times = [all_stim_times; cur_stim_times];
            all_rel_stime = [all_rel_stime; cur_stim_times - trial_start_times(tt)];
            all_rel_etime = [all_rel_etime; trial_end_times(tt) - cur_stim_times];
            all_Op = [all_Op; cur_Op];
            all_phase = [all_phase; cur_phase'];
        end
        
        all_trial_start_times = [all_trial_start_times; trial_start_times'];
        all_trial_end_times = [all_trial_end_times; trial_end_times'];
        all_trial_durs = [all_trial_durs; trial_durs'];
        all_trial_expt = [all_trial_expt; ee*ones(length(trial_durs),1)];
    end
    %%
    all_Vmat = [];
    all_tax = [];
    all_outgoing_sacs = [];
    all_returning_sacs = [];
    all_stim_on = [];
    all_stim_off = [];
    all_msacs = [];
    for ee = 1:length(cur_expt_set);
        cur_expt = cur_expt_set(ee);
        Vmat = [];
        for ll = 1:length(use_lfps)
            fprintf('Electrode %d of %d\n',ll,length(use_lfps));
            filename = sprintf('Expt%d.p%dFullV.mat',cur_expt,use_lfps(ll));
            load(filename);
            V = double(FullV.V);
            %     V = V + FullV.sumscale*sumv;
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
        
        %% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
        [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,t_ax(([1 end])));
        
        eye_dt = median(diff(eye_ts_interp));
        eye_fs = 1/eye_dt;
        lEyeXY = eye_vals_interp(:,1:2);
        clear sm_avg_eyepos eye_vel
        sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),3);
        sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),3);
        eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
        eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
        
        eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
        
        interp_eye_vals = interp1(eye_ts_interp,lEyeXY,t_ax);
        interp_eye_speed = interp1(eye_ts_interp,eye_speed,t_ax);
        
        orth_eye_proj = (cos(bar_oris(EE)/180*pi)*interp_eye_vals(:,1) - ...
            sin(bar_oris(EE)/180*pi)*interp_eye_vals(:,2));
        
        %find threshold crossings
        sac_thresh = 25;
        saccade_inds = 1 + find(interp_eye_speed((1:end-1)) < sac_thresh & interp_eye_speed((2:end)) > sac_thresh);
        saccade_inds = unique(saccade_inds);
        
        isi = [Inf diff(saccade_inds)]/Fsd;
        min_isi = 0.05;
        bad_sacs = find(isi < min_isi);
        saccade_inds(bad_sacs) = [];
        
        peri_thresh = 3;
        sac_starts = nan(size(saccade_inds));
        sac_stops = nan(size(saccade_inds));
        for i = 1:length(saccade_inds)
            cur_up = find(interp_eye_speed(1:saccade_inds(i)) < peri_thresh,1,'last');
            cur_down = find(interp_eye_speed(saccade_inds(i):end) < peri_thresh,1,'first');
            if ~isempty(cur_up) & ~isempty(cur_down)
                sac_starts(i) = cur_up;
                sac_stops(i) = saccade_inds(i) + cur_down;
            end
        end
        bad = find(isnan(sac_starts));
        saccade_inds(bad) = []; sac_starts(bad) = []; sac_stops(bad) = [];
        
        sac_starts = sac_starts - round(Fsd*0.01);
        sac_stops = sac_stops + round(Fsd*0.01);
        bad = find(sac_starts <= 0 | sac_starts > length(orth_eye_proj));
        saccade_inds(bad) = []; sac_starts(bad) = []; sac_stops(bad) = [];
        pre_pos = orth_eye_proj(sac_starts);
        post_pos = orth_eye_proj(sac_stops);
        sac_dist = abs(post_pos-pre_pos);
        
        cur_trial_set = find(all_trial_expt==ee);
        n_trials = length(cur_trial_set);
        in_trial = zeros(size(t_ax));
        for tt = 1:n_trials
            cur_set = find(t_ax >= all_trial_start_times(cur_trial_set(tt)) & t_ax < all_trial_end_times(cur_trial_set(tt)));
            in_trial(cur_set) = 1;
        end
        use_saccade_inds = find(ismember(saccade_inds,find(in_trial==1)));
        
        interp_rel_stime = interp1(all_stim_times,all_rel_stime,t_ax);
        interp_rel_etime = interp1(all_stim_times,all_rel_etime,t_ax);
        
%         outgoing_sacs = use_saccade_inds(pre_pos(use_saccade_inds) > -0.6 & post_pos(use_saccade_inds) < -1);
%         returning_sacs = use_saccade_inds(pre_pos(use_saccade_inds) < -1 & post_pos(use_saccade_inds) > -0.6);
%         outgoing_sacs = saccade_inds(outgoing_sacs);
%         returning_sacs = saccade_inds(returning_sacs);
%         msacs = saccade_inds(use_saccade_inds(sac_dist(use_saccade_inds)' < 1 & ~ismember(use_saccade_inds,outgoing_sacs) & ...
%             ~ismember(use_saccade_inds,returning_sacs)));
    if bar_oris(EE) == 0
        outgoing_sacs = saccade_inds(pre_pos > -0.9 & post_pos < -0.9);
        returning_sacs = saccade_inds(pre_pos < -1 & post_pos > -1);
    elseif bar_oris(EE) == 135
        outgoing_sacs = saccade_inds(pre_pos > -1 & post_pos < -1);
        returning_sacs = saccade_inds(pre_pos < -1.1 & post_pos > -1.2);
    elseif bar_oris(EE) == 45
        outgoing_sacs = saccade_inds(pre_pos > -1 & post_pos < -1);
        returning_sacs = saccade_inds(pre_pos < -1.1 & post_pos > -1);
    elseif bar_oris(EE) == 90
        outgoing_sacs = saccade_inds(pre_pos > -1 & post_pos < -1);
        returning_sacs = saccade_inds(pre_pos < -1 & post_pos > -1);
    end
    msacs = saccade_inds(sac_dist' < 1); msacs(ismember(msacs,outgoing_sacs)) = []; msacs(ismember(msacs,returning_sacs)) = [];
    
    use_out = find(interp_rel_stime(outgoing_sacs) > 0.8 & interp_rel_stime(outgoing_sacs) < 1.3);
    use_ret = find(interp_rel_stime(returning_sacs) > 1.4 & interp_rel_stime(returning_sacs) < 1.9);
    
    outgoing_sacs = outgoing_sacs(use_out);
    returning_sacs = returning_sacs(use_ret);
    
%     bad_sacs = find(all_rel_stimes(cur_tr_inds(saccade_inds)) < beg_sbuffer | all_rel_etimes(cur_tr_inds(saccade_inds)) < end_sbuffer);
%     outgoing_sacs(ismember(outgoing_sacs,saccade_inds(bad_sacs))) = [];
%     returning_sacs(ismember(returning_sacs,saccade_inds(bad_sacs))) = [];
%     msacs(ismember(msacs,saccade_inds(bad_sacs))) = [];
%     saccade_inds(bad_sacs) = [];
        
        outgoing_sacs(outgoing_sacs < backlag | outgoing_sacs > size(Vmat,1) - forwardlag) = [];
        returning_sacs(returning_sacs < backlag | returning_sacs > size(Vmat,1) - forwardlag) = [];
        msacs(msacs < backlag | msacs > size(Vmat,1) - forwardlag) = [];
        
        stim_on = 1+find(interp_rel_stime(1:end-1) < 0.7 & interp_rel_stime(2:end) > 0.7);
        stim_off = 1+find(interp_rel_stime(1:end-1) < 1.4 & interp_rel_stime(2:end) > 1.4);
        
        all_stim_on = [all_stim_on; stim_on' + length(all_tax)];
        all_stim_off = [all_stim_off; stim_off' + length(all_tax)];
        all_outgoing_sacs = [all_outgoing_sacs; outgoing_sacs'+length(all_tax)];
        all_returning_sacs = [all_returning_sacs; returning_sacs' + length(all_tax)];
        all_msacs = [all_msacs; msacs' + length(all_tax)];
        
        all_Vmat = [all_Vmat; Vmat];
        all_tax = [all_tax; t_ax'];
        
    end
    
    %     all_Vmat = zscore(all_Vmat);
    all_V_std = std(all_Vmat);
    %%
    interp_rel_stime = interp1(all_stim_times,all_rel_stime,all_tax);
    interp_rel_etime = interp1(all_stim_times,all_rel_etime,all_tax);
    
    trial_lfp_inds = round(interp1(all_tax,1:length(all_tax),all_trial_start_times));
    %     stim_lfp_inds = round(interp1(all_tax,1:length(all_tax),all_stim_times));
    n_trials = length(trial_lfp_inds);
    
    %     use_stims = find(all_rel_stime > backlag & all_rel_etime > forwardlag/Fsd);
    use_outgoing_inds = all_outgoing_sacs(interp_rel_stime(all_outgoing_sacs) > 0.25 & interp_rel_etime(all_outgoing_sacs) > 0.15);
    use_returning_inds = all_returning_sacs(interp_rel_stime(all_returning_sacs) > 0.25 & interp_rel_etime(all_returning_sacs) > 0.15);
    use_msac_inds = all_msacs(interp_rel_stime(all_msacs) > 0.25 & interp_rel_etime(all_msacs) > forwardlag/Fsd);
    
    use_stimon_inds = all_stim_on(interp_rel_stime(all_stim_on) > 0.25 & interp_rel_etime(all_stim_on) > 0.15);
    use_stimoff_inds = all_stim_off(interp_rel_stime(all_stim_off) > 0.25 & interp_rel_etime(all_stim_off) > 0.15);
    
    n_on_stims(EE) = length(use_stimon_inds);
    on_trig_mat = zeros(n_on_stims(EE),length(lags),length(use_lfps));
    for ii = 1:length(use_stimon_inds)
        cur_inds = (use_stimon_inds(ii) - backlag):(use_stimon_inds(ii) + forwardlag);
        on_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
    end
    on_trig_avgs(EE,:,:) = squeeze(mean(on_trig_mat));
    on_trig_sem(EE,:,:) = squeeze(std(on_trig_mat))/sqrt(length(use_stimon_inds));
    
    n_off_stims(EE) = length(use_stimoff_inds);
    off_trig_mat = zeros(n_off_stims(EE),length(lags),length(use_lfps));
    for ii = 1:length(use_stimoff_inds)
        cur_inds = (use_stimoff_inds(ii) - backlag):(use_stimoff_inds(ii) + forwardlag);
        off_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
    end
    off_trig_avgs(EE,:,:) = squeeze(mean(off_trig_mat));
    off_trig_sem(EE,:,:) = squeeze(std(off_trig_mat))/sqrt(length(use_stimoff_inds));
    
    n_outsacs(EE) = length(use_outgoing_inds);
    if n_outsacs(EE) > 10
    outsac_trig_mat = zeros(n_outsacs(EE),length(lags),length(use_lfps));
    for ii = 1:length(use_outgoing_inds)
        cur_inds = (use_outgoing_inds(ii) - backlag):(use_outgoing_inds(ii) + forwardlag);
        outsac_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
    end
    outsac_trig_avgs(EE,:,:) = squeeze(mean(outsac_trig_mat));
    outsac_trig_sem(EE,:,:) = squeeze(std(outsac_trig_mat))/sqrt(length(use_outgoing_inds));
    end
    
    n_retsacs(EE) = length(use_returning_inds);
    if n_retsacs(EE) > 10
    retsac_trig_mat = zeros(length(use_returning_inds),length(lags),length(use_lfps));
    for ii = 1:length(use_returning_inds)
        cur_inds = (use_returning_inds(ii) - backlag):(use_returning_inds(ii) + forwardlag);
        retsac_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
    end
    retsac_trig_avgs(EE,:,:) = squeeze(mean(retsac_trig_mat));
    retsac_trig_sem(EE,:,:) = squeeze(std(retsac_trig_mat))/sqrt(length(use_returning_inds));
    end
    
    n_msacs(EE) = length(use_msac_inds);
    msac_trig_mat = zeros(length(use_msac_inds),length(lags),length(use_lfps));
    for ii = 1:length(use_msac_inds)
        cur_inds = (use_msac_inds(ii) - backlag):(use_msac_inds(ii) + forwardlag);
        msac_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
    end
    msac_trig_avgs(EE,:,:) = squeeze(mean(msac_trig_mat));
    msac_trig_sem(EE,:,:) = squeeze(std(msac_trig_mat))/sqrt(length(use_msac_inds));
    
    nu_trials(EE) = n_trials;
    full_trig_mat = zeros(n_trials,length(flags),length(use_lfps));
    for ii = 1:n_trials
        cur_inds = (trial_lfp_inds(ii)):(trial_lfp_inds(ii) + fmaxlag);
        full_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
    end
    full_trig_avgs(EE,:,:) = squeeze(mean(full_trig_mat));
    full_trig_sem(EE,:,:) = squeeze(std(full_trig_mat))/sqrt(n_trials);
    
    
    sac_t_ax = linspace(0,2,200);
    all_sac_inds = [all_msacs; all_outgoing_sacs; all_returning_sacs];
    sac_hist(EE,:) = hist(interp_rel_stime(all_sac_inds),sac_t_ax);
    
    out_hist(EE,:) = hist(interp_rel_stime(use_outgoing_inds),sac_t_ax);
    ret_hist(EE,:) = hist(interp_rel_stime(use_returning_inds),sac_t_ax);
    msac_hist(EE,:) = hist(interp_rel_stime(use_msac_inds),sac_t_ax);
    
end

%%
save event_trg_avg_lfps_all_sim2 lags Fsd flags sac_t_ax *_hist *_trig_avgs* *_trig_sem n_* nu_* all_V_std
%%
close all
load ./event_trg_avg_lfps_all_gray2.mat
gray_out_avgs = outsac_trig_avgs;
gray_ret_avgs = retsac_trig_avgs;
gray_msac_avgs = msac_trig_avgs;
gray_on_avgs = on_trig_avgs;
gray_off_avgs = off_trig_avgs;
gray_out_sems = outsac_trig_sem;
gray_ret_sems = retsac_trig_sem;
gray_msac_sems = msac_trig_sem;
gray_on_sems = on_trig_sem;
gray_off_sems = off_trig_sem;
gray_nlfps = size(gray_out_avgs,3);
gray_lags = lags;

load ./event_trg_avg_lfps_all_im2.mat
im_out_avgs = outsac_trig_avgs;
im_ret_avgs = retsac_trig_avgs;
im_msac_avgs = msac_trig_avgs;
im_on_avgs = on_trig_avgs;
im_off_avgs = off_trig_avgs;
im_out_sems = outsac_trig_sem;
im_ret_sems = retsac_trig_sem;
im_msac_sems = msac_trig_sem;
im_on_sems = on_trig_sem;
im_off_sems = off_trig_sem;
im_nlfps = size(im_out_avgs,3);
im_lags = lags;

load ./event_trg_avg_lfps_all_sim2.mat
sim_out_avgs = on_trig_avgs;
sim_ret_avgs = off_trig_avgs;
sim_msac_avgs = msac_trig_avgs;
sim_out_sems = on_trig_sem;
sim_ret_sems = off_trig_sem;
sim_msac_sems = msac_trig_sem;
sim_nlfps = size(sim_out_avgs,3);
sim_lags = lags;

if (gray_nlfps ~= im_nlfps) | (im_nlfps ~= sim_nlfps)
    disp('Error different LFP sets')
end

% load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
% X_pos = ArrayConfig.X;
% Y_pos = ArrayConfig.Y;

for ll = 1:length(use_lfps)
    
    fprintf('Channel %d of %d\n',ll,length(use_lfps));
%     fprintf('X: %d  Y: %d\n',X_pos(use_lfps(ll)),Y_pos(use_lfps(ll)));

    subplot(3,4,1)
    shadedErrorBar(gray_lags/Fsd,squeeze(gray_out_avgs(1,:,ll)),squeeze(gray_out_sems(1,:,ll)),{'b'})
    hold on
    shadedErrorBar(gray_lags/Fsd,squeeze(gray_ret_avgs(1,:,ll)),squeeze(gray_ret_sems(1,:,ll)),{'r'})
%     errorbar(lags/Fsd,squeeze(gray_msac_avgs(1,:,ll)),squeeze(gray_msac_sems(1,:,ll)),'g')
    title('Gray back 0deg')
    axis tight
    yl1 = ylim();
    xl = xlim();
    line([0 0],yl1,'color','k')
    line(xl,[0 0],'color','k')
    subplot(3,4,5)
    shadedErrorBar(im_lags/Fsd,squeeze(im_out_avgs(1,:,ll)),squeeze(im_out_sems(1,:,ll)),{'b'})
    hold on
    shadedErrorBar(im_lags/Fsd,squeeze(im_ret_avgs(1,:,ll)),squeeze(im_ret_sems(1,:,ll)),{'r'})
%     errorbar(lags/Fsd,squeeze(im_msac_avgs(1,:,ll)),squeeze(im_msac_sems(1,:,ll)),'g')
    axis tight
    title('IM back 0deg')
    yl2 = ylim();
    line(xl,[0 0],'color','k')
     line([0 0],yl2,'color','k')
    subplot(3,4,9)
    shadedErrorBar(sim_lags/Fsd,squeeze(sim_out_avgs(1,:,ll)),squeeze(sim_out_sems(1,:,ll)),{'b'})
    hold on
    shadedErrorBar(sim_lags/Fsd,squeeze(sim_ret_avgs(1,:,ll)),squeeze(sim_ret_sems(1,:,ll)),{'r'})
%     errorbar(lags/Fsd,squeeze(sim_msac_avgs(1,:,ll)),squeeze(sim_msac_sems(1,:,ll)),'g')
    axis tight
    title('SIM IM 0deg')
    yl3 = ylim();
    line(xl,[0 0],'color','k')
    line([0 0],yl3,'color','k')
    lv = min([yl1(1),yl2(1),yl3(1)]);
    hv = max([yl1(2),yl2(2),yl3(2)]);
    
    
    subplot(3,4,2)
    shadedErrorBar(gray_lags/Fsd,squeeze(gray_out_avgs(2,:,ll)),squeeze(gray_out_sems(2,:,ll)),{'b'})
    hold on
    shadedErrorBar(gray_lags/Fsd,squeeze(gray_ret_avgs(2,:,ll)),squeeze(gray_ret_sems(2,:,ll)),{'r'})
%     errorbar(lags/Fsd,squeeze(gray_msac_avgs(2,:,ll)),squeeze(gray_msac_sems(2,:,ll)),'g')
    axis tight
    title('Gray back 45deg')
    yl1 = ylim();
    xl = xlim();
    line([0 0],yl1,'color','k')
    line(xl,[0 0],'color','k')
    subplot(3,4,6)
    shadedErrorBar(im_lags/Fsd,squeeze(im_out_avgs(2,:,ll)),squeeze(im_out_sems(2,:,ll)),{'b'})
    hold on
    shadedErrorBar(im_lags/Fsd,squeeze(im_ret_avgs(2,:,ll)),squeeze(im_ret_sems(2,:,ll)),{'r'})
%     plot(lags/Fsd,squeeze(gray_ret_avgs(2,:,ll))+squeeze(sim_ret_avgs(2,:,ll)),'k','linewidth',2)
%     errorbar(lags/Fsd,squeeze(im_msac_avgs(2,:,ll)),squeeze(im_msac_sems(2,:,ll)),'g')
    title('IM back 45deg')
    axis tight
    yl2 = ylim();
    line(xl,[0 0],'color','k')
    line([0 0],yl2,'color','k')
    subplot(3,4,10)
    shadedErrorBar(sim_lags/Fsd,squeeze(sim_out_avgs(2,:,ll)),squeeze(sim_out_sems(2,:,ll)),{'b'})
    hold on
    shadedErrorBar(sim_lags/Fsd,squeeze(sim_ret_avgs(2,:,ll)),squeeze(sim_ret_sems(2,:,ll)),{'r'})
%     errorbar(lags/Fsd,squeeze(sim_msac_avgs(2,:,ll)),squeeze(sim_msac_sems(2,:,ll)),'g')
    axis tight
    title('SIM IM 45deg')
    yl3 = ylim();
    line(xl,[0 0],'color','k')   
    line([0 0],yl3,'color','k')
    clv = min([yl1(1),yl2(1),yl3(1)]);
    chv = max([yl1(2),yl2(2),yl3(2)]);
    lv = min(lv,clv); hv = max(hv,chv);
 
    
    
    subplot(3,4,3)
    shadedErrorBar(gray_lags/Fsd,squeeze(gray_out_avgs(3,:,ll)),squeeze(gray_out_sems(3,:,ll)),{'b'})
    hold on
    shadedErrorBar(gray_lags/Fsd,squeeze(gray_ret_avgs(3,:,ll)),squeeze(gray_ret_sems(3,:,ll)),{'r'})
%     errorbar(lags/Fsd,squeeze(gray_msac_avgs(3,:,ll)),squeeze(gray_msac_sems(3,:,ll)),'g')
    axis tight
    title('Gray back 90deg')
    yl1 = ylim();
    xl = xlim();
    line([0 0],yl1,'color','k')
    line(xl,[0 0],'color','k')
    subplot(3,4,7)
    shadedErrorBar(im_lags/Fsd,squeeze(im_out_avgs(3,:,ll)),squeeze(im_out_sems(3,:,ll)),{'b'})
    hold on
    shadedErrorBar(im_lags/Fsd,squeeze(im_ret_avgs(3,:,ll)),squeeze(im_ret_sems(3,:,ll)),{'r'})
%     errorbar(lags/Fsd,squeeze(im_msac_avgs(3,:,ll)),squeeze(im_msac_sems(3,:,ll)),'g')
    title('IM back 90deg')
    axis tight
    yl2 = ylim();
    line(xl,[0 0],'color','k')
    line([0 0],yl2,'color','k')
    subplot(3,4,11)
    shadedErrorBar(sim_lags/Fsd,squeeze(sim_out_avgs(3,:,ll)),squeeze(sim_out_sems(3,:,ll)),{'b'})
    hold on
    shadedErrorBar(sim_lags/Fsd,squeeze(sim_ret_avgs(3,:,ll)),squeeze(sim_ret_sems(3,:,ll)),{'r'})
%     errorbar(lags/Fsd,squeeze(sim_msac_avgs(3,:,ll)),squeeze(sim_msac_sems(3,:,ll)),'g')
    title('SIM IM 90deg')
    axis tight
    yl3 = ylim();
    line(xl,[0 0],'color','k')   
    line([0 0],yl3,'color','k')
    clv = min([yl1(1),yl2(1),yl3(1)]);
    chv = max([yl1(2),yl2(2),yl3(2)]);
    lv = min(lv,clv); hv = max(hv,chv);

    
    subplot(3,4,4)
    shadedErrorBar(gray_lags/Fsd,squeeze(gray_out_avgs(4,:,ll)),squeeze(gray_out_sems(4,:,ll)),{'b'})
    hold on
    shadedErrorBar(gray_lags/Fsd,squeeze(gray_ret_avgs(4,:,ll)),squeeze(gray_ret_sems(4,:,ll)),{'r'})
%     errorbar(lags/Fsd,squeeze(gray_msac_avgs(4,:,ll)),squeeze(gray_msac_sems(4,:,ll)),'g')
    title('Gray back 135deg')
    axis tight
    yl1 = ylim();
    xl = xlim();
    line([0 0],yl1,'color','k')
    line(xl,[0 0],'color','k')
    subplot(3,4,8)
    shadedErrorBar(im_lags/Fsd,squeeze(im_out_avgs(4,:,ll)),squeeze(im_out_sems(4,:,ll)),{'b'})
    hold on
    shadedErrorBar(im_lags/Fsd,squeeze(im_ret_avgs(4,:,ll)),squeeze(im_ret_sems(4,:,ll)),{'r'})
%     errorbar(lags/Fsd,squeeze(im_msac_avgs(4,:,ll)),squeeze(im_msac_sems(4,:,ll)),'g')
    title('IM back 135deg')
    axis tight
    yl2 = ylim();
    line(xl,[0 0],'color','k')
    line([0 0],yl2,'color','k')
    subplot(3,4,12)
    shadedErrorBar(sim_lags/Fsd,squeeze(sim_out_avgs(4,:,ll)),squeeze(sim_out_sems(4,:,ll)),{'b'})
    hold on
    shadedErrorBar(sim_lags/Fsd,squeeze(sim_ret_avgs(4,:,ll)),squeeze(sim_ret_sems(4,:,ll)),{'r'})
%     errorbar(lags/Fsd,squeeze(sim_msac_avgs(4,:,ll)),squeeze(sim_msac_sems(4,:,ll)),'g')
    axis tight
    title('SIM IM 135deg')
    yl3 = ylim();
    line(xl,[0 0],'color','k')   
    line([0 0],yl3,'color','k')
    clv = min([yl1(1),yl2(1),yl3(1)]);
    chv = max([yl1(2),yl2(2),yl3(2)]);
    lv = min(lv,clv); hv = max(hv,chv);
    
    
    for jj = 1:12
        subplot(3,4,jj)
        ylim([lv hv])
    end
    pause
    clf
end

%%
    subplot(3,4,1)
    errorbar(lags/Fsd,squeeze(mean(gray_out_avgs(1,:,:),3)),squeeze(std(gray_out_avgs(1,:,:),[],3))/sqrt(96),'b')
    hold on
    errorbar(lags/Fsd,squeeze(mean(gray_ret_avgs(1,:,:),3)),squeeze(std(gray_ret_avgs(1,:,:),[],3))/sqrt(96),'r')
    %     errorbar(lags/Fsd,squeeze(gray_msac_avgs(1,:,ll)),squeeze(gray_msac_sems(1,:,ll)),'g')
    title('Gray back 0deg')
    axis tight
    yl1 = ylim();
    xl = xlim();
    line([0 0],yl1,'color','k')
    line(xl,[0 0],'color','k')
    subplot(3,4,5)
    errorbar(lags/Fsd,squeeze(mean(im_out_avgs(1,:,:),3)),squeeze(std(im_out_avgs(1,:,:),[],3))/sqrt(96),'b')
    hold on
    errorbar(lags/Fsd,squeeze(mean(im_ret_avgs(1,:,:),3)),squeeze(std(im_ret_avgs(1,:,:),[],3))/sqrt(96),'r')
    %     errorbar(lags/Fsd,squeeze(im_msac_avgs(1,:,ll)),squeeze(im_msac_sems(1,:,ll)),'g')
    axis tight
    title('IM back 0deg')
    yl2 = ylim();
    line(xl,[0 0],'color','k')
    line([0 0],yl2,'color','k')
    subplot(3,4,9)
    errorbar(lags/Fsd,squeeze(sim_out_avgs(1,:,ll)),squeeze(sim_out_sems(1,:,ll)),'b')
    hold on
    errorbar(lags/Fsd,squeeze(sim_ret_avgs(1,:,ll)),squeeze(sim_ret_sems(1,:,ll)),'r')
    %     errorbar(lags/Fsd,squeeze(sim_msac_avgs(1,:,ll)),squeeze(sim_msac_sems(1,:,ll)),'g')
    axis tight
    title('SIM IM 0deg')
    yl3 = ylim();
    line(xl,[0 0],'color','k')
    line([0 0],yl3,'color','k')
    lv = min([yl1(1),yl2(1),yl3(1)]);
    hv = max([yl1(2),yl2(2),yl3(2)]);
    
    
    subplot(3,4,2)
    errorbar(lags/Fsd,squeeze(gray_out_avgs(2,:,ll)),squeeze(gray_out_sems(2,:,ll)),'b')
    hold on
    errorbar(lags/Fsd,squeeze(gray_ret_avgs(2,:,ll)),squeeze(gray_ret_sems(2,:,ll)),'r')
    %     errorbar(lags/Fsd,squeeze(gray_msac_avgs(2,:,ll)),squeeze(gray_msac_sems(2,:,ll)),'g')
    axis tight
    title('Gray back 45deg')
    yl1 = ylim();
    xl = xlim();
    line([0 0],yl1,'color','k')
    line(xl,[0 0],'color','k')
    subplot(3,4,6)
    errorbar(lags/Fsd,squeeze(im_out_avgs(2,:,ll)),squeeze(im_out_sems(2,:,ll)),'b')
    hold on
    errorbar(lags/Fsd,squeeze(im_ret_avgs(2,:,ll)),squeeze(im_ret_sems(2,:,ll)),'r')
    %     plot(lags/Fsd,squeeze(gray_ret_avgs(2,:,ll))+squeeze(sim_ret_avgs(2,:,ll)),'k','linewidth',2)
    %     errorbar(lags/Fsd,squeeze(im_msac_avgs(2,:,ll)),squeeze(im_msac_sems(2,:,ll)),'g')
    title('IM back 45deg')
    axis tight
    yl2 = ylim();
    line(xl,[0 0],'color','k')
    line([0 0],yl2,'color','k')
    subplot(3,4,10)
    errorbar(lags/Fsd,squeeze(sim_out_avgs(2,:,ll)),squeeze(sim_out_sems(2,:,ll)),'b')
    hold on
    errorbar(lags/Fsd,squeeze(sim_ret_avgs(2,:,ll)),squeeze(sim_ret_sems(2,:,ll)),'r')
    %     errorbar(lags/Fsd,squeeze(sim_msac_avgs(2,:,ll)),squeeze(sim_msac_sems(2,:,ll)),'g')
    axis tight
    title('SIM IM 45deg')
    yl3 = ylim();
    line(xl,[0 0],'color','k')
    line([0 0],yl3,'color','k')
    clv = min([yl1(1),yl2(1),yl3(1)]);
    chv = max([yl1(2),yl2(2),yl3(2)]);
    lv = min(lv,clv); hv = max(hv,chv);
    
    
    
    subplot(3,4,3)
    errorbar(lags/Fsd,squeeze(gray_out_avgs(3,:,ll)),squeeze(gray_out_sems(3,:,ll)),'b')
    hold on
    errorbar(lags/Fsd,squeeze(gray_ret_avgs(3,:,ll)),squeeze(gray_ret_sems(3,:,ll)),'r')
    %     errorbar(lags/Fsd,squeeze(gray_msac_avgs(3,:,ll)),squeeze(gray_msac_sems(3,:,ll)),'g')
    axis tight
    title('Gray back 90deg')
    yl1 = ylim();
    xl = xlim();
    line([0 0],yl1,'color','k')
    line(xl,[0 0],'color','k')
    subplot(3,4,7)
    errorbar(lags/Fsd,squeeze(im_out_avgs(3,:,ll)),squeeze(im_out_sems(3,:,ll)),'b')
    hold on
    errorbar(lags/Fsd,squeeze(im_ret_avgs(3,:,ll)),squeeze(im_ret_sems(3,:,ll)),'r')
    %     errorbar(lags/Fsd,squeeze(im_msac_avgs(3,:,ll)),squeeze(im_msac_sems(3,:,ll)),'g')
    title('IM back 90deg')
    axis tight
    yl2 = ylim();
    line(xl,[0 0],'color','k')
    line([0 0],yl2,'color','k')
    subplot(3,4,11)
    errorbar(lags/Fsd,squeeze(sim_out_avgs(3,:,ll)),squeeze(sim_out_sems(3,:,ll)),'b')
    hold on
    errorbar(lags/Fsd,squeeze(sim_ret_avgs(3,:,ll)),squeeze(sim_ret_sems(3,:,ll)),'r')
    %     errorbar(lags/Fsd,squeeze(sim_msac_avgs(3,:,ll)),squeeze(sim_msac_sems(3,:,ll)),'g')
    title('SIM IM 90deg')
    axis tight
    yl3 = ylim();
    line(xl,[0 0],'color','k')
    line([0 0],yl3,'color','k')
    clv = min([yl1(1),yl2(1),yl3(1)]);
    chv = max([yl1(2),yl2(2),yl3(2)]);
    lv = min(lv,clv); hv = max(hv,chv);
    
    
    subplot(3,4,4)
    errorbar(lags/Fsd,squeeze(gray_out_avgs(4,:,ll)),squeeze(gray_out_sems(4,:,ll)),'b')
    hold on
    errorbar(lags/Fsd,squeeze(gray_ret_avgs(4,:,ll)),squeeze(gray_ret_sems(4,:,ll)),'r')
    %     errorbar(lags/Fsd,squeeze(gray_msac_avgs(4,:,ll)),squeeze(gray_msac_sems(4,:,ll)),'g')
    title('Gray back 135deg')
    axis tight
    yl1 = ylim();
    xl = xlim();
    line([0 0],yl1,'color','k')
    line(xl,[0 0],'color','k')
    subplot(3,4,8)
    errorbar(lags/Fsd,squeeze(im_out_avgs(4,:,ll)),squeeze(im_out_sems(4,:,ll)),'b')
    hold on
    errorbar(lags/Fsd,squeeze(im_ret_avgs(4,:,ll)),squeeze(im_ret_sems(4,:,ll)),'r')
    %     errorbar(lags/Fsd,squeeze(im_msac_avgs(4,:,ll)),squeeze(im_msac_sems(4,:,ll)),'g')
    title('IM back 135deg')
    axis tight
    yl2 = ylim();
    line(xl,[0 0],'color','k')
    line([0 0],yl2,'color','k')
    subplot(3,4,12)
    errorbar(lags/Fsd,squeeze(sim_out_avgs(4,:,ll)),squeeze(sim_out_sems(4,:,ll)),'b')
    hold on
    errorbar(lags/Fsd,squeeze(sim_ret_avgs(4,:,ll)),squeeze(sim_ret_sems(4,:,ll)),'r')
    %     errorbar(lags/Fsd,squeeze(sim_msac_avgs(4,:,ll)),squeeze(sim_msac_sems(4,:,ll)),'g')
    axis tight
    title('SIM IM 135deg')
    yl3 = ylim();
    line(xl,[0 0],'color','k')
    line([0 0],yl3,'color','k')
    clv = min([yl1(1),yl2(1),yl3(1)]);
    chv = max([yl1(2),yl2(2),yl3(2)]);
    lv = min(lv,clv); hv = max(hv,chv);
