clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
addpath('~/James_scripts/bruce/G081/');
%%
stim_fs = 100; %in Hz
Fs = 3e4;
dsf = 100;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[2 50]/niqf);
use_lfps = [1:12:96];

backlag = round(Fsd*0.1);
forwardlag = round(Fsd*0.4);
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
    
                cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE) & expt_sim_sacs == 0);
    %     cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 0 & expt_bar_ori == bar_oris(EE));
%     cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE) & expt_sim_sacs > 0);
    
    cur_expt_set(cur_expt_set<= 8) = []; %this expt had continuous bar luminance
    cur_expt_set(cur_expt_set==51) = []; %no SE
    cur_expt_set(cur_expt_set==50) = []; %no SE
    cur_expt_set(cur_expt_set==49) = []; %no SE
    cur_expt_set(cur_expt_set==48) = []; %no SE
    cur_expt_set(cur_expt_set==47) = []; %no SE
    cur_expt_set(cur_expt_set==46) = []; %no SE
    
    cur_expt_set(cur_expt_set <= 51) = []; %all with non-unique stim seqs
    cur_expt_set(cur_expt_set > 60) = []; %no rect
    
    %
    all_stim_times = [];
    all_rel_stime = [];
    all_rel_etime = [];
    all_Op = [];
    all_phase = [];
    all_se = [];
    all_trial_start_times = [];
    all_trial_end_times = [];
    all_trial_durs = [];
    all_trial_expt = [];
    all_trial_num = [];
    trial_cnt = 0;
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
            cur_se = Expts{cur_expt}.Trials(tt).se;
            
            all_stim_times = [all_stim_times; cur_stim_times];
            all_rel_stime = [all_rel_stime; cur_stim_times - trial_start_times(tt)];
            all_rel_etime = [all_rel_etime; trial_end_times(tt) - cur_stim_times];
            all_Op = [all_Op; cur_Op];
            all_se = [all_se; cur_se];
            all_phase = [all_phase; cur_phase'];
            all_trial_num = [all_trial_num; (trial_cnt + tt)*ones(length(cur_stim_times),1)];
        end
        
        trial_cnt = trial_cnt + n_trials;
        
        all_trial_start_times = [all_trial_start_times; trial_start_times'];
        all_trial_end_times = [all_trial_end_times; trial_end_times'];
        all_trial_durs = [all_trial_durs; trial_durs'];
        all_trial_expt = [all_trial_expt; ee*ones(length(trial_durs),1)];
    end
    un_im_set{EE} = unique(all_se);
    
    %%
    all_Vmat = [];
    all_tax = [];
    all_outgoing_sacs = [];
    all_returning_sacs = [];
    all_stim_on = [];
    all_stim_off = [];
    all_msacs = [];
    all_outgoing_sac_stimid = [];
    all_returning_sac_stimid = [];
    all_stim_on_stimid = [];
    all_stim_off_stimid = [];
    all_msac_stimid = [];
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
        
        isi = [0 diff(saccade_inds)]/Fsd;
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
        interp_trial_num = round(interp1(all_stim_times,all_trial_num,t_ax));
        
        outgoing_sacs = use_saccade_inds(pre_pos(use_saccade_inds) > -0.6 & post_pos(use_saccade_inds) < -1);
        returning_sacs = use_saccade_inds(pre_pos(use_saccade_inds) < -1 & post_pos(use_saccade_inds) > -0.6);
        outgoing_sacs = saccade_inds(outgoing_sacs);
        returning_sacs = saccade_inds(returning_sacs);
        msacs = saccade_inds(use_saccade_inds(sac_dist(use_saccade_inds)' < 1 & ~ismember(use_saccade_inds,outgoing_sacs) & ...
            ~ismember(use_saccade_inds,returning_sacs)));
        
        outgoing_sacs(outgoing_sacs < backlag | outgoing_sacs > size(Vmat,1) - forwardlag) = [];
        returning_sacs(returning_sacs < backlag | returning_sacs > size(Vmat,1) - forwardlag) = [];
        msacs(msacs < backlag | msacs > size(Vmat,1) - forwardlag) = [];
        
        stim_on = 1+find(interp_rel_stime(1:end-1) < 0.7 & interp_rel_stime(2:end) > 0.7);
        stim_off = 1+find(interp_rel_stime(1:end-1) < 1.4 & interp_rel_stime(2:end) > 1.4);
        
        cur_outgoing_stimid = all_se(interp_trial_num(outgoing_sacs));
        cur_returning_stimid = all_se(interp_trial_num(returning_sacs));
        cur_msac_stimid = all_se(interp_trial_num(msacs));
        cur_stimon_stimid = all_se(interp_trial_num(stim_on));
        cur_stimoff_stimid = all_se(interp_trial_num(stim_off));
        
        all_stim_on = [all_stim_on; stim_on' + length(all_tax)];
        all_stim_off = [all_stim_off; stim_off' + length(all_tax)];
        all_outgoing_sacs = [all_outgoing_sacs; outgoing_sacs'+length(all_tax)];
        all_returning_sacs = [all_returning_sacs; returning_sacs' + length(all_tax)];
        all_msacs = [all_msacs; msacs' + length(all_tax)];
        
        all_stim_on_stimid = [all_stim_on_stimid; cur_stimon_stimid];
        all_stim_off_stimid = [all_stim_off_stimid; cur_stimoff_stimid];
        all_outgoing_sac_stimid = [all_outgoing_sac_stimid; cur_outgoing_stimid];
        all_returning_sac_stimid = [all_returning_sac_stimid; cur_returning_stimid];
        all_msac_stimid = [all_msac_stimid; cur_msac_stimid];
        
        all_Vmat = [all_Vmat; Vmat];
        all_tax = [all_tax; t_ax'];
        
    end
    
    %%
    interp_rel_stime = interp1(all_stim_times,all_rel_stime,all_tax);
    interp_rel_etime = interp1(all_stim_times,all_rel_etime,all_tax);
    
    trial_lfp_inds = round(interp1(all_tax,1:length(all_tax),all_trial_start_times));
    n_trials = length(trial_lfp_inds);
    
    use_outgoing_inds = find(interp_rel_stime(all_outgoing_sacs) > 0.25 & interp_rel_etime(all_outgoing_sacs) > 0.15);
    use_returning_inds = find(interp_rel_stime(all_returning_sacs) > 0.25 & interp_rel_etime(all_returning_sacs) > 0.15);
    use_msac_inds = find(interp_rel_stime(all_msacs) > 0.25 & interp_rel_etime(all_msacs) > forwardlag/Fsd);
    use_stimon_inds = find(interp_rel_stime(all_stim_on) > 0.25 & interp_rel_etime(all_stim_on) > 0.15);
    use_stimoff_inds = find(interp_rel_stime(all_stim_off) > 0.25 & interp_rel_etime(all_stim_off) > 0.15);
    
    outgoing_stimids = all_outgoing_sac_stimid(use_outgoing_inds);
    returning_stimids = all_returning_sac_stimid(use_returning_inds);
    msac_stimids = all_msac_stimid(use_msac_inds);
    on_stimids = all_stim_on_stimid(use_stimon_inds);
    off_stimids = all_stim_off_stimid(use_stimoff_inds);
    
    for jj = 1:length(un_im_set{EE})
        cur_set = find(on_stimids==un_im_set{EE}(jj));
        n_on_stims(EE,jj) = length(cur_set);
        if n_on_stims(EE,jj) >= 5
            on_trig_mat = zeros(n_on_stims(EE,jj),length(lags),length(use_lfps));
            for ii = 1:n_on_stims(EE,jj)
                cur_inds = (all_stim_on(use_stimon_inds(cur_set(ii))) - backlag):(all_stim_on(use_stimon_inds(cur_set(ii))) + forwardlag);
                on_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
            end
            on_trig_avgs(EE,jj,:,:) = squeeze(mean(on_trig_mat));
            on_trig_sem(EE,jj,:,:) = squeeze(std(on_trig_mat))/sqrt(n_on_stims(EE,jj));
        else
            on_trig_avgs(EE,jj,:,:) = nan;
            on_trig_sem(EE,jj,:,:) = nan;
        end
    end
    
    for jj = 1:length(un_im_set{EE})
        cur_set = find(off_stimids==un_im_set{EE}(jj));
        n_off_stims(EE,jj) = length(cur_set);
        if n_off_stims(EE,jj) >= 5
            off_trig_mat = zeros(n_off_stims(EE,jj),length(lags),length(use_lfps));
            for ii = 1:n_off_stims(EE,jj)
                cur_inds = (all_stim_off(use_stimoff_inds(cur_set(ii))) - backlag):(all_stim_off(use_stimoff_inds(cur_set(ii))) + forwardlag);
                off_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
            end
            off_trig_avgs(EE,jj,:,:) = squeeze(mean(off_trig_mat));
            off_trig_sem(EE,jj,:,:) = squeeze(std(off_trig_mat))/sqrt(n_off_stims(EE,jj));
        else
            off_trig_avgs(EE,jj,:,:) = nan;
            off_trig_sem(EE,jj,:,:) = nan;
        end
    end
    
    for jj = 1:length(un_im_set{EE})
        cur_set = find(outgoing_stimids==un_im_set{EE}(jj));
        n_out_stims(EE,jj) = length(cur_set);
        if n_out_stims(EE,jj) >= 5
            out_trig_mat = zeros(n_out_stims(EE,jj),length(lags),length(use_lfps));
            for ii = 1:n_out_stims(EE,jj)
                cur_inds = (all_outgoing_sacs(use_outgoing_inds(cur_set(ii))) - backlag):(all_outgoing_sacs(use_outgoing_inds(cur_set(ii))) + forwardlag);
                out_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
            end
            out_trig_avgs(EE,jj,:,:) = squeeze(mean(out_trig_mat));
            out_trig_sem(EE,jj,:,:) = squeeze(std(out_trig_mat))/sqrt(n_out_stims(EE,jj));
        else
            out_trig_avgs(EE,jj,:,:) = nan;
            out_trig_sem(EE,jj,:,:) = nan;
        end
    end
    
    for jj = 1:length(un_im_set{EE})
        cur_set = find(returning_stimids==un_im_set{EE}(jj));
        n_ret_stims(EE,jj) = length(cur_set);
        if n_ret_stims(EE,jj) >= 5
            ret_trig_mat = zeros(n_ret_stims(EE,jj),length(lags),length(use_lfps));
            for ii = 1:n_ret_stims(EE,jj)
                cur_inds = (all_returning_sacs(use_returning_inds(cur_set(ii))) - backlag):(all_returning_sacs(use_returning_inds(cur_set(ii))) + forwardlag);
                ret_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
            end
            ret_trig_avgs(EE,jj,:,:) = squeeze(mean(ret_trig_mat));
            ret_trig_sem(EE,jj,:,:) = squeeze(std(ret_trig_mat))/sqrt(n_ret_stims(EE,jj));
        else
            ret_trig_avgs(EE,jj,:,:) = nan;
            ret_trig_sem(EE,jj,:,:) = nan;
        end
    end
    
    for jj = 1:length(un_im_set{EE})
        cur_set = find(msac_stimids==un_im_set{EE}(jj));
        n_msac_stims(EE,jj) = length(cur_set);
        if n_msac_stims(EE,jj) >= 5
            msac_trig_mat = zeros(n_msac_stims(EE,jj),length(lags),length(use_lfps));
            for ii = 1:n_msac_stims(EE,jj)
                cur_inds = (all_msacs(use_msac_inds(cur_set(ii))) - backlag):(all_msacs(use_msac_inds(cur_set(ii))) + forwardlag);
                msac_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
            end
            msac_trig_avgs(EE,jj,:,:) = squeeze(mean(msac_trig_mat));
            msac_trig_sem(EE,jj,:,:) = squeeze(std(msac_trig_mat))/sqrt(n_msac_stims(EE,jj));
        else
            msac_trig_avgs(EE,jj,:,:) = nan;
            msac_trig_sem(EE,jj,:,:) = nan;
        end
    end
    
    
end

%%
save stim_trg_avg_lfps_all_im lags Fsd *_trig_avgs* *_trig_sem n_* un_im_set
%%

close all
use_stims = [403004 419002 433008 434005];
% use_stims = [503002 572001];
cmap = jet(length(use_stims));
cmap(2,:) = [1 0 0];
for ll = 1:length(use_lfps)
    fprintf('Channel %d of %d\n',ll,length(use_lfps));
    
    for ii = 1:length(use_stims)
        cur_ss = find(un_im_set{1}==use_stims(ii));
        subplot(2,4,1)
        shadedErrorBar(lags/Fsd,squeeze(on_trig_avgs(1,cur_ss,:,ll)),squeeze(on_trig_sem(1,cur_ss,:,ll)),{'color',cmap(ii,:)});
        hold on
        axis tight
        subplot(2,4,5)
        shadedErrorBar(lags/Fsd,squeeze(off_trig_avgs(1,cur_ss,:,ll)),squeeze(off_trig_sem(1,cur_ss,:,ll)),{'color',cmap(ii,:)});
        hold on
        axis tight
        
        cur_ss = find(un_im_set{2}==use_stims(ii));
        subplot(2,4,2)
        shadedErrorBar(lags/Fsd,squeeze(on_trig_avgs(2,cur_ss,:,ll)),squeeze(on_trig_sem(2,cur_ss,:,ll)),{'color',cmap(ii,:)});
        hold on
        axis tight
        subplot(2,4,6)
        shadedErrorBar(lags/Fsd,squeeze(off_trig_avgs(2,cur_ss,:,ll)),squeeze(off_trig_sem(2,cur_ss,:,ll)),{'color',cmap(ii,:)});
        hold on
        axis tight
        
        cur_ss = find(un_im_set{3}==use_stims(ii));
        subplot(2,4,3)
        shadedErrorBar(lags/Fsd,squeeze(on_trig_avgs(3,cur_ss,:,ll)),squeeze(on_trig_sem(3,cur_ss,:,ll)),{'color',cmap(ii,:)});
        hold on
        axis tight
        subplot(2,4,7)
        shadedErrorBar(lags/Fsd,squeeze(off_trig_avgs(3,cur_ss,:,ll)),squeeze(off_trig_sem(3,cur_ss,:,ll)),{'color',cmap(ii,:)});
        hold on
        axis tight
        
%         cur_ss = find(un_im_set{4}==use_stims(ii));
%         subplot(2,4,4)
%         shadedErrorBar(lags/Fsd,squeeze(on_trig_avgs(4,cur_ss,:,ll)),squeeze(on_trig_sem(4,cur_ss,:,ll)),{'color',cmap(ii,:)});
%         hold on
%         axis tight
%         subplot(2,4,8)
%         shadedErrorBar(lags/Fsd,squeeze(off_trig_avgs(3,cur_ss,:,ll)),squeeze(off_trig_sem(3,cur_ss,:,ll)),{'color',cmap(ii,:)});
%         hold on
%         axis tight
    end
    
    pause
    clf
end