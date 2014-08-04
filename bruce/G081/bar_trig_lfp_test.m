clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
addpath('~/James_scripts/bruce/G081/');
%%
stim_fs = 100; %in Hz
Fs = 3e4;
dsf = 60;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[2]/niqf,'high');
use_lfps = [1:1:96];

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
for EE = 1:length(bar_oris);
    fprintf('Analyzing %d ori expts\n',bar_oris(EE));
    
    cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE) & expt_sim_sacs > 0);
    
%         %expts with X deg bars and gray back (sim sacs)
%         cur_expt_set = find(is_bar_expt==1 & expt_image_back == 0 & expt_bar_ori == bar_oris(EE));
% %         cur_expt_set(cur_expt_set==2) = []; %this expt had continuous bar luminance
%         cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?
    cur_expt_set(cur_expt_set<= 8) = []; %this expt had continuous bar luminance
    %
    for ee = 1:length(cur_expt_set);
        fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
        cur_expt = cur_expt_set(ee);
        
        n_trials = length(Expts{cur_expt}.Trials);
        trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
        trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
        trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
        
        expt_stim_times = [];
        expt_rel_stime = [];
        expt_rel_etime = [];
        expt_Op = [];
        expt_phase = [];
        for tt = 1:n_trials
            cur_stim_times = Expts{cur_expt}.Trials(tt).Start;
            cur_Op = Expts{cur_expt}.Trials(tt).Op;
            cur_phase = Expts{cur_expt}.Trials(tt).ph;
            
            expt_stim_times = [expt_stim_times; cur_stim_times];
            expt_rel_stime = [expt_rel_stime; cur_stim_times/1e4 - trial_start_times(tt)];
            expt_rel_etime = [expt_rel_etime; trial_end_times(tt) - cur_stim_times/1e4];
            expt_Op = [expt_Op; cur_Op];
            expt_phase = [expt_phase; cur_phase'];
        end
        expt_stim_times = expt_stim_times/1e4;
        %%
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
        
        in_trial = zeros(size(t_ax));
        for tt = 1:n_trials
            cur_set = find(t_ax >= trial_start_times(tt) & t_ax < trial_end_times(tt));
            in_trial(cur_set) = 1;
        end
        use_saccade_inds = find(ismember(saccade_inds,find(in_trial==1)));
%         
%         outgoing_sacs = use_saccade_inds(pre_pos(use_saccade_inds) > -0.6 & post_pos(use_saccade_inds) < -1);
%         returning_sacs = use_saccade_inds(pre_pos(use_saccade_inds) < -1 & post_pos(use_saccade_inds) > -0.6);
%         outgoing_sacs = saccade_inds(outgoing_sacs);
%         returning_sacs = saccade_inds(returning_sacs);
        
        
        interp_expt_rel_stime = interp1(expt_stim_times,expt_rel_stime,t_ax);
        interp_expt_rel_etime = interp1(expt_stim_times,expt_rel_etime,t_ax);
        
        outgoing_sacs = 1+find(interp_expt_rel_stime(1:end-1) < 0.7 & interp_expt_rel_stime(2:end) > 0.7);
        returning_sacs = 1+find(interp_expt_rel_stime(1:end-1) < 1.4 & interp_expt_rel_stime(2:end) > 1.4);
        msacs = saccade_inds(use_saccade_inds(sac_dist(use_saccade_inds)' < 1 & ~ismember(use_saccade_inds,outgoing_sacs) & ...
            ~ismember(use_saccade_inds,returning_sacs)));
        %%
        trial_lfp_inds = round(interp1(t_ax,1:length(t_ax),trial_start_times));
        
        expt_lfp_inds = round(interp1(t_ax,1:length(t_ax),expt_stim_times));
        unique_bar_pos = unique(expt_Op);
        n_bar_pos = length(unique_bar_pos);
        
        backlag = round(Fsd*0.1);
        forwardlag = round(Fsd*0.3);
        lags = -backlag:forwardlag;
        
        outgoing_sacs(outgoing_sacs < backlag | outgoing_sacs > size(Vmat,1) - forwardlag) = [];
        returning_sacs(returning_sacs < backlag | returning_sacs > size(Vmat,1) - forwardlag) = [];
        msacs(msacs < backlag | msacs > size(Vmat,1) - forwardlag) = [];
       
        use_stims = find(expt_rel_stime > backlag & expt_rel_etime > forwardlag/Fsd);
        use_outgoing_inds = outgoing_sacs(interp_expt_rel_stime(outgoing_sacs) > 0.25 & interp_expt_rel_etime(outgoing_sacs) > 0.15);
        use_returning_inds = returning_sacs(interp_expt_rel_stime(returning_sacs) > 0.25 & interp_expt_rel_etime(returning_sacs) > 0.15);
        use_msac_inds = msacs(interp_expt_rel_stime(msacs) > 0.25 & interp_expt_rel_etime(msacs) > forwardlag/Fsd);
        
        outsac_trig_avgs(ee,:,:) = zeros(length(lags),length(use_lfps));
        for ii = 1:length(use_outgoing_inds)
            cur_inds = (use_outgoing_inds(ii) - backlag):(use_outgoing_inds(ii) + forwardlag);
            outsac_trig_avgs(ee,:,:) = squeeze(outsac_trig_avgs(ee,:,:)) + Vmat(cur_inds,:);
        end
        outsac_trig_avgs(ee,:,:) = outsac_trig_avgs(ee,:,:)/length(use_outgoing_inds);
        
        retsac_trig_avgs(ee,:,:) = zeros(length(lags),length(use_lfps));
        for ii = 1:length(use_returning_inds)
            cur_inds = (use_returning_inds(ii) - backlag):(use_returning_inds(ii) + forwardlag);
            retsac_trig_avgs(ee,:,:) = squeeze(retsac_trig_avgs(ee,:,:)) + Vmat(cur_inds,:);
        end
        retsac_trig_avgs(ee,:,:) = retsac_trig_avgs(ee,:,:)/length(use_returning_inds);
        
        msac_trig_avgs(ee,:,:) = zeros(length(lags),length(use_lfps));
        for ii = 1:length(use_msac_inds)
            cur_inds = (use_msac_inds(ii) - backlag):(use_msac_inds(ii) + forwardlag);
            msac_trig_avgs(ee,:,:) = squeeze(msac_trig_avgs(ee,:,:)) + Vmat(cur_inds,:);
        end
        msac_trig_avgs(ee,:,:) = msac_trig_avgs(ee,:,:)/length(use_msac_inds);

        comb_trig_avgs(ee,:,:,:) = zeros(n_bar_pos,length(lags),length(use_lfps));
        for bb = 1:n_bar_pos
            cur_bar_inds = use_stims(expt_Op(use_stims)==unique_bar_pos(bb));
            for ii = 1:length(cur_bar_inds)
                cur_inds = (expt_lfp_inds(cur_bar_inds(ii)) - backlag):(expt_lfp_inds(cur_bar_inds(ii)) + forwardlag);
                comb_trig_avgs(ee,bb,:,:) = squeeze(comb_trig_avgs(ee,bb,:,:)) + Vmat(cur_inds,:);
            end
            comb_trig_avgs(ee,bb,:,:) = comb_trig_avgs(ee,bb,:,:)/length(cur_bar_inds);
        end
        
        black_trig_avgs(ee,:,:,:) = zeros(n_bar_pos,length(lags),length(use_lfps));
        for bb = 1:n_bar_pos
            cur_bar_inds = use_stims(expt_Op(use_stims)==unique_bar_pos(bb) & expt_phase(use_stims) == 0);
            for ii = 1:length(cur_bar_inds)
                cur_inds = (expt_lfp_inds(cur_bar_inds(ii)) - backlag):(expt_lfp_inds(cur_bar_inds(ii)) + forwardlag);
                black_trig_avgs(ee,bb,:,:) = squeeze(black_trig_avgs(ee,bb,:,:)) + Vmat(cur_inds,:);
            end
            black_trig_avgs(ee,bb,:,:) = black_trig_avgs(ee,bb,:,:)/length(cur_bar_inds);
        end
        
        white_trig_avgs(ee,:,:,:) = zeros(n_bar_pos,length(lags),length(use_lfps));
        for bb = 1:n_bar_pos
            cur_bar_inds = use_stims(expt_Op(use_stims)==unique_bar_pos(bb) & expt_phase(use_stims) == pi);
            for ii = 1:length(cur_bar_inds)
                cur_inds = (expt_lfp_inds(cur_bar_inds(ii)) - backlag):(expt_lfp_inds(cur_bar_inds(ii)) + forwardlag);
                white_trig_avgs(ee,bb,:,:) = squeeze(white_trig_avgs(ee,bb,:,:)) + Vmat(cur_inds,:);
            end
            white_trig_avgs(ee,bb,:,:) = white_trig_avgs(ee,bb,:,:)/length(cur_bar_inds);
        end
 
        fmaxlag = round(Fsd*2);
        flags = (0:round(2*Fsd))/Fsd;
        full_trig_avgs(ee,:,:) = zeros(length(flags),length(use_lfps));
        for ii = 1:n_trials
            cur_inds = (trial_lfp_inds(ii)):(trial_lfp_inds(ii) + fmaxlag);
            full_trig_avgs(ee,:,:) = squeeze(full_trig_avgs(ee,:,:)) + Vmat(cur_inds,:);
        end
        full_trig_avgs(ee,:,:) = full_trig_avgs(ee,:,:)/n_trials;
        
     use_n_trials(EE,ee) = n_trials;
   end
    
    %%
    avg_full_trig_avgs(EE,:,:) = squeeze(mean(full_trig_avgs));
    avg_black_trig_avgs(EE,:,:,:) = squeeze(mean(black_trig_avgs));
    avg_white_trig_avgs(EE,:,:,:) = squeeze(mean(white_trig_avgs));
    avg_comb_trig_avgs(EE,:,:,:) = squeeze(mean(comb_trig_avgs));
    avg_outsac_trig_avgs(EE,:,:) = squeeze(mean(outsac_trig_avgs));
    avg_retsac_trig_avgs(EE,:,:) = squeeze(mean(retsac_trig_avgs));
    avg_msac_trig_avgs(EE,:,:) = squeeze(mean(msac_trig_avgs));
    std_full_trig_avgs(EE,:,:,:) = squeeze(std(full_trig_avgs));
    std_black_trig_avgs(EE,:,:,:) = squeeze(std(black_trig_avgs));
    std_white_trig_avgs(EE,:,:,:) = squeeze(std(white_trig_avgs));
    std_comb_trig_avgs(EE,:,:,:) = squeeze(std(comb_trig_avgs));
    std_outsac_trig_avgs(EE,:,:) = squeeze(std(outsac_trig_avgs));
    std_retsac_trig_avgs(EE,:,:) = squeeze(std(retsac_trig_avgs));
    std_msac_trig_avgs(EE,:,:) = squeeze(std(msac_trig_avgs));
    
end

%%
ab_comb_avgs = (mean(avg_comb_trig_avgs,2));
ms_comb_trig_avgs = bsxfun(@minus,avg_comb_trig_avgs,ab_comb_avgs);

%%
save simsac_im_trig_avgs ms_comb_trig_avgs lags forwardlag backlag use_lfps unique_bar_pos avg_* std_*

%%
for i = 1:96
plot(flags,squeeze(avg_full_trig_avgs(:,:,i)))
yl = ylim();
line([0.7 0.7],yl,'color','k');
line([1.4 1.4],yl,'color','k')
pause
clf
end
%%
load ./images_im_trig_avgs
for i = 1:length(use_lfps)
    load ./images_im_trig_avgs
subplot(3,3,1)
    errorbar(lags/Fsd,squeeze(avg_outsac_trig_avgs(1,:,i)),squeeze(std_outsac_trig_avgs(1,:,i))/2)
    shg
    hold on
    errorbar(lags/Fsd,squeeze(avg_outsac_trig_avgs(2,:,i)),squeeze(std_outsac_trig_avgs(2,:,i))/2,'r')
    errorbar(lags/Fsd,squeeze(avg_outsac_trig_avgs(3,:,i)),squeeze(std_outsac_trig_avgs(3,:,i))/2,'k')
    errorbar(lags/Fsd,squeeze(avg_outsac_trig_avgs(4,:,i)),squeeze(std_outsac_trig_avgs(4,:,i))/2,'g')
    xlim([-0.1 0.3])
    subplot(3,3,4)
    errorbar(lags/Fsd,squeeze(avg_retsac_trig_avgs(1,:,i)),squeeze(std_retsac_trig_avgs(1,:,i))/2)
    shg
    hold on
    errorbar(lags/Fsd,squeeze(avg_retsac_trig_avgs(2,:,i)),squeeze(std_retsac_trig_avgs(2,:,i))/2,'r')
    errorbar(lags/Fsd,squeeze(avg_retsac_trig_avgs(3,:,i)),squeeze(std_retsac_trig_avgs(3,:,i))/2,'k')
    errorbar(lags/Fsd,squeeze(avg_retsac_trig_avgs(4,:,i)),squeeze(std_retsac_trig_avgs(4,:,i))/2,'g')
    xlim([-0.1 0.3])
    subplot(3,3,7)
    errorbar(lags/Fsd,squeeze(avg_msac_trig_avgs(1,:,i)),squeeze(std_msac_trig_avgs(1,:,i))/2)
    shg
    hold on
    errorbar(lags/Fsd,squeeze(avg_msac_trig_avgs(2,:,i)),squeeze(std_msac_trig_avgs(2,:,i))/2,'r')
    errorbar(lags/Fsd,squeeze(avg_msac_trig_avgs(3,:,i)),squeeze(std_msac_trig_avgs(3,:,i))/2,'k')
    errorbar(lags/Fsd,squeeze(avg_msac_trig_avgs(4,:,i)),squeeze(std_msac_trig_avgs(4,:,i))/2,'g')
    xlim([-0.1 0.3])
    
    
    load ./gray_im_trig_avgs
subplot(3,3,2)
    errorbar(lags/Fsd,squeeze(avg_outsac_trig_avgs(1,:,i)),squeeze(std_outsac_trig_avgs(1,:,i))/2)
    shg
    hold on
    errorbar(lags/Fsd,squeeze(avg_outsac_trig_avgs(2,:,i)),squeeze(std_outsac_trig_avgs(2,:,i))/2,'r')
    errorbar(lags/Fsd,squeeze(avg_outsac_trig_avgs(3,:,i)),squeeze(std_outsac_trig_avgs(3,:,i))/2,'k')
    errorbar(lags/Fsd,squeeze(avg_outsac_trig_avgs(4,:,i)),squeeze(std_outsac_trig_avgs(4,:,i))/2,'g')
    xlim([-0.1 0.3])
    subplot(3,3,5)
    errorbar(lags/Fsd,squeeze(avg_retsac_trig_avgs(1,:,i)),squeeze(std_retsac_trig_avgs(1,:,i))/2)
    shg
    hold on
    errorbar(lags/Fsd,squeeze(avg_retsac_trig_avgs(2,:,i)),squeeze(std_retsac_trig_avgs(2,:,i))/2,'r')
    errorbar(lags/Fsd,squeeze(avg_retsac_trig_avgs(3,:,i)),squeeze(std_retsac_trig_avgs(3,:,i))/2,'k')
    errorbar(lags/Fsd,squeeze(avg_retsac_trig_avgs(4,:,i)),squeeze(std_retsac_trig_avgs(4,:,i))/2,'g')
    xlim([-0.1 0.3])
    subplot(3,3,8)
    errorbar(lags/Fsd,squeeze(avg_msac_trig_avgs(1,:,i)),squeeze(std_msac_trig_avgs(1,:,i))/2)
    shg
    hold on
    errorbar(lags/Fsd,squeeze(avg_msac_trig_avgs(2,:,i)),squeeze(std_msac_trig_avgs(2,:,i))/2,'r')
    errorbar(lags/Fsd,squeeze(avg_msac_trig_avgs(3,:,i)),squeeze(std_msac_trig_avgs(3,:,i))/2,'k')
    errorbar(lags/Fsd,squeeze(avg_msac_trig_avgs(4,:,i)),squeeze(std_msac_trig_avgs(4,:,i))/2,'g')
    xlim([-0.1 0.3])

    load ./simsac_im_trig_avgs
subplot(3,3,3)
    errorbar(lags/Fsd,squeeze(avg_outsac_trig_avgs(1,:,i)),squeeze(std_outsac_trig_avgs(1,:,i))/2)
    shg
    hold on
    errorbar(lags/Fsd,squeeze(avg_outsac_trig_avgs(2,:,i)),squeeze(std_outsac_trig_avgs(2,:,i))/2,'r')
    errorbar(lags/Fsd,squeeze(avg_outsac_trig_avgs(3,:,i)),squeeze(std_outsac_trig_avgs(3,:,i))/2,'k')
    errorbar(lags/Fsd,squeeze(avg_outsac_trig_avgs(4,:,i)),squeeze(std_outsac_trig_avgs(4,:,i))/2,'g')
    xlim([-0.1 0.3])
    subplot(3,3,6)
    errorbar(lags/Fsd,squeeze(avg_retsac_trig_avgs(1,:,i)),squeeze(std_retsac_trig_avgs(1,:,i))/2)
    shg
    hold on
    errorbar(lags/Fsd,squeeze(avg_retsac_trig_avgs(2,:,i)),squeeze(std_retsac_trig_avgs(2,:,i))/2,'r')
    errorbar(lags/Fsd,squeeze(avg_retsac_trig_avgs(3,:,i)),squeeze(std_retsac_trig_avgs(3,:,i))/2,'k')
    errorbar(lags/Fsd,squeeze(avg_retsac_trig_avgs(4,:,i)),squeeze(std_retsac_trig_avgs(4,:,i))/2,'g')
    xlim([-0.1 0.3])
    subplot(3,3,9)
    errorbar(lags/Fsd,squeeze(avg_msac_trig_avgs(1,:,i)),squeeze(std_msac_trig_avgs(1,:,i))/2)
    shg
    hold on
    errorbar(lags/Fsd,squeeze(avg_msac_trig_avgs(2,:,i)),squeeze(std_msac_trig_avgs(2,:,i))/2,'r')
    errorbar(lags/Fsd,squeeze(avg_msac_trig_avgs(3,:,i)),squeeze(std_msac_trig_avgs(3,:,i))/2,'k')
    errorbar(lags/Fsd,squeeze(avg_msac_trig_avgs(4,:,i)),squeeze(std_msac_trig_avgs(4,:,i))/2,'g')
    xlim([-0.1 0.3])
    pause
    clf
    
end
%%
cx = [-4 4]*1e-6;
% for ll = 1:length(use_lfps)
%     ll
%     subplot(3,1,1)
%     pcolor(lags/Fsd,unique_bar_pos,squeeze(avg_white_trig_avgs(:,:,ll)));shading flat
%     caxis(cx)
%     subplot(3,1,2)
%     pcolor(lags/Fsd,unique_bar_pos,squeeze(avg_black_trig_avgs(:,:,ll)));shading flat
%     caxis(cx)
%     subplot(3,1,3)
%     pcolor(lags/Fsd,unique_bar_pos,squeeze(avg_comb_trig_avgs(:,:,ll)));shading flat
%     caxis(cx)
%     pause
%     clf
% end
load ./gray_im_trig_avgs
im_comb_trig_avgs = ms_comb_trig_avgs;
load ./gray_back_trig_avgs.mat
for ll = 1:length(use_lfps)
    ll
    subplot(4,2,1)
    pcolor(lags/Fsd,unique_bar_pos,squeeze(ms_comb_trig_avgs(1,:,:,ll)));shading flat
    caxis(cx)
    title('0 Deg');
    subplot(4,2,2)
    pcolor(lags/Fsd,unique_bar_pos,squeeze(ms_comb_trig_avgs(2,:,:,ll)));shading flat
    caxis(cx)
    title('45 Deg');
    subplot(4,2,3)
    pcolor(lags/Fsd,unique_bar_pos,squeeze(ms_comb_trig_avgs(3,:,:,ll)));shading flat
    caxis(cx)
    title('90 Deg');
    subplot(4,2,4)
    pcolor(lags/Fsd,unique_bar_pos,squeeze(ms_comb_trig_avgs(4,:,:,ll)));shading flat
    caxis(cx)
    title('135 Deg');
    
    subplot(4,2,5)
    pcolor(lags/Fsd,unique_bar_pos,squeeze(im_comb_trig_avgs(1,:,:,ll)));shading flat
    caxis(cx)
    title('0 Deg');
    subplot(4,2,6)
    pcolor(lags/Fsd,unique_bar_pos,squeeze(im_comb_trig_avgs(2,:,:,ll)));shading flat
    caxis(cx)
    title('45 Deg');
    subplot(4,2,7)
    pcolor(lags/Fsd,unique_bar_pos,squeeze(im_comb_trig_avgs(3,:,:,ll)));shading flat
    caxis(cx)
    title('90 Deg');
    subplot(4,2,8)
    pcolor(lags/Fsd,unique_bar_pos,squeeze(im_comb_trig_avgs(4,:,:,ll)));shading flat
    caxis(cx)
    title('135 Deg');
    pause
    clf
end

% for ll = 1:length(use_lfps)
%     ll
%     subplot(2,2,1)
%     pcolor(lags/Fsd,unique_bar_pos,squeeze(ms_comb_trig_avgs(1,:,:,ll)));shading flat
%     caxis(cx)
%     title('0 Deg');
%     subplot(2,2,2)
%     pcolor(lags/Fsd,unique_bar_pos,squeeze(ms_comb_trig_avgs(2,:,:,ll)));shading flat
%     caxis(cx)
%     title('45 Deg');
%     subplot(2,2,3)
%     pcolor(lags/Fsd,unique_bar_pos,squeeze(ms_comb_trig_avgs(3,:,:,ll)));shading flat
%     caxis(cx)
%     title('90 Deg');
%     subplot(2,2,4)
%     pcolor(lags/Fsd,unique_bar_pos,squeeze(ms_comb_trig_avgs(4,:,:,ll)));shading flat
%     caxis(cx)
%     title('135 Deg');
%     pause
%     clf
% end
%
%
% %%
% load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
% X_pos = ArrayConfig.X;
% Y_pos = ArrayConfig.Y;
% [~,y_ord] = sort(Y_pos(use_lfps));
% [~,x_ord] = sort(X_pos(use_lfps));
%
% for bb = 1:n_bar_pos
%     bb
%     subplot(2,1,1)
%     pcolor(lags/Fsd,1:length(use_lfps),squeeze(avg_comb_trig_avgs(bb,:,x_ord))');shading flat
%     caxis([-4 4]*1e-6)
%     xlim([0 0.15])
%
%     subplot(2,1,2)
%     pcolor(lags/Fsd,1:length(use_lfps),squeeze(avg_comb_trig_avgs(bb,:,y_ord))');shading flat
%     caxis([-4 4]*1e-6)
%     xlim([0 0.15])
%
%     pause
%     clf
% end

%%
cur_expt_set = find(is_bar_expt==1 & expt_image_back == 0 & expt_bar_ori == 0);
cur_expt_set(cur_expt_set==2) = []; %this expt had continuous bar luminance
cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?

cur_expt = cur_expt_set(1);

n_trials = length(Expts{cur_expt}.Trials);
trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;

expt_stim_times = [];
expt_rel_stime = [];
expt_rel_etime = [];
expt_Op = [];
expt_phase = [];
for tt = 1:n_trials
    cur_stim_times = Expts{cur_expt}.Trials(tt).Start;
    cur_Op = Expts{cur_expt}.Trials(tt).Op;
    cur_phase = Expts{cur_expt}.Trials(tt).ph;
    
    expt_stim_times = [expt_stim_times; cur_stim_times];
    expt_rel_stime = [expt_rel_stime; cur_stim_times/1e4 - trial_start_times(tt)];
    expt_rel_etime = [expt_rel_etime; trial_end_times(tt) - cur_stim_times/1e4];
    expt_Op = [expt_Op; cur_Op];
    expt_phase = [expt_phase; cur_phase'];
end
expt_stim_times = expt_stim_times/1e4;
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

expt_lfp_inds = round(interp1(t_ax,1:length(t_ax),expt_stim_times));
unique_bar_pos = unique(expt_Op);
n_bar_pos = length(unique_bar_pos);
use_stims = find(expt_rel_stime > 0.1 & expt_rel_etime > forwardlag/Fsd);

plags = sum(lags > 0);
nlags = sum(lags < 0);
zpad = plags - nlags;
pad_avgs = cat(3,zeros(4,n_bar_pos,zpad,length(use_lfps)),ms_comb_trig_avgs);
pad_avgs(:,:,1:plags,:) = 0;
for ll = 1:length(use_lfps)
    lin_pred(ll,:) = zeros(1,length(t_ax));
    for bb = 1:n_bar_pos
        cur_lin_filt = squeeze(pad_avgs(1,bb,:,ll));
        cur_lin_filt = cur_lin_filt - mean(cur_lin_filt);
        cur_bar_inds = use_stims(expt_Op(use_stims)==unique_bar_pos(bb));
        null_sig = zeros(size(t_ax));
        null_sig(expt_lfp_inds(cur_bar_inds)) = 1;
        cur_pred = conv(null_sig,cur_lin_filt,'same');
        lin_pred(ll,:) = lin_pred(ll,:) + cur_pred;
    end
end

trial_start_lfp_inds = round(interp1(t_ax,1:length(t_ax),trial_start_times));
trial_end_lfp_inds = round(interp1(t_ax,1:length(t_ax),trial_end_times));
use_lfp_inds = zeros(size(t_ax));
for i = 1:length(trial_start_lfp_inds)
    use_lfp_inds((trial_start_lfp_inds(i)+round(Fsd*0.25)):(trial_end_lfp_inds(i)-round(Fsd*0.1))) = 1;
end
use_lfp_inds = logical(use_lfp_inds);

