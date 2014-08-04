clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
addpath('~/James_scripts/bruce/G081/');
%%
stim_fs = 100; %in Hz
use_sus = 1:96;

%% PARSE TRIAL DATA STRUCTURES
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

load ./all_un_bar_pos
n_bar_pos = size(all_un_bar_pos,1);

%%
poss_bar_oris = [0 45 90 135];
for EE = 1:4
    %% USE ONLY GRAY BACKGROUND DATA
    flen = 12;
    beg_buffer = round(stim_fs*0.2);
    bar_oris = poss_bar_oris(EE);
    un_bar_pos = all_un_bar_pos(:,EE);
    
    fprintf('Analyzing %d ori expts\n',bar_oris);
    
    %expts with X deg bars and any back (including sim sacs)
    cur_expt_set = find(is_bar_expt==1 & expt_bar_ori == bar_oris);
    
    cur_expt_set(cur_expt_set<=8) = []; %this expt had continuous bar luminance
    cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?
    
    cur_expt_set(cur_expt_set >= 46 & cur_expt_set <= 51) = []; %problem with image background
    cur_expt_set(cur_expt_set > 60) = []; %no rect
    
    %% COMPUTE TRIAL DATA
    all_stim_times = [];
    all_rel_stimes = [];
    all_rel_etimes = [];
    all_phase = [];
    all_Op = [];
    all_bar_mat = [];
    all_used_inds = [];
    all_binned_spks = [];
    all_exptvec = [];
    all_trialvec = [];
    for ee = 1:length(cur_expt_set);
        fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
        cur_expt = cur_expt_set(ee);
        fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
        load(fname);
        
        n_trials = length(Expts{cur_expt}.Trials);
        trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
        trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
        trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
        
        for tt = 1:n_trials
            cur_stim_times = Expts{cur_expt}.Trials(tt).Start/1e4;
            cur_Op = Expts{cur_expt}.Trials(tt).Op;
            cur_phase = Expts{cur_expt}.Trials(tt).ph;
            
            cur_t_edges = [cur_stim_times; Expts{cur_expt}.Trials(tt).End(end)/1e4];
            cur_binned_spks = nan(length(cur_stim_times),length(use_sus));
            for cc = 1:length(use_sus)
                cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges);
                cur_binned_spks(:,cc) = cur_hist(1:end-1);
            end
            
            cur_bar_mat = zeros(length(cur_stim_times),n_bar_pos);
            for bb = 1:n_bar_pos
                cur_set = find(cur_Op==un_bar_pos(bb));
                pset = cur_phase(cur_set) == 0;
                nset = cur_phase(cur_set) == pi;
                cur_bar_mat(cur_set(pset),bb) = 1;
                cur_bar_mat(cur_set(nset),bb) = -1;
                %                 cur_bar_mat(cur_set,bb) = 1;
            end
            
            bar_Xmat = makeStimRows(cur_bar_mat,flen);
            cur_used_inds = ones(length(cur_stim_times),1);
            cur_used_inds(1:flen) = 0;
            cur_used_inds(1:beg_buffer) = 0;
            
            all_stim_times = [all_stim_times; cur_stim_times];
            all_rel_stimes = [all_rel_stimes; cur_stim_times- trial_start_times(tt)];
            all_rel_etimes = [all_rel_etimes; trial_end_times(tt) - cur_stim_times];
            all_Op = [all_Op; cur_Op];
            all_phase = [all_phase; cur_phase'];
            all_used_inds = [all_used_inds; cur_used_inds];
            all_bar_mat = [all_bar_mat; bar_Xmat];
            all_binned_spks = [all_binned_spks; cur_binned_spks];
            all_exptvec = [all_exptvec; ones(size(cur_stim_times))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_stim_times))*tt];
        end
        
    end
    
    %% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
    load ./jbeG081.em.mat
    all_eye_vals = [];
    all_eye_speed = [];
    all_eye_ts = [];
    for ee = 1:length(cur_expt_set);
        cur_set = find(all_exptvec==ee);
        [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_stim_times(cur_set([1 end])));
        
        eye_dt = median(diff(eye_ts_interp));
        eye_fs = 1/eye_dt;
        lEyeXY = eye_vals_interp(:,1:2);
        rEyeXY = eye_vals_interp(:,3:4);
        clear sm_avg_eyepos eye_vel
        sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),3);
        sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),3);
        eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
        eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
        
        eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
        
        all_eye_vals = [all_eye_vals; lEyeXY rEyeXY];
        all_eye_speed = [all_eye_speed; eye_speed];
        all_eye_ts = [all_eye_ts; eye_ts_interp'];
    end
    back_pts = 1 + find(diff(all_eye_ts) <= 0);
    double_samples = [];
    for i = 1:length(back_pts)
        next_forward = find(all_eye_ts > all_eye_ts(back_pts(i)-1),1,'first');
        double_samples = [double_samples back_pts(i):next_forward];
    end
    all_eye_ts(double_samples) = [];
    all_eye_speed(double_samples) = [];
    all_eye_vals(double_samples,:) = [];
    
    interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_stim_times);
    interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_stim_times);
    
    par_eye_pos = interp_eye_vals(:,1)*cos(bar_oris*pi/180) - interp_eye_vals(:,2)*sin(bar_oris*pi/180);
    % orth_eye_pos = interp_eye_vals(:,1)*sin(bar_oris*pi/180) - interp_eye_vals(:,2)*cos(bar_oris*pi/180);
    % orth_eye_pos = interp_eye_vals(:,2);
    % par_eye_pos = interp_eye_vals(:,1);
    % bad_pts = find(abs(orth_eye_pos) > 1);
    
    %% IF YOU WANT TO DO CROSS-VALIDATION ON INITIAL MODEL FITS, OTHERWISE,
    [c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
    n_trials = length(ia);
    % bad_trials = unique(ic(bad_pts)); %trials with putative blinks
    
    tr_set = 1:n_trials;
    tr_set = find(ismember(1:n_trials,tr_set));
    
    cur_tr_inds = find(ismember(ic,tr_set));
    
    cur_tr_inds(all_used_inds(cur_tr_inds) == 0) = [];
    
    Xmat = all_bar_mat;
    %% double sample Xmat
    bar_dx = 0.08;
    up_fac = 2;
    bar_dxd = bar_dx/up_fac;
    SDIM = n_bar_pos;
    NT = length(cur_tr_inds);
    resh_X = reshape(Xmat(cur_tr_inds,:)',[flen SDIM NT]);
    
    if up_fac > 1
        fprintf('Up-sampling X-matrix\n');
        dresh_X = [];
        for i = 1:SDIM-1
            dresh_X = cat(2,dresh_X,resh_X(:,i,:));
            dresh_X = cat(2,dresh_X,0.5*resh_X(:,i,:) + 0.5*resh_X(:,i+1,:));
        end
        dresh_X = cat(2,dresh_X,resh_X(:,end,:));
        
        up_SDIM = up_fac*SDIM-1;
    else
        up_SDIM = SDIM;
        dresh_X = resh_X;
    end
    X_z = reshape(dresh_X,flen*up_SDIM,NT)';
    
    %% LOAD AND INTERPOLATE MODELS
    load ./CellList.mat
    good_sus = find(all(CellList(:,:,1) > 0));
    
    fname = sprintf('./full_eye_correct_%ddeg.mat',bar_oris);
    load(fname);
    used_tr_inds = find(ismember(eye_times,all_stim_times(cur_tr_inds)));
    
    %% RECONSTRUCT MAP STIMULUS
    resh_X = reshape(X_z',[flen up_SDIM NT]);
    resh_X_sh = zeros(size(resh_X));
    for ii = 1:NT
        d2 = dist_shift2d(resh_X(:,:,ii), -it_shift_cor{end}(used_tr_inds(ii)),2,0);
        resh_X_sh(:,:,ii) = d2;
    end
    X_sh = reshape(resh_X_sh,flen*up_SDIM,NT)';
    
    %%
    use_eye_speed = interp_eye_speed(cur_tr_inds);
    
    trial_start_inds = find(diff(all_trialvec(cur_tr_inds)) ~= 0) + 1;
    end_sbuffer = 0.2;
    beg_sbuffer = 0.2;
    sac_thresh = 10;
    min_isi = 0.05;
    saccade_inds = 1 + find(use_eye_speed(1:end-1) < sac_thresh & use_eye_speed(2:end) > sac_thresh);
    saccade_inds = unique(saccade_inds);
    isis = [Inf; diff(saccade_inds)]/stim_fs;
    double_sacs = find(isis < min_isi);
    saccade_inds(double_sacs) = [];
    
    peri_thresh = 3;
    sac_starts = nan(size(saccade_inds));
    sac_stops = nan(size(saccade_inds));
    for i = 1:length(saccade_inds)
        cur_up = find(use_eye_speed(1:saccade_inds(i)) < peri_thresh,1,'last');
        cur_down = find(use_eye_speed(saccade_inds(i):end) < peri_thresh,1,'first');
        if ~isempty(cur_up) & ~isempty(cur_down)
            sac_starts(i) = cur_up;
            sac_stops(i) = saccade_inds(i) + cur_down;
        end
    end
    bad = find(isnan(sac_starts));
    saccade_inds(bad) = []; sac_starts(bad) = []; sac_stops(bad) = [];
    
    sac_starts = sac_starts - round(stim_fs*0.01);
    sac_stops = sac_stops + round(stim_fs*0.01);
    pre_pos = par_eye_pos(cur_tr_inds(sac_starts));
    post_pos = par_eye_pos(cur_tr_inds(sac_stops));
    
    
    sac_dist = abs(post_pos-pre_pos);
    if bar_oris == 0
        outgoing_sacs = saccade_inds(pre_pos > -0.9 & post_pos < -0.9);
        returning_sacs = saccade_inds(pre_pos < -1 & post_pos > -1);
    elseif bar_oris == 135
        outgoing_sacs = saccade_inds(pre_pos > -1 & post_pos < -1);
        returning_sacs = saccade_inds(pre_pos < -1.1 & post_pos > -1.2);
    elseif bar_oris == 45
        outgoing_sacs = saccade_inds(pre_pos > -1 & post_pos < -1);
        returning_sacs = saccade_inds(pre_pos < -1.1 & post_pos > -1);
    elseif bar_oris == 90
        outgoing_sacs = saccade_inds(pre_pos > -1 & post_pos < -1);
        returning_sacs = saccade_inds(pre_pos < -1 & post_pos > -1);
    end
    msacs = saccade_inds(sac_dist' < 1); msacs(ismember(msacs,outgoing_sacs)) = []; msacs(ismember(msacs,returning_sacs)) = [];
    
    use_out = find(all_rel_stimes(cur_tr_inds(outgoing_sacs)) > 0.8 & all_rel_stimes(cur_tr_inds(outgoing_sacs)) < 1.3);
    use_ret = find(all_rel_stimes(cur_tr_inds(returning_sacs)) > 1.4 & all_rel_stimes(cur_tr_inds(returning_sacs)) < 1.9);
    
    outgoing_sacs = outgoing_sacs(use_out);
    returning_sacs = returning_sacs(use_ret);
    
    bad_sacs = find(all_rel_stimes(cur_tr_inds(saccade_inds)) < beg_sbuffer | all_rel_etimes(cur_tr_inds(saccade_inds)) < end_sbuffer);
    outgoing_sacs(ismember(outgoing_sacs,saccade_inds(bad_sacs))) = [];
    returning_sacs(ismember(returning_sacs,saccade_inds(bad_sacs))) = [];
    msacs(ismember(msacs,saccade_inds(bad_sacs))) = [];
    saccade_inds(bad_sacs) = [];
    
    onstim_inds = 1+find(all_rel_stimes(cur_tr_inds(1:end-1)) < 0.7 & all_rel_stimes(cur_tr_inds(2:end)) > 0.7);
    offstim_inds = 1+find(all_rel_stimes(cur_tr_inds(1:end-1)) < 1.4 & all_rel_stimes(cur_tr_inds(2:end)) > 1.4);
    cstim_inds = unique([onstim_inds(:); offstim_inds(:)]);
    
    
    used_saccade_inds = [outgoing_sacs; returning_sacs];
    %%
    dt = 0.01;
    tent_centers = [0:dt:0.65];
    tent_centers = round(tent_centers/dt);
    tbmat = construct_tent_bases(tent_centers,1);
    [ntents,tblen] = size(tbmat);
    
    shift = round(0.2/dt);
    tbmat = [zeros(ntents,tblen-2*shift-1) tbmat];
    tent_centers = tent_centers-shift;
    
    sac_inds = zeros(NT,1);
    sac_inds(used_saccade_inds) = 1;
    sac_Tmat = zeros(NT,ntents);
    for i = 1:ntents
        sac_Tmat(:,i) = conv(sac_inds,tbmat(i,:),'same');
    end
    
    osac_inds = zeros(NT,1);
    osac_inds(outgoing_sacs) = 1;
    osac_Tmat = zeros(NT,ntents);
    for i = 1:ntents
        osac_Tmat(:,i) = conv(osac_inds,tbmat(i,:),'same');
    end
    
    rsac_inds = zeros(NT,1);
    rsac_inds(returning_sacs) = 1;
    rsac_Tmat = zeros(NT,ntents);
    for i = 1:ntents
        rsac_Tmat(:,i) = conv(rsac_inds,tbmat(i,:),'same');
    end
    
    msac_inds = zeros(NT,1);
    msac_inds(msacs) = 1;
    msac_Tmat = zeros(NT,ntents);
    for i = 1:ntents
        msac_Tmat(:,i) = conv(msac_inds,tbmat(i,:),'same');
    end
    
    stim_inds = zeros(NT,1);
    stim_inds(cstim_inds) = 1;
    stim_Tmat = zeros(NT,ntents);
    for i = 1:ntents
        stim_Tmat(:,i) = conv(stim_inds,tbmat(i,:),'same');
    end
    
    stim_inds = zeros(NT,1);
    stim_inds(onstim_inds) = 1;
    onstim_Tmat = zeros(NT,ntents);
    for i = 1:ntents
        onstim_Tmat(:,i) = conv(stim_inds,tbmat(i,:),'same');
    end
    stim_inds = zeros(NT,1);
    stim_inds(offstim_inds) = 1;
    offstim_Tmat = zeros(NT,ntents);
    for i = 1:ntents
        offstim_Tmat(:,i) = conv(stim_inds,tbmat(i,:),'same');
    end
    
    %%
    gray_back_real = find(expt_image_back(cur_expt_set) == 0);
    im_back_real = find(expt_image_back(cur_expt_set) == 1 & expt_sim_sacs(cur_expt_set)==0);
    im_back_sim = find(expt_image_back(cur_expt_set) == 1 & expt_sim_sacs(cur_expt_set)>0);
    
    gray_back_real_inds = find(ismember(all_exptvec(cur_tr_inds),gray_back_real));
    im_back_real_inds = find(ismember(all_exptvec(cur_tr_inds),im_back_real));
    im_back_sim_inds = find(ismember(all_exptvec(cur_tr_inds),im_back_sim));
    
    %%
    if EE == 4
        it_quad_fit(end) = [];
    end
    lambda = 50;
    % clear sac_kern gray_sac_kern im_sac_kern sim_sac_kern
    for cc = 1:96
        cc
        
        Robs = all_binned_spks(cur_tr_inds,cc);
        tr_spkbns = convert_to_spikebins(Robs);
        
        %     [nll, pnll, lpen, prate, g] = getLL_GNM(it_quad_fit{end}(cc),X_sh,tr_spkbns,'none');
        %     Xmat = [trial_Tmat g];
        %     klen = size(Xmat,2);
        %     K0 = zeros(klen+1,1);
        %     [fitp_sac,grad] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [], 0);
        %     sac_kern(cc,:) = fitp_sac.k(1:ntents);
        
        lamrange2 = [lambda 1 ntents 0; lambda ntents+1 2*ntents 0];
        cur_Robs = Robs(gray_back_real_inds);
        tr_spkbns = convert_to_spikebins(cur_Robs);
        [nll, pnll, lpen, prate, g] = getLL_GNM(it_quad_fit{end}(cc),X_sh(gray_back_real_inds,:),tr_spkbns,'none');
        Xmat = [sac_Tmat(gray_back_real_inds,:) msac_Tmat(gray_back_real_inds,:) g];
        klen = size(Xmat,2);
        K0 = zeros(klen+1,1);
        [fitp_gray(cc),grad,se] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [klen+1], 0);
        gray_sac_kern(cc,:) = fitp_gray(cc).k(1:ntents);
        gray_msac_kern(cc,:) = fitp_gray(cc).k(ntents+1:2*ntents);
        gray_sac_se(cc,:) = se(1:ntents);
        gray_msac_se(cc,:) = se(ntents+1:2*ntents);
        
        lamrange2 = [lambda 1 ntents 0; lambda ntents+1 2*ntents 0];
        cur_Robs = Robs(im_back_real_inds);
        tr_spkbns = convert_to_spikebins(cur_Robs);
        [nll, pnll, lpen, prate, g] = getLL_GNM(it_quad_fit{end}(cc),X_sh(im_back_real_inds,:),tr_spkbns,'none');
        Xmat = [sac_Tmat(im_back_real_inds,:) msac_Tmat(im_back_real_inds,:) g];
        K0 = zeros(klen+1,1);
        [fitp_im(cc),grad,se] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [klen+1], 0);
        im_sac_kern(cc,:) = fitp_im(cc).k(1:ntents);
        im_msac_kern(cc,:) = fitp_im(cc).k(ntents+1:2*ntents);
        im_sac_se(cc,:) = se(1:ntents);
        im_msac_se(cc,:) = se(ntents+1:2*ntents);
        
        lamrange2 = [lambda 1 ntents 0; lambda ntents+1 2*ntents 0];
        cur_Robs = Robs(im_back_sim_inds);
        tr_spkbns = convert_to_spikebins(cur_Robs);
        [nll, pnll, lpen, prate, g] = getLL_GNM(it_quad_fit{end}(cc),X_sh(im_back_sim_inds,:),tr_spkbns,'none');
        Xmat = [msac_Tmat(im_back_sim_inds,:) stim_Tmat(im_back_sim_inds,:) g];
        klen = size(Xmat,2);
        K0 = zeros(klen+1,1);
        [fitp_sim(cc),grad,se] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [klen+1], 0);
        sim_msac_kern(cc,:) = fitp_sim(cc).k(1:ntents);
        sim_stim_kern(cc,:) = fitp_sim(cc).k(ntents+1:2*ntents);
        sim_msac_se(cc,:) = se(1:ntents);
        sim_stim_se(cc,:) = se(ntents+1:2*ntents);
        
        lamrange2 = [lambda 1 ntents 0; lambda ntents+1 2*ntents 0; lambda 2*ntents+1 3*ntents 0];
        cur_Robs = Robs(gray_back_real_inds);
        tr_spkbns = convert_to_spikebins(cur_Robs);
        [nll, pnll, lpen, prate, g] = getLL_GNM(it_quad_fit{end}(cc),X_sh(gray_back_real_inds,:),tr_spkbns,'none');
        Xmat = [osac_Tmat(gray_back_real_inds,:) rsac_Tmat(gray_back_real_inds,:) msac_Tmat(gray_back_real_inds,:) g];
        klen = size(Xmat,2);
        K0 = zeros(klen+1,1);
        [fitp_gray(cc),grad,se] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [klen+1], 0);
        gray_osac_kern2(cc,:) = fitp_gray(cc).k(1:ntents);
        gray_rsac_kern2(cc,:) = fitp_gray(cc).k(ntents+1:2*ntents);
        gray_msac_kern2(cc,:) = fitp_gray(cc).k(2*ntents+1:3*ntents);
        gray_osac_se2(cc,:) = se(1:ntents);
        gray_rsac_se2(cc,:) = se(ntents+1:2*ntents);
        gray_msac_se2(cc,:) = se(2*ntents+1:3*ntents);
        
        lamrange2 = [lambda 1 ntents 0; lambda ntents+1 2*ntents 0; lambda 2*ntents+1 3*ntents 0];
        cur_Robs = Robs(im_back_real_inds);
        tr_spkbns = convert_to_spikebins(cur_Robs);
        [nll, pnll, lpen, prate, g] = getLL_GNM(it_quad_fit{end}(cc),X_sh(im_back_real_inds,:),tr_spkbns,'none');
        Xmat = [osac_Tmat(im_back_real_inds,:) rsac_Tmat(im_back_real_inds,:) msac_Tmat(im_back_real_inds,:) g];
        K0 = zeros(klen+1,1);
        [fitp_im(cc),grad,se] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [klen+1], 0);
        im_rsac_kern2(cc,:) = fitp_im(cc).k(1:ntents);
        im_osac_kern2(cc,:) = fitp_im(cc).k(ntents+1:2*ntents);
        im_msac_kern2(cc,:) = fitp_im(cc).k(2*ntents+1:3*ntents);
        im_osac_se2(cc,:) = se(1:ntents);
        im_rsac_se2(cc,:) = se(ntents+1:2*ntents);
        im_msac_se2(cc,:) = se(2*ntents+1:3*ntents);
        
        
        lamrange2 = [lambda 1 ntents 0; lambda ntents+1 2*ntents 0; lambda 2*ntents+1 3*ntents 0];
        cur_Robs = Robs(im_back_sim_inds);
        tr_spkbns = convert_to_spikebins(cur_Robs);
        [nll, pnll, lpen, prate, g] = getLL_GNM(it_quad_fit{end}(cc),X_sh(im_back_sim_inds,:),tr_spkbns,'none');
        Xmat = [onstim_Tmat(im_back_sim_inds,:) offstim_Tmat(im_back_sim_inds,:) msac_Tmat(im_back_sim_inds,:) g];
        klen = size(Xmat,2);
        K0 = zeros(klen+1,1);
        [fitp_sim(cc),grad,se] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [klen+1], 0);
        sim_onstim_kern2(cc,:) = fitp_sim(cc).k(1:ntents);
        sim_offstim_kern2(cc,:) = fitp_sim(cc).k(ntents+1:2*ntents);
        sim_msac_kern2(cc,:) = fitp_sim(cc).k(2*ntents+1:3*ntents);
        sim_onstim_se2(cc,:) = se(1:ntents);
        sim_offstim_se2(cc,:) = se(ntents+1:2*ntents);
        sim_msac_se2(cc,:) = se(2*ntents+1:3*ntents);
        
    end
    
    %%
    sname = sprintf('unit_sac_kern_mods_unsm_%ddeg',bar_oris);
    save(sname,'*_kern*','*_se*','ntents','tent_centers','dt','up_SDIM','flen')
end

%%
for EE = 1:4
    sname = sprintf('unit_sac_kern_mods_%ddeg',poss_bar_oris(EE));
    load(sname);
    
    gray_sac_avg(EE,:,:) = gray_sac_kern;
    gray_sac_sem(EE,:,:) = gray_sac_se;
    gray_osac_avg(EE,:,:) = gray_osac_kern2;
    gray_osac_sem(EE,:,:) = gray_osac_se2;
    gray_rsac_avg(EE,:,:) = gray_rsac_kern2;
    gray_rsac_sem(EE,:,:) = gray_rsac_se2;
    gray_msac_avg(EE,:,:) = gray_msac_kern2;
    gray_msac_sem(EE,:,:) = gray_msac_se2;
    
    im_sac_avg(EE,:,:) = im_sac_kern;
    im_sac_sem(EE,:,:) = im_sac_se;
    im_osac_avg(EE,:,:) = im_rsac_kern2; %o and r sacs were switched
    im_osac_sem(EE,:,:) = im_rsac_se2;
    im_rsac_avg(EE,:,:) = im_osac_kern2;
    im_rsac_sem(EE,:,:) = im_osac_se2;
    im_msac_avg(EE,:,:) = im_msac_kern2;
    im_msac_sem(EE,:,:) = im_msac_se2;
    
    sim_sac_avg(EE,:,:) = sim_stim_kern;
    sim_sac_sem(EE,:,:) = sim_stim_se;
    sim_osac_avg(EE,:,:) = sim_onstim_kern2;
    sim_osac_sem(EE,:,:) = sim_onstim_se2;
    sim_rsac_avg(EE,:,:) = sim_offstim_kern2;
    sim_rsac_sem(EE,:,:) = sim_offstim_se2;
    sim_msac_avg(EE,:,:) = sim_msac_kern2;
    sim_msac_sem(EE,:,:) = sim_msac_se2;
end

avg_gray_osac = squeeze(mean(gray_osac_avg(1:3,:,:),1));
avg_gray_rsac = squeeze(mean(gray_rsac_avg(1:3,:,:),1));
avg_gray_msac = squeeze(mean(gray_msac_avg(1:3,:,:),1));

avg_im_osac = squeeze(mean(im_osac_avg(1:3,:,:),1));
avg_im_rsac = squeeze(mean(im_rsac_avg(1:3,:,:),1));
avg_im_msac = squeeze(mean(im_msac_avg(1:3,:,:),1));

avg_sim_osac = squeeze(mean(sim_osac_avg(1:3,:,:),1));
avg_sim_rsac = squeeze(mean(sim_rsac_avg(1:3,:,:),1));
avg_sim_msac = squeeze(mean(sim_msac_avg(1:3,:,:),1));

avg_gray_sac = 0.5*avg_gray_osac + 0.5*avg_gray_rsac;
avg_im_sac = 0.5*avg_im_osac + 0.5*avg_im_rsac;

%%
load ./unit_ori_models
close all
for cc = 1:96
    cc
    f1 = figure(1);
    subplot(2,4,1);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(1,cc,:)),squeeze(gray_sac_sem(1,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(1,cc,:)),squeeze(im_sac_sem(1,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(2,4,2);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(2,cc,:)),squeeze(gray_sac_sem(2,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(2,cc,:)),squeeze(im_sac_sem(2,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(2,4,3);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(3,cc,:)),squeeze(gray_sac_sem(3,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(3,cc,:)),squeeze(im_sac_sem(3,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(2,4,4);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(4,cc,:)),squeeze(gray_sac_sem(4,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(4,cc,:)),squeeze(im_sac_sem(4,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    subplot(2,4,5);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_msac_avg(1,cc,:)),squeeze(gray_msac_sem(1,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_msac_avg(1,cc,:)),squeeze(im_msac_sem(1,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(2,4,6);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_msac_avg(2,cc,:)),squeeze(gray_msac_sem(2,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_msac_avg(2,cc,:)),squeeze(im_msac_sem(2,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(2,4,7);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_msac_avg(3,cc,:)),squeeze(gray_msac_sem(3,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_msac_avg(3,cc,:)),squeeze(im_msac_sem(3,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(2,4,8);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_msac_avg(4,cc,:)),squeeze(gray_msac_sem(4,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_msac_avg(4,cc,:)),squeeze(im_msac_sem(4,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    small = 0;
    big = 0;
    for i = 1:8
        subplot(2,4,i)
        cur_yl = ylim();
        if cur_yl(1) < small
            small = cur_yl(1);
        end
        if cur_yl(2) > big
            big = cur_yl(2);
        end
    end
    for i = 1:8
        subplot(2,4,i)
        ylim([small big])
    end
    
    set(f1,'Position',[100 650 1300 600])
    
    f2 = figure(2);
    cur_k = get_k_mat(quad_fit(cc));
    cur_k = reshape(cur_k,flen,SDIM);
    imagesc(un_bar_oris,1:flen,cur_k);
    mm = max(abs(cur_k(:)));
    caxis([-0.9*mm 0.9*mm])
    set(f2,'Position',[1300 650 600 400])
    
    pause
    close all
end

%%
close all
for cc = 1:96
    cc
    subplot(6,4,1);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(1,cc,:)),squeeze(gray_sac_sem(1,cc,:)),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,2);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(2,cc,:)),squeeze(gray_sac_sem(2,cc,:)),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,3);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(3,cc,:)),squeeze(gray_sac_sem(3,cc,:)),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,4);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(4,cc,:)),squeeze(gray_sac_sem(4,cc,:)),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    subplot(6,4,5);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_osac_avg(1,cc,:)),squeeze(gray_osac_sem(1,cc,:)),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(gray_rsac_avg(1,cc,:)),squeeze(gray_rsac_sem(1,cc,:)),{'ko-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,6);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_osac_avg(2,cc,:)),squeeze(gray_osac_sem(2,cc,:)),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(gray_rsac_avg(2,cc,:)),squeeze(gray_rsac_sem(2,cc,:)),{'ko-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,7);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_osac_avg(3,cc,:)),squeeze(gray_osac_sem(3,cc,:)),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(gray_rsac_avg(3,cc,:)),squeeze(gray_rsac_sem(3,cc,:)),{'ko-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,8);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_osac_avg(4,cc,:)),squeeze(gray_osac_sem(4,cc,:)),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(gray_rsac_avg(4,cc,:)),squeeze(gray_rsac_sem(4,cc,:)),{'ko-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    
    
    subplot(6,4,9);hold on
    shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(1,cc,:)),squeeze(im_sac_sem(1,cc,:)),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,10);hold on
    shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(2,cc,:)),squeeze(im_sac_sem(2,cc,:)),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,11);hold on
    shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(3,cc,:)),squeeze(im_sac_sem(3,cc,:)),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,12);hold on
    shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(4,cc,:)),squeeze(im_sac_sem(4,cc,:)),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    subplot(6,4,13);hold on
    shadedErrorBar(tent_centers*dt,squeeze(im_osac_avg(1,cc,:)),squeeze(im_osac_sem(1,cc,:)),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_rsac_avg(1,cc,:)),squeeze(im_rsac_sem(1,cc,:)),{'ko-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,14);hold on
    shadedErrorBar(tent_centers*dt,squeeze(im_osac_avg(2,cc,:)),squeeze(im_osac_sem(2,cc,:)),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_rsac_avg(2,cc,:)),squeeze(im_rsac_sem(2,cc,:)),{'ko-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,15);hold on
    shadedErrorBar(tent_centers*dt,squeeze(im_osac_avg(3,cc,:)),squeeze(im_osac_sem(3,cc,:)),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_rsac_avg(3,cc,:)),squeeze(im_rsac_sem(3,cc,:)),{'ko-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,16);hold on
    shadedErrorBar(tent_centers*dt,squeeze(im_osac_avg(4,cc,:)),squeeze(im_osac_sem(4,cc,:)),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_rsac_avg(4,cc,:)),squeeze(im_rsac_sem(4,cc,:)),{'ko-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    
    subplot(6,4,17);hold on
    shadedErrorBar(tent_centers*dt,squeeze(sim_sac_avg(1,cc,:)),squeeze(sim_sac_sem(1,cc,:)),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,18);hold on
    shadedErrorBar(tent_centers*dt,squeeze(sim_sac_avg(2,cc,:)),squeeze(sim_sac_sem(2,cc,:)),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,19);hold on
    shadedErrorBar(tent_centers*dt,squeeze(sim_sac_avg(3,cc,:)),squeeze(sim_sac_sem(3,cc,:)),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,20);hold on
    shadedErrorBar(tent_centers*dt,squeeze(sim_sac_avg(4,cc,:)),squeeze(sim_sac_sem(4,cc,:)),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    subplot(6,4,21);hold on
    shadedErrorBar(tent_centers*dt,squeeze(sim_osac_avg(1,cc,:)),squeeze(sim_osac_sem(1,cc,:)),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(sim_rsac_avg(1,cc,:)),squeeze(sim_rsac_sem(1,cc,:)),{'ko-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,22);hold on
    shadedErrorBar(tent_centers*dt,squeeze(sim_osac_avg(2,cc,:)),squeeze(sim_osac_sem(2,cc,:)),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(sim_rsac_avg(2,cc,:)),squeeze(sim_rsac_sem(2,cc,:)),{'ko-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,23);hold on
    shadedErrorBar(tent_centers*dt,squeeze(sim_osac_avg(3,cc,:)),squeeze(sim_osac_sem(3,cc,:)),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(sim_rsac_avg(3,cc,:)),squeeze(sim_rsac_sem(3,cc,:)),{'ko-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    subplot(6,4,24);hold on
    shadedErrorBar(tent_centers*dt,squeeze(sim_osac_avg(4,cc,:)),squeeze(sim_osac_sem(4,cc,:)),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(sim_rsac_avg(4,cc,:)),squeeze(sim_rsac_sem(4,cc,:)),{'ko-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    pause
    clf
end

%%
close all
for cc = 1:96
    cc
    subplot(3,4,1);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(1,cc,:)),squeeze(gray_sac_sem(1,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(gray_msac_avg(1,cc,:)),squeeze(gray_msac_sem(1,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    title('Gray back 0 deg')
    subplot(3,4,2);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(2,cc,:)),squeeze(gray_sac_sem(2,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(gray_msac_avg(2,cc,:)),squeeze(gray_msac_sem(2,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    title('Gray back 45 deg')
    subplot(3,4,3);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(3,cc,:)),squeeze(gray_sac_sem(3,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(gray_msac_avg(3,cc,:)),squeeze(gray_msac_sem(3,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    title('Gray back 90 deg')
    subplot(3,4,4);hold on
    shadedErrorBar(tent_centers*dt,squeeze(gray_sac_avg(4,cc,:)),squeeze(gray_sac_sem(4,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(gray_msac_avg(4,cc,:)),squeeze(gray_msac_sem(4,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    title('Gray back 135 deg')
    
    subplot(3,4,5);hold on
    shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(1,cc,:)),squeeze(im_sac_sem(1,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_msac_avg(1,cc,:)),squeeze(im_msac_sem(1,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    title('Image back 0 deg')
    subplot(3,4,6);hold on
    shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(2,cc,:)),squeeze(im_sac_sem(2,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_msac_avg(2,cc,:)),squeeze(im_msac_sem(2,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    title('Image back 45 deg')
    subplot(3,4,7);hold on
    shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(3,cc,:)),squeeze(im_sac_sem(3,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_msac_avg(3,cc,:)),squeeze(im_msac_sem(3,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    title('Image back 90 deg')
    subplot(3,4,8);hold on
    shadedErrorBar(tent_centers*dt,squeeze(im_sac_avg(4,cc,:)),squeeze(im_sac_sem(4,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(im_msac_avg(4,cc,:)),squeeze(im_msac_sem(4,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    title('Image back 135 deg')
    
    subplot(3,4,9);hold on
    shadedErrorBar(tent_centers*dt,squeeze(sim_sac_avg(1,cc,:)),squeeze(sim_sac_sem(1,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(sim_msac_avg(1,cc,:)),squeeze(sim_msac_sem(1,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    title('Sim sac 0 deg')
    subplot(3,4,10);hold on
    shadedErrorBar(tent_centers*dt,squeeze(sim_sac_avg(2,cc,:)),squeeze(sim_sac_sem(2,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(sim_msac_avg(2,cc,:)),squeeze(sim_msac_sem(2,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    title('Sim sac 45 deg')
    subplot(3,4,11);hold on
    shadedErrorBar(tent_centers*dt,squeeze(sim_sac_avg(3,cc,:)),squeeze(sim_sac_sem(3,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(sim_msac_avg(3,cc,:)),squeeze(sim_msac_sem(3,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    title('Sim sac 90 deg')
    subplot(3,4,12);hold on
    shadedErrorBar(tent_centers*dt,squeeze(sim_sac_avg(4,cc,:)),squeeze(sim_sac_sem(4,cc,:)),{'b.-'},1);
    shadedErrorBar(tent_centers*dt,squeeze(sim_msac_avg(4,cc,:)),squeeze(sim_msac_sem(4,cc,:)),{'r.-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    title('Sim sac 135 deg')
    
    big = 0; small = 0;
    for i = 1:12
        subplot(3,4,i)
        yl = ylim();
        if yl(1) < small
            small = yl(1);
        end
        if yl(2) > big
            big = yl(2);
        end
    end
    for i = 1:12
        subplot(3,4,i)
        ylim([small big])
    end
    
    pause
    clf
end

%%
figure
subplot(3,4,1); hold on
shadedErrorBar(tent_centers*dt,squeeze(mean(gray_osac_avg(1,:,:),2)),squeeze(std(gray_osac_avg(1,:,:),[],2))/sqrt(96),{'r.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(gray_rsac_avg(1,:,:),2)),squeeze(std(gray_rsac_avg(1,:,:),[],2))/sqrt(96),{'k.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(gray_msac_avg(1,:,:),2)),squeeze(std(gray_msac_avg(1,:,:),[],2))/sqrt(96),{'b.-'},1);
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
subplot(3,4,2); hold on
shadedErrorBar(tent_centers*dt,squeeze(mean(gray_osac_avg(2,:,:),2)),squeeze(std(gray_osac_avg(2,:,:),[],2))/sqrt(96),{'r.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(gray_rsac_avg(2,:,:),2)),squeeze(std(gray_rsac_avg(2,:,:),[],2))/sqrt(96),{'k.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(gray_msac_avg(2,:,:),2)),squeeze(std(gray_msac_avg(2,:,:),[],2))/sqrt(96),{'.-'},1);
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
subplot(3,4,3); hold on
shadedErrorBar(tent_centers*dt,squeeze(mean(gray_osac_avg(3,:,:),2)),squeeze(std(gray_osac_avg(3,:,:),[],2))/sqrt(96),{'r.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(gray_rsac_avg(3,:,:),2)),squeeze(std(gray_rsac_avg(3,:,:),[],2))/sqrt(96),{'k.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(gray_msac_avg(3,:,:),2)),squeeze(std(gray_msac_avg(3,:,:),[],2))/sqrt(96),{'.-'},1);
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
subplot(3,4,4); hold on
shadedErrorBar(tent_centers*dt,squeeze(mean(gray_osac_avg(4,:,:),2)),squeeze(std(gray_osac_avg(4,:,:),[],2))/sqrt(96),{'r.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(gray_rsac_avg(4,:,:),2)),squeeze(std(gray_rsac_avg(4,:,:),[],2))/sqrt(96),{'k.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(gray_msac_avg(4,:,:),2)),squeeze(std(gray_msac_avg(4,:,:),[],2))/sqrt(96),{'.-'},1);
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
subplot(3,4,5); hold on
shadedErrorBar(tent_centers*dt,squeeze(mean(im_osac_avg(1,:,:),2)),squeeze(std(im_osac_avg(1,:,:),[],2))/sqrt(96),{'r.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(im_rsac_avg(1,:,:),2)),squeeze(std(im_rsac_avg(1,:,:),[],2))/sqrt(96),{'k.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(im_msac_avg(1,:,:),2)),squeeze(std(im_msac_avg(1,:,:),[],2))/sqrt(96),{'.-'},1);
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
subplot(3,4,6); hold on
shadedErrorBar(tent_centers*dt,squeeze(mean(im_osac_avg(2,:,:),2)),squeeze(std(im_osac_avg(2,:,:),[],2))/sqrt(96),{'r.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(im_rsac_avg(2,:,:),2)),squeeze(std(im_rsac_avg(2,:,:),[],2))/sqrt(96),{'k.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(im_msac_avg(2,:,:),2)),squeeze(std(im_msac_avg(2,:,:),[],2))/sqrt(96),{'.-'},1);
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
subplot(3,4,7); hold on
shadedErrorBar(tent_centers*dt,squeeze(mean(im_osac_avg(3,:,:),2)),squeeze(std(im_osac_avg(3,:,:),[],2))/sqrt(96),{'r.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(im_rsac_avg(3,:,:),2)),squeeze(std(im_rsac_avg(3,:,:),[],2))/sqrt(96),{'k.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(im_msac_avg(3,:,:),2)),squeeze(std(im_msac_avg(3,:,:),[],2))/sqrt(96),{'.-'},1);
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
subplot(3,4,8); hold on
shadedErrorBar(tent_centers*dt,squeeze(mean(im_osac_avg(4,:,:),2)),squeeze(std(im_osac_avg(4,:,:),[],2))/sqrt(96),{'r.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(im_rsac_avg(4,:,:),2)),squeeze(std(im_rsac_avg(4,:,:),[],2))/sqrt(96),{'k.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(im_msac_avg(4,:,:),2)),squeeze(std(im_msac_avg(4,:,:),[],2))/sqrt(96),{'.-'},1);
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
subplot(3,4,9); hold on
shadedErrorBar(tent_centers*dt,squeeze(mean(sim_osac_avg(1,:,:),2)),squeeze(std(sim_osac_avg(1,:,:),[],2))/sqrt(96),{'r.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(sim_rsac_avg(1,:,:),2)),squeeze(std(sim_rsac_avg(1,:,:),[],2))/sqrt(96),{'k.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(sim_msac_avg(1,:,:),2)),squeeze(std(sim_msac_avg(1,:,:),[],2))/sqrt(96),{'.-'},1);
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
subplot(3,4,10); hold on
shadedErrorBar(tent_centers*dt,squeeze(mean(sim_osac_avg(2,:,:),2)),squeeze(std(sim_osac_avg(2,:,:),[],2))/sqrt(96),{'r.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(sim_rsac_avg(2,:,:),2)),squeeze(std(sim_rsac_avg(2,:,:),[],2))/sqrt(96),{'k.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(sim_msac_avg(2,:,:),2)),squeeze(std(sim_msac_avg(2,:,:),[],2))/sqrt(96),{'.-'},1);
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
subplot(3,4,11); hold on
shadedErrorBar(tent_centers*dt,squeeze(mean(sim_osac_avg(3,:,:),2)),squeeze(std(sim_osac_avg(3,:,:),[],2))/sqrt(96),{'r.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(sim_rsac_avg(3,:,:),2)),squeeze(std(sim_rsac_avg(3,:,:),[],2))/sqrt(96),{'k.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(sim_msac_avg(3,:,:),2)),squeeze(std(sim_msac_avg(3,:,:),[],2))/sqrt(96),{'.-'},1);
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
subplot(3,4,12); hold on
shadedErrorBar(tent_centers*dt,squeeze(mean(sim_osac_avg(4,:,:),2)),squeeze(std(sim_osac_avg(4,:,:),[],2))/sqrt(96),{'r.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(sim_rsac_avg(4,:,:),2)),squeeze(std(sim_rsac_avg(4,:,:),[],2))/sqrt(96),{'k.-'},1);
shadedErrorBar(tent_centers*dt,squeeze(mean(sim_msac_avg(4,:,:),2)),squeeze(std(sim_msac_avg(4,:,:),[],2))/sqrt(96),{'.-'},1);
axis tight
xl = xlim();
line(xl,[0 0],'color','k')

%%

figure;
subplot(3,1,1)
hold on
shadedErrorBar(tent_centers*dt,mean(avg_gray_osac),std(avg_gray_osac)/sqrt(96),{'r.-'})
shadedErrorBar(tent_centers*dt,mean(avg_gray_rsac),std(avg_gray_rsac)/sqrt(96),{'k.-'})
shadedErrorBar(tent_centers*dt,mean(avg_gray_msac),std(avg_gray_msac)/sqrt(96),{'.-'})
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
title('Gray background')
xlabel('Time (s)')
subplot(3,1,2)
hold on
shadedErrorBar(tent_centers*dt,mean(avg_im_osac),std(avg_im_osac)/sqrt(96),{'r.-'})
shadedErrorBar(tent_centers*dt,mean(avg_im_rsac),std(avg_im_rsac)/sqrt(96),{'k.-'})
shadedErrorBar(tent_centers*dt,mean(avg_im_msac),std(avg_im_msac)/sqrt(96),{'.-'})
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
title('Image background')
xlabel('Time (s)')
subplot(3,1,3)
hold on
shadedErrorBar(tent_centers*dt,mean(avg_sim_osac),std(avg_sim_osac)/sqrt(96),{'r.-'})
shadedErrorBar(tent_centers*dt,mean(avg_sim_rsac),std(avg_sim_rsac)/sqrt(96),{'k.-'})
shadedErrorBar(tent_centers*dt,mean(avg_sim_msac),std(avg_sim_msac)/sqrt(96),{'.-'})
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
title('Sim saccades')
xlabel('Time (s)')


figure;
subplot(3,1,1)
hold on
shadedErrorBar(tent_centers*dt,mean(avg_gray_osac),std(avg_gray_osac)/sqrt(96),{'r.-'})
shadedErrorBar(tent_centers*dt,mean(avg_im_osac),std(avg_im_osac)/sqrt(96),{'b.-'})
shadedErrorBar(tent_centers*dt,mean(avg_sim_osac),std(avg_sim_osac)/sqrt(96),{'k.-'})
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
title('Outgoing saccades')
xlabel('Time (s)')
subplot(3,1,2)
hold on
shadedErrorBar(tent_centers*dt,mean(avg_gray_rsac),std(avg_gray_rsac)/sqrt(96),{'r.-'})
shadedErrorBar(tent_centers*dt,mean(avg_im_rsac),std(avg_im_rsac)/sqrt(96),{'b.-'})
shadedErrorBar(tent_centers*dt,mean(avg_sim_rsac),std(avg_sim_rsac)/sqrt(96),{'k.-'})
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
title('Returning saccades')
xlabel('Time (s)')
subplot(3,1,3)
hold on
shadedErrorBar(tent_centers*dt,mean(avg_gray_msac),std(avg_gray_msac)/sqrt(96),{'r.-'})
shadedErrorBar(tent_centers*dt,mean(avg_im_msac),std(avg_im_msac)/sqrt(96),{'b.-'})
shadedErrorBar(tent_centers*dt,mean(avg_sim_msac),std(avg_sim_msac)/sqrt(96),{'k.-'})
axis tight
xl = xlim();
line(xl,[0 0],'color','k')
title('Microsaccades')
xlabel('Time (s)')

%%
up_bar_pos = (-floor(up_SDIM/2):floor(up_SDIM/2))*0.04;
close all
% interesting_us = [2 8 18 30 36 37 65 67 68 78 80 81 84 88 91 95 96];

%%
for cc = 1:96
    % for cc = good_sus
    % for cc = interesting_us
    cc
    subplot(2,3,1)
    hold on
    shadedErrorBar(tent_centers*dt,gray_sac_kern(cc,:),gray_sac_se(cc,:),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    subplot(2,3,4)
    shadedErrorBar(tent_centers*dt,gray_osac_kern2(cc,:),gray_osac_se2(cc,:),{'ro-'},1);
    hold on
    shadedErrorBar(tent_centers*dt,gray_rsac_kern2(cc,:),gray_rsac_se2(cc,:),{'ko-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    subplot(2,3,2)
    hold on
    shadedErrorBar(tent_centers*dt,im_sac_kern(cc,:),im_sac_se(cc,:),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    subplot(2,3,5)
    shadedErrorBar(tent_centers*dt,im_osac_kern2(cc,:),im_osac_se2(cc,:),{'ko-'},1);
    hold on
    shadedErrorBar(tent_centers*dt,im_rsac_kern2(cc,:),im_rsac_se2(cc,:),{'ro-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    
    subplot(2,3,3)
    hold on
    shadedErrorBar(tent_centers*dt,sim_stim_kern(cc,:),sim_stim_se(cc,:),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    subplot(2,3,6)
    shadedErrorBar(tent_centers*dt,sim_onstim_kern2(cc,:),sim_onstim_se2(cc,:),{'ko-'},1);
    hold on
    shadedErrorBar(tent_centers*dt,sim_offstim_kern2(cc,:),sim_offstim_se2(cc,:),{'ro-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    pause
    clf
    
    % sname = ['~/Desktop/' sprintf('Unit%d_sackerns',cc)];
    % fillPage(gcf,'papersize',[10 10]);
    % print(sname,'-dpng');
    % close
    
end

%%
up_bar_pos = (-floor(up_SDIM/2):floor(up_SDIM/2))*0.04;
close all
interesting_us = [2 8 18 30 36 37 65 67 78 80 81 84 88 91 95 96];
for cc = 1:96
    % for cc = good_sus
    % for cc = interesting_us
    cc
    k = get_k_mat(it_quad_fit{end}(cc));
    %     subplot(2,2,1)
    %     imagesc(up_bar_pos,up_bar_pos,reshape(k(:,1),flen,up_SDIM));
    %     ca = max(abs(caxis()));
    %     caxis([-ca ca]);%colorbar
    %     title('Linear')
    %     subplot(2,2,3)
    %     imagesc(up_bar_pos,up_bar_pos,reshape(k(:,2),flen,up_SDIM))
    %     ca = max(abs(caxis()));
    %     caxis([-ca ca]);%colorbar
    %     title('Quad')
    subplot(2,2,1)
    %     plot(tent_centers*dt,sac_kern(cc,:),'o-')
    hold on
    plot(tent_centers*dt,gray_sac_kern(cc,:),'ro-','linewidth',1);
    plot(tent_centers*dt,im_sac_kern(cc,:),'ko-','linewidth',1);
    plot(tent_centers*dt,sim_stim_kern(cc,:),'bo-','linewidth',1);
    legend('Gray back','Im back','Sim sac','Location','Southwest')
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    xlabel('Time lag (s)','fontsize',16)
    title('Saccade kernels')
    
    subplot(2,2,2)
    %     plot(tent_centers*dt,sac_kern(cc,:),'o-')
    hold on
    shadedErrorBar(tent_centers*dt,gray_sac_kern(cc,:),gray_sac_se(cc,:),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,im_sac_kern(cc,:),im_sac_se(cc,:),{'ko-'},1);
    shadedErrorBar(tent_centers*dt,sim_stim_kern(cc,:),sim_stim_se(cc,:),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    xlabel('Time lag (s)','fontsize',16)
    title('Saccade kernels')
    
    subplot(2,2,3)
    hold on
    plot(tent_centers*dt,gray_msac_kern(cc,:),'ro-','linewidth',1);
    plot(tent_centers*dt,im_msac_kern(cc,:),'ko-','linewidth',1);
    plot(tent_centers*dt,sim_msac_kern(cc,:),'bo-','linewidth',1);
    legend('Gray back','Im back','Sim sac','Location','Southwest')
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    xlabel('Time lag (s)','fontsize',16)
    title('Microaccade kernels')
    
    subplot(2,2,4)
    shadedErrorBar(tent_centers*dt,gray_msac_kern(cc,:),gray_msac_se(cc,:),{'ro-'},1);
    hold on
    shadedErrorBar(tent_centers*dt,im_msac_kern(cc,:),im_msac_se(cc,:),{'ko-'},1);
    shadedErrorBar(tent_centers*dt,sim_msac_kern(cc,:),sim_msac_se(cc,:),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    xlabel('Time lag (s)','fontsize',16)
    title('Microsaccade kernels')
    
    pause
    clf
    % sname = ['~/Desktop/' sprintf('Unit%d_sackerns',cc)];
    % fillPage(gcf,'papersize',[10 10]);
    % print(sname,'-dpng');
    % close
end

%%
for cc = 1:96
    % for cc = good_sus
    cc
    subplot(3,1,1)
    hold on
    shadedErrorBar(tent_centers*dt,gray_osac_kern2(cc,:),gray_osac_se2(cc,:),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,gray_rsac_kern2(cc,:),gray_rsac_se2(cc,:),{'ko-'},1);
    shadedErrorBar(tent_centers*dt,gray_msac_kern2(cc,:),gray_msac_se2(cc,:),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    subplot(3,1,2)
    hold on
    shadedErrorBar(tent_centers*dt,im_osac_kern2(cc,:),im_osac_se2(cc,:),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,im_rsac_kern2(cc,:),im_rsac_se2(cc,:),{'ko-'},1);
    shadedErrorBar(tent_centers*dt,im_msac_kern2(cc,:),im_msac_se2(cc,:),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    subplot(3,1,3)
    hold on
    shadedErrorBar(tent_centers*dt,sim_onstim_kern2(cc,:),sim_onstim_se2(cc,:),{'ro-'},1);
    shadedErrorBar(tent_centers*dt,sim_offstim_kern2(cc,:),sim_offstim_se2(cc,:),{'ko-'},1);
    shadedErrorBar(tent_centers*dt,sim_msac_kern2(cc,:),sim_msac_se2(cc,:),{'bo-'},1);
    axis tight
    xl = xlim();
    line(xl,[0 0],'color','k')
    
    pause
    clf
end
