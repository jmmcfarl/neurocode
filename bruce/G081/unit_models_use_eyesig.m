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

%% USE ONLY GRAY BACKGROUND DATA
flen = 15;
beg_buffer = round(stim_fs*0.15);
bar_oris = [0];
un_bar_pos = all_un_bar_pos(:,1);

fprintf('Analyzing %d ori expts\n',bar_oris);

%     cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE));

%expts with X deg bars and gray back (sim sacs)
cur_expt_set = find(is_bar_expt==1 & expt_image_back == 0 & expt_bar_ori == bar_oris);
cur_expt_set(cur_expt_set<=8) = []; %this expt had continuous bar luminance
cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?

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
trial_cnt = 0;
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
        all_trialvec = [all_trialvec; trial_cnt + ones(size(cur_stim_times))*tt];
    end
    trial_cnt = trial_cnt + n_trials;
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
    clear sm_avg_eyepos eye_vel
    sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),3);
    sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),3);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
    
    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
    
    all_eye_vals = [all_eye_vals; lEyeXY];
    all_eye_speed = [all_eye_speed; eye_speed];
    all_eye_ts = [all_eye_ts; eye_ts_interp'];
end
interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_stim_times);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_stim_times);

%%
orth_eye_pos = interp_eye_vals(:,2);
bad_pts = find(abs(orth_eye_pos) > 1);
bad_trials = unique(all_trialvec(bad_pts));
corrected_Op = all_Op - orth_eye_pos;
cor_bar_pos = linspace(-1.4,0.4,23);
cor_nbar_pos = length(cor_bar_pos);

bar_dists = nan(length(all_trialvec),cor_nbar_pos);
for bb = 1:cor_nbar_pos
    bar_dists(:,bb) = abs(corrected_Op - cor_bar_pos(bb));
end
[~,bar_inds] = min(bar_dists,[],2);

trial_flips = [1; 1+find(diff(all_trialvec) ~= 0); length(all_trialvec)+1];
all_cor_barmat = [];
for nn = 1:length(trial_flips)-1
    cur_inds = trial_flips(nn):(trial_flips(nn+1)-1);
    cur_bar_mat = zeros(length(cur_inds),cor_nbar_pos);
    for bb = 1:cor_nbar_pos
        cur_set = find(bar_inds(cur_inds)==(bb));
        pset = all_phase(cur_inds(cur_set)) == 0;
        nset = all_phase(cur_inds(cur_set)) == pi;
        cur_bar_mat(cur_set(pset),bb) = 1;
        cur_bar_mat(cur_set(nset),bb) = -1;      
    end
    bar_Xmat = makeStimRows(cur_bar_mat,flen);
    all_cor_barmat = [all_cor_barmat; bar_Xmat];
end

%% IF YOU WANT TO DO CROSS-VALIDATION ON INITIAL MODEL FITS, OTHERWISE, 
tr_inds = find(all_used_inds==1 & ~ismember(all_trialvec,bad_trials));
NT = length(tr_inds);
SDIM = n_bar_pos;
cSDIM = cor_nbar_pos;

%% FIT models assuming perfect fixation
load ./CellList.mat
good_sus = find(all(CellList(:,:,1) > 0));

klen = flen*SDIM;
stim_params.spatial_dims = 1;
stim_params.sdim = SDIM;
stim_params.flen = flen;
cklen = flen*cSDIM;
cstim_params.spatial_dims = 1;
cstim_params.sdim = cSDIM;
cstim_params.flen = flen;

for cc = 1:length(good_sus)
    fprintf('Fitting cell %d of %d\n',cc,length(good_sus));
    tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,good_sus(cc)));
    
    clear defmod
    defmod.lambda_L1x = 10;
    defmod.lambda_d2XT = 100;
    defmod.lambda_d2X = 250;
    defmod.lambda_d2T = 250;
    npq = 1;
    nnq = 0;
    clear kern_types
    kern_types{1} = 'lin';
    for i = 2:(1+npq+nnq)
        kern_types{i} = 'quad';
    end
    init_kerns = 0.01*randn(klen,1+npq+nnq);
    init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
    quad = createGNM(init_kerns,[1 ones(1,npq) -1*ones(1,nnq)],kern_types,defmod,stim_params);
    for i = 2:(1+npq+nnq)
        quad.mods(i).lambda_L1x = 1;
        quad.mods(i).lambda_d2XT = 10;
        quad.mods(i).lambda_d2X = 25;
        quad.mods(i).lambda_d2T = 25;
    end
    quad_fit(cc) = fitGNM_filters(quad,all_bar_mat(tr_inds,:),tr_spkbns,'none',[],1e-4,1e-6,0);
    orig_LL(cc) = quad_fit(cc).LL;
    
    Robs = all_binned_spks(tr_inds,good_sus(cc));
    avg_rate = mean(Robs);
    null_prate = ones(length(tr_inds),1)*avg_rate;  
    null_LL(cc) = -sum(Robs.*log(null_prate) - null_prate)/sum(Robs);
    
    
    init_kerns = 0.01*randn(cklen,1+npq+nnq);
    init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
    quad = createGNM(init_kerns,[1 ones(1,npq) -1*ones(1,nnq)],kern_types,defmod,cstim_params);
    for i = 2:(1+npq+nnq)
        quad.mods(i).lambda_L1x = 1;
        quad.mods(i).lambda_d2XT = 10;
        quad.mods(i).lambda_d2X = 25;
        quad.mods(i).lambda_d2T = 25;
    end
    cquad_fit(cc) = fitGNM_filters(quad,all_cor_barmat(tr_inds,:),tr_spkbns,'none',[],1e-4,1e-6,0);
    corig_LL(cc) = cquad_fit(cc).LL;
    
end
