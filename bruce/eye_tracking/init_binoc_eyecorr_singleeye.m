clear all
% close all
addpath('~/James_scripts/bruce/G081/');
addpath('~/James_scripts/bruce/7_15_scripts/pre_nrsa/')
addpath('~/James_scripts/bruce/initial_analyses')

Expt_name = 'G086';
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

anal_name = 'binoc_eyecorr';
refit_mods = 0;

%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 10;
nPix = 30;
stim_params = NIMcreate_stim_params([flen 2*nPix],dt);
Fr = 1;
bar_ori = 90;

%%
exclude_expts = {'rls.orXme'};
for ii = 1:length(Expts)
    expt_names{ii} = Expts{ii}.Header.expname;
    expt_dds(ii) = Expts{ii}.Stimvals.dd;
    expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
    expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
    expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
    excluded_type(ii) = strcmp(expt_names{ii},exclude_expts);
end

% cur_expt_set = find(expt_binoc(:) == 1 & expt_bar_ori(:) == 90 & expt_sac_dir(:) == 90 & ...
%     expt_dds(:) == 67 & expt_Fr(:) == 1 & expt_npix(:) == 36 & ~excluded_type(:));
cur_expt_set = find(expt_binoc(:) == 1 & expt_bar_ori(:) == bar_ori & expt_sac_dir(:) == 90 & ...
    expt_dds(:) == 12 & expt_Fr(:) == 1 & ~excluded_type(:));

%% COMPUTE TRIAL DATA
all_stim_times = [];
all_used_inds = false(0);
Xmat = [];
all_binned_spks = [];
all_exptvec = [];
all_trialvec = [];
all_tsince_start = [];
all_trial_start_times = [];
all_trial_end_times = [];
trial_cnt = 0;
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_expt);
    load(fname);
    
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_expt}.Trials(:).id];
    [un_ids,id_inds] = unique(trial_ids);
    
    buffer_pix = floor((expt_npix(cur_expt) - nPix)/2);
    cur_use_pix = (1:nPix) + buffer_pix;
    
    trial_durs = trial_durs(id_inds);
    trial_start_times = trial_start_times(id_inds);
    trial_end_times = trial_end_times(id_inds);
    
    all_trial_start_times = [all_trial_start_times trial_start_times];
    all_trial_end_times = [all_trial_end_times trial_end_times];
    
    use_trials = find(trial_durs >= 0.5);
    extra_trials = find(use_trials > length(left_stim_mats));
    if ~isempty(extra_trials)
        fprintf('Dropping %d trials with missing stim data\n',length(extra_trials));
        use_trials(extra_trials) = [];
    end
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_expt}.Trials(use_trials(tt)).Start/1e4;
        
        if length(cur_stim_times) == 1 %for blocks where the trial stimulus time stamps arent recorded
            n_frames = size(left_stim_mats{use_trials(tt)},1);
            cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
            cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
        else %for trials where we do have the stim time stamps recorded
            cur_t_edges = [cur_stim_times; Expts{cur_expt}.Trials(use_trials(tt)).End(end)/1e4];
        end
        
        cur_stim_mat = zeros(length(cur_stim_times),nPix*2);
        cur_stim_mat(:,1:nPix) = left_stim_mats{use_trials(tt)}(:,cur_use_pix);
        cur_stim_mat(:,nPix+1:end) = right_stim_mats{use_trials(tt)}(:,cur_use_pix);
        cur_lstim_mat = left_stim_mats{use_trials(tt)}(:,cur_use_pix);
        cur_rstim_mat = right_stim_mats{use_trials(tt)}(:,cur_use_pix);
        
        cur_binned_spks = nan(length(cur_stim_times),96);
        for cc = 1:96
            cur_hist = histc(Clusters{cc}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
            %             temp = convert_to_spikebins(cur_hist(1:end-1));
        end
        
        bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
        
        cur_used_inds = true(length(cur_stim_times),1);
        cur_used_inds(1:flen) = false;
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_tsince_start = [all_tsince_start; cur_stim_times - trial_start_times(tt)];
        all_used_inds = [all_used_inds; cur_used_inds];
        Xmat = [Xmat; bar_Xmat];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(size(cur_stim_times))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_stim_times))*(tt + trial_cnt)];
    end
    trial_cnt = trial_cnt + 1;
end
all_t_axis = all_stim_times;

%% normalize stimulus variance
Xmat = Xmat/std(reshape(Xmat(all_used_inds,:),1,[]));

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
emfile = [data_dir '/jbe' Expt_name '.em.mat'];
load(emfile);

all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
all_eye_exptvec = [];
eye_smooth = 3;
for ee = 1:length(cur_expt_set);
    fprintf('Loading eye data for expt %d of %d\n',ee,length(cur_expt_set));
    cur_set = find(all_exptvec==ee);
    [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_t_axis(cur_set([1 end])));
    
    eye_dt = Expt.Header.CRrates(1);
    eye_fs = 1/eye_dt;
    lEyeXY = eye_vals_interp(:,1:2);
    rEyeXY = eye_vals_interp(:,3:4);
    
    %slight smoothing before computing speed
    sm_avg_eyepos = lEyeXY; eye_vel = lEyeXY; %initialization
    sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),eye_smooth);
    sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),eye_smooth);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
    
    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
    
    all_eye_vals = [all_eye_vals; lEyeXY rEyeXY];
    all_eye_speed = [all_eye_speed; eye_speed];
    all_eye_ts = [all_eye_ts; eye_ts_interp'];
    all_eye_exptvec = [all_eye_exptvec; ee*ones(size(eye_speed))];
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

interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_t_axis);

orth_eye_pos = interp_eye_vals(:,1)*cos(bar_ori*pi/180+pi/2) + interp_eye_vals(:,2)*sin(bar_ori*pi/180 + pi/2);
par_eye_pos = interp_eye_vals(:,1)*cos(bar_ori*pi/180) + interp_eye_vals(:,2)*sin(bar_ori*pi/180);
bad_pts = find(abs(orth_eye_pos) > 1); %cheap way of detecting blinks because the right eye signal was unreliable even for that at times

%% detect saccades
sac_thresh = 10;
peak_sig = [0; diff(sign(diff(all_eye_speed))); 0];
saccade_inds = find(peak_sig == -2 & all_eye_speed > sac_thresh);

peri_thresh = 3; %threshold eye speed for defining saccade boundary inds
thresh_cross_up = 1 + find(all_eye_speed(1:end-1) < peri_thresh & all_eye_speed(2:end) >= peri_thresh);
thresh_cross_down = 1 + find(all_eye_speed(1:end-1) >= peri_thresh & all_eye_speed(2:end) < peri_thresh);
sac_start_inds = nan(size(saccade_inds));
sac_stop_inds = nan(size(saccade_inds));
for ii = 1:length(saccade_inds)
    next_tc = find(thresh_cross_down > saccade_inds(ii),1,'first');
    if ~isempty(next_tc)
        sac_stop_inds(ii) = thresh_cross_down(next_tc);
    end
    prev_tc = find(thresh_cross_up < saccade_inds(ii),1,'last');
    if ~isempty(prev_tc)
        sac_start_inds(ii) = thresh_cross_up(prev_tc);
    end
    
end

%get rid of double-peaks
min_isi = 0.05; max_isi = Inf;
isis = [Inf; diff(sac_start_inds)]/eye_fs;
bad_isis = (isis < min_isi | isis > max_isi);
bad_sacs = find(isnan(sac_stop_inds) | isnan(sac_start_inds) | bad_isis);
saccade_inds(bad_sacs) = []; isis(bad_sacs) = []; sac_start_inds(bad_sacs) = []; sac_stop_inds(bad_sacs) = [];

saccade_times = all_eye_ts(saccade_inds);
sac_start_times = all_eye_ts(sac_start_inds);
sac_stop_times = all_eye_ts(sac_stop_inds);
sac_durs = sac_stop_times - sac_start_times;

sac_dbuff = round(0.005/eye_dt);
pre_inds = saccade_inds - sac_dbuff;
pre_inds(pre_inds < 1) = 1;
sac_pre_pos = all_eye_vals(pre_inds,:);
post_inds = saccade_inds + sac_dbuff;
post_inds(post_inds > length(all_eye_ts)) = length(all_eye_ts);
sac_post_pos = all_eye_vals(post_inds,:);

%use only left-eye signal here
sac_delta_pos = sac_post_pos(:,1:2) - sac_pre_pos(:,1:2);
sac_amps = sqrt(sum(sac_delta_pos.^2,2));
sac_dirs = atan2(sac_delta_pos(:,2),sac_delta_pos(:,1));

temp = ones(length(saccade_times),1);
saccades = struct('peak_time',mat2cell(saccade_times,temp),'start_time',mat2cell(sac_start_times,temp),...
    'stop_time',mat2cell(sac_stop_times,temp),'isi',mat2cell(isis,temp),...
    'duration',mat2cell(sac_durs,temp),'amplitude',mat2cell(sac_amps,temp),'direction',mat2cell(sac_dirs,temp),...
    'pre_Lx',mat2cell(sac_pre_pos(:,1),temp),'post_Lx',mat2cell(sac_post_pos(:,1),temp),...
    'pre_Ly',mat2cell(sac_pre_pos(:,2),temp),'post_Ly',mat2cell(sac_post_pos(:,2),temp),...
    'pre_Rx',mat2cell(sac_pre_pos(:,3),temp),'post_Rx',mat2cell(sac_post_pos(:,3),temp),...
    'pre_Ry',mat2cell(sac_pre_pos(:,4),temp),'post_Ry',mat2cell(sac_post_pos(:,4),temp));

% is_blink = find(sac_durs > 0.1);
% saccades(is_blink) = [];

sac_start_times = [saccades(:).start_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
bad = find(isnan(interp_sac_start_inds));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = [];
saccades(bad) = []; sac_start_times(bad) = [];

%throw out sacs where the interpolated time is too far from the true time
interp_sac_start_times = all_t_axis(interp_sac_start_inds);
time_diff = abs(interp_sac_start_times(:) - sac_start_times(:));
bad = find(time_diff > dt);
saccades(bad) = [];
interp_sac_start_inds(bad) = [];

%% SET UP USEABLE TRIALS
% [c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
[c,ia,ic] = unique(all_trialvec);
n_trials = length(ia);
% bad_trials = unique(ic(bad_pts)); %trials with putative blinks

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(ic,xv_set));
tr_inds = find(~ismember(ic,xv_set))';

tr_inds(~all_used_inds(tr_inds)) = [];
NT = length(tr_inds);

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
Xexpt = zeros(length(all_stim_times),length(cur_expt_set)-1);
for i = 1:length(cur_expt_set)-1
    cur_set = find(all_exptvec==i);
    Xexpt(cur_set,i) = 1;
end

%% FIT OR LOAD INITIAL MODELS
save_dir = [anal_dir '/' anal_name];
if ~exist(save_dir,'dir')
    system(['mkdir ' save_dir]);
end
init_mods_fname = [save_dir '/init_mods_sepeyef10_d22' '.mat'];

klen = flen*nPix;
stim_params = NIMcreate_stim_params([flen nPix],dt,1,1,length(cur_expt_set)-1);
optim_params.progTol = 1e-7;
optim_params.optTol = 1e-4;
silent = 1;
NL_types = {'lin', 'quad'};

XmatL = Xmat(:,1:klen);
XmatR = Xmat(:,(klen+1:end));

l_d2XT = 15000;
l_L1 = 150;
ql_L1 = 40;
ql_d2XT = 400;
XTmix = [1 1];

reg_params = NIMcreate_reg_params('lambda_d2XT',l_d2XT,'lambda_L1',l_L1,'temporal_boundaries','zero');
null_params = NIMcreate_stim_params([1 size(Xexpt,2)],dt);

if refit_mods == 1 || ~exist(init_mods_fname,'file')
    fprintf('Computing initial models\n');
    for cc = 1:96
        Robs = all_binned_spks(tr_inds,cc);
        null_mod(cc) = NIMinitialize_model(null_params,1,{'lin'},reg_params); %initialize NIM
        null_mod(cc) = NIMfit_filters(null_mod(cc),Robs,Xexpt(tr_inds,:)); %fit stimulus filters
        
        if xv_frac > 0
            Robsxv = all_binned_spks(xv_inds,cc);
            null_xvLL(cc) = NIMmodel_eval(null_mod(cc),Robsxv,Xexpt(xv_inds,:));
        end
        
        mod_signs = 1;
        fitLeft(cc,1) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
        fitLeft(cc,1) = NIMfit_filters(fitLeft(cc,1),Robs,XmatL(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
        fitRight(cc,1) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
        fitRight(cc,1) = NIMfit_filters(fitRight(cc,1),Robs,XmatR(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
        if xv_frac > 0
            [xvLL_left(cc,1),~,xv_predrate] = NIMmodel_eval(fitLeft(cc,1),Robsxv,XmatL(xv_inds,:),Xexpt(xv_inds,:));
            xv_imp_left(cc,1) = (xvLL_left(cc,1)-null_xvLL(cc))/log(2);
            [xvLL_right(cc,1),~,xv_predrate] = NIMmodel_eval(fitRight(cc,1),Robsxv,XmatR(xv_inds,:),Xexpt(xv_inds,:));
            xv_imp_right(cc,1) = (xvLL_right(cc,1)-null_xvLL(cc))/log(2);
            %         fprintf('LL imp: %.3f\n',xv_imp(cc,1));
        end
        
        mod_signs = [1 1];
        %add a filter
        rand_filt = randn(prod(stim_params.stim_dims),1)/prod(stim_params.stim_dims) * 1;
        fitLeft(cc,2) = NIMadd_NLinput(fitLeft(cc,1),'quad',1,rand_filt);
        %adjust regularization of new filter
        fitLeft(cc,2) = NIMadjust_regularization(fitLeft(cc,2),[length(mod_signs)],'lambda_d2XT',ql_d2XT);
        fitLeft(cc,2) = NIMadjust_regularization(fitLeft(cc,2),[length(mod_signs)],'lambda_L1',0);
        fitRight(cc,2) = NIMadd_NLinput(fitRight(cc,1),'quad',1,rand_filt);
        %adjust regularization of new filter
        fitRight(cc,2) = NIMadjust_regularization(fitRight(cc,2),[length(mod_signs)],'lambda_d2XT',ql_d2XT);
        fitRight(cc,2) = NIMadjust_regularization(fitRight(cc,2),[length(mod_signs)],'lambda_L1',0);
        
        fitLeft(cc,2) = NIMfit_filters(fitLeft(cc,2),Robs,XmatL(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
        fitRight(cc,2) = NIMfit_filters(fitRight(cc,2),Robs,XmatR(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
        
        %add in L1 for new filter and re-fit
        fitLeft(cc,2) = NIMadjust_regularization(fitLeft(cc,2),[length(mod_signs)],'lambda_L1',ql_L1);
        fitLeft(cc,2) = NIMfit_filters(fitLeft(cc,2),Robs,XmatL(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
        fitRight(cc,2) = NIMadjust_regularization(fitRight(cc,2),[length(mod_signs)],'lambda_L1',ql_L1);
        fitRight(cc,2) = NIMfit_filters(fitRight(cc,2),Robs,XmatR(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
        
        if xv_frac > 0
            [xvLL_left(cc,2),~,xv_predrate] = NIMmodel_eval(fitLeft(cc,2),Robsxv,XmatL(xv_inds,:),Xexpt(xv_inds,:));
            xv_imp_left(cc,2) = (xvLL_left(cc,2)-null_xvLL(cc))/log(2);
            [xvLL_right(cc,2),~,xv_predrate] = NIMmodel_eval(fitRight(cc,2),Robsxv,XmatR(xv_inds,:),Xexpt(xv_inds,:));
            xv_imp_right(cc,2) = (xvLL_right(cc,2)-null_xvLL(cc))/log(2);
            fprintf('LL imp 1filt: L: %.4f R: %.4f  2filt: L: %.4f R: %.4f\n',xv_imp_left(cc,1),xv_imp_right(cc,1),xv_imp_left(cc,2),xv_imp_right(cc,2));
            %     [fh] = NIMcheck_nonstat(fit0(cc),Robsxv,all_bar_mat(xv_inds,:),Xexpt(xv_inds,:),25);
        end
        save(init_mods_fname,'fitLeft','fitRight','xv_imp*','null_xvLL','xvLL*');
    end
else
    fprintf('Loading precomputed models\n');
    load(init_mods_fname);
end
%% INITIALIZE MODELS TO USE
thresh_xvimp = 1e-3;
for cc = 1:96
    if xv_imp_left(cc,2) > xv_imp_left(cc,1)
        init_mods_left(cc) = fitLeft(cc,2);
        init_xvimp_left(cc) = xv_imp_left(cc,2);
        use_quad_left(cc) = 1;
    else
        init_mods_left(cc) = fitLeft(cc,1);
        init_xvimp_left(cc) = xv_imp_left(cc,1);
        use_quad_left(cc) = 0;
    end
    if init_xvimp_left(cc) > thresh_xvimp
        useable_unit_left(cc) = 1;
    else
        useable_unit_left(cc) = 0;
    end
    
    if xv_imp_right(cc,2) > xv_imp_right(cc,1)
        init_mods_right(cc) = fitRight(cc,2);
        init_xvimp_right(cc) = xv_imp_right(cc,2);
        use_quad_right(cc) = 1;
    else
        init_mods_right(cc) = fitRight(cc,1);
        init_xvimp_right(cc) = xv_imp_right(cc,1);
        use_quad_right(cc) = 0;
    end
    if init_xvimp_right(cc) > thresh_xvimp
        useable_unit_right(cc) = 1;
    else
        useable_unit_right(cc) = 0;
    end
end

%% REDEFINE USED TIME POINTS AND RE-TUNE MODELS TO FULL TR SET
xv_frac = 0;
xv_set = [];
xv_inds = [];
tr_inds = find(~ismember(ic,xv_set))';

tr_inds(~all_used_inds(tr_inds)) = [];
NT = length(tr_inds);

use_XmatL = XmatL(tr_inds,:);
use_XmatR = XmatR(tr_inds,:);

silent = 1;
% if refit_mods == 1
for cc = 1:96
    fprintf('Unit %d of %d\n',cc,96);
    init_mods_left(cc) = NIMfit_filters(init_mods_left(cc),all_binned_spks(tr_inds,cc),use_XmatL,Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
    init_mod_LL_left(cc) = init_mods_left(cc).LL_seq(end);
    
    init_mods_right(cc) = NIMfit_filters(init_mods_right(cc),all_binned_spks(tr_inds,cc),use_XmatR,Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
    init_mod_LL_right(cc) = init_mods_right(cc).LL_seq(end);

    null_mod(cc) = NIMfit_filters(null_mod(cc),all_binned_spks(tr_inds,cc),Xexpt(tr_inds,:)); %fit stimulus filters
    null_LL(cc) = null_mod(cc).LL_seq(end);
end
save(init_mods_fname,'fit*','xv_imp*','null_xvLL','xvLL*','init_mod*','null_mod');
% else
%     for cc = 1:96
%         fprintf('Unit %d of %d\n',cc,96);
%         cur_mod_LL(cc) =  NIMmodel_eval(init_mod(cc),all_binned_spks(tr_inds,cc),Xmat(tr_inds,:),Xexpt(tr_inds,:));
%         cur_null_LL(cc) =  NIMmodel_eval(null_mod(cc),all_binned_spks(tr_inds,cc),Xexpt(tr_inds,:));
%         init_mod_LL(cc) = init_mod(cc).LL_seq(end);
%         null_LL(cc) = null_mod(cc).LL_seq(end);
%     end
%     if any(cur_mod_LL(:) ~= init_mod_LL(:))
%         fprintf('Warning: Mismatch with precomputed model LLs!\n');
%     end
%      if any(cur_null_LL(:) ~= null_LL(:))
%         fprintf('Warning: Mismatch with precomputed null LLs!\n');
%      end
% end
init_mod_LL_imp_left = (init_mod_LL_left - null_LL)/log(2);
init_mod_LL_imp_right = (init_mod_LL_right - null_LL)/log(2);

%% DEFINE POINTS TO ALLOW RAPID CHANGES (TRIAL STARTS AFTER SACCADES)
trial_start_inds = [1; find(diff(all_trialvec(tr_inds)) ~= 0) + 1];

used_saccade_inds = find(ismember(tr_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,tr_inds));

%push the effects of saccades forward in time
sac_shift = round(0.05/dt);

%define times when to resort to independent prior
use_prior = false(size(tr_inds));
use_prior(trial_start_inds) = true;
for i = 1:length(used_saccade_inds)
    next_trial = trial_start_inds(find(trial_start_inds > used_saccade_inds(i),1,'first'));
    %mark a point (forward in time) as a saccade if it occurs within the
    %same trial as the saccade
    if next_trial > used_saccade_inds(i) + sac_shift
        use_prior(used_saccade_inds(i) + sac_shift) = 1;
    end
end

%% INITIALIZE TRANSITION PRIORS FOR HMM
NT = length(tr_inds);
sp_dx = 0.05;
max_shift = 12;
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);

%overall prior on shifts
eps_prior_sigma = 0.125; %0.125 start
leps_prior = -(shifts*sp_dx).^2/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior,n_shifts,1);

cdist = squareform(pdist(shifts'*sp_dx));
deps_sigma = 0.015; %.015
lA = -cdist.^2/(2*deps_sigma^2);
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

%% SELECT WHICH EYE TO USE
use_Xmat = use_XmatL;
init_mod = init_mods_left;
useable_unit = useable_unit_left;
use_quad = use_quad_left;
init_xvimp = init_xvimp_left;
init_mod_LL_imp = init_mod_LL_imp_left;

% use_Xmat = use_XmatR;
% init_mod = init_mods_right;
% useable_unit = useable_unit_right;
% use_quad = use_quad_right;
% init_xvimp = init_xvimp_right;
% init_mod_LL_imp = init_mod_LL_imp_right;

%% ESTIMATE LL for each shift in each stimulus frame

tr_set = find(useable_unit);
n_tr_chs = length(tr_set);

Robs = all_binned_spks(tr_inds,tr_set);

%create filter banks for all units
lin_filt_bank = zeros(96,klen);
quad_filt_bank = zeros(96,klen);
lin_kerns = zeros(96,stim_params.lin_dims);
mod_offsets = nan(1,96);
for cc = 1:96
    cur_k = [init_mod(cc).mods(:).filtK];
    lin_filt_bank(cc,:) = cur_k(:,1);
    if use_quad(cc)
        quad_filt_bank(cc,:) = cur_k(:,2);
    end
    mod_offsets(cc) = init_mod(cc).spk_NL_params(1);
    lin_kerns(cc,:) = init_mod(cc).kLin';
end
lin_filt_bank = reshape(lin_filt_bank',[stim_params.stim_dims(1:2) 96]);
quad_filt_bank = reshape(quad_filt_bank',[stim_params.stim_dims(1:2) 96]);
lin_filt_bank = lin_filt_bank(:,:,tr_set);
quad_filt_bank = quad_filt_bank(:,:,tr_set);

%indicator predictions
expt_out = Xexpt(tr_inds,:)*lin_kerns(tr_set,:)';

%precompute LL at all shifts for all units
LLs = nan(NT,length(tr_set),n_shifts);
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',xx,n_shifts);
    cur_lin_shift = shift_matrix_Nd(lin_filt_bank,shifts(xx),2);
    cur_lin_shift = reshape(cur_lin_shift,klen,n_tr_chs);
    
    cur_quad_shift = shift_matrix_Nd(quad_filt_bank,shifts(xx),2);
    cur_quad_shift = reshape(cur_quad_shift,klen,n_tr_chs);
    
    %outputs of stimulus models at current X-matrix shift
    lin_out = use_Xmat*cur_lin_shift;
    quad_out = (use_Xmat*cur_quad_shift).^2;
    gfun = lin_out + quad_out;
    gfun = bsxfun(@plus,gfun,mod_offsets(tr_set));
    
    %add contributions from extra lin kernels
    gfun = gfun + expt_out;
    
    too_large = gfun > 50;
    pred_rate = log(1+exp(gfun));
    pred_rate(too_large) = gfun(too_large);
    
    pred_rate(pred_rate < 1e-20) = 1e-20;
    
    LLs(:,:,xx) = Robs.*log(pred_rate) - pred_rate;
end

%% DO INITIAL EYE INFERENCE AND MODEL RE-FITTING FOR XV SET

poss_xv_set = tr_set(init_xvimp(tr_set) > 10e-3);

resh_X = reshape(use_Xmat',[flen nPix NT]);
resh_X_sh = zeros(size(resh_X));
for XV = 1:length(poss_xv_set)
    xv_cell = poss_xv_set(XV);
    
    cur_tr_set = find(tr_set~=xv_cell);
    
    fprintf('Using XV cell %d,  %d of %d\n',xv_cell,XV,length(poss_xv_set));
    frame_LLs = squeeze(sum(LLs(:,cur_tr_set,:),2));
    
    %%
    lB = frame_LLs;
    lalpha=zeros(NT,n_shifts);
    lbeta = zeros(NT,n_shifts);
    lscale=zeros(NT,1); %initialize rescaling parameters
    %compute rescaled forward messages
    lalpha(1,:) = leps_prior + lB(1,:);
    lscale(1)=logsumexp(lalpha(1,:));
    lalpha(1,:) = lalpha(1,:) - lscale(1);
    for t=2:NT
        %         if mod(t,1000)==0
        %             fprintf('%d of %d\n',t,NT);
        %         end
        if use_prior(t)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
        lscale(t) = logsumexp(lalpha(t,:));
        lalpha(t,:)= lalpha(t,:) - lscale(t);
    end
    
    %compute rescaled backward messages
    lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
    for t=NT-1:-1:1
        %         if mod(t,1000)==0
        %             fprintf('%d of %d\n',t,NT);
        %         end
        if use_prior(t+1)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
    end
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    
    gamma = exp(lgamma);
    post_mean(:,XV) = sum(bsxfun(@times,gamma,shifts),2);
    post_std(:,XV) = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - post_mean(:,XV).^2);
    
    %% RECONSTRUCT MAP STIMULUS (OLD METHOD)
    [max_post,max_loc] = max(lgamma,[],2);
    
    shift_cor(:,XV) = shifts(max_loc);
    % shift_cor(:,XV) = post_mean;
    for ii = 1:NT
        d2 = shift_matrix_Nd(resh_X(:,:,ii), -shift_cor(ii,XV),2);
        resh_X_sh(:,:,ii) = d2;
    end
    X_sh = reshape(resh_X_sh,klen,NT)';
    
    %% REFIT MODELS
    %adjust regularization of new filter
    ref_mod{1}(xv_cell) = init_mod(xv_cell);
    ref_mod{1}(xv_cell) = NIMfit_filters(ref_mod{1}(xv_cell),all_binned_spks(tr_inds,xv_cell),X_sh,Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
    
    ref_LL_imp(1,xv_cell) = (ref_mod{1}(xv_cell).LL_seq(end) - null_LL(xv_cell))/log(2);
    
    fprintf('Original: %.4f  New: %.4f\n',init_mod_LL_imp(xv_cell),ref_LL_imp(1,xv_cell));
    
end

%% if there are any units in the training set that weren't allowed to be XV units, then refit those, using a all-cell eye position
if any(~ismember(tr_set,poss_xv_set))
    
    fprintf('Using XV cell %d,  %d of %d\n',xv_cell,XV,length(poss_xv_set));
    frame_LLs = squeeze(sum(LLs,2));
    
    %%
    lB = frame_LLs;
    lalpha=zeros(NT,n_shifts);
    lbeta = zeros(NT,n_shifts);
    lscale=zeros(NT,1); %initialize rescaling parameters
    %compute rescaled forward messages
    lalpha(1,:) = leps_prior + lB(1,:);
    lscale(1)=logsumexp(lalpha(1,:));
    lalpha(1,:) = lalpha(1,:) - lscale(1);
    for t=2:NT
        if use_prior(t)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
        lscale(t) = logsumexp(lalpha(t,:));
        lalpha(t,:)= lalpha(t,:) - lscale(t);
    end
    
    %compute rescaled backward messages
    lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
    for t=NT-1:-1:1
        if use_prior(t+1)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
    end
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    
    gamma = exp(lgamma);
    fpost_mean = sum(bsxfun(@times,gamma,shifts),2);
    fpost_std = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - post_mean(:,XV).^2);
    
    %% RECONSTRUCT MAP STIMULUS (OLD METHOD)
    [max_post,max_loc] = max(lgamma,[],2);
    
    fshift_cor = shifts(max_loc);
    for ii = 1:NT
        %         d2 = dist_shift2d(resh_X(:,:,ii), -shift_cor(ii,XV),2,0);
        d2 = shift_matrix_Nd(resh_X(:,:,ii), -shift_cor(ii,XV),2);
        resh_X_sh(:,:,ii) = d2;
    end
    X_sh = reshape(resh_X_sh,klen,NT)';
    
    %% REFIT MODELS
    cur_fit_set = setdiff(tr_set,poss_xv_set);
    for cc = 1:length(cur_fit_set)
        fprintf('Refitting unit %d of %d in tr-only set\n',cc,length(cur_fit_set));
        cur_unit = cur_fit_set(cc);
        %adjust regularization of new filter
        ref_mod{1}(cur_unit) = init_mod(cur_unit);
        ref_mod{1}(cur_unit) = NIMfit_filters(ref_mod{1}(cur_unit),all_binned_spks(tr_inds,cur_unit),X_sh,Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
        
        ref_LL_imp(1,cur_unit) = (ref_mod{1}(cur_unit).LL_seq(end) - null_LL(cur_unit))/log(2);
        
        fprintf('Original: %.4f  New: %.4f\n',init_mod_LL_imp(cur_unit),ref_LL_imp(1,cur_unit));
    end
else
    fshift_cor = mean(shift_cor);
    %jacknife error est
    fpost_std = (length(poss_xv_set)-1)*mean(var(post_std,[],2));
    fpost_mean = mean(post_mean,2);
end


%% PLOT IMPROVEMENT
mdl = LinearModel.fit(init_mod_LL_imp(tr_set)',ref_LL_imp(1,tr_set)');
xx = linspace(0,0.2,100);
[ypred,pred_errs] = predict(mdl,xx');
figure;hold on
plot(xx,ypred,'k','linewidth',2);plot(xx,pred_errs,'r--')
plot(init_mod_LL_imp(tr_set),ref_LL_imp(1,tr_set),'o')
line([0 0.5],[0 0.5])
% xlim([0 0.5]); ylim([0 0.5])

beta = mdl.Coefficients.Estimate;
fprintf('%.3f-fold improvement on xval\n',beta(2));

%%
sname = 'lefteye_f10_d12_mods';
save(sname,'shift_cor','post_mean','post_std','fpost_mean','fpost_std','fshift_cor','ref*','*left')
% sname = 'righteye_f10_d12_mods';
% save(sname,'shift_cor','post_mean','post_std','fpost_mean','fpost_std','fshift_cor','ref*','*right')

%% INITIALIZE FOR ITERATIVE EYE_CORRECTIONS
it_shift_cor{1} = shift_cor;
it_post_mean{1} = post_mean;
it_post_std{1} = post_std;
it_fpost_mean{1} = fpost_mean;
it_fpost_std{1} = fpost_std;
it_fshift_cor{1} = fshift_cor;

%% NOW ITERATE

n_ITER = 2;
for nn = 2:n_ITER+1
    
    fprintf('Eye correction iteration %d of %d\n',nn-1,n_ITER);
    
    %create filter banks for all units
    lin_filt_bank = zeros(length(tr_set),full_klen);
    quad_filt_bank = zeros(length(tr_set),full_klen);
    lin_kerns = zeros(length(tr_set),stim_params.lin_dims);
    mod_offsets = nan(1,length(tr_set));
    for cc = 1:length(tr_set)
        cur_k = [ref_mod{nn-1}(tr_set(cc)).mods(:).filtK];
        lin_filt_bank(cc,:) = cur_k(:,1);
        if use_quad(tr_set(cc))
            quad_filt_bank(cc,:) = cur_k(:,2);
        end
        mod_offsets(cc) = ref_mod{nn-1}(tr_set(cc)).spk_NL_params(1);
        lin_kerns(cc,:) = ref_mod{nn-1}(tr_set(cc)).kLin';
    end
    lin_filt_bank = reshape(lin_filt_bank',[stim_params.stim_dims(1:2) length(tr_set)]);
    quad_filt_bank = reshape(quad_filt_bank',[stim_params.stim_dims(1:2) length(tr_set)]);
    
    lin_Lfilts = lin_filt_bank(:,left_eye_inds,:);
    lin_Rfilts = lin_filt_bank(:,~left_eye_inds,:);
    quad_Lfilts = quad_filt_bank(:,left_eye_inds,:);
    quad_Rfilts = quad_filt_bank(:,~left_eye_inds,:);
    
    %indicator predictions
    expt_out = Xexpt(tr_inds,:)*lin_kerns';
    
    %precompute LL at all shifts for all units
    LLs = nan(NT,length(tr_set),n_shifts);
    for xx = 1:length(shifts)
        fprintf('Shift %d of %d\n',xx,n_shifts);
        cur_lin_shift_L = shift_matrix_Nd(lin_Lfilts,shifts(xx),2);
        cur_lin_shift_R = shift_matrix_Nd(lin_Rfilts,shifts(xx),2);
        cur_lin_shift = cat(2,cur_lin_shift_L,cur_lin_shift_R);
        cur_lin_shift = reshape(cur_lin_shift,full_klen,n_tr_chs);
        
        cur_quad_shift_L = shift_matrix_Nd(quad_Lfilts,shifts(xx),2);
        cur_quad_shift_R = shift_matrix_Nd(quad_Rfilts,shifts(xx),2);
        cur_quad_shift = cat(2,cur_quad_shift_L,cur_quad_shift_R);
        cur_quad_shift = reshape(cur_quad_shift,full_klen,n_tr_chs);
        
        %outputs of stimulus models at current X-matrix shift
        lin_out = use_Xmat*cur_lin_shift;
        quad_out = (use_Xmat*cur_quad_shift).^2;
        gfun = lin_out + quad_out;
        gfun = bsxfun(@plus,gfun,mod_offsets);
        
        %add contributions from extra lin kernels
        gfun = gfun + expt_out;
        
        too_large = gfun > 50;
        pred_rate = log(1+exp(gfun));
        pred_rate(too_large) = gfun(too_large);
        
        pred_rate(pred_rate < 1e-20) = 1e-20;
        
        LLs(:,:,xx) = Robs.*log(pred_rate) - pred_rate;
    end
    
    for XV = 1:length(tr_set)
        xv_cell = tr_set(XV);
        
        cur_tr_inds = find(tr_set~=xv_cell);
        
        fprintf('Using XV cell %d,  %d of %d\n',xv_cell,XV,length(tr_set));
        frame_LLs = squeeze(sum(LLs(:,cur_tr_inds,:),2));
        
        %%
        lB = frame_LLs;
        lalpha=zeros(NT,n_shifts);
        lbeta = zeros(NT,n_shifts);
        lscale=zeros(NT,1); %initialize rescaling parameters
        %compute rescaled forward messages
        lalpha(1,:) = leps_prior + lB(1,:);
        lscale(1)=logsumexp(lalpha(1,:));
        lalpha(1,:) = lalpha(1,:) - lscale(1);
        for t=2:NT
            if use_prior(t)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
            lscale(t) = logsumexp(lalpha(t,:));
            lalpha(t,:)= lalpha(t,:) - lscale(t);
        end
        
        %compute rescaled backward messages
        lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
        for t=NT-1:-1:1
            if use_prior(t+1)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lf1 = lbeta(t+1,:) + lB(t+1,:);
            lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
        end
        lgamma= lalpha + lbeta;
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
        
        gamma = exp(lgamma);
        post_mean(:,XV) = sum(bsxfun(@times,gamma,shifts),2);
        post_std(:,XV) = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - post_mean(:,XV).^2);
        
        %% RECONSTRUCT MAP STIMULUS (OLD METHOD)
        [max_post,max_loc] = max(lgamma,[],2);
        
        shift_cor(:,XV) = shifts(max_loc);
        for ii = 1:NT
            %         d2 = dist_shift2d(resh_X(:,:,ii), -shift_cor(ii,XV),2,0);
            d2_left = shift_matrix_Nd(resh_X(:,left_eye_inds,ii), -shift_cor(ii,XV),2);
            d2_right = shift_matrix_Nd(resh_X(:,~left_eye_inds,ii), -shift_cor(ii,XV),2);
            resh_X_sh(:,:,ii) = cat(2,d2_left,d2_right);
        end
        X_sh = reshape(resh_X_sh,full_klen,NT)';
        
        %% REFIT MODELS
        %adjust regularization of new filter
        ref_mod{nn}(xv_cell) = ref_mod{nn-1}(xv_cell);
        %     ref_mod{1}(xv_cell) = NIMadjust_regularization(ref_mod{1}(xv_cell),[1 2],'lambda_L1',0);
        ref_mod{nn}(xv_cell) = NIMfit_filters(ref_mod{nn}(xv_cell),all_binned_spks(tr_inds,xv_cell),X_sh,Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
        
        ref_LL_imp(nn,xv_cell) = (ref_mod{nn}(xv_cell).LL_seq(end) - null_LL(xv_cell))/log(2);
        
        fprintf('Original: %.4f  Previous: %.4f  New: %.4f\n',init_mod_LL_imp(xv_cell),ref_LL_imp(nn-1,xv_cell),ref_LL_imp(nn,xv_cell));
        
    end
    
    %% if there are any units in the training set that weren't allowed to be XV units, then refit those, using a all-cell eye position
    if any(~ismember(tr_set,poss_xv_set))
        
        fprintf('Using XV cell %d,  %d of %d\n',xv_cell,XV,length(poss_xv_set));
        frame_LLs = squeeze(sum(LLs,2));
        
        %%
        lB = frame_LLs;
        lalpha=zeros(NT,n_shifts);
        lbeta = zeros(NT,n_shifts);
        lscale=zeros(NT,1); %initialize rescaling parameters
        %compute rescaled forward messages
        lalpha(1,:) = leps_prior + lB(1,:);
        lscale(1)=logsumexp(lalpha(1,:));
        lalpha(1,:) = lalpha(1,:) - lscale(1);
        for t=2:NT
            if use_prior(t)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
            lscale(t) = logsumexp(lalpha(t,:));
            lalpha(t,:)= lalpha(t,:) - lscale(t);
        end
        
        %compute rescaled backward messages
        lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
        for t=NT-1:-1:1
            if use_prior(t+1)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lf1 = lbeta(t+1,:) + lB(t+1,:);
            lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
        end
        lgamma= lalpha + lbeta;
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
        
        gamma = exp(lgamma);
        fpost_mean = sum(bsxfun(@times,gamma,shifts),2);
        fpost_std = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - post_mean(:,XV).^2);
        
        %% RECONSTRUCT MAP STIMULUS (OLD METHOD)
        [max_post,max_loc] = max(lgamma,[],2);
        
        fshift_cor = shifts(max_loc);
        for ii = 1:NT
            %         d2 = dist_shift2d(resh_X(:,:,ii), -shift_cor(ii,XV),2,0);
            d2_left = shift_matrix_Nd(resh_X(:,left_eye_inds,ii), -shift_cor(ii,XV),2);
            d2_right = shift_matrix_Nd(resh_X(:,~left_eye_inds,ii), -shift_cor(ii,XV),2);
            resh_X_sh(:,:,ii) = cat(2,d2_left,d2_right);
        end
        X_sh = reshape(resh_X_sh,full_klen,NT)';
        
        %% REFIT MODELS
        cur_fit_set = setdiff(tr_set,poss_xv_set);
        for cc = 1:length(cur_fit_set)
            fprintf('Refitting unit %d of %d in tr-only set\n',cc,length(cur_fit_set));
            cur_unit = cur_fit_set(cc);
            %adjust regularization of new filter
            ref_mod{nn}(cur_unit) = ref_mod{nn-1}(cur_unit);
            ref_mod{nn}(cur_unit) = NIMfit_filters(ref_mod{nn}(cur_unit),all_binned_spks(tr_inds,cur_unit),X_sh,Xexpt(tr_inds,:),[],silent,optim_params,full_L2_mats); %fit stimulus filters
            ref_LL_imp(nn,cur_unit) = (ref_mod{nn}(cur_unit).LL_seq(end) - null_LL(cur_unit))/log(2);
            fprintf('Original: %.4f  Previous: %.4f  New: %.4f\n',init_mod_LL_imp(cur_unit),ref_LL_imp(nn-1,cur_unit),ref_LL_imp(nn,cur_unit));
        end
    else
        fshift_cor = mean(shift_cor);
        %jacknife error est
        fpost_std = (length(poss_xv_set)-1)*mean(var(post_std,[],2));
        fpost_mean = mean(post_mean,2);
    end
    
    %%
    it_post_mean{nn} = post_mean;
    it_post_std{nn} = post_std;
    it_shift_cor{nn} = shift_cor;
    it_fpost_mean{nn} = fpost_mean;
    it_fpost_std{nn} = fpost_std;
    it_fshift_cor{nn} = fshift_cor;
    
    %% PLOT PROGRESS
    mdl = LinearModel.fit(init_mod_LL_imp(tr_set)',ref_LL_imp(nn,tr_set)');
    xx = linspace(0,0.2,100);
    [ypred,pred_errs] = predict(mdl,xx');
    figure;hold on
    plot(xx,ypred,'k','linewidth',2);plot(xx,pred_errs,'r--')
    plot(init_mod_LL_imp(tr_set),ref_LL_imp(nn,tr_set),'o')
    plot(init_mod_LL_imp(tr_set),ref_LL_imp(nn-1,tr_set),'r*')
    line([0 0.5],[0 0.5])
    % xlim([0 0.5]); ylim([0 0.5])
    
    beta = mdl.Coefficients.Estimate;
    fprintf('%.3f-fold improvement on xval\n',beta(2));
    
    mdl = LinearModel.fit(ref_LL_imp(nn-1,tr_set)',ref_LL_imp(nn,tr_set)');
    beta = mdl.Coefficients.Estimate;
    fprintf('%.3f-fold previous improvement on xval\n',beta(2));
    
    %     xx = linspace(0,0.2,100);
    %     [ypred,pred_errs] = predict(mdl,xx');
    %     subplot(2,1,2)
    %     hold on
    %     plot(xx,ypred,'k','linewidth',2);plot(xx,pred_errs,'r--')
    %     plot(prev_imp,new_imp,'o')
    %     plot(prev_imp_nonxv,new_imp_nonxv,'k*')
    %     line([0 0.2],[0 0.2])
    %     xlim([0 0.25]); ylim([0 0.25])
    
    
    %     %% SAVE PROGRESS
    %     eye_times = all_stim_times(tr_inds);
    %     save full_eye_correct_v6_90deg *_LL it* tr_set non_xv_cells tr_inds eye_times
    
end


%%
measured_pre_Lx = [saccades(:).pre_Lx];
measured_post_Lx = [saccades(:).post_Lx];
measured_delta_pos = measured_post_Lx(used_saccade_set) - measured_pre_Lx(used_saccade_set);

for ii = 1:length(it_fpost_mean)
    inferred_pos = mean(it_fpost_mean{ii},2);
    inferred_pre_pos = inferred_pos(used_saccade_inds-1);
    inferred_post_pos = inferred_pos(used_saccade_inds + 6);
    inferred_delta_pos = (inferred_post_pos - inferred_pre_pos)*sp_dx;
    
    sac_delta_corr(ii) =  corr(measured_delta_pos(:),inferred_delta_pos(:),'type','spearman');
end