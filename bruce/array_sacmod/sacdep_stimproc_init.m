clear all

Expt_name = 'G086';
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

anal_name = 'monoc_eyecorr_hbar';

%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 12;
nPix = 30;
stim_params = NIMcreate_stim_params([flen nPix],dt);
Fr = 1;
bar_ori = 0;
beg_buffer = round(0.15/dt);
end_buffer = round(0.05/dt);

%%
exclude_expts = {'rls.orXme'};
for ii = 1:length(Expts)
    expt_names{ii} = Expts{ii}.Header.expname;
    expt_dds(ii) = Expts{ii}.Stimvals.dd;
    expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
    expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
    expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
    expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
    expt_ntrials(ii) = length(Expts{ii}.Trials);
    excluded_type(ii) = strcmp(expt_names{ii},exclude_expts);
end

cur_expt_set = find(expt_binoc(:) == 0 & expt_bar_ori(:) == bar_ori & expt_sac_dir(:) == 0 & ...
    expt_dds(:) == 12 & expt_Fr(:) == 1 & ~excluded_type(:));

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');

all_stim_times = [];
all_used_inds = false(0);
Xmat = [];
all_binned_spks = [];
all_exptvec = [];
all_trialvec = [];
all_tsince_start = [];
all_ttill_end = [];
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
    
    if length(un_ids) < length(trial_ids)
        fprintf('Warning: repeat trial ids\n');
    end
    
    buffer_pix = floor((expt_npix(cur_expt) - nPix)/2);
    cur_use_pix = (1:nPix) + buffer_pix;
    
    trial_durs = trial_durs(id_inds);
    trial_start_times = trial_start_times(id_inds);
    trial_end_times = trial_end_times(id_inds);
    
    
    use_trials = find(trial_durs >= 0.5);
    extra_trials = find(use_trials > length(left_stim_mats));
    if ~isempty(extra_trials)
        fprintf('Dropping %d trials with missing stim data\n',length(extra_trials));
        use_trials(extra_trials) = [];
    end
    
    all_trial_start_times = [all_trial_start_times trial_start_times(use_trials)];
    all_trial_end_times = [all_trial_end_times trial_end_times(use_trials)];
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_expt}.Trials(use_trials(tt)).Start/1e4;
        
        if length(cur_stim_times) == 1 %for blocks where the trial stimulus time stamps arent recorded
            n_frames = size(left_stim_mats{use_trials(tt)},1);
            cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
        end
        
        %eliminate any stim onset time stamps occuring after trial end
        cur_stim_times(cur_stim_times >= trial_end_times(use_trials(tt))) = [];
        
        cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
        
        cur_binned_spks = nan(length(cur_stim_times),96);
        for cc = 1:96
            [cur_hist,cur_bins] = histc(Clusters{cc}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
            cur_spk_used{cc}(cur_bins(1:end-1) > 0) = true;
        end
        
        nframes = size(left_stim_mats{use_trials(tt)},1);
        use_frames = min(length(cur_stim_times),nframes);
        cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
        bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
        
        cur_used_inds = true(length(cur_stim_times),1);
        cur_used_inds([1:beg_buffer (end-end_buffer+1):end]) = false;
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_tsince_start = [all_tsince_start; cur_stim_times - trial_start_times(tt)];
        all_ttill_end = [all_ttill_end; cur_stim_times(end)-cur_stim_times];
        all_used_inds = [all_used_inds; cur_used_inds];
        Xmat = [Xmat; bar_Xmat];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(size(cur_stim_times))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_stim_times))*(tt + trial_cnt)];
    end
    trial_cnt = trial_cnt + n_trials;
end
all_t_axis = all_stim_times;

%% normalize stimulus variance
Xmat = Xmat/std(reshape(Xmat(all_used_inds,:),1,[]));

%% SET UP USEABLE TRIALS
% [c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
[c,ia,ic] = unique(all_trialvec);
n_trials = length(ia);
% bad_trials = unique(ic(bad_pts)); %trials with putative blinks

xv_frac = 0.4;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(ic,xv_set));
tr_inds = find(~ismember(ic,xv_set));

tr_inds(~all_used_inds(tr_inds)) = [];
xv_inds(~all_used_inds(xv_inds)) = [];
NT = length(tr_inds);
trxv_inds = sort([xv_inds; tr_inds]);

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
Xexpt = zeros(length(all_stim_times),length(cur_expt_set)-1);
for i = 1:length(cur_expt_set)-1
    cur_set = find(all_exptvec==i);
    Xexpt(cur_set,i) = 1;
end

%% STA/STC test
fname = ['~/Analysis/bruce/' Expt_name '/hbar_stc'];

nneg = 3;
npos = 3;
ndims = size(Xmat,2);
unit_sta = nan(96,ndims);
unit_stcs = nan(96,ndims,npos+nneg);
for cc = 1:96
    spikebins = convert_to_spikebins(all_binned_spks(trxv_inds,cc));
    spike_cond_stim = Xmat(trxv_inds(spikebins),:);
    sta      = mean(spike_cond_stim) - mean(Xmat(trxv_inds,:));
    sta = sta/norm(sta);
    proj_mat = sta'/(sta*sta')*sta;
    stim_proj = Xmat(trxv_inds,:) - Xmat(trxv_inds,:)*proj_mat;
    % stim_proj = stim_emb;
    stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
    [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
    stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

    figure('visible','off');
    subplot(3,3,1)
    imagesc(reshape(sta,flen,nPix));
    ca = max(abs(sta)); caxis([-ca ca]);
    for ii = 1:3
        subplot(3,3,3+ii)
        imagesc(reshape(stcs(:,ii),flen,nPix));
        ca = max(abs(stcs(:,ii))); caxis([-ca ca]);
    end
    for ii = 1:3
        subplot(3,3,6+ii)
        imagesc(reshape(stcs(:,end-ii+1),flen,nPix));
        ca = max(abs(stcs(:,end-ii+1))); caxis([-ca ca]);
    end
    fillPage(gcf,'papersize',[10 7.5]);
    print(fname,'-dpsc','-append');
    close all
% pause
% close all
    unit_sta(cc,:) = sta;
    unit_stcs(cc,:,:) = stcs;
end


%% FIT OR LOAD INITIAL MODELS
refit_mods = 1;
save_dir = [anal_dir '/' anal_name];
if ~exist(save_dir,'dir')
    system(['mkdir ' save_dir]);
end
init_mods_fname = [save_dir '/init_mods_hbar' '.mat'];

klen = flen*nPix; 
stim_params = NIMcreate_stim_params([flen nPix],dt,1,1,length(cur_expt_set)-1);
optim_params.progTol = 1e-7;
optim_params.optTol = 1e-4;
silent = 1;
NL_types = {'lin', 'quad', 'quad','quad','quad'};

l_d2XT = 5000;
l_L1 = 100;
ql_L1 = 20;
ql_d2XT = 500;
XTmix = [1 1];

reg_params = NIMcreate_reg_params('lambda_d2XT',l_d2XT,'lambda_L1',l_L1);
null_params = NIMcreate_stim_params([1 size(Xexpt,2)],dt);

if refit_mods == 1 || ~exist(init_mods_fname,'file')
    fprintf('Computing initial models\n');
    for cc = [42]
        fprintf('Cell %d of %d\n',cc,96);
        Robs = all_binned_spks(tr_inds,cc);
        null_mod(cc) = NIMinitialize_model(null_params,1,{'lin'},reg_params); %initialize NIM
        null_mod(cc) = NIMfit_filters(null_mod(cc),Robs,Xexpt(tr_inds,:)); %fit stimulus filters
        
        if xv_frac > 0
            Robsxv = all_binned_spks(xv_inds,cc);
            null_xvLL(cc) = NIMmodel_eval(null_mod(cc),Robsxv,Xexpt(xv_inds,:));
        end
        
        NL_types = {'lin', 'quad', 'quad','quad','quad'};
        mod_signs = [1 1 1];
        fit0(cc) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
        %adjust regularization of new filter
        fit0(cc) = NIMadjust_regularization(fit0(cc),[2:length(mod_signs)],'lambda_d2XT',ql_d2XT);
        fit0(cc) = NIMadjust_regularization(fit0(cc),[2:length(mod_signs)],'lambda_L1',0);
        fit0(cc) = NIMfit_filters(fit0(cc),Robs,Xmat(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
        %add in L1 for new filter and re-fit
        fit0(cc) = NIMadjust_regularization(fit0(cc),[2:length(mod_signs)],'lambda_L1',ql_L1);
        fit0(cc) = NIMfit_filters(fit0(cc),Robs,Xmat(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
        if xv_frac > 0
            [xvLL(cc),~,xv_predrate] = NIMmodel_eval(fit0(cc),Robsxv,Xmat(xv_inds,:),Xexpt(xv_inds,:));
            xv_imp(cc) = (xvLL(cc)-null_xvLL(cc))/log(2);
            fprintf('LL imp: %.3f\n',xv_imp(cc));
        end
        
        fit1(cc) = NIMadd_NLinput(fit0(cc),'quad',-1);
        fit1(cc) = NIMadjust_regularization(fit1(cc),[length(mod_signs)],'lambda_d2XT',ql_d2XT);
        fit1(cc) = NIMadjust_regularization(fit1(cc),[length(mod_signs)],'lambda_L1',0);
        fit1(cc) = NIMfit_filters(fit1(cc),Robs,Xmat(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
        
        
        NL_types = {'threshlin', 'threshlin', 'quad','quad'};
        mod_signs = [1 -1 1 1];
        nfit0(cc) = NIMinitialize_model(stim_params,mod_signs,NL_types,reg_params); %initialize NIM
        %adjust regularization of new filter
        nfit0(cc) = NIMadjust_regularization(nfit0(cc),[1:length(mod_signs)],'lambda_d2XT',2000);
        nfit0(cc) = NIMadjust_regularization(nfit0(cc),[1:length(mod_signs)],'lambda_L1',50);
        nfit0(cc) = NIMadjust_regularization(nfit0(cc),[3:length(mod_signs)],'lambda_d2XT',ql_d2XT);
        nfit0(cc) = NIMadjust_regularization(nfit0(cc),[3:length(mod_signs)],'lambda_L1',0);
        nfit0(cc) = NIMfit_filters(nfit0(cc),Robs,Xmat(tr_inds,:),Xexpt(tr_inds,:),[],silent,optim_params); %fit stimulus filters
        if xv_frac > 0
            [xvLL_n(cc),~,xv_predrate] = NIMmodel_eval(nfit0(cc),Robsxv,Xmat(xv_inds,:),Xexpt(xv_inds,:));
            xv_imp_n(cc) = (xvLL_n(cc)-null_xvLL(cc))/log(2);
            fprintf('LL imp: %.3f\n',xv_imp_n(cc));
        end        
        
%         save(init_mods_fname,'fit0','xv_imp','null_xvLL','xvLL');
    end
else
%     fprintf('Loading precomputed models\n');
%     load(init_mods_fname);
end
%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)

emfile = ['jbe' Expt_name '.em.mat'];
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

%% IDENTIFY SACCADES FOR ANALYSIS
sac_start_times = [saccades(:).start_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > dt);
sac_start_times(bad_sacs) = [];
saccades(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];

max_dur = 0.1;
sac_durs = [saccades(:).duration];
is_blink = sac_durs > max_dur;

isis = [saccades(:).isi];
is_sacburst = (isis(1:end-1) < 0.15 | isis(2:end) < 0.15);

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps < 1;

sac_post_Lx = [saccades(:).post_Lx];
sac_post_Ly = [saccades(:).post_Ly];
sac_pre_Lx = [saccades(:).pre_Lx];
sac_pre_Ly = [saccades(:).pre_Ly];
sac_deltaX = sac_post_Lx - sac_pre_Lx;
sac_deltaY = sac_post_Ly - sac_pre_Ly;

sac_dirs = [0 90];
sac_expts = all_exptvec(interp_sac_start_inds);
delta_sacpar = nan(size(saccades));
for bb = 1:2    
    cur_bar_expts = find(expt_bar_ori(cur_expt_set) == sac_dirs(bb));
    cur_sac_set = find(ismember(sac_expts,cur_bar_expts));
    
    cur_delta_sac_par = abs(sac_deltaX(cur_sac_set)*cos(sac_dirs(bb)*pi/180) + sac_deltaY(cur_sac_set)*sin(sac_dirs(bb)*pi/180));    
    delta_sacpar(cur_sac_set) = cur_delta_sac_par;
end
is_gsac = delta_sacpar' > 1;

sim_sac_expts = find(~expt_has_ds(cur_expt_set));

sim_sac_times = [0.7 1.4 2.1 2.8 3.5];
all_sim_sacs = [];
sim_expt_inds = find(ismember(all_exptvec,sim_sac_expts));
sim_sacs = cell(length(sim_sac_times),1);
for ii = 1:length(sim_sac_times)
    sim_sacs{ii} = sim_expt_inds(all_tsince_start(sim_expt_inds(1:end-1)) < sim_sac_times(ii) & ...
        all_tsince_start(sim_expt_inds(2:end)) >= sim_sac_times(ii));
    all_sim_sacs = [all_sim_sacs; sim_sacs{ii}];
end

%define which saccades to use
used_msacs = find(is_micro & ~is_blink);
used_gsacs = find(is_gsac & ~is_blink);

gsac_start_inds = interp_sac_start_inds(used_gsacs);
msac_start_inds = interp_sac_start_inds(used_msacs);

%%
sac_backlag = round(0.1/dt);
sac_forlag = round(0.35/dt);
lags = -sac_backlag:sac_forlag;

trial_start_inds = 1 + find(diff(all_trialvec) ~= 0);
trial_stop_inds = trial_start_inds - 1;
trial_indices = [[1; trial_start_inds(:)] [trial_stop_inds(:); length(all_t_axis)]];

gsac_indicator = zeros(size(all_t_axis));
gsac_indicator(gsac_start_inds) = 1;

gsac_Xmat = zeros(length(all_t_axis),length(lags));
for ii = 1:size(trial_indices,1)
    cur_inds = trial_indices(ii,1):trial_indices(ii,2);
    cur_sig = gsac_indicator(cur_inds);
    
    fc = [cur_sig((sac_backlag+1):end); zeros(sac_backlag,1)];
    fr = [cur_sig((sac_backlag+1):-1:1); zeros(length(lags) - (sac_backlag+1),1)];
    cur_emb_mat = toeplitz(fc,fr);
    
    gsac_Xmat(cur_inds,:) = cur_emb_mat;
end

msac_indicator = zeros(size(all_t_axis));
msac_indicator(msac_start_inds) = 1;

msac_Xmat = zeros(length(all_t_axis),length(lags));
for ii = 1:size(trial_indices,1)
    cur_inds = trial_indices(ii,1):trial_indices(ii,2);
    cur_sig = msac_indicator(cur_inds);
    
    fc = [cur_sig((sac_backlag+1):end); zeros(sac_backlag,1)];
    fr = [cur_sig((sac_backlag+1):-1:1); zeros(length(lags) - (sac_backlag+1),1)];
    cur_emb_mat = toeplitz(fc,fr);
    
    msac_Xmat(cur_inds,:) = cur_emb_mat;
end

%%
% rel_gsac_inds = find(gsac_indicator(trxv_inds)==1);
% for cc = 1:96
% 
% Robs = all_binned_spks(trxv_inds,cc);
% 
% [ev_avg(cc,:),lags] = get_event_trig_avg(Robs,rel_gsac_inds,sac_backlag,sac_forlag,[],all_trialvec(trxv_inds));
% end


cc=42;
cur_mod = nfit0(cc);
nmods = length(cur_mod.mods);
Robs = all_binned_spks(tr_inds,cc);
[LL, penLL, pred_rate, G, gint, fgint] = NIMmodel_eval(cur_mod,Robs,Xmat(tr_inds,:),Xexpt(tr_inds,:));
pred_mat = [];
for ii = 1:nmods
    pred_mat = [pred_mat fgint(:,ii)];
end
[fitmod] = regGLM_fit(pred_mat,Robs,[],[],[],[],0);
mod_weights = fitmod.K;

[LL, penLL, pred_rate, G, gint, fgint] = NIMmodel_eval(cur_mod,all_binned_spks(:,cc),Xmat,Xexpt);
fgint = bsxfun(@times,fgint,mod_weights');

nmods = size(fgint,2);
pred_mat = [];
for ii = 1:nmods
    pred_mat = [pred_mat bsxfun(@times,msac_Xmat,nfit0(cc).mods(ii).sign*fgint(:,ii))];
end
for ii = 1:nmods 
    pred_mat = [pred_mat fgint(:,ii)];
end

L2_params = [];
for ii = 1:nmods
    cur_krange = (ii-1)*length(lags) + [1 length(lags)];
    L2_params = create_L2_params(L2_params,cur_krange,length(lags),2,1,0);
end

[fitmod] = regGLM_fit(pred_mat(tr_inds,:),Robs,L2_params,[ones(length(L2_params),1)*5],[],[],0);
kfit = fitmod.K(1:length(lags)*nmods);
kmat = reshape(kfit,length(lags),nmods);

Robsxv = all_binned_spks(xv_inds,cc);
[xvLL, penLL, pred_rate, G] = regGLM_eval(fitmod,Robsxv,pred_mat(xv_inds,:));
q

nmods = size(fgint,2);
pred_mat = [msac_Xmat];
for ii = 1:nmods
    pred_mat = [pred_mat fgint(:,ii)];
end

L2_params = [];
cur_krange = [1 length(lags)];
L2_params = create_L2_params(L2_params,cur_krange,length(lags),2,1,0);

[fitmod] = regGLM_fit(pred_mat(tr_inds,:),Robs,L2_params,[1],[],[],0);

[xvLL2, penLL, pred_rate, G] = regGLM_eval(fitmod,Robsxv,pred_mat(xv_inds,:));



net_out = sum(fgint,2);
pred_mat = [msac_Xmat bsxfun(@times,msac_Xmat,net_out)];
for ii = 1:nmods
    pred_mat = [pred_mat fgint(:,ii)];
end

L2_params = [];
cur_krange = [1 length(lags)];
L2_params = create_L2_params(L2_params,cur_krange,length(lags),2,1,0);
cur_krange = [length(lags)+1 2*length(lags)];
L2_params = create_L2_params(L2_params,cur_krange,length(lags),2,1,0);

[fitmod] = regGLM_fit(pred_mat(tr_inds,:),Robs,L2_params,[10 50],[],[],0);

[xvLL3, penLL, pred_rate, G] = regGLM_eval(fitmod,Robsxv,pred_mat(xv_inds,:));



%%
cur_stim_params = NIMcreate_stim_params([length(lags)],1,1,1,nmods);
cur_reg_params = NIMcreate_reg_params('lambda_d2T',200,'temporal_boundaries','zero','lambda_L1',5);
temp_nim = NIMinitialize_model(cur_stim_params,[1 -1],{'threshlin','threshlin'},cur_reg_params);
cur_nimfit = NIMfit_filters(temp_nim,Robs,gsac_Xmat(tr_inds,:),fgint(tr_inds,:),[],0);
 
