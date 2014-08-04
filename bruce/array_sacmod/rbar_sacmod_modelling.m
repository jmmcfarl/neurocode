clear all
% close all

Expt_name = 'G088';
data_dir = ['~/Data/bruce/' Expt_name];

cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 12;
nPix = 24;
stim_params = NIMcreate_stim_params([flen nPix],dt);
Fr = 1;

%%
include_expts = {'rls.Fa', 'rls.FaXimi'};
expt_varsac = false(length(Expts),1);
for ii = 1:length(Expts)
    expt_names{ii} = Expts{ii}.Header.expname;
    expt_dds(ii) = Expts{ii}.Stimvals.dd;
    expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
    if isfield(Expts{ii}.Trials(1),'Fa')
        if length(Expts{ii}.Trials(1).Fa) > 1
            expt_sac_dir(ii) = mod(Expts{ii}.Trials(1).Fa(1),180);
            expt_varsac(ii) = true;
        else
            expt_sac_dir(ii) = mod(Expts{ii}.Trials(1).Fa,180);
        end
    else
        expt_sac_dir(ii) = nan;
    end
    expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
    expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
    included_type(ii) = any(strcmp(expt_names{ii},include_expts));
    expt_ntrials(ii) = length(Expts{ii}.Trials);
end
expt_has_ds(isnan(expt_has_ds)) = 0;
expt_has_ds(expt_has_ds == -1) = 0;
expt_binoc(isnan(expt_binoc)) = 0;
expt_orthsac = expt_has_ds' & (expt_sac_dir ~= expt_bar_ori);

% cur_expt_set = find(included_type & expt_has_ds' == 1);
cur_expt_set = find(included_type & ~expt_binoc' & expt_npix' == 36 & expt_Fr == 1 & ~expt_orthsac & ~expt_varsac');

if strcmp(Expt_name,'G087')
    cur_expt_set(cur_expt_set == 15) = []; %only 6 trials and causes problems
end
%%
min_tdur = 0.5;

all_stim_times = [];
all_Xmat = [];
all_binned_spks = [];
all_exptvec = [];
all_trialvec = [];
all_tsince_start = [];
all_ttill_end = [];
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
    
    use_trials = find(trial_durs' >= min_tdur & ~isnan(expt_rc_matches{cur_expt}));
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
        
        nframes = size(left_stim_mats{use_trials(tt)},1);
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_tdur/dt
            use_frames = min(length(cur_stim_times),nframes);
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
            
            cur_binned_spks = nan(length(cur_stim_times),96);
            for cc = 1:96
                cur_hist = histc(Clusters{cc}.times,cur_t_edges);
                cur_binned_spks(:,cc) = cur_hist(1:end-1);
                %             temp = convert_to_spikebins(cur_hist(1:end-1));
            end
            
            all_tsince_start = [all_tsince_start; cur_stim_times - trial_start_times(use_trials(tt))];
            all_ttill_end = [all_ttill_end; trial_end_times(use_trials(tt)) - cur_stim_times];
            all_stim_times = [all_stim_times; cur_stim_times];
            all_Xmat = [all_Xmat; bar_Xmat];
            all_binned_spks = [all_binned_spks; cur_binned_spks];
            all_exptvec = [all_exptvec; ones(size(cur_stim_times))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_stim_times))*tt];
        end
    end
end

all_t_axis = all_stim_times;
%%
tr_inds = find(all_tsince_start > 0.4);
%%
to_print = 1;
stc_figname = ['~/Analysis/bruce/' Expt_name '/all_gs_stc'];

hori_expts = find(expt_bar_ori(cur_expt_set) == 0);
ver_expts = find(expt_bar_ori(cur_expt_set) == 90);

first_fig = 1;

nneg = 5;
npos = 5;
ndims = size(all_Xmat,2);
unit_sta = nan(96,ndims);
unit_stcs = nan(96,ndims,npos+nneg);
unit_evecs = nan(96,ndims);
for cc = 1:96
    
    %for horizontal bar ori
    cur_tr_inds = tr_inds(ismember(all_exptvec(tr_inds),hori_expts));
    
    spikebins = convert_to_spikebins(all_binned_spks(cur_tr_inds,cc));
    spike_cond_stim = all_Xmat(cur_tr_inds(spikebins),:);
    sta      = mean(spike_cond_stim) - mean(all_Xmat(cur_tr_inds,:));
    sta = sta/norm(sta);
    proj_mat = sta'/(sta*sta')*sta;
    stim_proj = all_Xmat(cur_tr_inds,:) - all_Xmat(cur_tr_inds,:)*proj_mat;
    % stim_proj = stim_emb;
    stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
    [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
    stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
    
    if to_print == 1
        figure('visible','off');
    else
        figure
    end
    subplot(3,npos,1)
    imagesc(reshape(sta,flen,nPix));
    ca = max(abs(sta)); caxis([-ca ca]);
    title(sprintf('Unit %d Hori',cc));
    subplot(3,npos,2)
    plot(diag(evals),'.')
    subplot(3,npos,3)
    plot(diag(evals),'.'); xlim([ndims-50 ndims+1])
    subplot(3,npos,4)
    plot(diag(evals),'.'); xlim([0 50])
    for ii = 1:npos
        subplot(3,npos,npos+ii)
        imagesc(reshape(stcs(:,ii),flen,nPix));
        ca = max(abs(stcs(:,ii))); caxis([-ca ca]);
    end
    for ii = 1:nneg
        subplot(3,npos,2*npos+ii)
        imagesc(reshape(stcs(:,end-nneg+ii),flen,nPix));
        ca = max(abs(stcs(:,end-nneg+ii))); caxis([-ca ca]);
    end
    if to_print == 1
        fillPage(gcf,'papersize',[10 7.5]);
        if first_fig == 1
            print(stc_figname,'-dpsc');
        else
            print(stc_figname,'-dpsc','-append');
        end
        close all
    else
        pause
        close all
    end
    unit_sta_hor(cc,:) = sta;
    unit_stcs_hor(cc,:,:) = stcs;
    unit_evecs_hor(cc,:) = diag(evals);
    
    
    %for vertical bar ori
    cur_tr_inds = tr_inds(ismember(all_exptvec(tr_inds),ver_expts));
    
    spikebins = convert_to_spikebins(all_binned_spks(cur_tr_inds,cc));
    spike_cond_stim = all_Xmat(cur_tr_inds(spikebins),:);
    sta      = mean(spike_cond_stim) - mean(all_Xmat(cur_tr_inds,:));
    sta = sta/norm(sta);
    proj_mat = sta'/(sta*sta')*sta;
    stim_proj = all_Xmat(cur_tr_inds,:) - all_Xmat(cur_tr_inds,:)*proj_mat;
    % stim_proj = stim_emb;
    stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
    [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
    stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
    
    if to_print == 1
        figure('visible','off');
    else
        figure
    end
    subplot(3,npos,1)
    imagesc(reshape(sta,flen,nPix));
    ca = max(abs(sta)); caxis([-ca ca]);
    title(sprintf('Unit %d Ver',cc));
    subplot(3,npos,2)
    plot(diag(evals),'.')
    subplot(3,npos,3)
    plot(diag(evals),'.'); xlim([ndims-50 ndims+1])
    subplot(3,npos,4)
    plot(diag(evals),'.'); xlim([0 50])
    for ii = 1:npos
        subplot(3,npos,npos+ii)
        imagesc(reshape(stcs(:,ii),flen,nPix));
        ca = max(abs(stcs(:,ii))); caxis([-ca ca]);
    end
    for ii = 1:nneg
        subplot(3,npos,2*npos+ii)
        imagesc(reshape(stcs(:,end-nneg+ii),flen,nPix));
        ca = max(abs(stcs(:,end-nneg+ii))); caxis([-ca ca]);
    end
    if to_print == 1
        fillPage(gcf,'papersize',[10 7.5]);
        %         if first_fig == 1
        %         print(stc_figname,'-dpsc');
        %         else
        print(stc_figname,'-dpsc','-append');
        %         end
        close all
    else
        pause
        close all
    end
    unit_sta_ver(cc,:) = sta;
    unit_stcs_ver(cc,:,:) = stcs;
    unit_evecs_ver(cc,:) = diag(evals);
    
    first_fig = 0;
end

sname = ['~/Analysis/bruce/' Expt_name '/all_gs_stc_data'];
save(sname,'unit_sta*','unit_stcs*','unit_evecs*');

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
cd(data_dir);

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
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%%
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

%% PICK OUT SACCADES FOR ANALYSIS
sac_start_times = [saccades(:).start_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > 1/dt | isnan(all_trialvec(interp_sac_start_inds))');
sac_start_times(bad_sacs) = [];
saccades(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];

sac_stop_times = [saccades(:).stop_time];
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));
interp_sac_stop_inds(isnan(interp_sac_stop_inds)) = 1;

sac_peak_times = [saccades(:).peak_time];
interp_sac_peak_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_peak_times));
interp_sac_peak_inds(isnan(interp_sac_peak_inds)) = 1;

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

if strcmp(Expt_name,'G081')
    sac_dirs = [0 45 90 135];
else
    sac_dirs = [0 90];
end
sac_start_expts = all_exptvec(interp_sac_start_inds);
delta_sacpar = nan(size(saccades));
for bb = 1:length(sac_dirs)
    cur_bar_expts = find(expt_bar_ori(cur_expt_set) == sac_dirs(bb));
    cur_sac_set = find(ismember(sac_start_expts,cur_bar_expts));
    
    cur_delta_sac_par = abs(sac_deltaX(cur_sac_set)*cos(sac_dirs(bb)*pi/180) + sac_deltaY(cur_sac_set)*sin(sac_dirs(bb)*pi/180));
    delta_sacpar(cur_sac_set) = cur_delta_sac_par;
end
is_gsac = delta_sacpar' > 1;

sim_sac_expts = find(~expt_has_ds(cur_expt_set));

if strcmp(Expt_name,'G081')
    sim_sac_times = [0.7 1.4];
else
    sim_sac_times = [0.7 1.4 2.1 2.8 3.5];
end
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

% all_sim_sacs(all_tsince_start(all_sim_sacs) <= 1) = [];
% used_msacs(all_tsince_start(interp_sac_start_inds(used_msacs)) <= 1) = [];
% used_gsacs(all_tsince_start(interp_sac_start_inds(used_gsacs)) <= 1) = [];


imback_gs_expts = find(expt_has_ds(cur_expt_set) & expt_imback(cur_expt_set)');
grayback_gs_expts = find(expt_has_ds(cur_expt_set) & ~expt_imback(cur_expt_set)' & ~expt_binoc(cur_expt_set));
hori_expts = find(expt_bar_ori(cur_expt_set) == 0);
ver_expts = find(expt_bar_ori(cur_expt_set) == 90);

trial_start_inds = find(all_tsince_start < 1/dt);
simsac_trial_starts = trial_start_inds(ismember(all_exptvec(trial_start_inds),sim_sac_expts));
grayback_trial_starts = trial_start_inds(ismember(all_exptvec(trial_start_inds),grayback_gs_expts));
imback_trial_starts = trial_start_inds(ismember(all_exptvec(trial_start_inds),imback_gs_expts));

used_msacs_grayback = used_msacs(ismember(all_exptvec(interp_sac_start_inds(used_msacs)),grayback_gs_expts));
used_msacs_imback = used_msacs(ismember(all_exptvec(interp_sac_start_inds(used_msacs)),imback_gs_expts));
used_gsacs_grayback = used_gsacs(ismember(all_exptvec(interp_sac_start_inds(used_gsacs)),grayback_gs_expts));
used_gsacs_imback = used_gsacs(ismember(all_exptvec(interp_sac_start_inds(used_gsacs)),imback_gs_expts));

used_msacs_hori = used_msacs(ismember(all_exptvec(interp_sac_start_inds(used_msacs)),hori_expts));
used_msacs_ver = used_msacs(ismember(all_exptvec(interp_sac_start_inds(used_msacs)),ver_expts));
used_gsacs_hori = used_gsacs(ismember(all_exptvec(interp_sac_start_inds(used_gsacs)),hori_expts));
used_gsacs_ver = used_gsacs(ismember(all_exptvec(interp_sac_start_inds(used_gsacs)),ver_expts));

simsacs_hori = all_sim_sacs(ismember(all_exptvec(all_sim_sacs),hori_expts));
simsacs_ver = all_sim_sacs(ismember(all_exptvec(all_sim_sacs),ver_expts));


%%
gsac_backlag = 0.2;
gsac_forlag = 0.5;
gsac_bin_width = 2;
gsac_binspace = gsac_bin_width*dt;
gsac_bin_edges = -(gsac_backlag-gsac_binspace/2):gsac_binspace:(gsac_forlag+gsac_binspace/2);
gsac_bin_cents = 0.5*gsac_bin_edges(1:end-1) + 0.5*gsac_bin_edges(2:end);
n_gsac_bins = length(gsac_bin_cents);

gsac_inds = interp_sac_start_inds(used_gsacs);
msac_inds = interp_sac_start_inds(used_msacs);

trial_sac_mat = zeros(length(all_t_axis),n_gsac_bins);
for i = 1:n_gsac_bins
    for cc = 1:gsac_bin_width
%                 cur_inds = gsac_inds + round(gsac_bin_cents(i)/dt) - floor(gsac_bin_width/2) + cc;
        cur_inds = msac_inds + round(gsac_bin_cents(i)/dt) - floor(gsac_bin_width/2) + cc;
        %                 cur_inds = all_sim_sacs + round(gsac_bin_cents(i)/dt) - floor(gsac_bin_width/2) + cc;
        cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
        trial_sac_mat(cur_inds,i) = 1;
    end
end

stc_thresh = -1.5e-3;
nbins = 25;
bin_edges = linspace(0,100,nbins+1);
bin_cents = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);

smooth_fac = 1;

cur_figname = ['~/Analysis/bruce/' Expt_name '/all_gs_modNL'];

%%
% sus = [12 42 43 48 49 91 94]; %for G086
% to_print = 0;
% % for cc = 1:96
% for cc = sus
% % cc = 94;
%
%
%     Robs = all_binned_spks(tr_inds,cc);
%
%     cur_evec_diff = diff(fliplr(unit_evecs(cc,:)));
%     npos_stc(cc) = find(cur_evec_diff > stc_thresh,1,'first')-1;
%     nneg_stc(cc) = find(cur_evec_diff(end:-1:1) > stc_thresh,1,'first')-1;
%     g_sta = all_Xmat(tr_inds,:)*unit_sta(cc,:)';
%     g_posstc = zeros(length(tr_inds),npos_stc(cc));
%     g_negstc = zeros(length(tr_inds),nneg_stc(cc));
%     for gg = 1:npos_stc(cc)
%         g_posstc(:,gg) = all_Xmat(tr_inds,:)*unit_stcs(cc,:,gg)';
%     end
%     for gg = 1:nneg_stc(cc)
%         g_negstc(:,gg) =  all_Xmat(tr_inds,:)*unit_stcs(cc,:,end-nneg+gg)';
%     end
%
%     cur_Xmat = [g_sta g_posstc.^2 g_negstc.^2];
%     [fitmod] = regGLM_fit(cur_Xmat,Robs,[],[],[],[],0);
%     g_full = cur_Xmat*fitmod.K;
%
%     pp_sta = prctile(g_sta,bin_edges);
%     pp_full = prctile(g_full,bin_edges);
%     gsta_bin_cents(cc,:) = 0.5*pp_sta(1:end-1) + 0.5*pp_sta(2:end);
%     gfull_bin_cents(cc,:) = 0.5*pp_full(1:end-1) + 0.5*pp_full(2:end);
%     stas_g(cc,:,:) = zeros(n_gsac_bins,nbins);
%     full_g(cc,:,:) = zeros(n_gsac_bins,nbins);
%     for tt = 1:n_gsac_bins
%         cur_inds = find(trial_sac_mat(tr_inds,tt) == 1);
%         for ii = 1:nbins
%             cur_set = cur_inds(g_sta(cur_inds) >= pp_sta(ii) & g_sta(cur_inds) < pp_sta(ii+1));
%             stas_g(cc,tt,ii) = mean(Robs(cur_set));
%             cur_set = cur_inds(g_full(cur_inds) >= pp_full(ii) & g_full(cur_inds) < pp_full(ii+1));
%             full_g(cc,tt,ii) = mean(Robs(cur_set));
%         end
%         stas_g(cc,tt,:) = jmm_smooth_1d_cor(squeeze(stas_g(cc,tt,:)),smooth_fac);
%         full_g(cc,tt,:) = jmm_smooth_1d_cor(squeeze(full_g(cc,tt,:)),smooth_fac);
%     end
%     ov_fullg(cc,:) = zeros(1,nbins);
%     ov_stag(cc,:) = zeros(1,nbins);
%     for ii = 1:nbins
%         cur_set = find(g_full >= pp_full(ii) & g_full < pp_full(ii+1));
%         ov_fullg(cc,ii) = mean(Robs(cur_set));
%         cur_set = find(g_sta >= pp_sta(ii) & g_sta < pp_sta(ii+1));
%         ov_stag(cc,ii) = mean(Robs(cur_set));
%     end
%     ov_fullg(cc,:) = jmm_smooth_1d_cor(ov_fullg(cc,:),smooth_fac);
%     ov_stag(cc,:) = jmm_smooth_1d_cor(ov_stag(cc,:),smooth_fac);
%
%     [full_g_dist(cc,:),full_g_x(cc,:)] = ksdensity(g_full);
%     [sta_g_dist(cc,:),sta_g_x(cc,:)] = ksdensity(g_sta);
%     avg_sacta = squeeze(mean(full_g(cc,:,:),3));
%     [~,max_point] = max(avg_sacta);
%     [~,min_point] = min(avg_sacta);
%
%     point_range = min_point:max_point;
%     cmap = jet(length(point_range));
%
%     if to_print == 1
%         figure('visible','off');
%     else
%         figure
%     end
%     subplot(1,3,1)
%     plot(gsac_bin_cents*dt,avg_sacta/dt,'ko-')
%     hold on
%     plot(gsac_bin_cents(max_point)*dt,avg_sacta(max_point)/dt,'r*','markersize',16)
%     plot(gsac_bin_cents(min_point)*dt,avg_sacta(min_point)/dt,'b*','markersize',16);
%     title(sprintf('Unit %d using %d Pos %d Neg',cc,npos_stc(cc),nneg_stc(cc)));
%     subplot(1,3,2)
%     % plot(gfull_bin_cents(cc,:),squeeze(full_g(cc,max_point,:))/dt,'ro-')
%     hold on
%     % plot(gfull_bin_cents(cc,:),squeeze(full_g(cc,min_point,:))/dt,'bo-')
%     plot(gfull_bin_cents(cc,:),ov_fullg(cc,:)/dt,'ko-')
%     for ii = 1:length(point_range)
%         plot(gfull_bin_cents(cc,:),squeeze(full_g(cc,point_range(ii),:))/dt,'color',cmap(ii,:));
%     end
%     xlim(gfull_bin_cents(cc,[1 end]))
%     yl = ylim();
%     subplot(1,3,3)
%     % plot(gfull_bin_cents(cc,:),squeeze(full_g(cc,max_point,:))/dt,'ro-')
%     hold on
%     % plot(gfull_bin_cents(cc,:),squeeze(full_g(cc,min_point,:))/dt,'bo-')
%     plot(gsta_bin_cents(cc,:),ov_stag(cc,:)/dt,'ko-')
%     for ii = 1:length(point_range)
%         plot(gsta_bin_cents(cc,:),squeeze(stas_g(cc,point_range(ii),:))/dt,'color',cmap(ii,:));
%     end
%     xlim(gsta_bin_cents(cc,[1 end]))
%     ylim(yl);
%     if to_print == 1
%         fillPage(gcf,'papersize',[12 5]);
%         if cc == 1
%             print(cur_figname,'-dpsc');
%         else
%             print(cur_figname,'-dpsc','-append');
%         end
%         close all
%     else
% %         pause
% %         close all
%     end
%
%
%     if to_print == 1
%         figure('visible','off');
%     else
%         figure
%     end
%     subplot(3,3,1)
%     imagesc(reshape(unit_sta(cc,:),flen,nPix));
%     ca = max(abs(unit_sta(cc,:))); caxis([-ca ca]);
%     title(sprintf('Unit %d',cc));
%     for ii = 1:min(npos_stc(cc),3)
%         subplot(3,3,3+ii)
%         imagesc(reshape(unit_stcs(cc,:,ii),flen,nPix));
%         ca = max(abs(unit_stcs(cc,:,ii))); caxis([-ca ca]);
%     end
%     for ii = 1:min(nneg_stc(cc),3)
%         subplot(3,3,6+ii)
%         imagesc(reshape(unit_stcs(cc,:,end-nneg+ii),flen,nPix));
%         ca = max(abs(unit_stcs(cc,:,end-nneg+ii))); caxis([-ca ca]);
%     end
%     if to_print == 1
%         fillPage(gcf,'papersize',[10 7.5]);
%         print(cur_figname,'-dpsc','-append');
%         close all
%     else
%         pause
%         close all
%     end
%
%
% end

%% USING MODELS
fullmod_figname = ['~/Analysis/bruce/' Expt_name '/all_fullqmod_msacmodNL'];

to_print = 1;
refit_models = 0;

sus = [12 42 43 48 49 91 94]; %for G086
use_data_range = 1:(nbins);
l_d2XT = 500;
l_L1 = 50;
ql_L1 = 1;
ql_d2XT = 50;

stim_params = NIMcreate_stim_params([flen nPix],dt);
reg_params = NIMcreate_reg_params('lambda_d2XT',l_d2XT,'lambda_L1',l_L1);
init_reg_params = NIMcreate_reg_params('lambda_d2XT',l_d2XT);
silent = 0;
init_optim_params.progTol = 1e-7;
init_optim_params.optTol = 1e-4;
init_optim_params.maxIter = 100;

optim_params.progTol = 1e-7;
optim_params.optTol = 1e-4;
NT = length(tr_inds);
sub_samp_fac = 10;

first_fig = 1;
for cc = 1:96
    % for cc = sus
    
    
    %% FOR HORI BARS
    cur_tr_inds = tr_inds(ismember(all_exptvec(tr_inds),hori_expts));
    NT = length(cur_tr_inds);
    Robs = all_binned_spks(cur_tr_inds,cc);
    
    cur_evec_diff = diff(fliplr(unit_evecs_hor(cc,:)));
    npos_stc(cc) = find(cur_evec_diff > stc_thresh,1,'first')-1;
    nneg_stc(cc) = find(cur_evec_diff(end:-1:1) > stc_thresh,1,'first')-1;
    
    if ~exist('fit1_hor','var') || cc > length(fit1_hor) || isempty(fit1_hor(cc).mods) || refit_models == 1
        fprintf('Fitting unit %d, with %d exc and %d sup squared\n',cc,npos_stc(cc),nneg_stc(cc));
        mod_signs = [1 ones(1,npos_stc(cc)) -1*ones(1,nneg_stc(cc))];
        NL_types = cat(2,'lin', repmat({'quad'},1,npos_stc(cc) + nneg_stc(cc)));
        %     NL_types = repmat({'threshlin'},1,1+npos_stc(cc) + nneg_stc(cc));
        
        fit0_hor(cc) = NIMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params); %initialize NIM
        %adjust regularization of new filter
        fit0_hor(cc) = NIMadjust_regularization(fit0_hor(cc),[2:length(mod_signs)],'lambda_d2XT',ql_d2XT);
        
        sub_sample = randperm(NT);
        sub_sample = sub_sample(1:round(NT/sub_samp_fac));
        
        fit0_hor(cc) = NIMfit_filters(fit0_hor(cc),Robs(sub_sample),all_Xmat(cur_tr_inds(sub_sample),:),[],[],silent,init_optim_params); %fit stimulus filters
        fit0_hor(cc) = NIMfit_filters(fit0_hor(cc),Robs,all_Xmat(cur_tr_inds,:),[],[],silent,init_optim_params); %fit stimulus filters
        
        fit1_hor(cc) = NIMadjust_regularization(fit0_hor(cc),1,'lambda_L1',l_L1);
        fit1_hor(cc) = NIMadjust_regularization(fit1_hor(cc),[2:length(mod_signs)],'lambda_L1',ql_L1);
        fit1_hor(cc) = NIMfit_filters(fit1_hor(cc),Robs,all_Xmat(cur_tr_inds,:),[],[],silent,optim_params); %fit stimulus filters
    end
    [LL, penLL, pred_rate, G, gint, fgint] = NIMmodel_eval(fit1_hor(cc),Robs,all_Xmat(cur_tr_inds,:));
    G = G - fit1_hor(cc).spk_NL_params(1);
    fit2_hor = NIMfit_logexp_spkNL_withoffset(fit1_hor(cc),Robs,all_Xmat(cur_tr_inds,:),[],0,[]);
    ov_spk_nl_params(cc,:) = fit2_hor.spk_NL_params;
    
    pp_nim = prctile(G,bin_edges);
    gnim_bin_cents(cc,:) = 0.5*pp_nim(1:end-1) + 0.5*pp_nim(2:end);
    gnim_bin_cents(cc,end) = 2*gnim_bin_cents(cc,end-1) - gnim_bin_cents(cc,end-2);
    gnim_bin_cents(cc,1) = 2*gnim_bin_cents(cc,2) - gnim_bin_cents(cc,3);
    ov_nimg(cc,:) = zeros(1,nbins);
    for ii = 1:nbins
        cur_set = find(G >= pp_nim(ii) & G < pp_nim(ii+1));
        ov_nimg(cc,ii) = mean(Robs(cur_set));
    end
    ov_nimg(cc,:) = jmm_smooth_1d_cor(ov_nimg(cc,:),smooth_fac);
    [nim_g_dist(cc,:),nim_g_x(cc,:)] = ksdensity(G);
    %     ov_pfit(cc,:) = polyfit(gnim_bin_cents(cc,use_data_range),ov_nimg(cc,use_data_range),1);
    ov_pred_resp(cc,:) = logexp_fun(gnim_bin_cents(cc,:),ov_spk_nl_params(cc,:));
    
    %     ov_proj_zero = -ov_pfit(cc,2)/ov_pfit(cc,1);
    
    avg_sacta = zeros(1,n_gsac_bins);
    nim_g(cc,:,:) = zeros(n_gsac_bins,nbins);
    for tt = 1:n_gsac_bins
        cur_inds = find(trial_sac_mat(cur_tr_inds,tt) == 1);
        avg_sacta(tt) = mean(Robs(cur_inds));
        for ii = 1:nbins
            cur_set = cur_inds(G(cur_inds) >= pp_nim(ii) & G(cur_inds) < pp_nim(ii+1));
            nim_g(cc,tt,ii) = mean(Robs(cur_set));
        end
        %         pfit(cc,tt,:) = polyfit(gnim_bin_cents(cc,use_data_range),squeeze(nim_g(cc,tt,use_data_range))',1);
        %         cur_offset(tt) = polyval(squeeze(pfit(cc,tt,:)),ov_proj_zero);
        cur_fit = NIMfit_logexp_spkNL_withoffset(fit2_hor,Robs(cur_inds),all_Xmat(cur_tr_inds(cur_inds),:),[],0,[1 2]);
        cur_spk_nl_params(tt,:) = cur_fit.spk_NL_params;
        pred_resp(tt,:) = logexp_fun(gnim_bin_cents(cc,:),cur_spk_nl_params(tt,:));
        nim_g(cc,tt,:) = jmm_smooth_1d_cor(squeeze(nim_g(cc,tt,:)),smooth_fac);
    end
    
        poss_spks = 0:10;
    n_symbs = length(poss_spks);

    avg_rate = mean(Robs);
    n_pts = 100;
    [G_py,G_px] = ksdensity(G);
    [G_hy,G_hx] = hist(G,n_pts);
    G_hy = G_hy/sum(G_hy);
    %evaluate conditional information values
    h_pred_resp = logexp_fun(G_hx,ov_spk_nl_params(cc,:));
    poiss_prob = poisspdf(repmat(poss_spks,n_pts,1),repmat(h_pred_resp',1,n_symbs));
    marg_prob = sum(bsxfun(@times,poiss_prob,G_hy'));
%     ov_info_hor(cc) = sum(G_hy.*h_pred_resp.*log2(h_pred_resp/avg_rate));
    ov_info_hor(cc) = sum(G_hy'.*sum(poiss_prob.*log2(bsxfun(@rdivide,poiss_prob,marg_prob)),2));
    for tt = 1:n_gsac_bins
        cur_pred_resp = logexp_fun(G_hx,cur_spk_nl_params(tt,:));
        poiss_prob = poisspdf(repmat(poss_spks,n_pts,1),repmat(cur_pred_resp',1,n_symbs));
        marg_prob = sum(bsxfun(@times,poiss_prob,G_hy'));
        sacdep_info_hor(cc,tt) = sum(G_hy'.*sum(poiss_prob.*log2(bsxfun(@rdivide,poiss_prob,marg_prob)),2));
        
%         %        cur_pred_resp = squeeze(nim_g(cc,tt,:));
%         %        cur_pred_resp = interp1(gnim_bin_cents(cc,:),cur_pred_resp,G_hx);
%         sacdep_info_hor(cc,tt) = sum(G_hy.*cur_pred_resp.*log2(cur_pred_resp/avg_sacta(tt)));
    end
    
    
    
    %linear term
    lin_fg = fgint(:,1);
    pp_nim_int(1,:) = prctile(lin_fg,bin_edges);
    lingnim_bin_cents(1,:) = 0.5*pp_nim_int(1,1:end-1) + 0.5*pp_nim_int(1,2:end);
    lingnim_bin_cents(1,end) = 2*lingnim_bin_cents(1,end-1) - lingnim_bin_cents(1,end-2);
    lingnim_bin_cents(1,1) = 2*lingnim_bin_cents(1,2) - lingnim_bin_cents(1,3);
    nim_ling = zeros(n_gsac_bins,nbins);
    for tt = 1:n_gsac_bins
        cur_inds = find(trial_sac_mat(cur_tr_inds,tt) == 1);
        for ii = 1:nbins
            cur_set = cur_inds(lin_fg(cur_inds) >= pp_nim_int(1,ii) & lin_fg(cur_inds) < pp_nim_int(1,ii+1));
            nim_ling(tt,ii) = mean(Robs(cur_set));
        end
        nim_ling(tt,:) = jmm_smooth_1d_cor(squeeze(nim_ling(tt,:)),smooth_fac);
    end
    ov_nim_ling(1,:) = zeros(1,nbins);
    for ii = 1:nbins
        cur_set = find(lin_fg >= pp_nim_int(1,ii) & lin_fg < pp_nim_int(1,ii+1));
        ov_nim_ling(1,ii) = mean(Robs(cur_set));
    end
    ov_nim_ling(1,:) = jmm_smooth_1d_cor(ov_nim_ling(1,:),smooth_fac);
    
    %excitatory squared terms
    if npos_stc(cc) > 0
        exc_fg = sum(fgint(:,2:(1+npos_stc(cc))),2);
        pp_nim_int(1,:) = prctile(exc_fg,bin_edges);
        excgnim_bin_cents(1,:) = 0.5*pp_nim_int(1,1:end-1) + 0.5*pp_nim_int(1,2:end);
        excgnim_bin_cents(1,end) = 2*excgnim_bin_cents(1,end-1) - excgnim_bin_cents(1,end-2);
        excgnim_bin_cents(1,1) = 2*excgnim_bin_cents(1,2) - excgnim_bin_cents(1,3);
        nim_excg = zeros(n_gsac_bins,nbins);
        for tt = 1:n_gsac_bins
            cur_inds = find(trial_sac_mat(cur_tr_inds,tt) == 1);
            for ii = 1:nbins
                cur_set = cur_inds(exc_fg(cur_inds) >= pp_nim_int(1,ii) & exc_fg(cur_inds) < pp_nim_int(1,ii+1));
                nim_excg(tt,ii) = mean(Robs(cur_set));
            end
            nim_excg(tt,:) = jmm_smooth_1d_cor(squeeze(nim_excg(tt,:)),smooth_fac);
        end
        ov_nim_excg(1,:) = zeros(1,nbins);
        for ii = 1:nbins
            cur_set = find(exc_fg >= pp_nim_int(1,ii) & exc_fg < pp_nim_int(1,ii+1));
            ov_nim_excg(1,ii) = mean(Robs(cur_set));
        end
        ov_nim_excg(1,:) = jmm_smooth_1d_cor(ov_nim_excg(1,:),smooth_fac);
    else
        exc_fg = [];
    end
    %suppressive squared terms
    if nneg_stc(cc) > 0
        sup_fg = sum(fgint(:,(2+npos_stc(cc)):end),2);
        pp_nim_int(1,:) = prctile(sup_fg,bin_edges);
        supgnim_bin_cents(1,:) = 0.5*pp_nim_int(1,1:end-1) + 0.5*pp_nim_int(1,2:end);
        supgnim_bin_cents(1,end) = 2*supgnim_bin_cents(1,end-1) - supgnim_bin_cents(1,end-2);
        supgnim_bin_cents(1,1) = 2*supgnim_bin_cents(1,2) - supgnim_bin_cents(1,3);
        nim_supg = zeros(n_gsac_bins,nbins);
        for tt = 1:n_gsac_bins
            cur_inds = find(trial_sac_mat(cur_tr_inds,tt) == 1);
            for ii = 1:nbins
                cur_set = cur_inds(sup_fg(cur_inds) >= pp_nim_int(1,ii) & sup_fg(cur_inds) < pp_nim_int(1,ii+1));
                nim_supg(tt,ii) = mean(Robs(cur_set));
            end
            nim_supg(tt,:) = jmm_smooth_1d_cor(squeeze(nim_supg(tt,:)),smooth_fac);
        end
        ov_nim_supg(1,:) = zeros(1,nbins);
        for ii = 1:nbins
            cur_set = find(sup_fg >= pp_nim_int(1,ii) & sup_fg < pp_nim_int(1,ii+1));
            ov_nim_supg(1,ii) = mean(Robs(cur_set));
        end
        ov_nim_supg(1,:) = jmm_smooth_1d_cor(ov_nim_supg(1,:),smooth_fac);
    else
        sub_fg = [];
    end
    
    %     cnt = 0;
    %     L2_params = create_L2_params([],[1 n_gsac_bins]+cnt,n_gsac_bins);
    %     cnt = cnt + n_gsac_bins;
    %     L2_params = create_L2_params([],[1 n_gsac_bins]+cnt,n_gsac_bins);
    %     cnt = cnt + n_gsac_bins;
    %     if npos_stc(cc) > 0
    %         L2_params = create_L2_params(L2_params,[1 n_gsac_bins]+cnt,n_gsac_bins);
    %         cnt = cnt + n_gsac_bins;
    %     end
    %     if nneg_stc(cc) > 0
    %         L2_params = create_L2_params(L2_params,[1 n_gsac_bins]+cnt,n_gsac_bins);
    %         cnt = cnt + n_gsac_bins;
    %     end
    %     cur_Xmat = [trial_sac_mat(cur_tr_inds,:) bsxfun(@times,trial_sac_mat(cur_tr_inds,:),lin_fg) bsxfun(@times,trial_sac_mat(cur_tr_inds,:),exc_fg) bsxfun(@times,trial_sac_mat(cur_tr_inds,:),sup_fg) lin_fg exc_fg sup_fg];
    %     [fitmod] = regGLM_fit(cur_Xmat,Robs,L2_params,ones(1,length(L2_params))*200,[],[],0);
    
    
    
    %     avg_sacta = squeeze(mean(nim_g(cc,:,:),3));
    [~,min_point] = min(avg_sacta);
    [~,max_point] = max(avg_sacta(min_point:end));
    max_point = max_point + min_point-1;
    
    %     point_range = (min_point-2):(max_point+2);
    point_range = (min_point):(max_point);
    cmap = jet(length(point_range));
    
    if to_print == 1
        figure('visible','off');
    else
        figure
    end
    subplot(1,3,1)
    plot(gsac_bin_cents,avg_sacta/dt,'ko-','linewidth',2)
    hold on
    for ii = 1:length(point_range)
        plot(gsac_bin_cents(point_range(ii)),avg_sacta(point_range(ii))/dt,'*','markersize',16,'color',cmap(ii,:),'linewidth',2)
    end
    title(sprintf('HORI Unit %d using %d Pos %d Neg',cc,npos_stc(cc),nneg_stc(cc)));
    xlabel('Time (s)'); ylabel('Firing rate (Hz)');
    subplot(1,3,2)
    hold on
    %     plot(gnim_bin_cents(cc,use_data_range),ov_nimg(cc,use_data_range)/dt,'ko','markersize',8,'linewidth',2)
    plot(gnim_bin_cents(cc,use_data_range),ov_pred_resp(cc,use_data_range)/dt,'k','linewidth',1)
    for ii = 1:length(point_range)
        %         plot(gnim_bin_cents(cc,use_data_range),squeeze(nim_g(cc,point_range(ii),use_data_range))/dt,'o','color',cmap(ii,:),'linewidth',1);
        plot(gnim_bin_cents(cc,use_data_range),pred_resp(point_range(ii),:)/dt,'color',cmap(ii,:),'linewidth',1);
    end
    xlim_jmm(gnim_bin_cents(cc,use_data_range([1 end])))
    yl = ylim();
    plot(G_px,G_py/max(G_py)*yl(2)*0.8,'k--');
    xlabel('Stim output'); ylabel('Firing rate (Hz)');
    subplot(1,3,3)
    plot(gsac_bin_cents,sacdep_info_hor(cc,:)/dt,'o-')
    xlim_jmm(gsac_bin_cents([1 end]));
    xl = xlim();
    line(xl,ov_info_hor([cc cc])/dt,'color','k','linestyle','--');
    xlabel('Time (s)'); ylabel('Information (bits/sec)');
    if to_print == 1
        fillPage(gcf,'papersize',[18 6]);
        if first_fig == 1
            print(fullmod_figname,'-dpsc');
        else
            print(fullmod_figname,'-dpsc','-append');
        end
        close all
    else
        %                 pause
        %                 close all
    end
    
    
    if to_print == 1
        figure('visible','off');
    else
        figure
    end
    filtK = [fit1_hor(cc).mods(:).filtK];
    filtK = bsxfun(@rdivide,filtK,max(abs(filtK)));
    max_modnum = min(3,max(npos_stc(cc),nneg_stc(cc)));
    subplot(3,4,1)
    imagesc(reshape(filtK(:,1),flen,nPix)); caxis([-1 1]);
    subplot(3,4,4);
    hold on
    plot(lingnim_bin_cents(use_data_range),ov_nim_ling(use_data_range)/dt,'k','linewidth',2)
    for ii = 1:length(point_range)
        plot(lingnim_bin_cents(use_data_range),squeeze(nim_ling(point_range(ii),use_data_range))/dt,'.-','color',cmap(ii,:),'linewidth',1);
    end
    xlim_jmm(lingnim_bin_cents(use_data_range([1 end])));
    for mm = 1:min(3,npos_stc(cc))
        subplot(3,4,mm+4)
        imagesc(reshape(filtK(:,1+mm),flen,nPix));caxis([-1 1])
    end
    if npos_stc(cc) > 0
        subplot(3,4,8);
        hold on
        plot(excgnim_bin_cents(use_data_range),ov_nim_excg(use_data_range)/dt,'k','linewidth',2)
        for ii = 1:length(point_range)
            plot(excgnim_bin_cents(use_data_range),squeeze(nim_excg(point_range(ii),use_data_range))/dt,'.-','color',cmap(ii,:),'linewidth',1);
        end
        xlim_jmm(excgnim_bin_cents(use_data_range([1 end])));
    end
    for mm = 1:min(3,nneg_stc(cc))
        subplot(3,4,8+mm)
        imagesc(reshape(filtK(:,npos_stc(cc)+1+mm),flen,nPix));caxis([-1 1])
    end
    if nneg_stc(cc) > 0
        subplot(3,4,12);
        hold on
        plot(supgnim_bin_cents(use_data_range),ov_nim_supg(use_data_range)/dt,'k','linewidth',2)
        for ii = 1:length(point_range)
            plot(supgnim_bin_cents(use_data_range),squeeze(nim_supg(point_range(ii),use_data_range))/dt,'.-','color',cmap(ii,:),'linewidth',1);
        end
        xlim_jmm(supgnim_bin_cents(use_data_range([1 end])));
    end
    if to_print == 1
        fillPage(gcf,'papersize',[15 12]);
        print(fullmod_figname,'-dpsc','-append');
        close all
    else
        %         pause
        %         close all
    end
    
    
    
    %% FOR VER BARS
    cur_evec_diff = diff(fliplr(unit_evecs_ver(cc,:)));
    npos_stc(cc) = find(cur_evec_diff > stc_thresh,1,'first')-1;
    nneg_stc(cc) = find(cur_evec_diff(end:-1:1) > stc_thresh,1,'first')-1;
    
    cur_tr_inds = tr_inds(ismember(all_exptvec(tr_inds),ver_expts));
    NT = length(cur_tr_inds);
    
    Robs = all_binned_spks(cur_tr_inds,cc);
    
    if  ~exist('fit1_ver','var') || cc > length(fit1_ver) || isempty(fit1_ver(cc).mods) || refit_models == 1
        fprintf('Fitting unit %d, with %d exc and %d sup squared\n',cc,npos_stc(cc),nneg_stc(cc));
        mod_signs = [1 ones(1,npos_stc(cc)) -1*ones(1,nneg_stc(cc))];
        NL_types = cat(2,'lin', repmat({'quad'},1,npos_stc(cc) + nneg_stc(cc)));
        %     NL_types = repmat({'threshlin'},1,1+npos_stc(cc) + nneg_stc(cc));
        
        fit0_ver(cc) = NIMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params); %initialize NIM
        %adjust regularization of new filter
        fit0_ver(cc) = NIMadjust_regularization(fit0_ver(cc),[2:length(mod_signs)],'lambda_d2XT',ql_d2XT);
        
        sub_sample = randperm(NT);
        sub_sample = sub_sample(1:round(NT/sub_samp_fac));
        
        fit0_ver(cc) = NIMfit_filters(fit0_ver(cc),Robs(sub_sample),all_Xmat(cur_tr_inds(sub_sample),:),[],[],silent,init_optim_params); %fit stimulus filters
        fit0_ver(cc) = NIMfit_filters(fit0_ver(cc),Robs,all_Xmat(cur_tr_inds,:),[],[],silent,init_optim_params); %fit stimulus filters
        
        fit1_ver(cc) = NIMadjust_regularization(fit0_ver(cc),1,'lambda_L1',l_L1);
        fit1_ver(cc) = NIMadjust_regularization(fit1_ver(cc),[2:length(mod_signs)],'lambda_L1',ql_L1);
        fit1_ver(cc) = NIMfit_filters(fit1_ver(cc),Robs,all_Xmat(cur_tr_inds,:),[],[],silent,optim_params); %fit stimulus filters
    end
    [LL, penLL, pred_rate, G, gint, fgint] = NIMmodel_eval(fit1_ver(cc),Robs,all_Xmat(cur_tr_inds,:));
    G = G - fit1_ver(cc).spk_NL_params(1);
    fit2_ver = NIMfit_logexp_spkNL_withoffset(fit1_ver(cc),Robs,all_Xmat(cur_tr_inds,:),[],0,[]);
    ov_spk_nl_params(cc,:) = fit2_ver.spk_NL_params;
    
    pp_nim = prctile(G,bin_edges);
    gnim_bin_cents(cc,:) = 0.5*pp_nim(1:end-1) + 0.5*pp_nim(2:end);
    gnim_bin_cents(cc,end) = 2*gnim_bin_cents(cc,end-1) - gnim_bin_cents(cc,end-2);
    gnim_bin_cents(cc,1) = 2*gnim_bin_cents(cc,2) - gnim_bin_cents(cc,3);
    ov_nimg(cc,:) = zeros(1,nbins);
    for ii = 1:nbins
        cur_set = find(G >= pp_nim(ii) & G < pp_nim(ii+1));
        ov_nimg(cc,ii) = mean(Robs(cur_set));
    end
    ov_nimg(cc,:) = jmm_smooth_1d_cor(ov_nimg(cc,:),smooth_fac);
    [nim_g_dist(cc,:),nim_g_x(cc,:)] = ksdensity(G);
    %     ov_pfit(cc,:) = polyfit(gnim_bin_cents(cc,use_data_range),ov_nimg(cc,use_data_range),1);
    ov_pred_resp(cc,:) = logexp_fun(gnim_bin_cents(cc,:),ov_spk_nl_params(cc,:));
    
    %     ov_proj_zero = -ov_pfit(cc,2)/ov_pfit(cc,1);
    
    avg_sacta = zeros(1,n_gsac_bins);
    nim_g(cc,:,:) = zeros(n_gsac_bins,nbins);
    for tt = 1:n_gsac_bins
        cur_inds = find(trial_sac_mat(cur_tr_inds,tt) == 1);
        avg_sacta(tt) = mean(Robs(cur_inds));
        for ii = 1:nbins
            cur_set = cur_inds(G(cur_inds) >= pp_nim(ii) & G(cur_inds) < pp_nim(ii+1));
            nim_g(cc,tt,ii) = mean(Robs(cur_set));
        end
        %         pfit(cc,tt,:) = polyfit(gnim_bin_cents(cc,use_data_range),squeeze(nim_g(cc,tt,use_data_range))',1);
        %         cur_offset(tt) = polyval(squeeze(pfit(cc,tt,:)),ov_proj_zero);
        cur_fit = NIMfit_logexp_spkNL_withoffset(fit2_ver,Robs(cur_inds),all_Xmat(cur_tr_inds(cur_inds),:),[],0,[1 2]);
        cur_spk_nl_params(tt,:) = cur_fit.spk_NL_params;
        pred_resp(tt,:) = logexp_fun(gnim_bin_cents(cc,:),cur_spk_nl_params(tt,:));
        nim_g(cc,tt,:) = jmm_smooth_1d_cor(squeeze(nim_g(cc,tt,:)),smooth_fac);
    end
    
    poss_spks = 0:10;
    n_symbs = length(poss_spks);
    
    avg_rate = mean(Robs);
    n_pts = 100;
    [G_py,G_px] = ksdensity(G);
    [G_hy,G_hx] = hist(G,n_pts);
    G_hy = G_hy/sum(G_hy);
    %evaluate conditional information values
%     h_pred_resp = logexp_fun(G_hx,ov_spk_nl_params(cc,:));
%     ov_info_ver(cc) = sum(G_hy.*h_pred_resp.*log2(h_pred_resp/avg_rate));
    h_pred_resp = logexp_fun(G_hx,ov_spk_nl_params(cc,:));
    poiss_prob = poisspdf(repmat(poss_spks,n_pts,1),repmat(h_pred_resp',1,n_symbs));
    marg_prob = sum(bsxfun(@times,poiss_prob,G_hy'));
    ov_info_ver(cc) = sum(G_hy'.*sum(poiss_prob.*log2(bsxfun(@rdivide,poiss_prob,marg_prob)),2));
    for tt = 1:n_gsac_bins
        cur_pred_resp = logexp_fun(G_hx,cur_spk_nl_params(tt,:));
        poiss_prob = poisspdf(repmat(poss_spks,n_pts,1),repmat(cur_pred_resp',1,n_symbs));
        marg_prob = sum(bsxfun(@times,poiss_prob,G_hy'));
        sacdep_info_ver(cc,tt) = sum(G_hy'.*sum(poiss_prob.*log2(bsxfun(@rdivide,poiss_prob,marg_prob)),2));
        
%         %        cur_pred_resp = squeeze(nim_g(cc,tt,:));
%         %        cur_pred_resp = interp1(gnim_bin_cents(cc,:),cur_pred_resp,G_hx);
%         sacdep_info_ver(cc,tt) = sum(G_hy.*cur_pred_resp.*log2(cur_pred_resp/avg_sacta(tt)));
    end
    
    
    %linear term
    lin_fg = fgint(:,1);
    pp_nim_int(1,:) = prctile(lin_fg,bin_edges);
    lingnim_bin_cents(1,:) = 0.5*pp_nim_int(1,1:end-1) + 0.5*pp_nim_int(1,2:end);
    lingnim_bin_cents(1,end) = 2*lingnim_bin_cents(1,end-1) - lingnim_bin_cents(1,end-2);
    lingnim_bin_cents(1,1) = 2*lingnim_bin_cents(1,2) - lingnim_bin_cents(1,3);
    nim_ling = zeros(n_gsac_bins,nbins);
    for tt = 1:n_gsac_bins
        cur_inds = find(trial_sac_mat(cur_tr_inds,tt) == 1);
        for ii = 1:nbins
            cur_set = cur_inds(lin_fg(cur_inds) >= pp_nim_int(1,ii) & lin_fg(cur_inds) < pp_nim_int(1,ii+1));
            nim_ling(tt,ii) = mean(Robs(cur_set));
        end
        nim_ling(tt,:) = jmm_smooth_1d_cor(squeeze(nim_ling(tt,:)),smooth_fac);
    end
    ov_nim_ling(1,:) = zeros(1,nbins);
    for ii = 1:nbins
        cur_set = find(lin_fg >= pp_nim_int(1,ii) & lin_fg < pp_nim_int(1,ii+1));
        ov_nim_ling(1,ii) = mean(Robs(cur_set));
    end
    ov_nim_ling(1,:) = jmm_smooth_1d_cor(ov_nim_ling(1,:),smooth_fac);
    
    %excitatory squared terms
    if npos_stc(cc) > 0
        exc_fg = sum(fgint(:,2:(1+npos_stc(cc))),2);
        pp_nim_int(1,:) = prctile(exc_fg,bin_edges);
        excgnim_bin_cents(1,:) = 0.5*pp_nim_int(1,1:end-1) + 0.5*pp_nim_int(1,2:end);
        excgnim_bin_cents(1,end) = 2*excgnim_bin_cents(1,end-1) - excgnim_bin_cents(1,end-2);
        excgnim_bin_cents(1,1) = 2*excgnim_bin_cents(1,2) - excgnim_bin_cents(1,3);
        nim_excg = zeros(n_gsac_bins,nbins);
        for tt = 1:n_gsac_bins
            cur_inds = find(trial_sac_mat(cur_tr_inds,tt) == 1);
            for ii = 1:nbins
                cur_set = cur_inds(exc_fg(cur_inds) >= pp_nim_int(1,ii) & exc_fg(cur_inds) < pp_nim_int(1,ii+1));
                nim_excg(tt,ii) = mean(Robs(cur_set));
            end
            nim_excg(tt,:) = jmm_smooth_1d_cor(squeeze(nim_excg(tt,:)),smooth_fac);
        end
        ov_nim_excg(1,:) = zeros(1,nbins);
        for ii = 1:nbins
            cur_set = find(exc_fg >= pp_nim_int(1,ii) & exc_fg < pp_nim_int(1,ii+1));
            ov_nim_excg(1,ii) = mean(Robs(cur_set));
        end
        ov_nim_excg(1,:) = jmm_smooth_1d_cor(ov_nim_excg(1,:),smooth_fac);
    end
    %suppressive squared terms
    if nneg_stc(cc) > 0
        sup_fg = sum(fgint(:,(2+npos_stc(cc)):end),2);
        pp_nim_int(1,:) = prctile(sup_fg,bin_edges);
        supgnim_bin_cents(1,:) = 0.5*pp_nim_int(1,1:end-1) + 0.5*pp_nim_int(1,2:end);
        supgnim_bin_cents(1,end) = 2*supgnim_bin_cents(1,end-1) - supgnim_bin_cents(1,end-2);
        supgnim_bin_cents(1,1) = 2*supgnim_bin_cents(1,2) - supgnim_bin_cents(1,3);
        nim_supg = zeros(n_gsac_bins,nbins);
        for tt = 1:n_gsac_bins
            cur_inds = find(trial_sac_mat(cur_tr_inds,tt) == 1);
            for ii = 1:nbins
                cur_set = cur_inds(sup_fg(cur_inds) >= pp_nim_int(1,ii) & sup_fg(cur_inds) < pp_nim_int(1,ii+1));
                nim_supg(tt,ii) = mean(Robs(cur_set));
            end
            nim_supg(tt,:) = jmm_smooth_1d_cor(squeeze(nim_supg(tt,:)),smooth_fac);
        end
        ov_nim_supg(1,:) = zeros(1,nbins);
        for ii = 1:nbins
            cur_set = find(sup_fg >= pp_nim_int(1,ii) & sup_fg < pp_nim_int(1,ii+1));
            ov_nim_supg(1,ii) = mean(Robs(cur_set));
        end
        ov_nim_supg(1,:) = jmm_smooth_1d_cor(ov_nim_supg(1,:),smooth_fac);
    end
    
    
    %     avg_sacta = squeeze(mean(nim_g(cc,:,:),3));
    [~,min_point] = min(avg_sacta);
    [~,max_point] = max(avg_sacta(min_point:end));
    max_point = max_point + min_point-1;
    
    %     point_range = (min_point-2):(max_point+2);
    point_range = (min_point):(max_point);
    cmap = jet(length(point_range));
    
    if to_print == 1
        figure('visible','off');
    else
        figure
    end
    subplot(1,3,1)
    plot(gsac_bin_cents,avg_sacta/dt,'ko-','linewidth',2)
    hold on
    for ii = 1:length(point_range)
        plot(gsac_bin_cents(point_range(ii)),avg_sacta(point_range(ii))/dt,'*','markersize',16,'color',cmap(ii,:),'linewidth',2)
    end
    title(sprintf('VERT Unit %d using %d Pos %d Neg',cc,npos_stc(cc),nneg_stc(cc)));
    xlabel('Time (s)'); ylabel('Firing rate (Hz)');
    subplot(1,3,2)
    hold on
    %     plot(gnim_bin_cents(cc,use_data_range),ov_nimg(cc,use_data_range)/dt,'ko','markersize',8,'linewidth',2)
    plot(gnim_bin_cents(cc,use_data_range),ov_pred_resp(cc,use_data_range)/dt,'k','linewidth',1)
    for ii = 1:length(point_range)
        %         plot(gnim_bin_cents(cc,use_data_range),squeeze(nim_g(cc,point_range(ii),use_data_range))/dt,'o','color',cmap(ii,:),'linewidth',1);
        plot(gnim_bin_cents(cc,use_data_range),pred_resp(point_range(ii),:)/dt,'color',cmap(ii,:),'linewidth',1);
    end
    xlim_jmm(gnim_bin_cents(cc,use_data_range([1 end])))
    yl = ylim();
    plot(G_px,G_py/max(G_py)*yl(2)*0.8,'k--');
    xlabel('Stim output'); ylabel('Firing rate (Hz)');
    subplot(1,3,3)
    plot(gsac_bin_cents,sacdep_info_ver(cc,:)/dt,'o-')
    xlabel('Time (s)'); ylabel('Information (bits/sec)');
    xlim_jmm(gsac_bin_cents([1 end]));
    xl = xlim();
    line(xl,ov_info_ver([cc cc])/dt,'color','k','linestyle','--');
    if to_print == 1
        fillPage(gcf,'papersize',[18 6]);
        %         if first_fig == 1
        %             print(fullmod_figname,'-dpsc');
        %         else
        print(fullmod_figname,'-dpsc','-append');
        %         end
        close all
    else
        %                 pause
        %                 close all
    end
    
    if to_print == 1
        figure('visible','off');
    else
        figure
    end
    filtK = [fit1_ver(cc).mods(:).filtK];
    filtK = bsxfun(@rdivide,filtK,max(abs(filtK)));
    max_modnum = min(3,max(npos_stc(cc),nneg_stc(cc)));
    subplot(3,4,1)
    imagesc(reshape(filtK(:,1),flen,nPix)); caxis([-1 1]);
    subplot(3,4,4);
    hold on
    plot(lingnim_bin_cents(use_data_range),ov_nim_ling(use_data_range)/dt,'k','linewidth',2)
    for ii = 1:length(point_range)
        plot(lingnim_bin_cents(use_data_range),squeeze(nim_ling(point_range(ii),use_data_range))/dt,'.-','color',cmap(ii,:),'linewidth',1);
    end
    xlim_jmm(lingnim_bin_cents(use_data_range([1 end])));
    for mm = 1:min(3,npos_stc(cc))
        subplot(3,4,mm+4)
        imagesc(reshape(filtK(:,1+mm),flen,nPix));caxis([-1 1])
    end
    if npos_stc(cc) > 0
        subplot(3,4,8);
        hold on
        plot(excgnim_bin_cents(use_data_range),ov_nim_excg(use_data_range)/dt,'k','linewidth',2)
        for ii = 1:length(point_range)
            plot(excgnim_bin_cents(use_data_range),squeeze(nim_excg(point_range(ii),use_data_range))/dt,'.-','color',cmap(ii,:),'linewidth',1);
        end
        xlim_jmm(excgnim_bin_cents(use_data_range([1 end])));
    end
    for mm = 1:min(3,nneg_stc(cc))
        subplot(3,4,8+mm)
        imagesc(reshape(filtK(:,npos_stc(cc)+1+mm),flen,nPix));caxis([-1 1])
    end
    if nneg_stc(cc) > 0
        subplot(3,4,12);
        hold on
        plot(supgnim_bin_cents(use_data_range),ov_nim_supg(use_data_range)/dt,'k','linewidth',2)
        for ii = 1:length(point_range)
            plot(supgnim_bin_cents(use_data_range),squeeze(nim_supg(point_range(ii),use_data_range))/dt,'.-','color',cmap(ii,:),'linewidth',1);
        end
        xlim(supgnim_bin_cents(use_data_range([1 end])));
    end
    if to_print == 1
        fillPage(gcf,'papersize',[15 12]);
        print(fullmod_figname,'-dpsc','-append');
        close all
    else
        %         pause
        %         close all
    end
    
    first_fig = 0;
end

%%
% clear sac_sta sac_trig_avg
% cc = 42;
% spikebins = convert_to_spikebins(all_binned_spks(tr_inds,cc));
% spike_cond_stim = all_Xmat(tr_inds(spikebins),:);
% sta      = mean(spike_cond_stim) - mean(all_Xmat(tr_inds,:));
% ov_sta = sta/norm(sta);
%
% proj_mat = ov_sta'/(ov_sta*ov_sta')*ov_sta;
% stim_proj = all_Xmat(tr_inds,:) - all_Xmat(tr_inds,:)*proj_mat;
% % stim_proj = stim_emb;
% stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
% [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
% ov_stcs  = evecs(:,end);
%
% sac_trig_avg(cc,:) = zeros(1,n_gsac_bins);
% for tt = 1:n_gsac_bins
%     tt
%
%     cur_inds = tr_inds(trial_sac_mat(tr_inds,tt) == 1);
%     sac_trig_avg(cc,tt) = mean(all_binned_spks(cur_inds,cc));
%
%     spikebins = convert_to_spikebins(all_binned_spks(cur_inds,cc));
%     spike_cond_stim = all_Xmat(cur_inds(spikebins),:);
%     sta      = mean(spike_cond_stim) - mean(all_Xmat(cur_inds,:));
%     sac_sta(tt,:) = sta/norm(sta);
%
%     proj_mat = sac_sta(tt,:)'/(sac_sta(tt,:)*sac_sta(tt,:)')*sac_sta(tt,:);
%     stim_proj = all_Xmat(cur_inds,:) - all_Xmat(cur_inds,:)*proj_mat;
%     % stim_proj = stim_emb;
%     stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
%     [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
%     sac_stcs(tt,:)  = evecs(:,end);
%
%
% end
%
% %%
% for tt = 1:n_gsac_bins
% subplot(3,1,1)
% % imagesc(reshape(ov_sta,flen,nPix))
% % ca = max(abs(ov_sta)); caxis([-ca ca]);
% %  subplot(3,1,2)
% %    imagesc(reshape(sac_sta(tt,:),flen,nPix))
% %     caxis([-ca ca]);
%
% imagesc(reshape(ov_stcs,flen,nPix))
% ca = max(abs(ov_stcs)); caxis([-ca ca]);
%  subplot(3,1,2)
%    imagesc(reshape(sac_stcs(tt,:),flen,nPix))
%     caxis([-ca ca]);
%
% subplot(3,1,3)
%  plot(gsac_bin_cents*dt,sac_trig_avg(cc,:),'.-')
%  hold on
%  plot(gsac_bin_cents(tt)*dt,sac_trig_avg(cc,tt),'ro')
%
%  pause
%  clf
% end