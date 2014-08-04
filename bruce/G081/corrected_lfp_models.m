clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
addpath('~/James_scripts/bruce/G081/');
%%
stim_fs = 100; %in Hz
use_sus = 1:96;
Fs = 3e4;
dsf = 60;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1 60]/niqf);
use_lfps = [1:1:96];

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
flen = 20;
beg_buffer = round(stim_fs*0.15);
bar_oris = [90];
un_bar_pos = all_un_bar_pos(:,3);

fprintf('Analyzing %d ori expts\n',bar_oris);

%expts with X deg bars and any back (including sim sacs)
cur_expt_set = find(is_bar_expt==1 & expt_bar_ori == bar_oris);

% %expts with X deg bars and any back (including sim sacs)
% cur_expt_set = find(is_bar_expt==1 & expt_bar_ori == bar_oris & expt_image_back == 0);

% %expts with X deg bars and any back (including sim sacs)
% cur_expt_set = find(is_bar_expt==1 & expt_bar_ori == bar_oris & expt_image_back == 1);

cur_expt_set(cur_expt_set<=8) = []; %this expt had continuous bar luminance
cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?
cur_expt_set(cur_expt_set > 60) = [];
%% LOAD AND INTERPOLATE MODELS
load ./CellList.mat
good_sus = find(all(CellList(:,:,1) > 0));

load ./full_eye_correct_90deg
eye_times = eye_times - 0.05;
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
all_trial_start_times = [];
all_trial_end_times = [];
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    n_trials = length(Expts{cur_expt}.Trials);
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    
    all_trial_start_times = [all_trial_start_times; trial_start_times(:)];
    all_trial_end_times = [all_trial_end_times; trial_end_times(:)];
    
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
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_rel_stimes = [all_rel_stimes; cur_stim_times- trial_start_times(tt)];
        all_rel_etimes = [all_rel_etimes; trial_end_times(tt) - cur_stim_times];
        all_Op = [all_Op; cur_Op];
        all_phase = [all_phase; cur_phase'];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(size(cur_stim_times))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_stim_times))*tt];
    end
    
end

%%
all_Vmat = [];
all_tax = [];
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
    
    all_Vmat = [all_Vmat; Vmat];
    all_tax = [all_tax; t_ax'];
end

%%
in_trial = zeros(size(all_tax));
interp_trial_start_inds = round(interp1(all_tax,1:length(all_tax),all_trial_start_times));
interp_trial_end_inds = round(interp1(all_tax,1:length(all_tax),all_trial_end_times));
n_trials = length(all_trial_start_times);
for i = 1:n_trials
    cur_inds = interp_trial_start_inds(i):interp_trial_end_inds(i);
    in_trial(cur_inds) = 1;
end
all_Vmat = bsxfun(@minus,all_Vmat,mean(all_Vmat(in_trial==1,:)));
all_Vmat = bsxfun(@rdivide,all_Vmat,std(all_Vmat(in_trial==1,:)));
%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
% load ./jbeG081.em.mat
% all_eye_vals = [];
% all_eye_speed = [];
% all_eye_ts = [];
% for ee = 1:length(cur_expt_set);
%     cur_set = find(all_exptvec==ee);
%     [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_stim_times(cur_set([1 end])));
%     
%     eye_dt = median(diff(eye_ts_interp));
%     eye_fs = 1/eye_dt;
%     lEyeXY = eye_vals_interp(:,1:2);
%     rEyeXY = eye_vals_interp(:,3:4);
%     clear sm_avg_eyepos eye_vel
%     sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),3);
%     sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),3);
%     eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
%     eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
%     
%     eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
%     
%     all_eye_vals = [all_eye_vals; lEyeXY rEyeXY];
%     all_eye_speed = [all_eye_speed; eye_speed];
%     all_eye_ts = [all_eye_ts; eye_ts_interp'];
% end
% back_pts = 1 + find(diff(all_eye_ts) <= 0);
% double_samples = [];
% for i = 1:length(back_pts)
%     next_forward = find(all_eye_ts > all_eye_ts(back_pts(i)-1),1,'first');
%     double_samples = [double_samples back_pts(i):next_forward];
% end
% all_eye_ts(double_samples) = [];
% all_eye_speed(double_samples) = [];
% all_eye_vals(double_samples,:) = [];
% 
% interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_stim_times);
% interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_stim_times);
% 
% orth_eye_pos = interp_eye_vals(:,2);
% bad_pts = find(abs(orth_eye_pos) > 1);

%%
interp_eye_corrections = interp1(eye_times,it_shift_cor{end}*0.04,all_stim_times);
all_cor_Op = all_Op - interp_eye_corrections;

unique_bar_pos = unique(all_Op);
n_bar_pos = length(unique_bar_pos);
up_bar_ax = (unique_bar_pos(1)-0.08):0.04:(unique_bar_pos(end)+0.08);
up_n_bar_pos = length(up_bar_ax);

all_up_bar_dists = nan(length(all_Op),up_n_bar_pos);
for i = 1:up_n_bar_pos
    all_up_bar_dists(:,i) = abs(all_cor_Op - up_bar_ax(i));
end
[~,all_up_bar_inds] = min(all_up_bar_dists,[],2);

% trial_lfp_inds = round(interp1(all_tax,1:length(all_tax),trial_start_times));

expt_lfp_inds = round(interp1(all_tax,1:length(all_tax),all_stim_times));

backlag = round(Fsd*0.05);
forwardlag = round(Fsd*0.25);
lags = -backlag:forwardlag;

%%
back_buffer = 0.2;
use_stims = find(all_rel_stimes > backlag/Fsd & all_rel_stimes > back_buffer & all_rel_etimes > forwardlag/Fsd);
comb_trig_avgs = zeros(n_bar_pos,length(lags),length(use_lfps));
for bb = 1:n_bar_pos
    cur_bar_inds = use_stims(all_Op(use_stims)==unique_bar_pos(bb));
    for ii = 1:length(cur_bar_inds)
        cur_inds = (expt_lfp_inds(cur_bar_inds(ii)) - backlag):(expt_lfp_inds(cur_bar_inds(ii)) + forwardlag);
        comb_trig_avgs(bb,:,:) = squeeze(comb_trig_avgs(bb,:,:)) + all_Vmat(cur_inds,:);
    end
    comb_trig_avgs(bb,:,:) = comb_trig_avgs(bb,:,:)/length(cur_bar_inds);
    cur_n(bb) = length(cur_bar_inds);
end

comb_trig_avgs_cor = zeros(up_n_bar_pos,length(lags),length(use_lfps));
for bb = 1:up_n_bar_pos
    cur_bar_inds = use_stims(all_up_bar_inds(use_stims)==bb);
    for ii = 1:length(cur_bar_inds)
        cur_inds = (expt_lfp_inds(cur_bar_inds(ii)) - backlag):(expt_lfp_inds(cur_bar_inds(ii)) + forwardlag);
        comb_trig_avgs_cor(bb,:,:) = squeeze(comb_trig_avgs_cor(bb,:,:)) + all_Vmat(cur_inds,:);
    end
    comb_trig_avgs_cor(bb,:,:) = comb_trig_avgs_cor(bb,:,:)/length(cur_bar_inds);
    cur_n_cor(bb) = length(cur_bar_inds);
end

black_trig_avgs_cor = zeros(up_n_bar_pos,length(lags),length(use_lfps));
for bb = 1:up_n_bar_pos
    cur_bar_inds = use_stims(all_up_bar_inds(use_stims)==bb & all_phase(use_stims) == 0);
    for ii = 1:length(cur_bar_inds)
        cur_inds = (expt_lfp_inds(cur_bar_inds(ii)) - backlag):(expt_lfp_inds(cur_bar_inds(ii)) + forwardlag);
        black_trig_avgs_cor(bb,:,:) = squeeze(black_trig_avgs_cor(bb,:,:)) + all_Vmat(cur_inds,:);
    end
    black_trig_avgs_cor(bb,:,:) = black_trig_avgs_cor(bb,:,:)/length(cur_bar_inds);
end

white_trig_avgs_cor = zeros(up_n_bar_pos,length(lags),length(use_lfps));
for bb = 1:up_n_bar_pos
    cur_bar_inds = use_stims(all_up_bar_inds(use_stims)==bb & all_phase(use_stims) == pi);
    for ii = 1:length(cur_bar_inds)
        cur_inds = (expt_lfp_inds(cur_bar_inds(ii)) - backlag):(expt_lfp_inds(cur_bar_inds(ii)) + forwardlag);
        white_trig_avgs_cor(bb,:,:) = squeeze(white_trig_avgs_cor(bb,:,:)) + all_Vmat(cur_inds,:);
    end
    white_trig_avgs_cor(bb,:,:) = white_trig_avgs_cor(bb,:,:)/length(cur_bar_inds);
end

%%
save cor_trig_avg_lfps_90deg_all_tshift_zsc lags Fsd up_n_bar_pos n_bar_pos un_bar_pos *trig_avgs* use_lfps unique_bar_pos up_bar_ax


%%
load ./cor_trig_avg_lfps_0deg_all_tshift.mat     

close all
for ll = 1:length(use_lfps)
    subplot(2,2,1)
    pcolor(lags/Fsd,1:up_n_bar_pos,squeeze(comb_trig_avgs_cor(:,:,ll)));shading flat
    ca = max(abs(caxis()));
    caxis([-ca ca])
    
    subplot(2,2,2)
    pcolor(lags/Fsd,1:n_bar_pos,squeeze(comb_trig_avgs(:,:,ll)));shading flat
    caxis([-ca ca])

        subplot(2,2,3)
    pcolor(lags/Fsd,1:up_n_bar_pos,squeeze(white_trig_avgs_cor(:,:,ll)));shading flat
    ca = max(abs(caxis()));
    caxis([-ca ca])
    
    subplot(2,2,4)
    pcolor(lags/Fsd,1:up_n_bar_pos,squeeze(black_trig_avgs_cor(:,:,ll)));shading flat
    caxis([-ca ca])
    
    ll
    pause
    clf
end
%%
load ./cor_trig_avg_lfps_0deg_all_tshift.mat     
comb_0 = comb_trig_avgs_cor;
up_bar_ax_0 = up_bar_ax;
load ./cor_trig_avg_lfps_90deg_all_tshift.mat     
comb_90 = comb_trig_avgs_cor;
up_bar_ax_90 = up_bar_ax;

close all
for ll = 1:length(use_lfps)
    subplot(2,1,1)
    pcolor(lags/Fsd,up_bar_ax_0,squeeze(comb_0(:,:,ll)));shading flat
    ca = max(abs(caxis()));
    caxis([-ca ca])
    yl = ylim();
    line([0 0],yl,'color','k')
    line([0.05 0.05],yl,'color','k')
    line([0.1 0.1],yl,'color','k')
    
    subplot(2,1,2)
    pcolor(lags/Fsd,up_bar_ax_90,squeeze(comb_90(:,:,ll)));shading flat
     ca = max(abs(caxis()));
   caxis([-ca ca])
    yl = ylim();
    line([0 0],yl,'color','k')
   line([0.05 0.05],yl,'color','k')
    line([0.1 0.1],yl,'color','k')

    ll
    pause
    clf
end

%%
% load ./cor_trig_avg_lfps_0deg_all_tshift.mat 

load ./cor_trig_avg_lfps_90deg_all_tshift.mat     
up_bar_ax = -up_bar_ax;

% load ./cor_trig_avg_lfps_45deg_all_tshift.mat     

cd ~/Data/bruce/7_15_12/G034/
load ./gabor_tracking_varmeans
gabor_params_used = gabor_params_f{2};
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

cd ~/Data/bruce/G081

%fit smoothed retinotopic surface
interp_x = nan(96,1);
interp_y = nan(96,1);
tempinds = zeros(10,10);
for i = 1:96
    tempx(Y_pos(i),X_pos(i)) = gabor_params_used(i,1);
    tempy(Y_pos(i),X_pos(i)) = gabor_params_used(i,2);
    tempinds(Y_pos(i),X_pos(i)) = i;
end
weights = ones(10,10);
weights(tempx==0) = 0;
xpos_interp = smoothn(tempx,weights,'robust');
ypos_interp = smoothn(tempy,weights,'robust');
used_inds = find(weights == 1);
tempinds = tempinds(used_inds);
interp_x(tempinds) = xpos_interp(used_inds);
interp_y(tempinds) = ypos_interp(used_inds);
xi = linspace(min(interp_x),max(interp_x),50);
yi = linspace(min(interp_y),max(interp_y),50);
[Xi,Yi] = meshgrid(xi,yi);

id_mat = nan(10,10);
for i = 1:10
    for j = 1:10
        cur = find(X_pos==j&Y_pos==i);
        if ~isempty(cur)
            id_mat(i,j) = cur;
        end
    end
end
use_ids = find(~isnan(id_mat));


% plot(interp_x,interp_y,'ro')

%%
cax = [-12 12]*1e-6;
for cc = 10:length(up_bar_ax)-5
cur_bpos = cc
response_mat = squeeze(comb_trig_avgs_cor(cur_bpos,:,:))';
% figure
% pcolor(lags/Fsd,1:96,response_mat);shading flat

subplot(2,3,1)
target_lag = 0.04;
lag_ind = find(lags/Fsd > target_lag,1,'first');
cur_set = response_mat(:,lag_ind);
F = TriScatteredInterp(interp_x,interp_y,cur_set);
Vq = F(Xi,Yi);
pcolor(xi,yi,Vq); shading flat; colorbar;
caxis(cax)
% xl = xlim();
% line(xl,up_bar_ax([cc cc]),'color','k')
yl = ylim();
line(up_bar_ax([cc cc]),yl,'color','k')
title('40 ms lag')

subplot(2,3,2)
target_lag = 0.05;
lag_ind = find(lags/Fsd > target_lag,1,'first');
cur_set = response_mat(:,lag_ind);
F = TriScatteredInterp(interp_x,interp_y,cur_set);
Vq = F(Xi,Yi);
pcolor(xi,yi,Vq); shading flat; colorbar;
caxis(cax)
% line(xl,up_bar_ax([cc cc]),'color','k')
line(up_bar_ax([cc cc]),yl,'color','k')
title('50 ms lag')

subplot(2,3,3)
target_lag = 0.075;
lag_ind = find(lags/Fsd > target_lag,1,'first');
cur_set = response_mat(:,lag_ind);
F = TriScatteredInterp(interp_x,interp_y,cur_set);
Vq = F(Xi,Yi);
pcolor(xi,yi,Vq); shading flat; colorbar;
caxis(cax)
% line(xl,up_bar_ax([cc cc]),'color','k')
line(up_bar_ax([cc cc]),yl,'color','k')
title('75 ms lag')

subplot(2,3,4)
target_lag = 0.1;
lag_ind = find(lags/Fsd > target_lag,1,'first');
cur_set = response_mat(:,lag_ind);
F = TriScatteredInterp(interp_x,interp_y,cur_set);
Vq = F(Xi,Yi);
pcolor(xi,yi,Vq); shading flat; colorbar;
caxis(cax)
% line(xl,up_bar_ax([cc cc]),'color','k')
line(up_bar_ax([cc cc]),yl,'color','k')
title('100 ms lag')

subplot(2,3,5)
target_lag = 0.125;
lag_ind = find(lags/Fsd > target_lag,1,'first');
cur_set = response_mat(:,lag_ind);
F = TriScatteredInterp(interp_x,interp_y,cur_set);
Vq = F(Xi,Yi);
pcolor(xi,yi,Vq); shading flat; colorbar;
caxis(cax)
% line(xl,up_bar_ax([cc cc]),'color','k')
line(up_bar_ax([cc cc]),yl,'color','k')
title('125 ms lag')

subplot(2,3,6)
target_lag = 0.15;
lag_ind = find(lags/Fsd > target_lag,1,'first');
cur_set = response_mat(:,lag_ind);
F = TriScatteredInterp(interp_x,interp_y,cur_set);
Vq = F(Xi,Yi);
pcolor(xi,yi,Vq); shading flat; colorbar;
caxis(cax)
% line(xl,up_bar_ax([cc cc]),'color','k')
line(up_bar_ax([cc cc]),yl,'color','k')
title('150 ms lag')

pause
clf
end

%%
close all
cc = 18;
up_bar_ax(cc)
cur_bpos = cc;
response_mat = squeeze(comb_trig_avgs_cor(cur_bpos,:,:))';
% response_matb = squeeze(black_trig_avgs_cor(cur_bpos,:,:))';
% response_matw = squeeze(white_trig_avgs_cor(cur_bpos,:,:))';

cax = [-8 8]*1e-6;
dt = 0.001;
cur_lag = 0;
while cur_lag < 0.25
    lag_ind = find(lags/Fsd > cur_lag,1,'first');
    
    subplot(2,1,1)
    cur_set = response_mat(:,lag_ind);
    F = TriScatteredInterp(interp_x,interp_y,cur_set);
    Vq = F(Xi,Yi);
    pcolor(xi,yi,Vq); shading flat; colorbar;
    caxis(cax)
    xl = xlim(); yl = ylim();
    line(xl,up_bar_ax([cc cc])-0.04,'color','k')
    line(xl,up_bar_ax([cc cc])+0.04,'color','k')
%     line(up_bar_ax([cc cc])-0.04,yl,'color','k')
%     line(up_bar_ax([cc cc])+0.04,yl,'color','k')
    title(sprintf('T: %.3f',cur_lag));
    cur_lag = cur_lag + dt;
   
    subplot(2,1,2)
    cur_mat = nan(10,10);
    cur_mat(use_ids) = cur_set(id_mat(use_ids));
     imagesc(cur_mat); shading flat; colorbar;
    caxis(cax)
    xl = xlim(); yl = ylim();
    line(xl,up_bar_ax([cc cc])-0.04,'color','k')
    line(xl,up_bar_ax([cc cc])+0.04,'color','k')
%     line(up_bar_ax([cc cc])-0.04,yl,'color','k')
%     line(up_bar_ax([cc cc])+0.04,yl,'color','k')
    title(sprintf('T: %.3f',cur_lag));
    cur_lag = cur_lag + dt;
   
    
    
    
    pause(0.075)
%     clf


%     cur_set = response_matb(:,lag_ind);
%     F = TriScatteredInterp(interp_x,interp_y,cur_set);
%     Vq = F(Xi,Yi);
%     subplot(2,1,1)
%     pcolor(xi,yi,Vq); shading flat; colorbar;
%     caxis([-1 1]*1e-5)
%     xl = xlim();
%     line(xl,up_bar_ax([cc cc])-0.04,'color','k')
%     line(xl,up_bar_ax([cc cc])+0.04,'color','k')
% %     line(up_bar_ax([cc cc])-0.04,yl,'color','k')
% %     line(up_bar_ax([cc cc])+0.04,yl,'color','k')
%     title(sprintf('T: %.3f',cur_lag));
% 
%     cur_set = response_matw(:,lag_ind);
%     F = TriScatteredInterp(interp_x,interp_y,cur_set);
%     Vq = F(Xi,Yi);
%         subplot(2,1,2)
%     pcolor(xi,yi,Vq); shading flat; colorbar;
%     caxis([-1 1]*1e-5)
%     line(xl,up_bar_ax([cc cc])-0.04,'color','k')
%     line(xl,up_bar_ax([cc cc])+0.04,'color','k')
% %     line(up_bar_ax([cc cc])-0.04,yl,'color','k')
% %     line(up_bar_ax([cc cc])+0.04,yl,'color','k')
%     title(sprintf('T: %.3f',cur_lag));

%     cur_lag = cur_lag + dt;
%     pause(0.1)
%     clf

end

%%
load ./cor_trig_avg_lfps_0deg_all_tshift.mat 
comb_trig_avgs_cor0 = comb_trig_avgs_cor;
up_bar_ax0 = up_bar_ax;

load ./cor_trig_avg_lfps_90deg_all_tshift.mat     
comb_trig_avgs_cor90 = comb_trig_avgs_cor;
up_bar_ax90 = -up_bar_ax;

load ./cor_trig_avg_lfps_45deg_all_tshift.mat     
comb_trig_avgs_cor45 = comb_trig_avgs_cor;
up_bar_ax45 = -up_bar_ax;

clear pc1* pc2*
for cc = 1:length(up_bar_ax);
cur_bpos = cc;

response_mat = squeeze(comb_trig_avgs_cor0(cur_bpos,:,:))';

[coeff,scores,latent] = princomp(response_mat');
cur_set = coeff(:,1);
F = TriScatteredInterp(interp_x,interp_y,cur_set);
pc10(cc,:,:) = F(Xi,Yi);
sc10(cc,:) = scores(:,1);
cur_set = coeff(:,2);
F = TriScatteredInterp(interp_x,interp_y,cur_set);
pc20(cc,:,:) = F(Xi,Yi);
sc20(cc,:) = scores(:,2);

response_mat = squeeze(comb_trig_avgs_cor90(cur_bpos,:,:))';

[coeff,scores,latent] = princomp(response_mat');
cur_set = coeff(:,1);
F = TriScatteredInterp(interp_x,interp_y,cur_set);
pc190(cc,:,:) = F(Xi,Yi);
sc190(cc,:) = scores(:,1);
cur_set = coeff(:,2);
F = TriScatteredInterp(interp_x,interp_y,cur_set);
pc290(cc,:,:) = F(Xi,Yi);
sc290(cc,:) = scores(:,2);

response_mat = squeeze(comb_trig_avgs_cor45(cur_bpos,:,:))';

[coeff,scores,latent] = princomp(response_mat');
cur_set = coeff(:,1);
F = TriScatteredInterp(interp_x,interp_y,cur_set);
pc145(cc,:,:) = F(Xi,Yi);
sc145(cc,:) = scores(:,1);
cur_set = coeff(:,2);
F = TriScatteredInterp(interp_x,interp_y,cur_set);
pc245(cc,:,:) = F(Xi,Yi);
sc245(cc,:) = scores(:,2);

end
%% 
close all
for cc = 1:length(up_bar_ax)
    subplot(2,3,1)
    pcolor(xi,yi,squeeze(pc10(cc,:,:)));shading flat
     xl = xlim(); yl = ylim();
    line(xl,up_bar_ax0([cc cc])-0.04,'color','k')
    line(xl,up_bar_ax0([cc cc])+0.04,'color','k')
    subplot(2,3,4)
    pcolor(xi,yi,squeeze(pc20(cc,:,:)));shading flat
     xl = xlim(); yl = ylim();
    line(xl,up_bar_ax0([cc cc])-0.04,'color','k')
    line(xl,up_bar_ax0([cc cc])+0.04,'color','k')

    subplot(2,3,2)
    pcolor(xi,yi,squeeze(pc145(cc,:,:)));shading flat
     xl = xlim(); yl = ylim();
%     line(up_bar_ax90([cc cc])-0.04,yl,'color','k')
%     line(up_bar_ax90([cc cc])+0.04,yl,'color','k')
    subplot(2,3,5)
    pcolor(xi,yi,squeeze(pc245(cc,:,:)));shading flat
     xl = xlim(); yl = ylim();
%     line(up_bar_ax90([cc cc])-0.04,yl,'color','k')
%     line(up_bar_ax90([cc cc])+0.04,yl,'color','k')

    subplot(2,3,3)
    pcolor(xi,yi,squeeze(pc190(cc,:,:)));shading flat
     xl = xlim(); yl = ylim();
    line(up_bar_ax90([cc cc])-0.04,yl,'color','k')
    line(up_bar_ax90([cc cc])+0.04,yl,'color','k')
    subplot(2,3,6)
    pcolor(xi,yi,squeeze(pc290(cc,:,:)));shading flat
     xl = xlim(); yl = ylim();
    line(up_bar_ax90([cc cc])-0.04,yl,'color','k')
    line(up_bar_ax90([cc cc])+0.04,yl,'color','k')
pause
   clf
end

%%
% load ./cor_trig_avg_lfps_0deg_im.mat
% im_white_avgs = white_trig_avgs_cor;
% im_black_avgs = black_trig_avgs_cor;
% 
% load ./cor_trig_avg_lfps_0deg_gray.mat
% gray_white_avgs = white_trig_avgs_cor;
% gray_black_avgs = black_trig_avgs_cor;
% 
% %%
% close all
% for ll = 1:length(use_lfps)
%     subplot(2,2,1)
%     pcolor(lags/Fsd,1:up_n_bar_pos,squeeze(im_white_avgs(:,:,ll)));shading flat
%     ca = max(abs(caxis()));
%     caxis([-ca ca])
%     
%     subplot(2,2,2)
%     pcolor(lags/Fsd,1:up_n_bar_pos,squeeze(im_black_avgs(:,:,ll)));shading flat
%     caxis([-ca ca])
% 
%         subplot(2,2,3)
%     pcolor(lags/Fsd,1:up_n_bar_pos,squeeze(gray_white_avgs(:,:,ll)));shading flat
%     ca = max(abs(caxis()));
%     caxis([-ca ca])
%     
%     subplot(2,2,4)
%     pcolor(lags/Fsd,1:up_n_bar_pos,squeeze(gray_black_avgs(:,:,ll)));shading flat
%     caxis([-ca ca])
%     
%     ll
%     pause
%     clf
% end