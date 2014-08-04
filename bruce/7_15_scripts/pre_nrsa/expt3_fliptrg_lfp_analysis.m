%% Load Data
clear all;
cd ~/Data/bruce/7_15_12/G034/

obj_info_dir = '~/James_scripts/data_processing/Images/object_im_info/';
nat_set = 1:685;
white_set = 686:913;
sim_set = 914:1141;
obj_set = 1142:1369;
noise_set = [white_set sim_set];

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt3_fixbased_data.mat

Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
use_win = [-7 7;-7 7];

min_fix_dur = 0.15;
use_lfps = [1:2:96];

Fs = 3e4;
dsf = 100;Fsd = Fs/dsf;
[b,a] = butter(2,[1 100]/(Fsd/2));
% [b2,a2] = butter(2,[6 12]/(Fsd/2));
[b2,a2] = butter(2,[4 15]/(Fsd/2));
[b3,a3] = butter(2,[30 70]/(Fsd/2));
gam_sm_win = 10;

alpha_dp = 0.15;
spk_smwin = 6;

forwardlag = round(Fsd*0.6);
backlag = round(Fsd*0.4);
lags = -backlag:forwardlag;

forwardlag_p = round(5*2*pi/alpha_dp);
backlag_p = round(2*2*pi/alpha_dp);
lags_phase = -backlag_p:forwardlag_p;

x_pix = xax(xpatch_inds);
y_pix = yax(ypatch_inds);
[X_pix,Y_pix] = meshgrid(x_pix,y_pix);

% load ./expt2_lfp_amp_models_full_v2
% % fft_pred_out = fft_pred_out(1:2:end,:,:);
% use_gamma_freq = 20;
% gamma_mod_out = squeeze(fft_pred_out(:,use_gamma_freq,:));
% gamma_mod_out = zscore(gamma_mod_out')';
% gamma_prctile_mat = prctile(gamma_mod_out',[25 75]);

%%
n_trials = length(trial_start_inds);
trial_stop_sacinds = trial_stop_inds;
for i = 1:n_trials
    cur_inds = trial_start_inds(i):trial_stop_inds(i);
    cur_sac_stops = find(full_insac(cur_inds) == 1,1,'first');
    if ~isempty(cur_sac_stops)
        cur_sac_stops = cur_inds(cur_sac_stops);
        trial_stop_sacinds(i) = cur_sac_stops;
    end
end
trial_stop_inds = trial_stop_sacinds;
trial_start_times = full_t(trial_start_inds);
trial_stop_times = full_t(trial_stop_inds);
trial_durs = trial_stop_times - trial_start_times;

min_trial_dur = 0.4;
used_trials = find(trial_durs >= min_trial_dur);
trial_start_inds = trial_start_inds(used_trials);
trial_stop_inds = trial_stop_inds(used_trials);
trial_imnum = trial_imnum(used_trials);
trial_stimnum = trial_stimnum(used_trials);
trial_start_times = trial_start_times(used_trials);
trial_stop_times = trial_stop_times(used_trials);

% resh_all_stims = resh_all_stims(used_trials,:);
% resh_all_obj = resh_all_obj(used_trials,:);

trial_stop_winds = trial_start_inds + round(min_trial_dur/dt);
trial_stop_wtimes = full_t(trial_stop_winds);

n_trials = length(trial_start_inds);
trial_expt_num = nan(n_trials,1);
full_intrial = zeros(size(full_t));
for i = 1:n_trials
    cur_inds = trial_start_inds(i):trial_stop_winds(i);
    trial_expt_num(i) = unique(full_expt_vec(cur_inds));
end

%%
% load ./gabor_initfits_d1p25_varrate
load ./gabor_tracking_varmeans
gabor_params_used = gabor_params_f{2};
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

%fit smoothed retinotopic surface
orientations = linspace(0,pi-pi/12,12);
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

%% find which LFP RFs are in object for each object image
is_in_object = nan(n_trials,length(use_lfps));
obj_images = find(ismember(trial_imnum,obj_set));
obj_stim_ids = trial_stimnum(obj_images);

fov_inds = find(sqrt((Y_pix).^2+(X_pix).^2) < 0.1);
clear rf_inds
for i = 1:length(use_lfps)
    cur_inds = find(sqrt((Y_pix-interp_y(use_lfps(i))).^2 + (X_pix-interp_x(use_lfps(i))).^2) < 0.075);
    rf_inds{i}= cur_inds;
end


n_obj_stims = length(obj_stim_ids);
fov_obj_ids = nan(n_obj_stims,1); %0=backgrnd, 1=texture_patch, 2=object
rf_obj_ids = nan(length(use_lfps),n_obj_stims); %0=backgrnd, 1=texture_patch, 2=object

cd(obj_info_dir)
for i = 1:n_obj_stims
    cur_obj_id = find(obj_set == trial_imnum(obj_images(i)));
    cur_fname = sprintf('info4%.4d.mat',cur_obj_id);
    load(cur_fname);
    
    cur_fov_id = round(mean(resh_all_obj(obj_stim_ids(i),fov_inds),2));
    if cur_fov_id == 0
        fov_obj_ids(i) = 0;
    else
        if strcmp(item_info(cur_fov_id).item_type,'texture_patch')
            fov_obj_ids(i) = 1;
        else
            fov_obj_ids(i) = 2;
        end
    end
    
    for cc = 1:length(use_lfps)
        cur_rf_id = round(mean(resh_all_obj(obj_stim_ids(i),rf_inds{cc}),2));
        if cur_rf_id == 0
            rf_obj_ids(cc,i) = 0;
        else
            if strcmp(item_info(cur_rf_id).item_type,'texture_patch')
                rf_obj_ids(cc,i) = 1;
            else
                rf_obj_ids(cc,i) = 2;
            end
        end
    end
end
%%
all_trg_avg_alphaf = [];
all_trg_Vbb = [];
all_trg_Valpha = [];
all_trg_ValphaA = [];
all_trg_Vgamma = [];
all_trg_Valphaf = [];
all_trg_Vgammaf = [];
all_trg_spkrates = [];
all_trgp_Vbb = [];
all_trgp_Valpha = [];
all_trgp_Vgamma = [];
all_trgp_Valphaf = [];
all_trgp_Vgammaf = [];
all_trgp_spkrates = [];

cd ~/Data/bruce/7_15_12/G034/
% Expt_nu = [3 8 15 19 26 29];
% Expt_nu = [13 14 15 16 25 28 29];
Expt_nu = [7:12]; %expt 3 34
for ee = 1:length(Expt_nu)
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    
    filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(1));
    load(filename);
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    %%
    all_Vbb = nan(length(t_ax),length(use_lfps));
    all_Valpha = nan(length(t_ax),length(use_lfps));
    all_Vgamma = nan(length(t_ax),length(use_lfps));
    all_alphaphase = nan(length(t_ax),length(use_lfps));
    all_alphafreq = nan(length(t_ax),length(use_lfps));
    all_Valphaamp = nan(length(t_ax),length(use_lfps));
    all_gammafreq = nan(length(t_ax),length(use_lfps));
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        V_alpha = filtfilt(b2,a2,V);
        V_gamma = filtfilt(b3,a3,V);
        V_bb = filtfilt(b,a,V);
        all_Vbb(:,ll) = V_bb;
        all_Valpha(:,ll) = V_alpha;
%         all_Vgamma(:,ll) = abs(hilbert(V_gamma));
        all_Vgamma(:,ll) = smooth(abs(hilbert(V_gamma)),gam_sm_win);

                all_Valphaamp(:,ll) = abs(hilbert(V_alpha));

        alpha_phase = angle(hilbert(V_alpha));
        alpha_phase = unwrap_phase_monotonic(alpha_phase);
        all_alphaphase(:,ll) = alpha_phase;
        alpha_freq = smooth(alpha_phase,25,'lowess')';
        alpha_freq = [nan diff(alpha_freq)]/(2*pi)*Fsd;
        all_alphafreq(:,ll) = alpha_freq;
        gamma_phase = angle(hilbert(V_gamma));
        gamma_phase = unwrap_phase_monotonic(gamma_phase);
        %         gamma_freq = sgolayfilt(gamma_phase,sg_ord,sg_frame);
        gamma_freq = smooth(gamma_phase,25,'lowess')';
        gamma_freq = [nan diff(gamma_freq)]/(2*pi)*Fsd;
        all_gammafreq(:,ll) = gamma_freq;
    end
    avg_alpha = mean(all_Valpha,2);
    avg_alpha_phase = angle(hilbert(avg_alpha));
    avg_alpha_phase = unwrap_phase_monotonic(avg_alpha_phase);
    avg_alpha_phase = smooth(avg_alpha_phase,25,'lowess');
    avg_alpha_freq = [nan diff(avg_alpha_phase')]/(2*pi)*Fsd;
    
    uniform_phase_grid = avg_alpha_phase(1):alpha_dp:avg_alpha_phase(end);
    
    all_Vbb = zscore(all_Vbb);
    all_Valpha = zscore(all_Valpha);
    all_Valphaamp = zscore(all_Valphaamp);
    all_Vgamma = zscore(all_Vgamma);
    
    binned_spks = nan(length(t_ax),length(use_lfps));
    smoothed_spks = nan(length(t_ax),length(use_lfps));
    for cc = 1:length(use_lfps)
        binned_spks(:,cc) = hist(Clusters{use_lfps(cc)}.times,t_ax);
        smoothed_spks(:,cc) = smooth(binned_spks(:,cc),spk_smwin);
    end
    smoothed_spks = zscore(smoothed_spks);
    
    Vbb_interp = interp1(avg_alpha_phase,all_Vbb,uniform_phase_grid);
    Valpha_interp = interp1(avg_alpha_phase,all_Valpha,uniform_phase_grid);
    ValphaA_interp = interp1(avg_alpha_phase,all_Valphaamp,uniform_phase_grid);
    Vgamma_interp = interp1(avg_alpha_phase,all_Vgamma,uniform_phase_grid);
    Vgammaf_interp = interp1(avg_alpha_phase,all_gammafreq,uniform_phase_grid);
    Valphaf_interp = interp1(avg_alpha_phase,all_alphafreq,uniform_phase_grid);
    avg_alphaf_interp = interp1(avg_alpha_phase,avg_alpha_freq,uniform_phase_grid);
    spks_interp = interp1(avg_alpha_phase,smoothed_spks,uniform_phase_grid);
    t_ax_interp = interp1(avg_alpha_phase,t_ax,uniform_phase_grid);
    
    %%
    cur_trials = find(trial_expt_num==Expt_nu(ee));
    cur_trial_start_inds = round(interp1(t_ax,1:length(t_ax),trial_start_times(cur_trials)));
    cur_trial_stop_inds = round(interp1(t_ax,1:length(t_ax),trial_stop_times(cur_trials)));
    
    cur_trial_start_pinds = round(interp1(t_ax_interp,1:length(t_ax_interp),trial_start_times(cur_trials)));
    cur_trial_stop_pinds = round(interp1(t_ax_interp,1:length(t_ax_interp),trial_stop_times(cur_trials)));
    
    cur_n_trials = length(cur_trials);
    %%
    trg_Valpha = nan(cur_n_trials,length(lags),length(use_lfps));
    trg_ValphaA = nan(cur_n_trials,length(lags),length(use_lfps));
    trg_Vgamma = nan(cur_n_trials,length(lags),length(use_lfps));
    trg_Vbb = nan(cur_n_trials,length(lags),length(use_lfps));
    trg_alphaf = nan(cur_n_trials,length(lags),length(use_lfps));
    trg_gammaf = nan(cur_n_trials,length(lags),length(use_lfps));
    trg_spkrates = nan(cur_n_trials,length(lags),length(use_lfps));
    trg_avg_alphaf = nan(cur_n_trials,length(lags));
    for i = 1:cur_n_trials
        cur_inds = (cur_trial_start_inds(i)-backlag):(cur_trial_start_inds(i)+forwardlag);
%         cur_inds(cur_inds > cur_trial_stop_inds(i)) = [];
        cl = length(cur_inds);
        trg_Vbb(i,1:cl,:) = all_Vbb(cur_inds,:);
        trg_Valpha(i,1:cl,:) = all_Valpha(cur_inds,:);
        trg_ValphaA(i,1:cl,:) = all_Valphaamp(cur_inds,:);
        trg_Vgamma(i,1:cl,:) = all_Vgamma(cur_inds,:);
        trg_alphaf(i,1:cl,:) = all_alphafreq(cur_inds,:);
        trg_gammaf(i,1:cl,:) = all_gammafreq(cur_inds,:);
        trg_spkrates(i,1:cl,:) = smoothed_spks(cur_inds,:);
        trg_avg_alphaf(i,1:cl) = avg_alpha_freq(cur_inds);
    end
    trgp_Valpha = nan(cur_n_trials,length(lags_phase),length(use_lfps));
    trgp_Vgamma = nan(cur_n_trials,length(lags_phase),length(use_lfps));
    trgp_Vbb = nan(cur_n_trials,length(lags_phase),length(use_lfps));
    trgp_alphaf = nan(cur_n_trials,length(lags_phase),length(use_lfps));
    trgp_gammaf = nan(cur_n_trials,length(lags_phase),length(use_lfps));
    trgp_spkrates = nan(cur_n_trials,length(lags_phase),length(use_lfps));
    for i = 1:cur_n_trials
        cur_inds = (cur_trial_start_pinds(i)-backlag_p):(cur_trial_start_pinds(i)+forwardlag_p);
%         cur_inds(cur_inds > cur_trial_stop_pinds(i)) = [];
        cl = length(cur_inds);
        trgp_Vbb(i,1:cl,:) = Vbb_interp(cur_inds,:);
        trgp_Valpha(i,1:cl,:) = Valpha_interp(cur_inds,:);
        trgp_Vgamma(i,1:cl,:) = Vgamma_interp(cur_inds,:);
        trgp_alphaf(i,1:cl,:) = Valphaf_interp(cur_inds,:);
        trgp_gammaf(i,1:cl,:) = Vgammaf_interp(cur_inds,:);
        trgp_spkrates(i,1:cl,:) = spks_interp(cur_inds,:);
    end
    
    all_trg_Vbb = cat(1,all_trg_Vbb,trg_Vbb);
    all_trg_Valpha = cat(1,all_trg_Valpha,trg_Valpha);
    all_trg_ValphaA = cat(1,all_trg_ValphaA,trg_ValphaA);
    all_trg_Vgamma = cat(1,all_trg_Vgamma,trg_Vgamma);
    all_trg_Valphaf = cat(1,all_trg_Valphaf,trg_alphaf);
    all_trg_Vgammaf = cat(1,all_trg_Vgammaf,trg_gammaf);
    all_trg_spkrates = cat(1,all_trg_spkrates,trg_spkrates);
    all_trg_avg_alphaf = [all_trg_avg_alphaf; trg_avg_alphaf];
    
    all_trgp_Vbb = cat(1,all_trgp_Vbb,trgp_Vbb);
    all_trgp_Valpha = cat(1,all_trgp_Valpha,trgp_Valpha);
    all_trgp_Vgamma = cat(1,all_trgp_Vgamma,trgp_Vgamma);
    all_trgp_Valphaf = cat(1,all_trgp_Valphaf,trgp_alphaf);
    all_trgp_Vgammaf = cat(1,all_trgp_Vgammaf,trgp_gammaf);
    all_trgp_spkrates = cat(1,all_trgp_spkrates,trgp_spkrates);
    
end

%%
cd ~/Data/bruce/7_15_12/G034
save expt3_fliptrg_analysis_v2 all_* lags* fov_* rf_* 

%%
close all
n_trials = size(all_trgp_Vbb,1);
for cur_ch = 1:24;
    subplot(2,1,1)
shadedErrorBar(lags/Fsd,squeeze(mean(all_trg_Vbb(:,:,cur_ch))),squeeze(std(all_trg_Vbb(:,:,cur_ch)))/sqrt(n_trials))
xlim([-0.1 0.45])
    subplot(2,1,2)
shadedErrorBar(lags_phase/(2*pi)*alpha_dp,squeeze(mean(all_trgp_Vbb(:,:,cur_ch))),squeeze(std(all_trgp_Vbb(:,:,cur_ch)))/sqrt(n_trials))
xlim([-1 4])

pause
clf
end

%%
noise_trials = find(ismember(trial_imnum,noise_set));
nat_trials = find(ismember(trial_imnum,nat_set));
obj_trials = find(ismember(trial_imnum,obj_set));

fov_inobj_trials = obj_trials(fov_obj_ids ~= 0);
fov_outobj_trials = obj_trials(fov_obj_ids == 0);
for i = 1:length(use_lfps)
%     good_stims{i} = sac_inds(find(gamma_mod_out(i,cur_fixs(sac_inds)) >= gamma_prctile_mat(2,i)));
%     bad_stims{i} = sac_inds(find(gamma_mod_out(i,cur_fixs(sac_inds)) <= gamma_prctile_mat(1,i)));
    rf_infore_trials{i} = obj_trials(rf_obj_ids(i,:) ~= 0);
    rf_inback_trials{i} = obj_trials(rf_obj_ids(i,:) == 0);
    rf_inobj_trials{i} = obj_trials(rf_obj_ids(i,:) == 2);
    rf_inpatch_trials{i} = obj_trials(rf_obj_ids(i,:) == 1);
    
%    uset = randperm(length(rf_inback_trials{i}));
%    uset(length(rf_infore_trials{i})+1:end) = [];
%    rf_inback_trials{i} = rf_inback_trials{i}(uset);
end

%%
nboot = 100;
for i = 1:length(use_lfps)
boot_error_plot(squeeze(all_trg_spkrates(nat_trials,:,i)),lags/Fsd,nboot,'k')
hold on
boot_error_plot(squeeze(all_trg_spkrates(noise_trials,:,i)),lags/Fsd,nboot,'r')
boot_error_plot(squeeze(all_trg_spkrates(obj_trials,:,i)),lags/Fsd,nboot,'b')
pause
    clf
end
%%
nboot = 100;
for i = 1:length(use_lfps)
boot_error_plot(squeeze(all_trgp_Vgammaf(nat_trials,:,i)),lags_phase,nboot,'k')
hold on
boot_error_plot(squeeze(all_trgp_Vgammaf(noise_trials,:,i)),lags_phase,nboot,'r')
boot_error_plot(squeeze(all_trgp_Vgammaf(obj_trials,:,i)),lags_phase,nboot,'b')
pause
    clf
end

%%
nboot = 100;
for i = 1:length(use_lfps)
    % boot_error_plot(squeeze(all_trg_Vgamma(nat_trials,:,i)),lags/Fsd,nboot,'k')
    % hold on
    % boot_error_plot(squeeze(all_trg_Vgamma(noise_trials,:,i)),lags/Fsd,nboot,'r')
    % boot_error_plot(squeeze(all_trg_Vgamma(fov_inobj_trials,:,i)),lags/Fsd,nboot,'b')
    % boot_error_plot(squeeze(all_trg_Vgamma(fov_outobj_trials,:,i)),lags/Fsd,nboot,'g')
    
    boot_error_plot(squeeze(all_trg_spkrates(nat_trials,:,i)),lags/Fsd,nboot,'k')
    hold on
    boot_error_plot(squeeze(all_trg_spkrates(noise_trials,:,i)),lags/Fsd,nboot,'r')
    boot_error_plot(squeeze(all_trg_spkrates(fov_inobj_trials,:,i)),lags/Fsd,nboot,'b')
    boot_error_plot(squeeze(all_trg_spkrates(fov_outobj_trials,:,i)),lags/Fsd,nboot,'g')
    
    xlim([-0.15 0.45])
    pause
    clf
end
%%
nboot = 100;
for i = 1:length(use_lfps)
%     boot_error_plot(squeeze(all_trg_Valpha(nat_trials,:,i)),lags/Fsd,nboot,'k')
%     hold on
%     boot_error_plot(squeeze(all_trg_Valpha(noise_trials,:,i)),lags/Fsd,nboot,'r')
%     boot_error_plot(squeeze(all_trg_Valpha(fov_inobj_trials,:,i)),lags/Fsd,nboot,'b')
%     boot_error_plot(squeeze(all_trg_Valpha(fov_outobj_trials,:,i)),lags/Fsd,nboot,'g')
 
%     boot_error_plot(squeeze(all_trg_Vbb(nat_trials,:,i)),lags/Fsd,nboot,'k')
%     hold on
%     boot_error_plot(squeeze(all_trg_Vbb(noise_trials,:,i)),lags/Fsd,nboot,'r')
%     boot_error_plot(squeeze(all_trg_Vbb(rf_inobj_trials{i},:,i)),lags/Fsd,nboot,'b')
%     boot_error_plot(squeeze(all_trg_Vbb(rf_outobj_trials{i},:,i)),lags/Fsd,nboot,'g')
  
    boot_error_plot(squeeze(all_trg_Vgamma(:,:,i)),lags/Fsd,nboot,'k')
    hold on
    boot_error_plot(squeeze(all_trg_Vgamma(rf_inobj_trials{i},:,i)),lags/Fsd,nboot,'r')
    boot_error_plot(squeeze(all_trg_Vgamma(rf_infore_trials{i},:,i)),lags/Fsd,nboot,'b')
    boot_error_plot(squeeze(all_trg_Vgamma(rf_inback_trials{i},:,i)),lags/Fsd,nboot,'g')

% n_for = length(rf_infore_trials{i});
% n_back = length(rf_inback_trials{i});
% temp = randperm(n_back);
% temp(n_for+1:end) = [];
% 
%     boot_error_plot(squeeze(all_trgp_Vgamma(:,:,i)),lags_phase/(2*pi)*alpha_dp,nboot,'k')
%     hold on
%     boot_error_plot(squeeze(all_trgp_Vgammaf(noise_trials,:,i)),lags_phase/(2*pi)*alpha_dp,nboot,'r')
%     boot_error_plot(squeeze(all_trgp_Vgammaf(nat_trials,:,i)),lags_phase/(2*pi)*alpha_dp,nboot,'g')

%     boot_error_plot(squeeze(all_trgp_Vgammaf(:,:,i)),lags_phase/(2*pi)*alpha_dp,nboot,'k')
%     hold on
% %     boot_error_plot(squeeze(all_trgp_Valpha(rf_inobj_trials{i},:,i)),lags_phase/(2*pi)*alpha_dp,nboot,'r')
%     boot_error_plot(squeeze(all_trgp_Vgammaf(rf_infore_trials{i},:,i)),lags_phase/(2*pi)*alpha_dp,nboot,'b')
%     boot_error_plot(squeeze(all_trgp_Vgammaf(rf_inback_trials{i},:,i)),lags_phase/(2*pi)*alpha_dp,nboot,'g')
%     boot_error_plot(squeeze(all_trgp_Valpha(rf_inback_trials{i}(temp),:,i)),lags_phase/(2*pi)*alpha_dp,nboot,'g')
%     boot_error_plot(squeeze(all_trg_Valpha(rf_inpatch_trials{i},:,i)),lags/Fsd,nboot,'c')
%     boot_error_plot(squeeze(all_trg_spkrates(nat_trials,:,i)),lags/Fsd,nboot,'k')
% %     hold on
%     boot_error_plot(squeeze(all_trg_spkrates(noise_trials,:,i)),lags/Fsd,nboot,'r')
%     boot_error_plot(squeeze(all_trg_spkrates(fov_inobj_trials,:,i)),lags/Fsd,nboot,'b')
%     boot_error_plot(squeeze(all_trg_spkrates(fov_outobj_trials,:,i)),lags/Fsd,nboot,'g')
    
%     xlim([-0.15 0.45])
    pause
    clf
end

%%
use = find(diff(trial_imnum) ~= 0)+1;
use2 = find(trial_imnum(2:end)-trial_imnum(1:end-1) == 0)+1;
use2(~ismember(use2,use+1)) = [];
use3 = find(trial_imnum(2:end)-trial_imnum(1:end-1) == 0)+1;
use3(~ismember(use3,use+2)) = [];
use4 = find(trial_imnum(2:end)-trial_imnum(1:end-1) == 0)+1;
use4(~ismember(use4,use+3)) = [];
use5 = [use2; use3; use4];
clear cur_avg_*
for i = 1:length(use_lfps)
    cur_avg_nat(i,:) = squeeze(mean(all_trg_Vbb(nat_trials,:,i)));
    cur_avg_noise(i,:) = squeeze(mean(all_trg_Vbb(noise_trials,:,i)));
    cur_avg_sig(i,:) = squeeze(mean(all_trg_Vbb(:,:,i)));
    cur_avg_sig2(i,:) = squeeze(mean(all_trg_Vbb(use,:,i)));
    cur_avg_sig3(i,:) = squeeze(mean(all_trg_Vbb(use2,:,i)));
    cur_avg_sig4(i,:) = squeeze(mean(all_trg_Vbb(use3,:,i)));
    cur_avg_sig5(i,:) = squeeze(mean(all_trg_Vbb(use4,:,i)));
    cur_avg_sig6(i,:) = squeeze(mean(all_trg_Vbb(use5,:,i)));
    cur_avg_sig_fore(i,:) = squeeze(mean(all_trg_Vbb(rf_infore_trials{i},:,i)));
    cur_avg_sig_back(i,:) = squeeze(mean(all_trg_Vbb(rf_inback_trials{i},:,i)));
    cur_avg_sig_obj(i,:) = squeeze(mean(all_trg_Vbb(rf_inobj_trials{i},:,i)));
end

close all
figure
hold on
shadedErrorBar(lags/Fsd,mean(cur_avg_nat),std(cur_avg_nat)/sqrt(length(use_lfps)),{'color','k'})
shadedErrorBar(lags/Fsd,mean(cur_avg_noise),std(cur_avg_noise)/sqrt(length(use_lfps)),{'color','g'})
% shadedErrorBar(lags/Fsd,mean(cur_avg_sig2),std(cur_avg_sig2)/sqrt(length(use_lfps)),{'color','r'})
% shadedErrorBar(lags/Fsd,mean(cur_avg_sig3),std(cur_avg_sig3)/sqrt(length(use_lfps)),{'color','b'})
% shadedErrorBar(lags/Fsd,mean(cur_avg_sig4),std(cur_avg_sig4)/sqrt(length(use_lfps)),{'color','k'})
% shadedErrorBar(lags/Fsd,mean(cur_avg_sig5),std(cur_avg_sig5)/sqrt(length(use_lfps)),{'color','g'})
% shadedErrorBar(lags/Fsd,mean(cur_avg_sig6),std(cur_avg_sig6)/sqrt(length(use_lfps)),{'color','k'})
shadedErrorBar(lags/Fsd,mean(cur_avg_sig_fore),std(cur_avg_sig_fore)/sqrt(length(use_lfps)),{'color','r'})
shadedErrorBar(lags/Fsd,mean(cur_avg_sig_back),std(cur_avg_sig_back)/sqrt(length(use_lfps)),{'color','b'})
xlabel('Time (s)','fontsize',16)

%%
for i = 1:length(use_lfps)
    cur_avg_sig(i,:) = squeeze(mean(all_trg_Valpha(:,:,i)));
    cur_avg_sig2(i,:) = squeeze(mean(all_trg_Valpha(use,:,i)));
    cur_avg_sig3(i,:) = squeeze(mean(all_trg_Valpha(use2,:,i)));
    cur_avg_sig4(i,:) = squeeze(mean(all_trg_Valpha(use3,:,i)));
    cur_avg_sig5(i,:) = squeeze(mean(all_trg_Valpha(use4,:,i)));
    cur_avg_sig6(i,:) = squeeze(mean(all_trg_Valpha(use5,:,i)));
    cur_avg_sig_fore(i,:) = squeeze(mean(all_trg_Valpha(rf_infore_trials{i},:,i)));
    cur_avg_sig_back(i,:) = squeeze(mean(all_trg_Valpha(rf_inback_trials{i},:,i)));
    cur_avg_sig_obj(i,:) = squeeze(mean(all_trg_Valpha(rf_inobj_trials{i},:,i)));
end

% close all
figure
hold on
shadedErrorBar(lags/Fsd,mean(cur_avg_sig2),std(cur_avg_sig2)/sqrt(length(use_lfps)),{'color','r'})
% shadedErrorBar(lags/Fsd,mean(cur_avg_sig3),std(cur_avg_sig3)/sqrt(length(use_lfps)),{'color','b'})
% shadedErrorBar(lags/Fsd,mean(cur_avg_sig4),std(cur_avg_sig4)/sqrt(length(use_lfps)),{'color','k'})
% shadedErrorBar(lags/Fsd,mean(cur_avg_sig5),std(cur_avg_sig5)/sqrt(length(use_lfps)),{'color','g'})
shadedErrorBar(lags/Fsd,mean(cur_avg_sig6),std(cur_avg_sig6)/sqrt(length(use_lfps)),{'color','k'})
% shadedErrorBar(lags/Fsd,mean(cur_avg_sig_fore),std(cur_avg_sig_fore)/sqrt(length(use_lfps)),{'color','r'})
% shadedErrorBar(lags/Fsd,mean(cur_avg_sig_back),std(cur_avg_sig_back)/sqrt(length(use_lfps)),{'color','b'})
xlabel('Time (s)','fontsize',16)

%%
use = find(diff(trial_imnum) ~= 0)+1;
use2 = find(trial_imnum(2:end)-trial_imnum(1:end-1) == 0)+1;
use2(~ismember(use2,use+1)) = [];
use3 = find(trial_imnum(2:end)-trial_imnum(1:end-1) == 0)+1;
use3(~ismember(use3,use+2)) = [];
use4 = find(trial_imnum(2:end)-trial_imnum(1:end-1) == 0)+1;
use4(~ismember(use4,use+3)) = [];
use5 = [use2; use3; use4];
clear cur_avg_*
for i = 1:length(use_lfps)
    cur_avg_sig(i,:) = squeeze(nanmean(all_trgp_Vbb(:,:,i)));
    cur_avg_sig2(i,:) = squeeze(nanmean(all_trgp_Vbb(use,:,i)));
    cur_avg_sig3(i,:) = squeeze(nanmean(all_trgp_Vbb(use2,:,i)));
    cur_avg_sig4(i,:) = squeeze(nanmean(all_trgp_Vbb(use3,:,i)));
    cur_avg_sig5(i,:) = squeeze(nanmean(all_trgp_Vbb(use4,:,i)));
    cur_avg_sig6(i,:) = squeeze(nanmean(all_trgp_Vbb(use5,:,i)));
    cur_avg_sig_fore(i,:) = squeeze(nanmean(all_trgp_Vbb(rf_infore_trials{i},:,i)));
    cur_avg_sig_back(i,:) = squeeze(nanmean(all_trgp_Vbb(rf_inback_trials{i},:,i)));
    cur_avg_sig_obj(i,:) = squeeze(nanmean(all_trgp_Vbb(rf_inobj_trials{i},:,i)));
end

close all
figure
hold on
% shadedErrorBar(lags_phase/(2*pi)*alpha_dp,mean(cur_avg_sig2),std(cur_avg_sig2)/sqrt(length(use_lfps)),{'color','r'})
% shadedErrorBar(lags_phase/(2*pi)*alpha_dp,mean(cur_avg_sig3),std(cur_avg_sig3)/sqrt(length(use_lfps)),{'color','b'})
% shadedErrorBar(lags_phase/(2*pi)*alpha_dp,mean(cur_avg_sig4),std(cur_avg_sig4)/sqrt(length(use_lfps)),{'color','k'})
% shadedErrorBar(lags_phase/(2*pi)*alpha_dp,mean(cur_avg_sig5),std(cur_avg_sig5)/sqrt(length(use_lfps)),{'color','g'})
% shadedErrorBar(lags_phase/(2*pi)*alpha_dp,mean(cur_avg_sig6),std(cur_avg_sig6)/sqrt(length(use_lfps)),{'color','k'})
shadedErrorBar(lags_phase/(2*pi)*alpha_dp,nanmean(cur_avg_sig_fore),nanstd(cur_avg_sig_fore)/sqrt(length(use_lfps)),{'color','r'})
shadedErrorBar(lags_phase/(2*pi)*alpha_dp,nanmean(cur_avg_sig_back),nanstd(cur_avg_sig_back)/sqrt(length(use_lfps)),{'color','b'})
xlabel('Time (s)','fontsize',16)

%%
clear all
cd ~/Data/bruce/7_15_12/G034/
load ./sac_trg_avgs_v6 all_*Vbb lags Fsd all_*Vgamma
% fix_trg_Vbb = squeeze(nanmean(all_strg_Vbb,1));
fix_trg_Vbb = squeeze(nanmean(all_strg_Vgamma,1));
mac_trg_Vbb = squeeze(nanmean(all_mactrg_Vbb,1));
mic_trg_Vbb = squeeze(nanmean(all_mictrg_Vbb,1));
f_lags = lags/Fsd;

load ./expt3_fliptrg_analysis_v2 all_trg_Vbb all_trg_Vgamma
% stim_trg_Vbb = squeeze(nanmean(all_trg_Vbb));
stim_trg_Vbb = squeeze(nanmean(all_trg_Vgamma));
st_lags = lags/Fsd;

% cd ~/Data/bruce/7_15_12/G035/
% load ./expt4_flip_sac_trg_avgs avg_ftrg_Vbb lags Fsd
% flash_trg_Vbb = avg_ftrg_Vbb;
% fl_lags = lags/Fsd;

load ./expt1_trig_avgs.mat

figure
hold on
shadedErrorBar(f_lags,nanmean(fix_trg_Vbb,2),nanstd(fix_trg_Vbb,[],2))
% shadedErrorBar(f_lags,nanmean(mic_trg_Vbb,2),nanstd(mic_trg_Vbb,[],2),{'color','g'})
shadedErrorBar(f_lags,nanmean(stim_trg_Vbb,2),nanstd(stim_trg_Vbb,[],2),{'color','r'})
% shadedErrorBar(f_lags,nanmean(flash_trg_Vbb,2),nanstd(flash_trg_Vbb,[],2)/sqrt(48),{'color','b'})
% shadedErrorBar(expt1_lags,nanmean(all_mtrig_V,2),nanstd(all_mtrig_V,[],2),{'color','b'})
% shadedErrorBar(f_lags,nanmean(mac_trg_Vbb,2),nanstd(mac_trg_Vbb,[],2))

xlim([-0.1 0.5])

%%
for i = 1:48
hold on
plot(f_lags,mic_trg_Vbb(:,i),'g')
plot(f_lags,stim_trg_Vbb(:,i),'r')
plot(expt1_lags,all_mtrig_V(:,i),'b')
plot(f_lags,mac_trg_Vbb(:,i),'k')

xlim([-0.1 0.5])

pause
clf
end