%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12/G034/
obj_info_dir = '~/James_scripts/data_processing/Images/object_im_info';

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt2_compiled_windata_d1p5

Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
use_win = [-7 7;-7 7];
x_pix = xax(xpatch_inds);
y_pix = yax(ypatch_inds);
[X_pix,Y_pix] = meshgrid(x_pix,y_pix);

min_fix_dur = 0.15;
use_lfps = [1:96];
% use_lfps = [1:2:96];

Fs = 3e4;
dsf = 100;Fsd = Fs/dsf;
[b,a] = butter(2,[1 100]/(Fsd/2));
% [b2,a2] = butter(2,[6 12]/(Fsd/2));
[b2,a2] = butter(2,[5 14]/(Fsd/2));
[b3,a3] = butter(2,[35 70]/(Fsd/2));
sg_ord = 5;
sg_frame = 15;
gam_sm_win = 10;

alpha_dp = 0.15;
alpha_ch = 1;

forwardlag = round(Fsd*0.6);
backlag = round(Fsd*0.4);
lags = -backlag:forwardlag;

forwardlag_p = round(5*2*pi/alpha_dp);
backlag_p = round(2*2*pi/alpha_dp);
lags_phase = -backlag_p:forwardlag_p;

nat_set = 1:685;
white_set = 686:913;
sim_set = 914:1141;
obj_set = 1142:1369;
noise_set = [white_set sim_set];

load ./expt2_lfp_amp_models_full_v2
% fft_pred_out = fft_pred_out(1:2:end,:,:);
use_gamma_freq = 23;
gamma_mod_out = squeeze(fft_pred_out(:,use_gamma_freq,:));
gamma_mod_out = zscore(gamma_mod_out')';
gamma_prctile_mat = prctile(gamma_mod_out',[25 75]);

spk_smwin = round(0.015*Fsd);
%%
n_fixs = length(full_fix_wends);
fix_expt_num = nan(n_fixs,1);
fix_idnum = nan(n_fixs,1);
for i = 1:n_fixs
    cur_inds = full_fix_starts(i):full_fix_wends(i);
    fix_expt_num(i) = unique(full_expt_vec(cur_inds));
    fix_idnum(i) = unique(full_fix_ids(cur_inds));
end

%%
cd ~/Data/bruce/7_15_12/G034
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

% find which LFP RFs are in object for each object image
is_in_object = nan(n_fixs,length(use_lfps));
obj_images = find(ismember(fix_imnum,obj_set));
obj_images(fix_imnum(obj_images)==1167) = []; %don't have info for this image...
obj_stim_ids = fix_idnum(obj_images);
n_obj_stims = length(obj_stim_ids);

fov_inds = find(sqrt((Y_pix).^2+(X_pix).^2) < 0.1);
clear rf_inds
for i = 1:length(use_lfps)
    cur_inds = find(sqrt((Y_pix-interp_y(use_lfps(i))).^2 + (X_pix-interp_x(use_lfps(i))).^2) < 0.075);
    rf_inds{i}= cur_inds;
end


fov_obj_ids = nan(n_fixs,1); %0=backgrnd, 1=texture_patch, 2=object
rf_obj_ids = nan(length(use_lfps),n_fixs); %0=backgrnd, 1=texture_patch, 2=object

cd(obj_info_dir)
for i = 1:n_obj_stims
    cur_obj_id = find(obj_set == fix_imnum(obj_images(i)));
    cur_fname = sprintf('info4%.4d.mat',cur_obj_id);
    load(cur_fname);
    
    cur_fov_id = round(mean(resh_all_obj(obj_stim_ids(i),fov_inds),2));
    if cur_fov_id == 0
        fov_obj_ids(obj_stim_ids(i)) = 0;
    else
        if strcmp(item_info(cur_fov_id).item_type,'texture_patch')
            fov_obj_ids(obj_stim_ids(i)) = 1;
        else
            fov_obj_ids(obj_stim_ids(i)) = 2;
        end
    end
    
    for cc = 1:length(use_lfps)
        cur_rf_id = round(mean(resh_all_obj(obj_stim_ids(i),rf_inds{cc}),2));
        if cur_rf_id == 0
            rf_obj_ids(cc,obj_stim_ids(i)) = 0;
        else
            if strcmp(item_info(cur_rf_id).item_type,'texture_patch')
                rf_obj_ids(cc,obj_stim_ids(i)) = 1;
            else
                rf_obj_ids(cc,obj_stim_ids(i)) = 2;
            end
        end
    end
end

%%
 cd ~/Data/bruce/7_15_12/G034/
% Expt_nu = [3 8 15 19 26 29];
Expt_nu = [13 14 15 16 25 28 29];
n_allunits = 96;
all_all_strg_Vbb =[];
all_all_strgp_Vbb = [];
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
    all_ValphaA = nan(length(t_ax),length(use_lfps));
    all_Vgamma = nan(length(t_ax),length(use_lfps));
    all_alphaphase = nan(length(t_ax),length(use_lfps));
    all_alphafreq = nan(length(t_ax),length(use_lfps));
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
        all_ValphaA(:,ll) = abs(hilbert(V_alpha));
        all_Vgamma(:,ll) = abs(hilbert(V_gamma));
        all_Vgamma(:,ll) = smooth(all_Vgamma(:,ll),gam_sm_win);
        
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
    % avg_alpha = all_Valpha(:,alpha_ch);
    avg_alpha_phase = angle(hilbert(avg_alpha));
    avg_alpha_phase = unwrap_phase_monotonic(avg_alpha_phase);
    avg_alpha_phase = smooth(avg_alpha_phase,25,'lowess');
    
    uniform_phase_grid = avg_alpha_phase(1):alpha_dp:avg_alpha_phase(end);
    
    all_Vbb = zscore(all_Vbb);
    all_Valpha = zscore(all_Valpha);
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
    ValphaA_interp = interp1(avg_alpha_phase,all_ValphaA,uniform_phase_grid);
    Vgamma_interp = interp1(avg_alpha_phase,all_Vgamma,uniform_phase_grid);
    Vgammaf_interp = interp1(avg_alpha_phase,all_gammafreq,uniform_phase_grid);
    Valphaf_interp = interp1(avg_alpha_phase,all_alphafreq,uniform_phase_grid);
    spks_interp = interp1(avg_alpha_phase,smoothed_spks,uniform_phase_grid);
    t_ax_interp = interp1(avg_alpha_phase,t_ax,uniform_phase_grid);
    Valphaphase_interp = interp1(avg_alpha_phase,avg_alpha_phase,uniform_phase_grid);
    
%     sm_win = 5;
%     cur_binned_spks = nan(96,length(t_ax));
%     for j = 1:96
%         cur_binned_spks(j,:) = smooth(histc(Clusters{j}.times,t_ax),sm_win);
%     end

    %%
    cur_fixs = find(fix_expt_num == ee);
    cur_n_fixs = length(cur_fixs);
    
    fix_start_inds = round(interp1(t_ax,1:length(t_ax),full_t(full_fix_starts(cur_fixs))));
    fix_stop_inds = round(interp1(t_ax,1:length(t_ax),full_t(full_fix_ends(cur_fixs))));
    
    fix_start_pinds = round(interp1(t_ax_interp,1:length(t_ax_interp),full_t(full_fix_starts(cur_fixs))));
    fix_stop_pinds = round(interp1(t_ax_interp,1:length(t_ax_interp),full_t(full_fix_ends(cur_fixs))));
    
    big_sacs = find(fix_sac_amps(cur_fixs) > 1);
    micro_sacs = find(fix_sac_amps(cur_fixs) < 1);
    sac_inds = find(~isnan(fix_sac_amps(cur_fixs)));
    flip_inds = find(isnan(fix_sac_amps(cur_fixs)));
    noise_inds = sac_inds(ismember(fix_imnum(cur_fixs(sac_inds)),noise_set));
    wnoise_inds = sac_inds(ismember(fix_imnum(cur_fixs(sac_inds)),white_set));
    nat_inds = sac_inds(ismember(fix_imnum(cur_fixs(sac_inds)),nat_set));
    obj_inds = sac_inds(ismember(fix_imnum(cur_fixs(sac_inds)),obj_set));
    for cc =1:length(use_lfps)
        good_stims{cc} = sac_inds(find(gamma_mod_out(cc,cur_fixs(sac_inds)) >= gamma_prctile_mat(2,cc)));
        bad_stims{cc} = sac_inds(find(gamma_mod_out(cc,cur_fixs(sac_inds)) <= gamma_prctile_mat(1,cc)));
        rf_infor_trials{cc} = obj_inds(rf_obj_ids(cc,obj_inds) == 1 | rf_obj_ids(cc,obj_inds) == 2);
        rf_inback_trials{cc} = obj_inds(rf_obj_ids(cc,obj_inds) == 0);
        rf_inobj_trials{cc} = obj_inds(rf_obj_ids(cc,obj_inds) == 2);
        rf_inpatch_trials{cc} = obj_inds(rf_obj_ids(cc,obj_inds) == 1);
    end
    %%
    trg_Valpha = nan(cur_n_fixs,length(lags),length(use_lfps));
    trg_ValphaA = nan(cur_n_fixs,length(lags),length(use_lfps));
    trg_Vgamma = nan(cur_n_fixs,length(lags),length(use_lfps));
    trg_Vbb = nan(cur_n_fixs,length(lags),length(use_lfps));
    trg_alphaf = nan(cur_n_fixs,length(lags),length(use_lfps));
    trg_gammaf = nan(cur_n_fixs,length(lags),length(use_lfps));
    trg_spkrates = nan(cur_n_fixs,length(lags),length(use_lfps));
    for i = 1:cur_n_fixs
        cur_inds = (fix_start_inds(i)-backlag):(fix_start_inds(i)+forwardlag);
        cur_inds(cur_inds > fix_stop_inds(i)) = [];
        cl = length(cur_inds);
        trg_Vbb(i,1:cl,:) = all_Vbb(cur_inds,:);
        trg_Valpha(i,1:cl,:) = all_Valpha(cur_inds,:);
        trg_ValphaA(i,1:cl,:) = all_ValphaA(cur_inds,:);
        trg_Vgamma(i,1:cl,:) = all_Vgamma(cur_inds,:);
        trg_alphaf(i,1:cl,:) = all_alphafreq(cur_inds,:);
        trg_gammaf(i,1:cl,:) = all_gammafreq(cur_inds,:);
        trg_spkrates(i,1:cl,:) = smoothed_spks(cur_inds,:);
    end
    
    trgp_Valpha = nan(cur_n_fixs,length(lags_phase),length(use_lfps));
    trgp_Vgamma = nan(cur_n_fixs,length(lags_phase),length(use_lfps));
    trgp_Vbb = nan(cur_n_fixs,length(lags_phase),length(use_lfps));
    trgp_alphaf = nan(cur_n_fixs,length(lags_phase),length(use_lfps));
    trgp_gammaf = nan(cur_n_fixs,length(lags_phase),length(use_lfps));
    trgp_spkrates = nan(cur_n_fixs,length(lags_phase),length(use_lfps));
    for i = 1:cur_n_fixs
        cur_inds = (fix_start_pinds(i)-backlag_p):(fix_start_pinds(i)+forwardlag_p);
        cur_inds(cur_inds > fix_stop_pinds(i)) = [];
        cl = length(cur_inds);
        trgp_Vbb(i,1:cl,:) = Vbb_interp(cur_inds,:);
        trgp_Valpha(i,1:cl,:) = Valpha_interp(cur_inds,:);
        trgp_Vgamma(i,1:cl,:) = Vgamma_interp(cur_inds,:);
        trgp_alphaf(i,1:cl,:) = Valphaf_interp(cur_inds,:);
        trgp_gammaf(i,1:cl,:) = Vgammaf_interp(cur_inds,:);
        trgp_spkrates(i,1:cl,:) = spks_interp(cur_inds,:);
    end
    
    all_all_strg_Vbb = cat(1,all_all_strg_Vbb, trg_Vbb(sac_inds,:,:));
    all_all_strgp_Vbb = cat(1,all_all_strgp_Vbb, trgp_Vbb(sac_inds,:,:));
    
    all_ftrg_Vbb(ee,:,:) = nanmean(trg_Vbb(flip_inds,:,:));
    all_ftrg_Valpha(ee,:,:) = nanmean(trg_Valpha(flip_inds,:,:));
    all_ftrg_Vgamma(ee,:,:) = nanmean(trg_Vgamma(flip_inds,:,:));
    all_ftrg_Vgammaf(ee,:,:) = nanmean(trg_gammaf(flip_inds,:,:));
    all_ftrg_Valphaf(ee,:,:) = nanmean(trg_alphaf(flip_inds,:,:));
    all_ftrg_spk(ee,:,:) = nanmean(trg_spkrates(flip_inds,:,:));
    all_ftrgp_Vbb(ee,:,:) = nanmean(trgp_Vbb(flip_inds,:,:));
    all_ftrgp_Valpha(ee,:,:) = nanmean(trgp_Valpha(flip_inds,:,:));
    all_ftrgp_Vgamma(ee,:,:) = nanmean(trgp_Vgamma(flip_inds,:,:));
    all_ftrgp_Vgammaf(ee,:,:) = nanmean(trgp_gammaf(flip_inds,:,:));
    all_ftrgp_Valphaf(ee,:,:) = nanmean(trgp_alphaf(flip_inds,:,:));
    all_ftrgp_spk(ee,:,:) = nanmean(trgp_spkrates(flip_inds,:,:));
    trial_n_flips(ee) = length(flip_inds);
    
    all_strg_Vbb(ee,:,:) = nanmean(trg_Vbb(sac_inds,:,:));
    all_strg_Valpha(ee,:,:) = nanmean(trg_Valpha(sac_inds,:,:));
    all_strg_ValphaA(ee,:,:) = nanmean(trg_ValphaA(sac_inds,:,:));
    all_strg_Vgamma(ee,:,:) = nanmean(trg_Vgamma(sac_inds,:,:));
    all_strg_Vgammaf(ee,:,:) = nanmean(trg_gammaf(sac_inds,:,:));
    all_strg_Valphaf(ee,:,:) = nanmean(trg_alphaf(sac_inds,:,:));
    all_strg_spk(ee,:,:) = nanmean(trg_spkrates(sac_inds,:,:));
    all_strgp_Vbb(ee,:,:) = nanmean(trgp_Vbb(sac_inds,:,:));
    all_strgp_Valpha(ee,:,:) = nanmean(trgp_Valpha(sac_inds,:,:));
    all_strgp_Vgamma(ee,:,:) = nanmean(trgp_Vgamma(sac_inds,:,:));
    all_strgp_Vgammaf(ee,:,:) = nanmean(trgp_gammaf(sac_inds,:,:));
    all_strgp_Valphaf(ee,:,:) = nanmean(trgp_alphaf(sac_inds,:,:));
    all_strgp_spk(ee,:,:) = nanmean(trgp_spkrates(sac_inds,:,:));
    trial_n_sacs(ee) = length(sac_inds);
   
    all_mactrg_Vbb(ee,:,:) = nanmean(trg_Vbb(big_sacs,:,:));
    all_mactrg_Valpha(ee,:,:) = nanmean(trg_Valpha(big_sacs,:,:));
    all_mactrg_ValphaA(ee,:,:) = nanmean(trg_ValphaA(big_sacs,:,:));
    all_mactrg_Vgamma(ee,:,:) = nanmean(trg_Vgamma(big_sacs,:,:));
    all_mactrg_Vgammaf(ee,:,:) = nanmean(trg_gammaf(big_sacs,:,:));
    all_mactrg_Valphaf(ee,:,:) = nanmean(trg_alphaf(big_sacs,:,:));
    all_mactrg_spk(ee,:,:) = nanmean(trg_spkrates(big_sacs,:,:));
    all_mactrgp_Vbb(ee,:,:) = nanmean(trgp_Vbb(big_sacs,:,:));
    all_mactrgp_Valpha(ee,:,:) = nanmean(trgp_Valpha(big_sacs,:,:));
    all_mactrgp_Vgamma(ee,:,:) = nanmean(trgp_Vgamma(big_sacs,:,:));
    all_mactrgp_Vgammaf(ee,:,:) = nanmean(trgp_gammaf(big_sacs,:,:));
    all_mactrgp_Valphaf(ee,:,:) = nanmean(trgp_alphaf(big_sacs,:,:));
    all_mactrgp_spk(ee,:,:) = nanmean(trgp_spkrates(big_sacs,:,:));
        trial_n_macsacs(ee) = length(big_sacs);
        
    all_mictrg_Vbb(ee,:,:) = nanmean(trg_Vbb(micro_sacs,:,:));
    all_mictrg_Valpha(ee,:,:) = nanmean(trg_Valpha(micro_sacs,:,:));
    all_mictrg_ValphaA(ee,:,:) = nanmean(trg_ValphaA(micro_sacs,:,:));
    all_mictrg_Vgamma(ee,:,:) = nanmean(trg_Vgamma(micro_sacs,:,:));
    all_mictrg_Vgammaf(ee,:,:) = nanmean(trg_gammaf(micro_sacs,:,:));
    all_mictrg_Valphaf(ee,:,:) = nanmean(trg_alphaf(micro_sacs,:,:));
    all_mictrg_spk(ee,:,:) = nanmean(trg_spkrates(micro_sacs,:,:));
    all_mictrgp_Vbb(ee,:,:) = nanmean(trgp_Vbb(micro_sacs,:,:));
    all_mictrgp_Valpha(ee,:,:) = nanmean(trgp_Valpha(micro_sacs,:,:));
    all_mictrgp_Vgamma(ee,:,:) = nanmean(trgp_Vgamma(micro_sacs,:,:));
    all_mictrgp_Vgammaf(ee,:,:) = nanmean(trgp_gammaf(micro_sacs,:,:));
    all_mictrgp_Valphaf(ee,:,:) = nanmean(trgp_alphaf(micro_sacs,:,:));
    all_mictrgp_spk(ee,:,:) = nanmean(trgp_spkrates(micro_sacs,:,:));
            trial_n_micsacs(ee) = length(micro_sacs);

    all_nattrg_Vbb(ee,:,:) = nanmean(trg_Vbb(nat_inds,:,:));
    all_nattrg_Valpha(ee,:,:) = nanmean(trg_Valpha(nat_inds,:,:));
    all_nattrg_ValphaA(ee,:,:) = nanmean(trg_ValphaA(nat_inds,:,:));
    all_nattrg_Vgamma(ee,:,:) = nanmean(trg_Vgamma(nat_inds,:,:));
    all_nattrg_Vgammaf(ee,:,:) = nanmean(trg_gammaf(nat_inds,:,:));
    all_nattrg_Valphaf(ee,:,:) = nanmean(trg_alphaf(nat_inds,:,:));
    all_nattrg_spk(ee,:,:) = nanmean(trg_spkrates(nat_inds,:,:));
    all_nattrgp_Vbb(ee,:,:) = nanmean(trgp_Vbb(nat_inds,:,:));
    all_nattrgp_Valpha(ee,:,:) = nanmean(trgp_Valpha(nat_inds,:,:));
    all_nattrgp_Vgamma(ee,:,:) = nanmean(trgp_Vgamma(nat_inds,:,:));
    all_nattrgp_Vgammaf(ee,:,:) = nanmean(trgp_gammaf(nat_inds,:,:));
    all_nattrgp_Valphaf(ee,:,:) = nanmean(trgp_alphaf(nat_inds,:,:));
    all_nattrgp_spk(ee,:,:) = nanmean(trgp_spkrates(nat_inds,:,:));
             trial_n_nat(ee) = length(nat_inds);
   
    all_noisetrg_Vbb(ee,:,:) = nanmean(trg_Vbb(noise_inds,:,:));
    all_noisetrg_Valpha(ee,:,:) = nanmean(trg_Valpha(noise_inds,:,:));
    all_noisetrg_ValphaA(ee,:,:) = nanmean(trg_ValphaA(noise_inds,:,:));
    all_noisetrg_Vgamma(ee,:,:) = nanmean(trg_Vgamma(noise_inds,:,:));
    all_noisetrg_Vgammaf(ee,:,:) = nanmean(trg_gammaf(noise_inds,:,:));
    all_noisetrg_Valphaf(ee,:,:) = nanmean(trg_alphaf(noise_inds,:,:));
    all_noisetrg_spk(ee,:,:) = nanmean(trg_spkrates(noise_inds,:,:));
    all_noisetrgp_Vbb(ee,:,:) = nanmean(trgp_Vbb(noise_inds,:,:));
    all_noisetrgp_Valpha(ee,:,:) = nanmean(trgp_Valpha(noise_inds,:,:));
    all_noisetrgp_Vgamma(ee,:,:) = nanmean(trgp_Vgamma(noise_inds,:,:));
    all_noisetrgp_Vgammaf(ee,:,:) = nanmean(trgp_gammaf(noise_inds,:,:));
    all_noisetrgp_Valphaf(ee,:,:) = nanmean(trgp_alphaf(noise_inds,:,:));
     all_noisetrgp_spk(ee,:,:) = nanmean(trgp_spkrates(noise_inds,:,:));
             trial_n_noise(ee) = length(noise_inds);
     
     all_objtrg_Vbb(ee,:,:) = nanmean(trg_Vbb(obj_inds,:,:));
     all_objtrg_Valpha(ee,:,:) = nanmean(trg_Valpha(obj_inds,:,:));
     all_objtrg_ValphaA(ee,:,:) = nanmean(trg_ValphaA(obj_inds,:,:));
     all_objtrg_Vgamma(ee,:,:) = nanmean(trg_Vgamma(obj_inds,:,:));
     all_objtrg_Vgammaf(ee,:,:) = nanmean(trg_gammaf(obj_inds,:,:));
     all_objtrg_Valphaf(ee,:,:) = nanmean(trg_alphaf(obj_inds,:,:));
     all_objtrg_spk(ee,:,:) = nanmean(trg_spkrates(obj_inds,:,:));
     all_objtrgp_Vbb(ee,:,:) = nanmean(trgp_Vbb(obj_inds,:,:));
     all_objtrgp_Valpha(ee,:,:) = nanmean(trgp_Valpha(obj_inds,:,:));
     all_objtrgp_Vgamma(ee,:,:) = nanmean(trgp_Vgamma(obj_inds,:,:));
     all_objtrgp_Vgammaf(ee,:,:) = nanmean(trgp_gammaf(obj_inds,:,:));
     all_objtrgp_Valphaf(ee,:,:) = nanmean(trgp_alphaf(obj_inds,:,:));
     all_objtrgp_spk(ee,:,:) = nanmean(trgp_spkrates(obj_inds,:,:));
             trial_n_obj(ee) = length(obj_inds);
     
    for cc = 1:length(use_lfps)
        all_goodtrg_Vbb(ee,:,cc) = nanmean(trg_Vbb(good_stims{cc},:,cc));
        all_goodtrg_Valpha(ee,:,cc) = nanmean(trg_Valpha(good_stims{cc},:,cc));
        all_goodtrg_ValphaA(ee,:,cc) = nanmean(trg_ValphaA(good_stims{cc},:,cc));
        all_goodtrg_Vgamma(ee,:,cc) = nanmean(trg_Vgamma(good_stims{cc},:,cc));
        all_goodtrg_Vgammaf(ee,:,cc) = nanmean(trg_gammaf(good_stims{cc},:,cc));
        all_goodtrg_Valphaf(ee,:,cc) = nanmean(trg_alphaf(good_stims{cc},:,cc));
        all_goodtrg_spk(ee,:,cc) = nanmean(trg_spkrates(good_stims{cc},:,cc));
        all_goodtrgp_Vbb(ee,:,cc) = nanmean(trgp_Vbb(good_stims{cc},:,cc));
        all_goodtrgp_Valpha(ee,:,cc) = nanmean(trgp_Valpha(good_stims{cc},:,cc));
        all_goodtrgp_Vgamma(ee,:,cc) = nanmean(trgp_Vgamma(good_stims{cc},:,cc));
        all_goodtrgp_Vgammaf(ee,:,cc) = nanmean(trgp_gammaf(good_stims{cc},:,cc));
        all_goodtrgp_Valphaf(ee,:,cc) = nanmean(trgp_alphaf(good_stims{cc},:,cc));
        all_goodtrgp_spk(ee,:,cc) = nanmean(trgp_spkrates(good_stims{cc},:,cc));
              trial_n_good(ee,cc) = length(good_stims{cc});
       
        all_badtrg_Vbb(ee,:,cc) = nanmean(trg_Vbb(bad_stims{cc},:,cc));
        all_badtrg_Valpha(ee,:,cc) = nanmean(trg_Valpha(bad_stims{cc},:,cc));
        all_badtrg_ValphaA(ee,:,cc) = nanmean(trg_ValphaA(bad_stims{cc},:,cc));
        all_badtrg_Vgamma(ee,:,cc) = nanmean(trg_Vgamma(bad_stims{cc},:,cc));
        all_badtrg_Vgammaf(ee,:,cc) = nanmean(trg_gammaf(bad_stims{cc},:,cc));
        all_badtrg_Valphaf(ee,:,cc) = nanmean(trg_alphaf(bad_stims{cc},:,cc));
        all_badtrg_spk(ee,:,cc) = nanmean(trg_spkrates(bad_stims{cc},:,cc));
        all_badtrgp_Vbb(ee,:,cc) = nanmean(trgp_Vbb(bad_stims{cc},:,cc));
        all_badtrgp_Valpha(ee,:,cc) = nanmean(trgp_Valpha(bad_stims{cc},:,cc));
        all_badtrgp_Vgamma(ee,:,cc) = nanmean(trgp_Vgamma(bad_stims{cc},:,cc));
        all_badtrgp_Vgammaf(ee,:,cc) = nanmean(trgp_gammaf(bad_stims{cc},:,cc));
        all_badtrgp_Valphaf(ee,:,cc) = nanmean(trgp_alphaf(bad_stims{cc},:,cc));
        all_badtrgp_spk(ee,:,cc) = nanmean(trgp_spkrates(bad_stims{cc},:,cc));
        trial_n_bad(ee,cc) = length(bad_stims{cc});
        
        all_fortrg_Vbb(ee,:,cc) = nanmean(trg_Vbb(rf_infor_trials{cc},:,cc));
        all_fortrg_Valpha(ee,:,cc) = nanmean(trg_Valpha(rf_infor_trials{cc},:,cc));
        all_fortrg_ValphaA(ee,:,cc) = nanmean(trg_ValphaA(rf_infor_trials{cc},:,cc));
        all_fortrg_Vgamma(ee,:,cc) = nanmean(trg_Vgamma(rf_infor_trials{cc},:,cc));
        all_fortrg_Vgammaf(ee,:,cc) = nanmean(trg_gammaf(rf_infor_trials{cc},:,cc));
        all_fortrg_Valphaf(ee,:,cc) = nanmean(trg_alphaf(rf_infor_trials{cc},:,cc));
        all_fortrg_spk(ee,:,cc) = nanmean(trg_spkrates(rf_infor_trials{cc},:,cc));
        all_fortrgp_Vbb(ee,:,cc) = nanmean(trgp_Vbb(rf_infor_trials{cc},:,cc));
        all_fortrgp_Valpha(ee,:,cc) = nanmean(trgp_Valpha(rf_infor_trials{cc},:,cc));
        all_fortrgp_Vgamma(ee,:,cc) = nanmean(trgp_Vgamma(rf_infor_trials{cc},:,cc));
        all_fortrgp_Vgammaf(ee,:,cc) = nanmean(trgp_gammaf(rf_infor_trials{cc},:,cc));
        all_fortrgp_Valphaf(ee,:,cc) = nanmean(trgp_alphaf(rf_infor_trials{cc},:,cc));
        all_fortrgp_spk(ee,:,cc) = nanmean(trgp_spkrates(rf_infor_trials{cc},:,cc));
        trial_n_for(ee,cc) = length(rf_infor_trials{cc});

            all_backtrg_Vbb(ee,:,cc) = nanmean(trg_Vbb(rf_inback_trials{cc},:,cc));
        all_backtrg_Valpha(ee,:,cc) = nanmean(trg_Valpha(rf_inback_trials{cc},:,cc));
        all_backtrg_ValphaA(ee,:,cc) = nanmean(trg_ValphaA(rf_inback_trials{cc},:,cc));
        all_backtrg_Vgamma(ee,:,cc) = nanmean(trg_Vgamma(rf_inback_trials{cc},:,cc));
        all_backtrg_Vgammaf(ee,:,cc) = nanmean(trg_gammaf(rf_inback_trials{cc},:,cc));
        all_backtrg_Valphaf(ee,:,cc) = nanmean(trg_alphaf(rf_inback_trials{cc},:,cc));
        all_backtrg_spk(ee,:,cc) = nanmean(trg_spkrates(rf_inback_trials{cc},:,cc));
        all_backtrgp_Vbb(ee,:,cc) = nanmean(trgp_Vbb(rf_inback_trials{cc},:,cc));
        all_backtrgp_Valpha(ee,:,cc) = nanmean(trgp_Valpha(rf_inback_trials{cc},:,cc));
        all_backtrgp_Vgamma(ee,:,cc) = nanmean(trgp_Vgamma(rf_inback_trials{cc},:,cc));
        all_backtrgp_Vgammaf(ee,:,cc) = nanmean(trgp_gammaf(rf_inback_trials{cc},:,cc));
        all_backtrgp_Valphaf(ee,:,cc) = nanmean(trgp_alphaf(rf_inback_trials{cc},:,cc));
        all_backtrgp_spk(ee,:,cc) = nanmean(trgp_spkrates(rf_inback_trials{cc},:,cc));
        trial_n_back(ee,cc) = length(rf_inback_trials{cc});
        
        all_inobjtrg_Vbb(ee,:,cc) = nanmean(trg_Vbb(rf_inobj_trials{cc},:,cc));
        all_inobjtrg_Valpha(ee,:,cc) = nanmean(trg_Valpha(rf_inobj_trials{cc},:,cc));
        all_inobjtrg_ValphaA(ee,:,cc) = nanmean(trg_ValphaA(rf_inobj_trials{cc},:,cc));
        all_inobjtrg_Vgamma(ee,:,cc) = nanmean(trg_Vgamma(rf_inobj_trials{cc},:,cc));
        all_inobjtrg_Vgammaf(ee,:,cc) = nanmean(trg_gammaf(rf_inobj_trials{cc},:,cc));
        all_inobjtrg_Valphaf(ee,:,cc) = nanmean(trg_alphaf(rf_inobj_trials{cc},:,cc));
        all_inobjtrg_spk(ee,:,cc) = nanmean(trg_spkrates(rf_inobj_trials{cc},:,cc));
        all_inobjtrgp_Vbb(ee,:,cc) = nanmean(trgp_Vbb(rf_inobj_trials{cc},:,cc));
        all_inobjtrgp_Valpha(ee,:,cc) = nanmean(trgp_Valpha(rf_inobj_trials{cc},:,cc));
        all_inobjtrgp_Vgamma(ee,:,cc) = nanmean(trgp_Vgamma(rf_inobj_trials{cc},:,cc));
        all_inobjtrgp_Vgammaf(ee,:,cc) = nanmean(trgp_gammaf(rf_inobj_trials{cc},:,cc));
        all_inobjtrgp_Valphaf(ee,:,cc) = nanmean(trgp_alphaf(rf_inobj_trials{cc},:,cc));
        all_inobjtrgp_spk(ee,:,cc) = nanmean(trgp_spkrates(rf_inobj_trials{cc},:,cc));
        trial_n_inobj(ee,cc) = length(rf_inobj_trials{cc});
    end
    %%
    
end

%%
save sac_trg_avgs_v6 lags Fsd all_* lags_phase alpha_dp trial*

%%
% close all
% % load ./expt3_fliptrg_analysis_v2
% for i = 1:24
%     subplot(2,1,1)
%     matrix_errorbar_plot(squeeze(all_all_strg_Vbb(:,:,(i-1)*2+1)),lags/Fsd,'k')
%     hold on
% %     matrix_errorbar_plot(squeeze(all_trg_Vbb(:,:,i)),lags/Fsd,'r')
%     xlim([-0.1 0.4])
%     subplot(2,1,2)
%     matrix_errorbar_plot(squeeze(all_all_strgp_Vbb(:,:,(i-1)*2+1)),lags_phase/(2*pi)*alpha_dp,'k')
%     hold on
%     matrix_errorbar_plot(squeeze(all_trgp_Vbb(:,:,i)),lags_phase/(2*pi)*alpha_dp,'r')
%     xlim([-1 4])    
%     pause
%     clf
% end
% 
%%
% load ./expt3_fliptrg_analysis all_trg_Vbb
% all_ftrg_Vbb = all_trg_Vbb;
% load ./sac_trg_avgs_v3 all_strg_Vbb
close all
for cc = 1 :length(use_lfps)
    cc
        matrix_errorbar_plot(squeeze(all_strg_Vbb(:,:,cc)),lags/Fsd,'k');
        hold on
%         matrix_errorbar_plot(squeeze(all_ftrg_Vbb(:,:,cc)),lags/Fsd,'b');
%         matrix_errorbar_plot(squeeze(all_strg_Vbb(:,:,cc)),lags/Fsd,'k');
%         hold on
%         matrix_errorbar_plot(squeeze(all_strg_Valpha(:,:,cc)),lags/Fsd,'b');
%         matrix_errorbar_plot(squeeze(all_strg_Vgamma(:,:,cc)),lags/Fsd,'r');
        
        matrix_errorbar_plot(squeeze(all_mactrg_Vbb(:,:,cc)),lags/Fsd,'r');
        matrix_errorbar_plot(squeeze(all_mictrg_Vbb(:,:,cc)),lags/Fsd,'b');
%         hold on
%         matrix_errorbar_plot(squeeze(all_mactrg_Valpha(:,:,cc)),lags/Fsd,'b');
%         matrix_errorbar_plot(squeeze(all_mactrg_Vgamma(:,:,cc)),lags/Fsd,'r');
%         matrix_errorbar_plot(squeeze(all_mactrg_spk(:,:,cc)),lags/Fsd,'g');
%         matrix_errorbar_plot(squeeze(all_strgp_Vbb(:,:,cc)),lags_phase/(2*pi),'k');
%         hold on
%         matrix_errorbar_plot(squeeze(all_ftrgp_Vbb(:,:,cc)),lags_phase/(2*pi),'b');
%         matrix_errorbar_plot(squeeze(all_goodtrg_Vbb(:,:,cc)),lags/Fsd,'b');
%         matrix_errorbar_plot(squeeze(all_badtrg_Vbb(:,:,cc)),lags/Fsd,'r');
    
%     matrix_errorbar_plot(squeeze(all_strg_Valpha(:,:,cc)),lags/Fsd,'k');

% subplot(2,1,1)
%     hold on
%     matrix_errorbar_plot(squeeze(all_goodtrg_Vbb(:,:,cc)),lags/Fsd,'b');
%     matrix_errorbar_plot(squeeze(all_badtrg_Vbb(:,:,cc)),lags/Fsd,'r');
% subplot(2,1,2)
% ndims = length(fft_beta{cc,use_gamma_freq});
% pcolor(reshape(fft_beta{cc,use_gamma_freq},sqrt(ndims),sqrt(ndims)));shading flat

%     matrix_errorbar_plot(squeeze(all_goodtrg_V(:,:,cc)-all_badtrg_Vgamma(:,:,cc)),lags/Fsd,'r');
    
%     matrix_errorbar_plot(squeeze(all_strg_Valphaf(:,:,cc)),lags/Fsd,'k');
%     hold on
%     matrix_errorbar_plot(squeeze(all_goodtrg_Valphaf(:,:,cc)),lags/Fsd,'b');
%     matrix_errorbar_plot(squeeze(all_badtrg_Valphaf(:,:,cc)),lags/Fsd,'r');

%     matrix_errorbar_plot(squeeze(all_strg_Valphaf(:,:,cc)),lags/Fsd,'k');
%     hold on
%     matrix_errorbar_plot(squeeze(all_noisetrg_Valphaf(:,:,cc)),lags/Fsd,'b');
%     matrix_errorbar_plot(squeeze(all_nattrg_Valphaf(:,:,cc)),lags/Fsd,'r');

%     matrix_errorbar_plot(squeeze(all_strg_Vgammaf(:,:,cc)),lags/Fsd,'k');
%     hold on
%     matrix_errorbar_plot(squeeze(all_noisetrg_Valpha(:,:,cc)),lags/Fsd,'b');
%     matrix_errorbar_plot(squeeze(all_nattrg_Valpha(:,:,cc)),lags/Fsd,'r');
%     matrix_errorbar_plot(squeeze(all_fortrg_Valpha(:,:,cc)),lags/Fsd,'g');
%     matrix_errorbar_plot(squeeze(all_backtrg_Valpha(:,:,cc)),lags/Fsd,'k');
% % 
%   
%     plot(lags/Fsd,squeeze(nansum(bsxfun(@times,all_strg_Valpha(:,:,cc),trial_n_sacs')))/sum(trial_n_sacs),'k')
%     hold on
%     plot(lags/Fsd,squeeze(nansum(bsxfun(@times,all_backtrg_Valpha(:,:,cc),trial_n_back(:,cc))))/sum(trial_n_back(:,cc)),'g')
%     plot(lags/Fsd,squeeze(nansum(bsxfun(@times,all_fortrg_Valpha(:,:,cc),trial_n_for(:,cc))))/sum(trial_n_for(:,cc)),'b')
%     plot(lags_phase,squeeze(nansum(bsxfun(@times,all_strgp_Valpha(:,:,cc),trial_n_sacs')))/sum(trial_n_sacs),'k')
%     hold on
%     plot(lags_phase,squeeze(nansum(bsxfun(@times,all_backtrgp_Valpha(:,:,cc),trial_n_back(:,cc))))/sum(trial_n_back(:,cc)),'g')
%     plot(lags_phase,squeeze(nansum(bsxfun(@times,all_fortrgp_Valpha(:,:,cc),trial_n_for(:,cc))))/sum(trial_n_for(:,cc)),'b')
%     plot(lags_phase,squeeze(nansum(bsxfun(@times,all_inobjtrgp_Valpha(:,:,cc),trial_n_inobj(:,cc))))/sum(trial_n_inobj(:,cc)),'r')

    %         matrix_errorbar_plot(squeeze(all_strg_Valpha(:,:,cc)),lags/Fsd,'k');
%     hold on
% %     matrix_errorbar_plot(squeeze(all_inobjtrg_Valpha(:,:,cc)),lags/Fsd,'r');
%     matrix_errorbar_plot(squeeze(all_backtrg_Valpha(:,:,cc)),lags/Fsd,'g');
%     matrix_errorbar_plot(squeeze(all_fortrg_Valpha(:,:,cc)),lags/Fsd,'b');
% 
    %     matrix_errorbar_plot(squeeze(all_strg_spk(:,:,cc)),lags/Fsd,'k');
%     hold on
%     matrix_errorbar_plot(squeeze(all_mactrg_spk(:,:,cc)),lags/Fsd,'b');
%     matrix_errorbar_plot(squeeze(all_mictrg_spk(:,:,cc)),lags/Fsd,'r');
%         matrix_errorbar_plot(squeeze(all_strg_Vgammaf(:,:,cc)),lags/Fsd,'k');
%     hold on
%     matrix_errorbar_plot(squeeze(all_goodtrg_Vbb(:,:,cc)),lags/Fsd,'b');
%     matrix_errorbar_plot(squeeze(all_badtrg_Vbb(:,:,cc)),lags/Fsd,'r');

%     matrix_errorbar_plot(squeeze(all_strg_Vgammaf(:,:,cc)),lags/Fsd,'k');
%     hold on
%     matrix_errorbar_plot(squeeze(all_mactrg_Vgammaf(:,:,cc)),lags/Fsd,'b');
%     matrix_errorbar_plot(squeeze(all_mictrg_Vgammaf(:,:,cc)),lags/Fsd,'r');

%     matrix_errorbar_plot(squeeze(all_strg_Valphaf(:,:,cc)),lags/Fsd,'k');
%     hold on
%     matrix_errorbar_plot(squeeze(all_mactrg_Valphaf(:,:,cc)),lags/Fsd,'b');
%     matrix_errorbar_plot(squeeze(all_mictrg_Valphaf(:,:,cc)),lags/Fsd,'r');
% 
%     matrix_errorbar_plot(squeeze(all_strg_Vbb(:,:,cc)),lags/Fsd,'k');
%     hold on
%     matrix_errorbar_plot(squeeze(all_mactrg_Vbb(:,:,cc)),lags/Fsd,'b');
%     matrix_errorbar_plot(squeeze(all_mictrg_Vbb(:,:,cc)),lags/Fsd,'r');
% 
%     xlim([-0.15 0.5])
%     yl = ylim();
%     line([0 0],yl,'color','k')
    pause
    clf
end

%%
for cc = 1 :length(use_lfps)
    cc
    %         matrix_errorbar_plot(squeeze(all_strgp_Vbb(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'k');
    %         hold on
    %         matrix_errorbar_plot(squeeze(all_goodtrgp_Vbb(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'b');
    %         matrix_errorbar_plot(squeeze(all_badtrgp_Vbb(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
    %
        matrix_errorbar_plot(squeeze(all_strgp_Vbb(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'k');
%         hold on
%         matrix_errorbar_plot(squeeze(all_strgp_Valpha(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'b');
%         matrix_errorbar_plot(squeeze(all_strgp_Vgamma(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
%         hold on
%         matrix_errorbar_plot(squeeze(all_mactrgp_Vbb(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'b');
%         matrix_errorbar_plot(squeeze(all_mictrgp_Vbb(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
%         matrix_errorbar_plot(squeeze(all_goodtrgp_Vgamma(:,:,cc)-all_badtrgp_Vgamma(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
%     
%         matrix_errorbar_plot(squeeze(all_strgp_Vgammaf(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'k');
%         hold on
%         matrix_errorbar_plot(squeeze(all_goodtrgp_Vgammaf(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'b');
%         matrix_errorbar_plot(squeeze(all_badtrgp_Vgammaf(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
%         matrix_errorbar_plot(squeeze(all_goodtrgp_Vgammaf(:,:,cc)-all_badtrgp_Vgammaf(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
    
%         matrix_errorbar_plot(squeeze(all_strgp_spk(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'k');
%         hold on
%         matrix_errorbar_plot(squeeze(all_goodtrgp_spk(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'b');
%         matrix_errorbar_plot(squeeze(all_badtrgp_spk(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
    
    %     matrix_errorbar_plot(squeeze(all_strgp_Vgamma(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'k');
    %     hold on
    %     matrix_errorbar_plot(squeeze(all_nattrgp_Vgamma(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'b');
    %     matrix_errorbar_plot(squeeze(all_noisetrgp_Vgamma(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
 
%             matrix_errorbar_plot(squeeze(all_strgp_spk(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'k');
%     hold on
% %     matrix_errorbar_plot(squeeze(all_inobjtrgp_Valpha(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
%     matrix_errorbar_plot(squeeze(all_backtrgp_spk(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'g');
%     matrix_errorbar_plot(squeeze(all_fortrgp_spk(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'b');

%     matrix_errorbar_plot(squeeze(all_strgp_Vgammaf(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'k');
%     hold on
%     matrix_errorbar_plot(squeeze(all_nattrgp_Vgammaf(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'b');
%     matrix_errorbar_plot(squeeze(all_noisetrgp_Vgammaf(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
%     matrix_errorbar_plot(squeeze(all_objtrgp_Vgammaf(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'g');

%         matrix_errorbar_plot(squeeze(all_strgp_Valpha(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'k');
%     hold on
% %     matrix_errorbar_plot(squeeze(all_inobjtrgp_Valpha(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
%     matrix_errorbar_plot(squeeze(all_backtrgp_Valpha(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'g');
%     matrix_errorbar_plot(squeeze(all_fortrgp_Valpha(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'b');

%     matrix_errorbar_plot(squeeze(all_strgp_Vgamma(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'k');
%     hold on
%     matrix_errorbar_plot(squeeze(all_mactrgp_Vgamma(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'b');
%     matrix_errorbar_plot(squeeze(all_mictrgp_Vgamma(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
%     matrix_errorbar_plot(squeeze(all_strgp_spk(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'k');
%     hold on
%     matrix_errorbar_plot(squeeze(all_mactrgp_spk(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'b');
%     matrix_errorbar_plot(squeeze(all_mictrgp_spk(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
%     
%     matrix_errorbar_plot(squeeze(all_strgp_Vgammaf(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'k');
%     hold on
%     matrix_errorbar_plot(squeeze(all_mactrgp_Vgammaf(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'b');
%     matrix_errorbar_plot(squeeze(all_mictrgp_Vgammaf(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');

xlim([-2 5])
    pause
    clf
end


%%
% all_avg_mactrg_Vbb = squeeze(mean(all_mactrg_Vbb,1))';
% all_avg_mactrg_Valpha = squeeze(mean(all_mactrg_Valpha,1))';
% all_avg_mactrg_Vgamma = squeeze(mean(all_mactrg_Vgamma,1))';
% all_avg_mactrg_spk = squeeze(mean(all_mactrg_spk,1))';
all_avg_mactrg_Vbb = squeeze(mean(all_strg_Vbb,1))';
all_avg_mactrg_Valpha = squeeze(mean(all_strg_Valpha,1))';
all_avg_mactrg_Vgamma = squeeze(mean(all_strg_Vgamma,1))';
all_avg_mactrg_spk = squeeze(mean(all_strg_spk,1))';

subplot(3,1,1)
        matrix_errorbar_plot(squeeze(all_avg_mactrg_Vbb),lags/Fsd,'k');
        hold on
        matrix_errorbar_plot(squeeze(all_avg_mactrg_Valpha),lags/Fsd,'b');
    xlim([-0.15 0.4])
    xlabel('Time since fixation onset (s)','fontsize',16)
    yl = ylim();
    line([0 0],yl,'color','k')
subplot(3,1,2)
    matrix_errorbar_plot(squeeze(all_avg_mactrg_Vgamma),lags/Fsd,'r');
    xlim([-0.15 0.4])
    xlabel('Time since fixation onset (s)','fontsize',16)
    yl = ylim();
    line([0 0],yl,'color','k')
subplot(3,1,3)
matrix_errorbar_plot(squeeze(all_avg_mactrg_spk),lags/Fsd,'g');
    xlim([-0.15 0.4])
    xlabel('Time since fixation onset (s)','fontsize',16)
    yl = ylim();
    line([0 0],yl,'color','k')

    %%
    for cc = 1:length(use_lfps);
    subplot(3,1,1)
        matrix_errorbar_plot(squeeze(all_mactrg_Vbb(:,:,cc)),lags/Fsd,'k');
        hold on
        matrix_errorbar_plot(squeeze(all_mactrg_Valpha(:,:,cc)),lags/Fsd,'b');
    xlim([-0.15 0.4])
    xlabel('Time since fixation onset (s)','fontsize',16)
    yl = ylim();
    line([0 0],yl,'color','k')
subplot(3,1,2)
    matrix_errorbar_plot(squeeze(all_mactrg_Vgamma(:,:,cc)),lags/Fsd,'r');
    xlim([-0.15 0.4])
    xlabel('Time since fixation onset (s)','fontsize',16)
    yl = ylim();
    line([0 0],yl,'color','k')
subplot(3,1,3)
matrix_errorbar_plot(squeeze(all_mactrg_spk(:,:,cc)),lags/Fsd,'g');
    xlim([-0.15 0.4])
    xlabel('Time since fixation onset (s)','fontsize',16)
    yl = ylim();
    line([0 0],yl,'color','k')

    pause
    clf
    end
    
    %%
        for cc = 1:length(use_lfps);
    subplot(3,1,1)
        matrix_errorbar_plot(squeeze(all_strg_Vbb(:,:,cc)),lags/Fsd,'k');
        hold on
        matrix_errorbar_plot(squeeze(all_strg_Valpha(:,:,cc)),lags/Fsd,'b');
    xlim([-0.15 0.4])
    xlabel('Time since fixation onset (s)','fontsize',16)
    yl = ylim();
    line([0 0],yl,'color','k')
subplot(3,1,2)
    matrix_errorbar_plot(squeeze(all_strg_Vgamma(:,:,cc)),lags/Fsd,'r');
    xlim([-0.15 0.4])
    xlabel('Time since fixation onset (s)','fontsize',16)
    yl = ylim();
    line([0 0],yl,'color','k')
subplot(3,1,3)
matrix_errorbar_plot(squeeze(all_strg_spk(:,:,cc)),lags/Fsd,'g');
    xlim([-0.15 0.4])
    xlabel('Time since fixation onset (s)','fontsize',16)
    yl = ylim();
    line([0 0],yl,'color','k')

    pause
    clf
    end

    %%
        for cc = 1:length(use_lfps);
    subplot(3,1,1)
        matrix_errorbar_plot(squeeze(all_mactrgp_Vbb(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'k');
        hold on
        matrix_errorbar_plot(squeeze(all_mactrgp_Valpha(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'b');
%     xlim([-0.15 0.4])
    xlabel('Time since fixation onset (s)','fontsize',16)
    yl = ylim();
    line([0 0],yl,'color','k')
subplot(3,1,2)
    matrix_errorbar_plot(squeeze(all_mactrgp_Vgamma(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'r');
%     xlim([-0.15 0.4])
    xlabel('Time since fixation onset (s)','fontsize',16)
    yl = ylim();
    line([0 0],yl,'color','k')
subplot(3,1,3)
matrix_errorbar_plot(squeeze(all_mactrgp_spk(:,:,cc)),lags_phase/(2*pi)*alpha_dp,'g');
%     xlim([-0.15 0.4])
    xlabel('Time since fixation onset (s)','fontsize',16)
    yl = ylim();
    line([0 0],yl,'color','k')

    pause
    clf
    end

    %%
    % all_avg_strg_Vbb = squeeze(nanmean(all_strg_Vbb,1))';
    % all_avg_nattrg_Vbb = squeeze(nanmean(all_nattrg_Vbb,1))';
    % all_avg_fortrg_Vbb = squeeze(nanmean(all_fortrg_Vbb,1))';
    % all_avg_backtrg_Vbb = squeeze(nanmean(all_backtrg_Vbb,1))';
    % all_avg_goodtrg_Vbb = squeeze(nanmean(all_goodtrg_Vbb,1))';
    % all_avg_badtrg_Vbb = squeeze(nanmean(all_badtrg_Vbb,1))';
    % all_avg_objtrg_Vbb = squeeze(nanmean(all_inobjtrg_Vbb,1))';
    all_avg_strg_Vbb = squeeze(nanmean(all_strg_Vgamma,3));
    all_avg_mactrg_Vbb = squeeze(nanmean(all_mactrg_Vgamma,3));
    all_avg_mictrg_Vbb = squeeze(nanmean(all_mictrg_Vgamma,3));
    all_avg_nattrg_Vbb = squeeze(nanmean(all_nattrg_Vbb,3));
    all_avg_fortrg_Vbb = squeeze(nanmean(all_fortrg_Vbb,3));
    all_avg_backtrg_Vbb = squeeze(nanmean(all_backtrg_Vbb,3));
    all_avg_goodtrg_Vbb = squeeze(nanmean(all_goodtrg_Vbb,3));
    all_avg_badtrg_Vbb = squeeze(nanmean(all_badtrg_Vbb,3));
    all_avg_objtrg_Vbb = squeeze(nanmean(all_inobjtrg_Vbb,3));
    
    matrix_errorbar_plot(squeeze(all_avg_strg_Vbb),lags/Fsd,'r');
    hold on
            matrix_errorbar_plot(squeeze(all_avg_mictrg_Vbb),lags/Fsd,'g');
    matrix_errorbar_plot(squeeze(all_avg_mactrg_Vbb),lags/Fsd,'b');
    %         matrix_errorbar_plot(squeeze(all_avg_nattrg_Vbb),lags/Fsd,'g');
%     matrix_errorbar_plot(squeeze(all_avg_fortrg_Vbb),lags/Fsd,'b');
%     matrix_errorbar_plot(squeeze(all_avg_backtrg_Vbb),lags/Fsd,'r');
    %         matrix_errorbar_plot(squeeze(all_avg_objtrg_Vbb),lags/Fsd,'g');
%     matrix_errorbar_plot(squeeze(all_avg_goodtrg_Vbb),lags/Fsd,'k');
%     matrix_errorbar_plot(squeeze(all_avg_badtrg_Vbb),lags/Fsd,'r');
    
    %     xlim([-0.15 0.5])
    yl = ylim();
    line([0 0],yl,'color','k')
    xlim([-0.2 0.5])
    xlabel('Time (s)','fontsize',16)
    ylabel('Amp (z)','fontsize',16)
    

    %%
%     all_avg_strg_Valpha = squeeze(nanmean(all_strg_Valpha,1))';
% all_avg_nattrg_Valpha = squeeze(nanmean(all_nattrg_Valpha,1))';
% all_avg_fortrg_Valpha = squeeze(nanmean(all_fortrg_Valpha,1))';
% all_avg_backtrg_Valpha = squeeze(nanmean(all_backtrg_Valpha,1))';
% all_avg_objtrg_Valpha = squeeze(nanmean(all_inobjtrg_Valpha,1))';
% all_avg_goodtrg_Valpha = squeeze(nanmean(all_goodtrg_Valpha,1))';
% all_avg_badtrg_Valpha = squeeze(nanmean(all_badtrg_Valpha,1))';
    all_avg_strg_Valpha = squeeze(nanmean(all_strg_ValphaA,3));
all_avg_nattrg_Valpha = squeeze(nanmean(all_nattrg_ValphaA,3));
all_avg_fortrg_Valpha = squeeze(nanmean(all_fortrg_ValphaA,3));
all_avg_backtrg_Valpha = squeeze(nanmean(all_backtrg_ValphaA,3));
all_avg_objtrg_Valpha = squeeze(nanmean(all_inobjtrg_ValphaA,3));
all_avg_goodtrg_Valpha = squeeze(nanmean(all_goodtrg_ValphaA,3));
all_avg_badtrg_Valpha = squeeze(nanmean(all_badtrg_ValphaA,3));

        matrix_errorbar_plot(squeeze(all_avg_strg_Valpha),lags/Fsd,'k');
        hold on
%         matrix_errorbar_plot(squeeze(all_avg_nattrg_Valpha),lags/Fsd,'g');
        matrix_errorbar_plot(squeeze(all_avg_fortrg_Valpha),lags/Fsd,'b');
        matrix_errorbar_plot(squeeze(all_avg_backtrg_Valpha),lags/Fsd,'r');
%         matrix_errorbar_plot(squeeze(all_avg_objtrg_Vbalpha),lags/Fsd,'g');
        matrix_errorbar_plot(squeeze(all_avg_goodtrg_Valpha),lags/Fsd,'g');
        matrix_errorbar_plot(squeeze(all_avg_badtrg_Valpha),lags/Fsd,'c');

        %     xlim([-0.15 0.5])
    yl = ylim();
    line([0 0],yl,'color','k')
xlim([-0.2 0.5])
xlabel('Time (s)','fontsize',16)
ylabel('Amp (z)','fontsize',16)

