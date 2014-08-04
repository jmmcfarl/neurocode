clear all
% close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

%%
load ./Expt3_fixbased_data_all.mat
Pix2Deg = 0.018837;
NT = length(full_stim_ids);
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

norm_fac = std(resh_all_stims(:));
resh_all_stims = resh_all_stims/norm_fac;

load ./NS40_gabor_mods.mat gabor*
clear gabor_*filt
for t = 1:96
    
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params_f(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params_f(t,1:6),pi/2);
    
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
end

all_gabor_out1 = resh_all_stims*gabor_emp1_filt';
all_gabor_out2 = resh_all_stims*gabor_emp2_filt';
spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);
% spatial_mean = mean(spatial_mod_out);
% spatial_std = std(spatial_mod_out);


%%
sim_sac_blocks = [14 18 24 28 37 40 50] - 6;

% repeat_inds = [1e3 2e3 3e3 2.01e5 2.02e5 2.03e5 4.01e5 4.02e5 4.03e5];
repeat_inds = [1e3 2e3 3e3 4.01e5 4.02e5 4.03e5];

% used_inds = find(ismember(full_seof_vec,repeat_inds));
used_inds = find(ismember(full_seof_vec,repeat_inds) & full_result_vec == 1);
% used_inds = 1:length(full_seof_vec);

%%
stim_fs = 1e4/117.5;
Fs = 3e4;
% dsf = 120;Fsd = Fs/dsf;
dsf = 80;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1]/niqf,'high');
[filt_blf,filt_alf] = butter(2,[1 25]/niqf);
[filt_blf2,filt_alf2] = butter(2,[2.5 25]/niqf);
use_lfps = [1:1:96];
use_units = 1:96;
% use_lfps = [1 63];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);

scales = logspace(log10(10),log10(125),30);
scales = scales*60/dsf;
scales(1:3) = []; %only up to 40 Hz...
% scales = scales*2;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);



%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
nearest_lfps = nan(length(use_lfps),1);
for ll = 1:96
    all_dists = sqrt((X_pos-X_pos(ll)).^2 + (Y_pos-Y_pos(ll)).^2);
    %    all_dists(ll) = inf;
    [~,best_loc] = min(all_dists(use_lfps));
    nearest_lfps(ll) = best_loc;
end

%%
stim_dur = 0.47;
trial_cnt = 1;
all_interp_phasegrams = [];
all_interp_ampgrams = [];
all_interp_lfps = [];
all_interp_lfps_lf = [];
all_interp_lfps_lf2 = [];

for bb = 1:length(sim_sac_blocks)
    
    %%
    fprintf('Analyzing block %d of %d\n',bb,length(sim_sac_blocks));
    load(sprintf('Expt%dClusterTimes.mat',sim_sac_blocks(bb)));
    
    filename = sprintf('Expt%dFullVmean.mat',sim_sac_blocks(bb));
    load(filename);
    
    Vmat = [];
    Vmatlf = [];
    Vmatlf2 = [];
    %     Vmatf = [];
    phasegrams = [];
    ampgrams = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',sim_sac_blocks(bb),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = V + FullV.sumscale*sumv;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        %         dV = [];
        dVf = [];
        dVlf = [];
        dVlf2 = [];
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            dec_v = decimate(V(cur_range),dsf);
            %             dV = [dV decimate(V(cur_range),dsf)];
            dVf = [dVf filtfilt(filt_b,filt_a,dec_v)]; %do some high-pass filtering
            dVlf = [dVlf filtfilt(filt_blf,filt_alf,dec_v)]; %do some high-pass filtering
            dVlf2 = [dVlf2 filtfilt(filt_blf2,filt_alf2,dec_v)]; %do some high-pass filtering
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(ll,:) = dVf;
        Vmatlf(ll,:) = dVlf;
        Vmatlf2(ll,:) = dVlf2;
        temp = cwt(dVf,scales,'cmor1-1');
        phasegrams(:,:,ll) = angle(temp)';
        ampgrams(:,:,ll) = abs(temp)';
    end
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,2)+1:end) = [];
    
    Vmat = Vmat';
    Vmatlf = Vmatlf';
    Vmatlf2 = Vmatlf2';
    cur_set = used_inds(full_expt_vec(used_inds) == sim_sac_blocks(bb));
    
    unwr_phasegram = unwrap(phasegrams);
    interp_phasegrams = interp1(t_ax,unwr_phasegram,full_t(cur_set));
    interp_phasegrams = mod(interp_phasegrams+pi,2*pi)-pi;
    interp_ampgrams = interp1(t_ax,ampgrams,full_t(cur_set));
    
    interp_lfps = interp1(t_ax,Vmat,full_t(cur_set));
    interp_lfps_lf = interp1(t_ax,Vmatlf,full_t(cur_set));
    interp_lfps_lf2 = interp1(t_ax,Vmatlf2,full_t(cur_set));
    
    all_interp_phasegrams = cat(1,all_interp_phasegrams,interp_phasegrams);
    all_interp_ampgrams = cat(1,all_interp_ampgrams,interp_ampgrams);
    all_interp_lfps = cat(1,all_interp_lfps,interp_lfps);
    all_interp_lfps_lf = cat(1,all_interp_lfps_lf,interp_lfps_lf);
    all_interp_lfps_lf2 = cat(1,all_interp_lfps_lf2,interp_lfps_lf2);
    clear unwr_phasegram phasegrams interp_phasegrams interp_ampgrams ampgrams Vmat interp_lfps*
    
end

%%
all_interp_ampgrams = bsxfun(@rdivide,all_interp_ampgrams,nanstd(all_interp_ampgrams));

%%
trial_set = unique(full_trial_vec(used_inds));
tot_n_trials = length(trial_set);
xv_frac = 1;
n_xv_trials = round(tot_n_trials*xv_frac);
xv_set = randperm(tot_n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(full_trial_vec(used_inds),trial_set(xv_set)));
tr_inds = find(~ismember(full_trial_vec(used_inds),trial_set(xv_set)));

bad_inds = find(isnan(all_interp_lfps(:,1)));
xv_inds(ismember(xv_inds,bad_inds)) = [];
tr_inds(ismember(tr_inds,bad_inds)) = [];

%%
flen_t = stim_dur;
tent_centers = [0:.009:flen_t];

tent_centers = round(tent_centers/dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

tbmat = [zeros(ntents,tblen-1) tbmat];

trial_binds = zeros(size(full_t));
trial_binds(trial_start_inds) = 1;
trial_einds = zeros(size(full_t));
trial_einds(trial_stop_inds) = 1;
trial_Tmat = zeros(length(used_inds),ntents);
for i = 1:ntents
    trial_Tmat(:,i) = conv(trial_binds(used_inds),tbmat(i,:),'same');
end

beg_dur = 0.15;
% late_indicator = zeros(size(used_inds));
% late_indicator(full_t_since_start(used_inds) >= beg_dur) = 1;
late_indicator = zeros(size(full_t));
for i = 1:length(trial_start_inds)
    late_indicator((trial_start_inds(i)+round(beg_dur/dt)):trial_stop_inds(i)) = 1;
end
late_indicator = late_indicator(used_inds);

transient_indicator = zeros(size(used_inds));
transient_indicator(full_t_since_start(used_inds) < 0.47) = 1;

% xv_inds(transient_indicator(xv_inds)==1) = [];

xv_inds_late = xv_inds(late_indicator(xv_inds)==1);
%%
used_stim_starts = find(trial_binds(used_inds)==1);
used_stim_stops = find(trial_einds(used_inds)==1);
used_stim_durs = used_stim_stops-used_stim_starts;

% min_tdur = round(stim_dur/dt)-1;
% used_trials = find(used_trial_durs >= min_tdur);
% used_trial_starts = used_trial_starts(used_trials);
% used_trial_stops = used_trial_stops(used_trials);

used_trial_starts = 1+find(diff(full_trial_vec(used_inds))~=0);
used_trial_stops = find(diff(full_trial_vec(used_inds))~=0);
used_trial_starts = [1; used_trial_starts];
used_trial_stops = [used_trial_stops; length(used_inds)];
used_trial_durs = used_trial_stops-used_trial_starts;

used_stim_imnum = full_image_vec(used_inds(used_stim_starts));
used_stim_imnum = mod(used_stim_imnum,10);


NT = length(used_inds);
phase_set = [reshape(cos(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps)) reshape(sin(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps))];
ampphase_set = [reshape(all_interp_ampgrams,NT,length(wfreqs)*length(use_lfps)).*reshape(cos(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps)) ...
    reshape(all_interp_ampgrams,NT,length(wfreqs)*length(use_lfps)).*reshape(sin(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps))];
phase_elec_set = ones(length(wfreqs),1)*(1:length(use_lfps));
phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];

phase_set = bsxfun(@times,phase_set,late_indicator');
ampphase_set = bsxfun(@times,ampphase_set,late_indicator');

%%
load ./sim_sac_phase_mod_v4_late_train_nt.mat
% load ./sim_sac_phase_mod_v4_all_train_nt.mat
% load ./sim_sac_phase_mod_v4_late_train_nt.mat
load sim_sac_tr_mean_rate tr_mean_rate
all_spatial_mod_out = bsxfun(@minus,spatial_mod_out,spatial_mean);
all_spatial_mod_out = bsxfun(@rdivide,all_spatial_mod_out,spatial_std);

all_spatial_mod_out = all_spatial_mod_out(full_stim_ids(used_inds),:);

to_filt = [stim_ind_to_filt stim_dep_to_filt];
to_ns_filt = stim_ind_to_ns_filt;
phase_filt = [phase_cfilt phase_sfilt];
ampphase_filt = [ampphase_cfilt ampphase_sfilt];
phasef_filt = [phase_cfilt phase_sfilt phase_ind_tfilt phase_dep_tfilt];
ampphasef_filt = [ampphase_cfilt ampphase_sfilt ampphase_ind_tfilt ampphase_dep_tfilt];
ampphaseo_filt = [ampphase_ind_tfilt ampphase_dep_tfilt];

%%
trial_len = round(stim_dur/dt)*8;

fresolve_predictor = nan(96,length(used_inds),length(wfreqs));
phase_predictor = nan(96,length(used_inds));
ampphase_predictor = nan(96,length(used_inds));
phasef_predictor = nan(96,length(used_inds));
ampphasef_predictor = nan(96,length(used_inds));
to_predictor = nan(96,length(used_inds));
to_ns_predictor = nan(96,length(used_inds));
psth_predictor = nan(96,length(used_inds));

for cc = 1:96
    fprintf('Cell %d of %d\n',cc,96);
    cur_lfp = nearest_lfps(cc);
    
    Tmat = [trial_Tmat bsxfun(@times,trial_Tmat,all_spatial_mod_out(:,cc))];
    
    for ss = 1:length(repeat_inds)
        
        cur_inds = xv_inds(full_seof_vec(used_inds(xv_inds)) ==repeat_inds(ss));
        n_trials(ss) = length(cur_inds)/trial_len;
        if mod(n_trials(ss),1) ~= 0; disp('Trial length error!'); end;
        
        use_set = find(phase_elec_set==cur_lfp);
        
        cur_fresolve_out = bsxfun(@times,ampphase_set(cur_inds,use_set(1:length(wfreqs))),ampphase_cfilt(cc,:));
        cur_fresolve_out = cur_fresolve_out + bsxfun(@times,ampphase_set(cur_inds,use_set((length(wfreqs)+1:end))),ampphase_sfilt(cc,:));
        
        cur_phasemod_out = phase_set(cur_inds,use_set)*phase_filt(cc,:)';
        cur_ampphasemod_out = ampphase_set(cur_inds,use_set)*ampphase_filt(cc,:)';
        cur_to_ns_out = trial_Tmat(cur_inds,:)*stim_ind_to_ns_filt(cc,:)';
        cur_to_out = Tmat(cur_inds,:)*to_filt(cc,:)';
        cur_oampmod_out = Tmat(cur_inds,:)*ampphaseo_filt(cc,:)';
        cur_fphasemod_out = [phase_set(cur_inds,use_set) Tmat(cur_inds,:)]*phasef_filt(cc,:)';
        cur_fampphasemod_out = [ampphase_set(cur_inds,use_set) Tmat(cur_inds,:)]*ampphasef_filt(cc,:)';
        
        fresolve_out{ss,cc} = reshape(cur_fresolve_out,trial_len,n_trials(ss),length(wfreqs));
        phasemod_out{ss,cc} = reshape(cur_phasemod_out,trial_len,n_trials(ss));
        ampmod_out{ss,cc} = reshape(cur_ampphasemod_out,trial_len,n_trials(ss));
        to_ns_mod_out{ss,cc} = reshape(cur_to_ns_out,trial_len,n_trials(ss));
        to_mod_out{ss,cc} = reshape(cur_to_out,trial_len,n_trials(ss));
        oampmod_out{ss,cc} = reshape(cur_oampmod_out,trial_len,n_trials(ss));
        fphasemod_out{ss,cc} = reshape(cur_fphasemod_out,trial_len,n_trials(ss));
        fampmod_out{ss,cc} = reshape(cur_fampphasemod_out,trial_len,n_trials(ss));
        
        sm_binned_spk_mat{ss,cc} = reshape(full_smbinned_spks(used_inds(cur_inds),cc),trial_len,n_trials(ss));
        binned_spk_mat{ss,cc} = reshape(full_binned_spks(used_inds(cur_inds),cc),trial_len,n_trials(ss));
        lfp_mat{ss,cc} = reshape(all_interp_lfps(cur_inds,cur_lfp),trial_len,n_trials(ss));
        lfp_mat2{ss,cc} = reshape(all_interp_lfps_lf2(cur_inds,cur_lfp),trial_len,n_trials(ss));
        
        %         white_lfp = reshape(all_interp_ampgrams(cur_inds,:,cur_lfp).*cos(all_interp_phasegrams(cur_inds,:,cur_lfp)),[trial_len,n_trials(ss),length(wfreqs)]);
        %         white_lfp_mat{ss,cc} = nansum(white_lfp,3);
        
        avg_fresolve_out(ss,cc,:,:) = nanmean(fresolve_out{ss,cc},2);
        avg_phasemod_out(ss,cc,:) = nanmean(phasemod_out{ss,cc},2);
        avg_ampmod_out(ss,cc,:) = nanmean(ampmod_out{ss,cc},2);
        avg_fphasemod_out(ss,cc,:) = nanmean(fphasemod_out{ss,cc},2);
        avg_oampmod_out(ss,cc,:) = nanmean(oampmod_out{ss,cc},2);
        avg_fampmod_out(ss,cc,:) = nanmean(fampmod_out{ss,cc},2);
        avg_to_mod_out(ss,cc,:) = nanmean(to_mod_out{ss,cc},2);
        avg_to_ns_mod_out(ss,cc,:) = nanmean(to_ns_mod_out{ss,cc},2);
        avg_binned_spk(ss,cc,:) = nanmean(binned_spk_mat{ss,cc},2);
        avg_lfp_mat(ss,cc,:) = nanmean(lfp_mat{ss,cc},2);
        avg_lfp_mat2(ss,cc,:) = nanmean(lfp_mat2{ss,cc},2);
        %         avg_wlfp_mat(ss,cc,:) = nanmean(white_lfp_mat{ss,cc},2);
        
        to_ns_predictor(cc,cur_inds) = cur_to_ns_out;
        to_predictor(cc,cur_inds) = cur_to_out;
        fresolve_predictor(cc,cur_inds,:) = cur_fresolve_out;
        phase_predictor(cc,cur_inds) = repmat(squeeze(avg_fphasemod_out(ss,cc,:)),[n_trials(ss) 1]);
        ampphase_predictor(cc,cur_inds) = repmat(squeeze(avg_fampmod_out(ss,cc,:)),[n_trials(ss) 1]);
        psth_predictor(cc,cur_inds) = repmat(squeeze(avg_binned_spk(ss,cc,:)),[n_trials(ss) 1]);
        phasef_predictor(cc,cur_inds) = cur_fphasemod_out;
        ampphasef_predictor(cc,cur_inds) = cur_fampphasemod_out;
        
    end
end

%% NOW FIT GLMS to AVG predictors
% for cc = 1:96
%     Robs = full_binned_spks(used_inds(tr_inds),cc);
%     xvRobs = full_binned_spks(used_inds(xv_inds),cc);
%
%     X = [to_predictor(cc,tr_inds)' phase_predictor(cc,tr_inds)'];
%     [avgphase_beta(cc,:)] = glmfit(X,Robs,'poisson');
%     X = [to_predictor(cc,xv_inds)' phase_predictor(cc,xv_inds)'];
%     pred_rate = glmval(avgphase_beta(cc,:)',X,'log');
%     avgphase_xvLL(cc) = -sum(xvRobs.*log(pred_rate)-pred_rate)/sum(xvRobs);
%
%     X = [to_predictor(cc,tr_inds)' ampphase_predictor(cc,tr_inds)'];
%     [avgampphase_beta(cc,:)] = glmfit(X,Robs,'poisson');
%     X = [to_predictor(cc,xv_inds)' ampphase_predictor(cc,xv_inds)'];
%      pred_rate = glmval(avgampphase_beta(cc,:)',X,'log');
%     avgampphase_xvLL(cc) = -sum(xvRobs.*log(pred_rate)-pred_rate)/sum(xvRobs);
%
% %     obs_to_psth = [];
% %     obs_avgphase_pred = [];
% %     obs_avgampphase_pred = [];
% %     obs_psth = [];
% %     for ss = 1:length(repeat_inds)
% %     to_pred_mat = exp(to_mod_out{ss,cc} + to_const(cc));
% %     to_pred_psth = nanmean(to_pred_mat,2);
% %     obs_to_psth = [obs_to_psth; to_pred_psth(:)];
% %     obs_avgphasse_pred = [obs_avgphase_pred; squeeze(avg_phasemod_out(ss,cc,:))];
% %     obs_avgampphase_pred = [obs_avgampphase_pred; squeeze(avg_ampmod_out(ss,cc,:))];
% %     obs_psth = [obs_psth; squeeze(avg_binned_spk(ss,cc,:))];
% %     end
%
% end

%%
use_xv_inds = xv_inds;
use_xv_inds(transient_indicator(use_xv_inds)==1) = [];
use_xv_inds_late = xv_inds_late;
use_xv_inds_late(transient_indicator(use_xv_inds_late)==1) = [];

for cc = 1:96
    Robs = full_binned_spks(used_inds(use_xv_inds),cc)';
    
    pred_rate = exp(to_ns_predictor(cc,use_xv_inds) + to_ns_const(cc));
    to_ns_LL(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    pred_rate = exp(to_predictor(cc,use_xv_inds) + to_const(cc));
    to_LL(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    pred_rate = exp(phase_predictor(cc,use_xv_inds) + sinphase_const(cc));
    phase_LL(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    pred_rate = exp(ampphase_predictor(cc,use_xv_inds) + ampsinphase_const(cc));
    ampphase_LL(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    pred_rate = exp(phasef_predictor(cc,use_xv_inds) + sinphase_const(cc));
    phasef_LL(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    pred_rate = exp(ampphasef_predictor(cc,use_xv_inds) + ampsinphase_const(cc));
    ampphasef_LL(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    pred_rate = psth_predictor(cc,use_xv_inds);
    pred_rate(pred_rate < 1e-20) = 1e-20;
    psth_LL(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    avg_rate = tr_mean_rate(cc);
    pred_rate = ones(size(Robs))*avg_rate;
    null_LL(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    
    Robs = full_binned_spks(used_inds(use_xv_inds_late),cc)';
    
    pred_rate = exp(to_ns_predictor(cc,use_xv_inds_late) + to_ns_const(cc));
    to_ns_LL_late(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    pred_rate = exp(to_predictor(cc,use_xv_inds_late) + to_const(cc));
    to_LL_late(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    pred_rate = exp(phase_predictor(cc,use_xv_inds_late) + sinphase_const(cc));
    phase_LL_late(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    pred_rate = exp(ampphase_predictor(cc,use_xv_inds_late) + ampsinphase_const(cc));
    ampphase_LL_late(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    pred_rate = exp(phasef_predictor(cc,use_xv_inds_late) + sinphase_const(cc));
    phasef_LL_late(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    pred_rate = exp(ampphasef_predictor(cc,use_xv_inds_late) + ampsinphase_const(cc));
    ampphasef_LL_late(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    pred_rate = psth_predictor(cc,use_xv_inds_late);
    pred_rate(pred_rate < 1e-20) = 1e-20;
    psth_LL_late(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
    
    avg_rate = tr_mean_rate(cc);
    pred_rate = ones(size(Robs))*avg_rate;
    null_LL_late(cc) = -sum(Robs.*log(pred_rate)-pred_rate)/sum(Robs);
end

to_ns_imp = (null_LL-to_ns_LL)/log(2);
to_imp = (null_LL-to_LL)/log(2);
phase_imp = (null_LL-phase_LL)/log(2);
ampphase_imp = (null_LL-ampphase_LL)/log(2);
fphase_imp = (null_LL-phasef_LL)/log(2);
fampphase_imp = (null_LL-ampphasef_LL)/log(2);

ucells = 1:96;
bad_units = find(to_ns_imp(ucells) < 0 | to_imp(ucells) < 0);
ucells(bad_units) = [];
figure
boxplot([to_ns_imp(ucells); to_imp(ucells); ampphase_imp(ucells); fampphase_imp(ucells)]');


to_rel_imp = to_imp - to_ns_imp;
ampphase_rel_imp = ampphase_imp - to_imp;
fampphase_rel_imp = fampphase_imp - ampphase_imp;
figure
boxplot([to_rel_imp(ucells); ampphase_rel_imp(ucells); fampphase_rel_imp(ucells)]');

to_rel_imp = -(to_LL_late - to_ns_LL_late)/log(2);
ampphase_rel_imp = -(ampphase_LL_late - to_LL_late)/log(2);
fampphase_rel_imp = -(ampphasef_LL_late - ampphase_LL_late)/log(2);
figure
boxplot([to_rel_imp(ucells); ampphase_rel_imp(ucells); fampphase_rel_imp(ucells)]');

%%
for cc = 1:96
    fprintf('Cell %d of %d\n',cc,96);
    for ss = 1:length(repeat_inds)
        
        cur_inds = xv_inds(full_seof_vec(used_inds(xv_inds)) ==repeat_inds(ss));
        n_trials(ss) = length(cur_inds)/trial_len;
        if mod(n_trials(ss),1) ~= 0; disp('Trial length error!'); end;
        
        use_set = find(phase_elec_set==cur_lfp);
        
        phasemod_out{ss,cc} = reshape(phase_predictor(cc,cur_inds),trial_len,n_trials(ss));
        ampmod_out{ss,cc} = reshape(ampphase_predictor(cc,cur_inds),trial_len,n_trials(ss));
        
        to_ns_pred{cc,ss} = exp(to_ns_mod_out{ss,cc} + to_ns_const(cc));
        to_pred{cc,ss} = exp(to_mod_out{ss,cc} + to_const(cc));
        fphasemod_pred{cc,ss} = exp(fphasemod_out{ss,cc} + sinphase_const(cc));
        fampmod_pred{cc,ss} = exp(fampmod_out{ss,cc} + ampsinphase_const(cc));
        phasemod_pred{cc,ss} = exp(phasemod_out{ss,cc} + sinphase_const(cc));
        ampmod_pred{cc,ss} = exp(ampmod_out{ss,cc} + ampsinphase_const(cc));
        
        avg_to_pred(cc,ss,:) = nanmean(to_pred{cc,ss},2);
        avg_ampmod_pred(cc,ss,:) = nanmean(ampmod_pred{cc,ss},2);
        %         Robs_mat = binned_spk_mat{ss,cc};
        %         to_ns_LL_mat{cc,ss} = -bsxfun(@minus,sum(Robs_mat.*log(to_ns_pred) - to_ns_pred),sum(Robs_mat));
        %         to_LL_mat{cc,ss} = -bsxfun(@minus,sum(Robs_mat.*log(to_pred) - to_pred),sum(Robs_mat));
        %         phase_LL_mat{cc,ss} = -bsxfun(@minus,sum(Robs_mat.*log(phasemod_pred) - phasemod_pred),sum(Robs_mat));
        %         amp_LL_mat{cc,ss} = -bsxfun(@minus,sum(Robs_mat.*log(ampmod_pred) - ampmod_pred),sum(Robs_mat));
        %         fphase_LL_mat{cc,ss} = -bsxfun(@minus,sum(Robs_mat.*log(fphasemod_pred) - fphasemod_pred),sum(Robs_mat));
        %         famp_LL_mat{cc,ss} = -bsxfun(@minus,sum(Robs_mat.*log(fampmod_pred) - fampmod_pred),sum(Robs_mat));
        
    end
end

%%
close all
trial_t = (1:trial_len)*dt;
stim_times = (40:40:320)/stim_fs;

for unit_num = 1:50
    % unit_num = 6;
    stim_num = 6;
    % unit_num = 4;
    
    % % close all
    % figure
    shadedErrorBar(trial_t,nanmean(sm_binned_spk_mat{stim_num,unit_num},2),nanstd(sm_binned_spk_mat{stim_num,unit_num},[],2)/sqrt(n_trials(stim_num)),{'color','k'});
    hold on
    plot(trial_t,nanmean(to_ns_pred{unit_num,stim_num},2),'r','linewidth',2)
    % plot(trial_t,nanmean(to_pred{unit_num,stim_num},2),'b','linewidth',2)
    % plot(trial_t,nanmean(ampmod_pred{unit_num,stim_num},2),'g','linewidth',2)
    plot(trial_t,nanmean(fampmod_pred{unit_num,stim_num},2),'g','linewidth',2)
    yl = ylim();
    for i = 1:length(stim_times)
        line(stim_times([i i]),yl,'color','k')
    end
    
    pause
    clf
end


%%
Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
siz = [1280 1280];
Fs = 1/Pix2Deg;
xax = linspace(-Nxp/2,Nxp/2,Nxp)/Fs; yax = linspace(-Nyp/2,Nyp/2,Nyp)/Fs;
[XAX,YAX] = meshgrid(xax,yax);
x0_avg = 0.35;
y0_avg = -0.4;
xw = 1;yw=1;
cd ~/Data/bruce/7_15_12/G034/
load ./gabor_tracking_varmeans
gabor_params_used = gabor_params_f{2};
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
%fit smoothed retinotopic surface
orientations = linspace(0,pi-pi/12,12);
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
used_winds = find(weights == 1);
tempinds = tempinds(used_winds);
interp_x(tempinds) = xpos_interp(used_winds);
interp_y(tempinds) = ypos_interp(used_winds);
xi = linspace(min(interp_x),max(interp_x),50);
yi = linspace(min(interp_y),max(interp_y),50);
[Xi,Yi] = meshgrid(xi,yi);
xxi = linspace(1,10,50);
[XXi,YYi] = meshgrid(xxi,xxi);


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

%%
close all

[~,use_ordx] = sort(X_pos(use_lfps));
[~,use_ordy] = sort(Y_pos(use_lfps));

stim_num = 6;
ttt = 373;

figure
subplot(2,1,1)
% pcolor(trial_t,1:96,exp(squeeze(avg_ampmod_out(stim_num,use_ordx,:))));shading flat
pcolor(trial_t,1:96,squeeze(avg_ampmod_out(stim_num,use_ordx,:)));shading flat
% caxis([0.8 1.2]);colorbar
caxis([-0.2 0.2]);colorbar
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
xlim([0.47 2.82])
xlim([0.47*2 0.47*5])
line(trial_t([ttt ttt]),yl,'color','w','linewidth',2)

% subplot(3,1,2)
% pcolor(trial_t,1:96,(squeeze(avg_lfp_mat(stim_num,use_ordx,:))));shading flat
% caxis([-10 5]*1e-5)
% colorbar
% yl = ylim();
%     for i = 1:length(stim_times)
%         line(stim_times([i i]),yl,'color','k')
%     end
%
% xlim([0.47 2.82])

subplot(2,1,2)
% pcolor(trial_t,1:96,(squeeze(avg_wlfp_mat(stim_num,use_ordx,:))));shading flat
pcolor(trial_t,1:96,(squeeze(avg_lfp_mat(stim_num,use_ordx,:))));shading flat
% caxis([-30 20])
caxis([-5 5]*1e-5)
colorbar
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
xlim([0.47 2.82])
xlim([0.47*2 0.47*5])
line(trial_t([ttt ttt]),yl,'color','w','linewidth',2)

figure
subplot(2,2,1)
cur_set = squeeze(avg_lfp_mat(stim_num,:,ttt));
temp = nan(10,10);
temp(use_ids) = cur_set(id_mat(use_ids));
% F = TriScatteredInterp(X_pos(:),Y_pos(:),cur_set(:));
% Vq = F(XXi,YYi);
% imagesc(xxi,xxi,Vq);colorbar; set(gca,'ydir','normal')
imagesc(temp);colorbar; set(gca,'ydir','normal')
% caxis([-2.5 2.5]*1e-5)
set(gca,'ydir','normal')
subplot(2,2,3)
F = TriScatteredInterp(interp_x,interp_y,cur_set(:));
Vq = F(Xi,Yi);
pcolor(xi,yi,Vq); shading flat; colorbar;
% caxis([-2.5 2.5]*1e-5)
subplot(2,2,2)
cur_set = squeeze(avg_ampmod_out(stim_num,:,ttt));
temp = nan(10,10);
temp(use_ids) = cur_set(id_mat(use_ids));
% F = TriScatteredInterp(X_pos(:),Y_pos(:),cur_set(:));
% Vq = F(XXi,YYi);
% imagesc(xxi,xxi,Vq);colorbar; set(gca,'ydir','normal')
imagesc(temp);colorbar; set(gca,'ydir','normal')
% caxis([-2.5 2.5]*1e-5)
set(gca,'ydir','normal')
subplot(2,2,4)
F = TriScatteredInterp(interp_x,interp_y,cur_set(:));
Vq = F(Xi,Yi);
pcolor(xi,yi,Vq); shading flat; colorbar;
% caxis([-2.5 2.5]*1e-5)

%%
% close all
ttt = 285;
stim_num = 6;
tr = 17; %17
tr_famp_out = nan(trial_len,96);
tr_wlfp_out = nan(trial_len,96);
for i = 1:96
    %     tr_famp_out(:,i) = exp(squeeze(fampmod_out{stim_num,i}(:,tr)));
    tr_famp_out(:,i) = squeeze(fampmod_out{stim_num,i}(:,tr));
    tr_wlfp_out(:,i) = squeeze(white_lfp_mat{stim_num,i}(:,tr));
    tr_lfp_out(:,i) = squeeze(lfp_mat{stim_num,i}(:,tr));
end

figure
subplot(2,1,1)
pcolor(trial_t,1:96,tr_famp_out(:,use_ordx)');shading flat
caxis([-0.7 0.7]);colorbar
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
xlim([0.47 2.82])
xlim([0.47*2 0.47*5])
line(trial_t([ttt ttt]),yl,'color','w','linewidth',2)

subplot(2,1,2)
pcolor(trial_t,1:96,tr_lfp_out(:,use_ordx)');shading flat
caxis([-1.25 1.25]*1e-4);colorbar
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
xlim([0.47 2.82])
xlim([0.47*2 0.47*5])
line(trial_t([ttt ttt]),yl,'color','w','linewidth',2)

figure
subplot(2,2,1)
cur_set = tr_famp_out(ttt,:)';
temp = nan(10,10);
temp(use_ids) = cur_set(id_mat(use_ids));
% F = TriScatteredInterp(X_pos(:),Y_pos(:),cur_set(:));
% Vq = F(XXi,YYi);
% imagesc(xxi,xxi,Vq);colorbar; set(gca,'ydir','normal')
imagesc(temp);colorbar; set(gca,'ydir','normal')
% caxis([-2.5 2.5]*1e-5)
set(gca,'ydir','normal')
subplot(2,2,3)
F = TriScatteredInterp(interp_x,interp_y,cur_set(:));
Vq = F(Xi,Yi);
pcolor(xi,yi,Vq); shading flat; colorbar;
% caxis([-2.5 2.5]*1e-5)
subplot(2,2,2)
cur_set = tr_lfp_out(ttt,:)';
temp = nan(10,10);
temp(use_ids) = cur_set(id_mat(use_ids));
% F = TriScatteredInterp(X_pos(:),Y_pos(:),cur_set(:));
% Vq = F(XXi,YYi);
% imagesc(xxi,xxi,Vq);colorbar; set(gca,'ydir','normal')
imagesc(temp);colorbar; set(gca,'ydir','normal')
% caxis([-2.5 2.5]*1e-5)
set(gca,'ydir','normal')
subplot(2,2,4)
F = TriScatteredInterp(interp_x,interp_y,cur_set(:));
Vq = F(Xi,Yi);
pcolor(xi,yi,Vq); shading flat; colorbar;
% caxis([-2.5 2.5]*1e-5)

%%

%%
bad_elecs = [16 67 92];
used_elecs = setdiff(use_lfps,bad_elecs);

el_dist = squareform(pdist([X_pos(use_lfps(used_elecs))' Y_pos(use_lfps(used_elecs))']));
el_dist = el_dist(:);
[XX,YY] = meshgrid(X_pos(use_lfps(used_elecs)),Y_pos(use_lfps(used_elecs)));
el_dist(XX < YY) = 0;
cset = find(el_dist > 0);
for i = 1:length(wfreqs)
    cur_corrmat = corr(squeeze(fresolve_predictor(used_elecs,use_xv_inds_late,i))');
    avg_freq_cor(i) = mean(cur_corrmat(cset));
    sem_freq_cor(i) = std(cur_corrmat(cset))/sqrt(length(cset));
end

dset = diag(length(used_elecs));
dset = logical(dset);
for i = 1:length(wfreqs)
    i
    for ss = 1:6
        for kk = 1:6
            cur_corrmat = corr(squeeze(avg_fresolve_out(ss,used_elecs,:,i))',squeeze(avg_fresolve_out(kk,used_elecs,:,i))');
            cur_corrmat(dset) = nan;
            avg_sd_freq_cor(i,ss,kk) = nanmean(cur_corrmat(:));
            %         sem_sd_freq_cor(i,ss) = std(cur_corrmat(cset))/sqrt(length(cset));
        end
    end
end
%%
all_ms_ampmod_out = [];
for ss = 1:6
    ss
    cur_ntrials = size(ampmod_out{ss,cc},2);
    cur_ms_ampmod_out = nan(96,cur_ntrials*560);
    for cc = 1:96
        temp = bsxfun(@minus,ampmod_out{ss,cc},mean(ampmod_out{ss,cc}));
        temp(1:80,:) = [];
        cur_ms_ampmod_out(cc,:) = temp(:);
    end
    all_ms_ampmod_out = cat(2,all_ms_ampmod_out,cur_ms_ampmod_out);
end
%%
stim_len = stim_dur/dt;
resh_phase_out = avg_phasemod_out(:,:,stim_len+1:end);
resh_phase_out = permute(resh_phase_out,[3 1 2]);
resh_phase_out = reshape(resh_phase_out,[stim_len*7*6,96]);

cur_corrmat = corr(resh_phase_out(:,used_elecs));
cur_corrmat = cur_corrmat(:);
cur_corrmat = cur_corrmat(cset);

xx = linspace(0.3,5,100);
figure
plot(el_dist(cset)*0.4,cur_corrmat,'r.','markersize',6)
hold on
[~,ord] = sort(el_dist(cset));
sm_fun = smooth(cur_corrmat(ord),5000,'rlowess');
plot(el_dist(cset(ord))*0.4,sm_fun,'k','linewidth',2)


% for cc = 1:96
%     cur_lfp = nearest_lfps(cc);
%     
%     use_set = find(phase_elec_set==cur_lfp);
%     cur_ampphasemod_out = ampphase_set(use_xv_inds,use_set)*ampphase_filt(cc,:)';
%     all_ampphasemod_out(cc,:) = cur_ampphasemod_out;
% end
bad_elecs = [16 67 92];
used_elecs = setdiff(use_lfps,bad_elecs);

% all_ampphasemod_out = all_ampphasemod_out';
all_ms_ampmod_out = all_ms_ampmod_out';
cur_corrmat2 = corr(all_ms_ampmod_out(:,used_elecs));
cur_corrmat2 = cur_corrmat2(:);
cur_corrmat2 = cur_corrmat2(cset);

xx = linspace(0.3,5,100);
% figure
plot(el_dist(cset)*0.4,cur_corrmat2,'.','markersize',6)
hold on
[~,ord] = sort(el_dist(cset));
sm_fun = smooth(cur_corrmat2(ord),5000,'rlowess');
plot(el_dist(cset(ord))*0.4,sm_fun,'k','linewidth',2)

