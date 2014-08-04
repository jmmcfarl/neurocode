clear all
% close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

%%
load ./Expt3_fixbased_data.mat
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
spatial_mean = mean(spatial_mod_out);
spatial_std = std(spatial_mod_out);
spatial_mod_out = bsxfun(@minus,spatial_mod_out,spatial_mean);
spatial_mod_out = bsxfun(@rdivide,spatial_mod_out,spatial_std);


%%
sim_sac_blocks = [14 18 24 28 37 40 50] - 6;

repeat_inds = [1e3 2e3 3e3 2.01e5 2.02e5 2.03e5 4.01e5 4.02e5 4.03e5];

used_inds = find(~ismember(full_seof_vec,repeat_inds));
% used_inds = 1:length(full_seof_vec);

%%
stim_fs = 1e4/117.5;
Fs = 3e4;
% dsf = 120;Fsd = Fs/dsf;
dsf = 80;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1]/niqf,'high');
use_lfps = [1:96];
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
% desired_dt = 0.003;
% desired_dt = 0.006;
all_interp_phasegrams = [];
all_interp_ampgrams = [];
% all_binned_spks = [];
% all_expt_t = [];
% all_expt_seof = [];
% all_t_since_start = [];
% all_trial_id = [];
% all_tnum = [];
for bb = 1:length(sim_sac_blocks)
    
    %%
    fprintf('Analyzing block %d of %d\n',bb,length(sim_sac_blocks));
    load(sprintf('Expt%dClusterTimes.mat',sim_sac_blocks(bb)));
    
    filename = sprintf('Expt%dFullVmean.mat',sim_sac_blocks(bb));
    load(filename);
    
    Vmat = [];
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
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            %             dV = [dV decimate(V(cur_range),dsf)];
            dVf = [dVf filtfilt(filt_b,filt_a,decimate(V(cur_range),dsf))]; %do some high-pass filtering
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(ll,:) = dVf;
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
        
    cur_set = used_inds(full_expt_vec(used_inds) == sim_sac_blocks(bb));
    
    unwr_phasegram = unwrap(phasegrams);
    interp_phasegrams = interp1(t_ax,unwr_phasegram,full_t(cur_set));
    interp_phasegrams = mod(interp_phasegrams+pi,2*pi)-pi;
    interp_ampgrams = interp1(t_ax,ampgrams,full_t(cur_set));
    
    all_interp_phasegrams = cat(1,all_interp_phasegrams,interp_phasegrams);
    all_interp_ampgrams = cat(1,all_interp_ampgrams,interp_ampgrams);
    clear unwr_phasegram phasegrams interp_phasegrams interp_ampgrams ampgrams Vmat
    
end

%%
all_interp_ampgrams = bsxfun(@rdivide,all_interp_ampgrams,nanstd(all_interp_ampgrams));

%%
trial_set = unique(full_trial_vec(used_inds));
n_trials = length(trial_set);

xv_frac = 0;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(full_trial_vec(used_inds),trial_set(xv_set)));
tr_inds = find(~ismember(full_trial_vec(used_inds),trial_set(xv_set)))';

bad_inds = find(isnan(all_interp_phasegrams(:,1,1)));
xv_inds(ismember(xv_inds,bad_inds)) = [];
tr_inds(ismember(tr_inds,bad_inds)) = [];

% xv_inds_first = xv_inds(all_t_since_start(xv_inds) <= 0.2);
% tr_inds_first = tr_inds(all_t_since_start(tr_inds) <= 0.2);
% % xv_inds_second = xv_inds(all_t_since_start(xv_inds) > 0.2 & all_t_since_start(xv_inds)  <= 0.4);
% % tr_inds_second = tr_inds(all_t_since_start(tr_inds) > 0.2 & all_t_since_start(tr_inds)  <= 0.4);
% xv_inds_second = xv_inds(all_t_since_start(xv_inds) > 0.15);
% tr_inds_second = tr_inds(all_t_since_start(tr_inds) > 0.15);

%%
flen_t = stim_dur;
tent_centers = [0:.009:flen_t];

tent_centers = round(tent_centers/dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

tbmat = [zeros(ntents,tblen-1) tbmat];

trial_inds = zeros(size(full_t));
trial_inds(trial_start_inds) = 1;
trial_Tmat = zeros(length(used_inds),ntents);
for i = 1:ntents
    trial_Tmat(:,i) = conv(trial_inds(used_inds),tbmat(i,:),'same');
end

beg_dur = 0.15;
% late_indicator = zeros(size(used_inds));
% late_indicator(full_t_since_start(used_inds) >= beg_dur) = 1;
transient_indicator = zeros(size(used_inds));
transient_indicator(full_t_since_start(used_inds) < 0.4) = 1;
tr_inds(transient_indicator(tr_inds)==1) = [];
xv_inds(transient_indicator(xv_inds)==1) = [];

late_indicator = zeros(size(full_t));
for i = 1:length(trial_start_inds)
    late_indicator((trial_start_inds(i)+round(beg_dur/dt)):trial_stop_inds(i)) = 1;
end
late_indicator = late_indicator(used_inds);

%%
close all
NL_type = 1; %exp

reg_params.dl1_ind = 0;
reg_params.dl1_dep = 0;
reg_params.dl2_ind = 5000;
reg_params.dl2_dep = 5000;
reg_params.dl2_freq_ind = 0;
reg_params.dl2_freq_dep = 0;
reg_params.dl2_time_ind = 500;
reg_params.dl2_time_dep = 1000;
reg_params.l1_ind = 0;
reg_params.l1_dep = 0;
reg_params.l1t_ind = 0;
reg_params.l1t_dep = 0;
reg_params.is_phase = 0;

% reg_params2.dl1_ind = 0;
% reg_params2.dl1_dep = 0;
% reg_params2.dl2_ind = 10000;
% reg_params2.dl2_dep = 10000;
% reg_params2.dl2_freq_ind = 0;
% reg_params2.dl2_freq_dep = 0;
% reg_params2.dl2_time_ind = 100;
% reg_params2.dl2_time_dep = 100;
% reg_params2.l1_ind = 0;
% reg_params2.l1_dep = 0;
% reg_params2.l1t_ind = 0;
% reg_params2.l1t_dep = 0;
% reg_params2.is_phase = 0;

silent = 1;
NT = length(used_inds);
% used_inds = 1:NT;
nbins = 30;
pax = linspace(-pi,pi,nbins+1);

phase_set = [reshape(cos(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps)) reshape(sin(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps))];
ampphase_set = [reshape(all_interp_ampgrams,NT,length(wfreqs)*length(use_lfps)).*reshape(cos(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps)) ...
    reshape(all_interp_ampgrams,NT,length(wfreqs)*length(use_lfps)).*reshape(sin(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps))];
% phase_elec_set = [repmat(1:length(use_lfps),1,length(wfreqs)) repmat(1:length(use_lfps),1,length(wfreqs))];
% phase_elec_set = phase_elec_set(:);
phase_elec_set = ones(length(wfreqs),1)*(1:length(use_lfps));
phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];

phase_set = bsxfun(@times,phase_set,late_indicator);
ampphase_set = bsxfun(@times,ampphase_set,late_indicator);

%%
for ll = 1:length(use_lfps)
    fprintf('Electrode %d of %d\n',ll,length(use_lfps));
    cur_units = find(nearest_lfps==ll);
    use_elecs = ll;
    use_set = find(phase_elec_set==use_elecs);
    %     X_p = [phase_set(tr_inds,use_set) trial_Tmat(tr_inds,:)];
    %     X_a = [ampphase_set(tr_inds,use_set) trial_Tmat(tr_inds,:)];
    %     X_p_xv = [phase_set(xv_inds,use_set) trial_Tmat(xv_inds,:)];
    %     X_a_xv = [ampphase_set(xv_inds,use_set) trial_Tmat(xv_inds,:)];
    
    for cc = 1:length(cur_units)
        fprintf('Cell %d of %d\n',cc,length(cur_units));
        cur_spatial_mod_out = spatial_mod_out(full_stim_ids(used_inds),cur_units(cc));
%         cur_spatial_mod_out = interp1(full_t,cur_spatial_mod_out,all_expt_t);
%         cur_tr_inds = tr_inds;
%         cur_tr_inds(isnan(cur_spatial_mod_out(cur_tr_inds))) = [];
%         cur_xv_inds = xv_inds;
%         cur_xv_inds(isnan(cur_spatial_mod_out(cur_xv_inds))) = [];
        
        
        Robs = full_binned_spks(used_inds(tr_inds),cur_units(cc));
        tr_spkbns = convert_to_spikebins(Robs);
        
        tr_mean_rate(cur_units(cc)) = mean(Robs);
        
        stim_params = [0,0,ntents];
        Tmat = [trial_Tmat(tr_inds,:)];
        klen = size(Tmat,2);
        K0 = zeros(klen+1,1);
        [fitp_to_ns] = fit_GLM_phase_model(Tmat, Robs, K0, 1, stim_params, reg_params,1,NL_type);
        stim_ind_to_ns_filt(cur_units(cc),:) = fitp_to_ns.k((1):(ntents));
        to_ns_const(cur_units(cc)) = fitp_to_ns.k(end);
        
        %fit model with stim ind and dep time course
        stim_params = [0,0,ntents];
        Tmat = [trial_Tmat(tr_inds,:) bsxfun(@times,trial_Tmat(tr_inds,:),cur_spatial_mod_out(tr_inds))];
        klen = size(Tmat,2);
        K0 = zeros(klen+1,1);
        [fitp_to] = fit_GLM_phase_model(Tmat, Robs, K0, 1, stim_params, reg_params,1,NL_type);
        stim_ind_to_filt(cur_units(cc),:) = fitp_to.k((1:ntents));
        stim_dep_to_filt(cur_units(cc),:) = fitp_to.k((ntents+1):end-1);        
        to_const(cur_units(cc)) = fitp_to.k(end);
        
%         cur_phases = squeeze(all_interp_phasegrams(cur_tr_inds,:,nearest_lfps(cur_units(cc))));
%         Pmat = [];
%         for ww = 1:length(wfreqs)
%             cur_tb = tbrep(cur_phases(:,ww),pax);
%             cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
%             cur_tb(:,end) = [];
%             Pmat = [Pmat cur_tb];
%         end
%         tr_late_indicator = logical(late_indicator(cur_tr_inds));
%         Pmat = bsxfun(@times,Pmat,tr_late_indicator);
%         
% %         Tmat = [trial_Tmat(cur_tr_inds,:) bsxfun(@times,trial_Tmat(cur_tr_inds,:),cur_spatial_mod_out(cur_tr_inds))];
%         Xmat = [Pmat Tmat];
%         stim_params = [nbins,length(wfreqs),ntents,1];
%         klen = size(Xmat,2);
%         K0 = zeros(klen+1,1);
%         [fitp] = fit_GLM_phase_model(Xmat, Robs, K0,1, stim_params, reg_params,1,NL_type);
%         stim_ind_pfilt(cur_units(cc),:) = fitp.k(1:nbins*length(wfreqs));
%         stim_ind_tfilt(cur_units(cc),:) = fitp.k((nbins*length(wfreqs)+1):(nbins*length(wfreqs)+ntents));
%         stim_dep_tfilt(cur_units(cc),:) = fitp.k((nbins*length(wfreqs)+ntents+1):(nbins*length(wfreqs)+2*ntents));
%         phase_const(cur_units(cc)) = fitp.k(end);
        
%         stim_params = [length(wfreqs),1];
        stim_params = [length(wfreqs),1,ntents];
        X = [phase_set(tr_inds,use_set) Tmat];
        klen = size(X,2);
        K0 = zeros(klen+1,1);
        [fitp_phase] = fit_GLM_phase_model(X, Robs, K0,silent, stim_params, reg_params,1,NL_type);
        phase_cfilt(cur_units(cc),:) = fitp_phase.k(1:length(wfreqs));
        phase_sfilt(cur_units(cc),:) = fitp_phase.k((length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
        phase_ind_tfilt(cur_units(cc),:) = fitp_phase.k((length(wfreqs)*2+1):(length(wfreqs)*2+ntents));
        phase_dep_tfilt(cur_units(cc),:) = fitp_phase.k((length(wfreqs)*2+ntents+1):(length(wfreqs)*2+2*ntents));
        sinphase_const(cur_units(cc)) = fitp_phase.k(end);
        
        stim_params = [length(wfreqs),1,ntents];
        X = [ampphase_set(tr_inds,use_set) Tmat];
        klen = size(X,2);
        K0 = zeros(klen+1,1);
        [fitp_ampphase] = fit_GLM_phase_model(X, Robs, K0,silent, stim_params, reg_params,1,NL_type);
        ampphase_cfilt(cur_units(cc),:) = fitp_ampphase.k(1:length(wfreqs));
        ampphase_sfilt(cur_units(cc),:) = fitp_ampphase.k((length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
        ampphase_ind_tfilt(cur_units(cc),:) = fitp_ampphase.k((length(wfreqs)*2+1):(length(wfreqs)*2+ntents));
        ampphase_dep_tfilt(cur_units(cc),:) = fitp_ampphase.k((length(wfreqs)*2+ntents+1):(length(wfreqs)*2+2*ntents));
        ampsinphase_const(cur_units(cc)) = fitp_ampphase.k(end);
        
        if xv_frac > 0
        
        xv_Robs = full_binned_spks(used_inds(xv_inds),cur_units(cc));
        
        Tmat = [trial_Tmat(xv_inds,:)];
        xv_pred_rate = Tmat*fitp_to_ns.k(1:end-1) + fitp_to_ns.k(end);
        if NL_type == 0
            xv_pred_rate = log(1+exp(xv_pred_rate));
        else
            xv_pred_rate = exp(xv_pred_rate);
        end
        xv_to_ns_LL(cur_units(cc)) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);


        Tmat = [trial_Tmat(xv_inds,:) bsxfun(@times,trial_Tmat(xv_inds,:),cur_spatial_mod_out(xv_inds))];
        xv_pred_rate = Tmat*fitp_to.k(1:end-1) + fitp_to.k(end);
        if NL_type == 0
            xv_pred_rate = log(1+exp(xv_pred_rate));
        else
            xv_pred_rate = exp(xv_pred_rate);
        end
        xv_to_LL(cur_units(cc)) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);
        
%         cur_phases = squeeze(all_interp_phasegrams(cur_xv_inds,:,nearest_lfps(cur_units(cc))));
%         Pmat = [];
%         for ww = 1:length(wfreqs)
%             cur_tb = tbrep(cur_phases(:,ww),pax);
%             cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
%             cur_tb(:,end) = [];
%             Pmat = [Pmat cur_tb];
%         end
%         tr_late_indicator = logical(late_indicator(cur_xv_inds));
%         Pmat = bsxfun(@times,Pmat,tr_late_indicator);
%         
%         Xmat = [Pmat Tmat];
%         xv_pred_rate = Xmat*fitp.k(1:end-1) + fitp.k(end);
%         if NL_type == 0
%             xv_pred_rate = log(1+exp(xv_pred_rate));
%         else
%             xv_pred_rate = exp(xv_pred_rate);
%         end
%         xv_LL(cur_units(cc)) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);
        
        Xmat = [phase_set(xv_inds,use_set) Tmat];
        xv_pred_rate = Xmat*fitp_phase.k(1:end-1) + fitp_phase.k(end);
        if NL_type == 0
            xv_pred_rate = log(1+exp(xv_pred_rate));
        else
            xv_pred_rate = exp(xv_pred_rate);
        end
        xv_phase_LL(cur_units(cc)) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);

        Xmat = [ampphase_set(xv_inds,use_set) Tmat];
        xv_pred_rate = Xmat*fitp_ampphase.k(1:end-1) + fitp_ampphase.k(end);
        if NL_type == 0
            xv_pred_rate = log(1+exp(xv_pred_rate));
        else
            xv_pred_rate = exp(xv_pred_rate);
        end
        xv_ampphase_LL(cur_units(cc)) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);
        
        avg_rate = mean(Robs);
        null_pred = avg_rate*ones(size(xv_Robs));
        xv_null_LL(cur_units(cc)) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);
        end
    end
end

%%

%%
phase_ampkern = sqrt(phase_cfilt.^2+phase_sfilt.^2);
phase_phasekern = -atan2(phase_cfilt,phase_sfilt)+pi/2;
ampphase_ampkern = sqrt(ampphase_cfilt.^2+ampphase_sfilt.^2);
ampphase_phasekern = -atan2(ampphase_cfilt,ampphase_sfilt)+pi/2;

%%
cd ~/Data/bruce/G075
save sim_sac_phase_mod_v4_late_train_nt wfreqs *filt *LL *const tent_centers ntents pax nbins use_lfps tr_mean_rate spatial_mean spatial_std
% save sim_sac_phase_mod_to ntents tent_centers to_*

%%
to_ns_imp = (xv_null_LL - xv_to_ns_LL)/log(2);
to_imp = (xv_null_LL - xv_to_LL)/log(2);
phase_imp = (xv_null_LL - xv_phase_LL)/log(2);
ampphase_imp = (xv_null_LL - xv_ampphase_LL)/log(2);

% ucells = single_units;
ucells = 1:96;
bad_units = find(to_ns_imp(ucells) < 0);
ucells(bad_units) = [];
boxplot([to_ns_imp(ucells); to_imp(ucells); phase_imp(ucells); ampphase_imp(ucells)]');

%%
% use_elecs = 1:length(use_lfps);
% use_set = find(ismember(phase_elec_set,use_elecs));
% X_p = [phase_set(tr_inds,use_set) trial_Tmat(tr_inds,:)];
% X_a = [ampphase_set(tr_inds,use_set) trial_Tmat(tr_inds,:)];
% X_p_xv = [phase_set(xv_inds,use_set) trial_Tmat(xv_inds,:)];
% X_a_xv = [ampphase_set(xv_inds,use_set) trial_Tmat(xv_inds,:)];
%
% for cc = 1:96
%     fprintf('Cell %d of %d\n',cc,96);
%
%     Robs = all_binned_spks(tr_inds,cc);
%     tr_spkbns = convert_to_spikebins(Robs);
%
%     xv_Robs = all_binned_spks(xv_inds,cc);
%
%     %         stim_params = [length(wfreqs),1];
%     stim_params = [length(wfreqs),length(use_lfps),ntents];
%     klen = size(X_p,2);
%     K0 = zeros(klen+1,1);
%     [fitp_phase] = fit_GLM_phase_model(X_p, Robs, K0,silent, stim_params, reg_params,1,NL_type);
%     full_phase_cfilt(cc,:) = fitp_phase.k(1:length(use_elecs)*length(wfreqs));
%     full_phase_sfilt(cc,:) = fitp_phase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
%     full_phase_tfilt(cc,:) = fitp_phase.k((length(use_elecs)*length(wfreqs)*2+1):end-1);
%
%     [fitp_ampphase] = fit_GLM_phase_model(X_a, Robs, K0,silent, stim_params, reg_params2,1,NL_type);
%     ampphase_cfilt(cc,:) = fitp_ampphase.k(1:length(use_elecs)*length(wfreqs));
%     ampphase_sfilt(cc,:) = fitp_ampphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
%     ampphase_tfilt(cc,:) = fitp_phase.k((length(use_elecs)*length(wfreqs)*2+1):end-1);
%
%     xv_phase_pred_rate = X_p_xv*fitp_phase.k(1:end-1) + fitp_phase.k(end);
%     if NL_type == 0
%         xv_phase_pred_rate = log(1+exp(xv_phase_pred_rate));
%     else
%         xv_phase_pred_rate = exp(xv_phase_pred_rate);
%     end
%     xv_phase_LL(cc) = -sum(xv_Robs.*log(xv_phase_pred_rate)-xv_phase_pred_rate)/sum(xv_Robs);
%
%     xv_ampphase_pred_rate = X_a_xv*fitp_ampphase.k(1:end-1) + fitp_ampphase.k(end);
%     if NL_type == 0
%         xv_ampphase_pred_rate = log(1+exp(xv_ampphase_pred_rate));
%     else
%         xv_ampphase_pred_rate = exp(xv_ampphase_pred_rate);
%     end
%     xv_ampphase_LL(cc) = -sum(xv_Robs.*log(xv_ampphase_pred_rate)-xv_ampphase_pred_rate)/sum(xv_Robs);
%
%     avg_rate = mean(Robs);
%     null_pred = avg_rate*ones(size(xv_Robs));
%     xv_null_LL(cc) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);
%
% end
%
