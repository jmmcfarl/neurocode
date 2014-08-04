clear all
% close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

stim_fs = 1e4/117.5;
Fs = 3e4;
% dsf = 120;Fsd = Fs/dsf;
dsf = 80;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1]/niqf,'high');
use_lfps = [1:8:96];
use_units = 1:96;
% use_lfps = [1 63];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);

scales = logspace(log10(8),log10(75),20);
scales = [scales 85 100 115 130 150 200 250];
scales = scales*60/dsf;
% scales = scales*2;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

%%
% loc_glob_blocks = [8 15 20 27 33 39 49] - 6;
stim_fs = 1e4/117.5;
grating_blocks = [12 34 52] - 6;
all_expt_id = [];
trial_start_times = [];
trial_stop_times = [];
trial_completed = [];
trial_expt_id = [];
for bb = 1:length(grating_blocks)
    n_trials(bb) = length(Expts{grating_blocks(bb)}.Trials);
    for i = 1:n_trials(bb)
        trial_start_times = [trial_start_times; Expts{grating_blocks(bb)}.Trials(i).Start(1)/1e4];
        trial_stop_times = [trial_stop_times; Expts{grating_blocks(bb)}.Trials(i).End(1)/1e4];
        trial_completed = [trial_completed; Expts{grating_blocks(bb)}.Trials(i).Result];
        trial_expt_id = [trial_expt_id; grating_blocks(bb)];
    end
end
trial_durs = trial_stop_times-trial_start_times;
min_trial_dur = 1;
%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
nearest_lfps = nan(length(use_lfps),1);
for ll = 1:96
   all_dists = sqrt((X_pos-X_pos(ll)).^2 + (Y_pos-Y_pos(ll)).^2); 
   all_dists(ll) = inf;
   [~,best_loc] = min(all_dists(use_lfps));
   nearest_lfps(ll) = best_loc;
end

%%
trial_cnt = 1;
desired_dt = 0.005;
all_interp_phasegrams = [];
all_interp_ampgrams = [];
all_binned_spks = [];
all_expt_t = [];
all_t_since_start = [];
all_trial_id = [];
for bb = 1:length(grating_blocks)
    
    %%
    fprintf('Analyzing block %d of %d\n',bb,length(grating_blocks));
    load(sprintf('Expt%dClusterTimes.mat',grating_blocks(bb)));
    
    filename = sprintf('Expt%dFullVmean.mat',grating_blocks(bb));
    load(filename);
    
    Vmat = [];
    %     Vmatf = [];
    phasegrams = [];
    ampgrams = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',grating_blocks(bb),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = V + FullV.sumscale*sumv;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dV = [];
        dVf = [];
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            dV = [dV decimate(V(cur_range),dsf)];
            dVf = [dVf filtfilt(filt_b,filt_a,decimate(V(cur_range),dsf))]; %do some high-pass filtering
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(ll,:) = dV;
        %         Vmatf(ll,:) = dVf;
        temp = cwt(dV,scales,'cmor1-1');
        phasegrams(:,:,ll) = angle(temp)';
        ampgrams(:,:,ll) = abs(temp)';
    end
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,2)+1:end) = [];
    
    cur_trials = find(trial_expt_id==grating_blocks(bb) & trial_durs>min_trial_dur);
    
    expt_t_axis = [];
    t_since_start = [];
    for tt = 1:length(cur_trials)
        cur_t_edges = trial_start_times(cur_trials(tt)):desired_dt:trial_stop_times(cur_trials(tt));
        binned_spks = nan(length(use_units),length(cur_t_edges)-1);
        for cc = 1:length(use_units)
            temp = histc(Clusters{use_units(cc)}.times,cur_t_edges);
            binned_spks(cc,:) = temp(1:end-1);
        end
        all_binned_spks = cat(1,all_binned_spks,binned_spks');
        all_trial_id = [all_trial_id; ones(length(cur_t_edges)-1,1)*trial_cnt];        
        trial_cnt = trial_cnt + 1;
        
        expt_t_axis = [expt_t_axis cur_t_edges(1:end-1)+desired_dt/2];
        t_since_start = [t_since_start cur_t_edges(1:end-1)-cur_t_edges(1)+desired_dt/2];
    end
    
    unwr_phasegram = unwrap(phasegrams);
    interp_phasegrams = interp1(t_ax,unwr_phasegram,expt_t_axis);
    interp_phasegrams = mod(interp_phasegrams+pi,2*pi)-pi;
    interp_ampgrams = interp1(t_ax,ampgrams,expt_t_axis);
    
    all_interp_phasegrams = cat(1,all_interp_phasegrams,interp_phasegrams);
    all_interp_ampgrams = cat(1,all_interp_ampgrams,interp_ampgrams);
    all_expt_t = cat(1,all_expt_t,expt_t_axis');
    all_t_since_start = cat(1,all_t_since_start,t_since_start');
end
%%
trial_set = unique(all_trial_id);
n_trials = length(trial_set);

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(all_trial_id,xv_set));
tr_inds = find(~ismember(all_trial_id,xv_set))';

%%
close all
NL_type = 0; %exp

reg_params.dl1_ind = 3000;
reg_params.dl1_dep = 0;
reg_params.dl2_ind = 30000;
reg_params.dl2_dep = 0;
reg_params.dl2_freq_ind = 300;
reg_params.dl2_freq_dep = 0;
reg_params.dl2_time_ind = 0;
reg_params.dl2_time_dep = 0;
reg_params.l1_ind = 0;
reg_params.l1_dep = 0;
reg_params.is_phase = 1;

% reg_params2.dl1_ind = 0;
% reg_params2.dl1_dep = 0;
% reg_params2.dl2_ind = 40000;
% reg_params2.dl2_dep = 40000;
% reg_params2.dl2_freq_ind = 20000;
% reg_params2.dl2_freq_dep = 20000;
% reg_params2.dl2_time_ind = 0;
% reg_params2.dl2_time_dep = 0;
% reg_params2.l1_ind = 0;
% reg_params2.l1_dep = 0;
% reg_params2.l1t_ind = 0;
% reg_params2.l1t_dep = 0;
% reg_params2.is_phase = 0;

silent = 1;
NT = length(all_expt_t);
used_inds = 1:NT;

nbins = 30;
pax = linspace(-pi,pi,nbins+1);
doub_pax = [pax(1:end-1) 2*pi+pax(1:end-1)];

new_phase_set = [reshape(cos(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps)) reshape(sin(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps))];
phase_elec_set = [repmat(1:length(use_lfps),1,length(wfreqs)) repmat(1:length(use_lfps),1,length(wfreqs))];
phase_elec_set = phase_elec_set(:);


for ll = 1:length(use_lfps)
    fprintf('Electrode %d of %d\n',ll,length(use_lfps));
    cur_units = find(nearest_lfps==ll);
    
    cur_alpha_phase = squeeze(all_interp_phasegrams(tr_inds,:,ll));
    Pmat = nan(length(tr_inds),length(wfreqs)*nbins);
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        Pmat(:,cur_set) = cur_tb;
    end
    cur_alpha_phase = squeeze(all_interp_phasegrams(xv_inds,:,ll));
    xvPmat = nan(length(xv_inds),length(wfreqs)*nbins);
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        xvPmat(:,cur_set) = cur_tb;
    end

    for cc = 1:length(cur_units)
        fprintf('Cell %d of %d\n',cc,length(cur_units));
        
        Robs = all_binned_spks(tr_inds,cur_units(cc));
        tr_spkbns = convert_to_spikebins(Robs);
        
        xv_Robs = all_binned_spks(xv_inds,cur_units(cc));

        stim_params = [nbins,length(wfreqs)];
        klen = size(Pmat,2);
        K0 = zeros(klen+1,1);
        [fitp_phase] = fit_GLM_phase_model(Pmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
        stim_ind_phase_pfilt(cur_units(cc),:) = fitp_phase.k(1:nbins*length(wfreqs));
        
        xv_phase_pred_rate = xvPmat*fitp_phase.k(1:end-1) + fitp_phase.k(end);
        if NL_type == 0
            xv_phase_pred_rate = log(1+exp(xv_phase_pred_rate));
        else
            xv_phase_pred_rate = exp(xv_phase_pred_rate);
        end
        xv_phase_LL(cur_units(cc)) = -sum(xv_Robs.*log(xv_phase_pred_rate)-xv_phase_pred_rate)/sum(xv_Robs);

        avg_rate = mean(Robs);
        null_pred = avg_rate*ones(size(xv_Robs));
        xv_null_LL(cur_units(cc)) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);       
        
    end
end

%%
load ./grating_40_phase_mod_data
grating_mod = stim_ind_phase_pfilt;
figure
hist(xv_phase_LL-xv_null_LL,100)

load ./noise_40_phase_mod_data
noise_mod = stim_ind_phase_pfilt;
figure
hist(xv_phase_LL-xv_null_LL,100)
%%
figure
for cc = 1:length(use_units)
    subplot(2,1,1)
    pcolor(pax(1:end-1),wfreqs,reshape(grating_mod(cc,:),nbins,length(wfreqs))');shading flat;colorbar
    subplot(2,1,2)
    pcolor(pax(1:end-1),wfreqs,reshape(noise_mod(cc,:),nbins,length(wfreqs))');shading flat;colorbar
    pause
    clf
end

%%
for cc = 1:length(use_units)
    pcolor(pax(1:end-1),wfreqs,reshape(stim_ind_phase_pfilt(cc,:),nbins,length(wfreqs))');shading flat
    pause
    clf
end

%%
for cc = 1:length(use_units)
kern = reshape(stim_ind_phase_pfilt(cc,:),nbins,length(wfreqs))';

    dom(cc,:) = (max(kern,[],2)-min(kern,[],2));

end

%%
save grating_40_phase_mod_data pax wfreqs nbins stim_ind* dom xv_*