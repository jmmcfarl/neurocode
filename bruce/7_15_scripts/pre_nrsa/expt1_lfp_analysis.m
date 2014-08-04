
clear all
% close all
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/7_15_12

cd G034/
load ./CellList.mat
load ./G034Expts.mat

load ./expt1_eyedata_34_new.mat


samp_fac = 0.5;
dt = samp_fac*118/1e4;
fst = 1/dt;
frames_per_jump = 44/samp_fac;

Pix2Deg = 0.018837;
dsfrac = 1.25;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);

rf_cent = [0.34 -0.43];
min_trial_dur = 1;
max_sac_amp = 0.5;
single_units = find(CellList(1,:,1) > 0);
n_sus = length(single_units);

Fs = 3e4;
dsf = 60;Fsd = Fs/dsf;
scales = logspace(log10(2.5),log10(85),40);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);
use_lfps = [1:8:96];
% use_lfps = [1];
% [b,a] = butter(2,[5 12]/(Fsd/2));
[b,a] = butter(2,[1 100]/(Fsd/2));
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

    forwardlag = round(0.5*Fsd);
    backlag = round(0.1*Fsd);
    lags = -backlag:forwardlag;

%%
% Expt_nu = [3 8 15 19 26 29];
% Expt_nu = [13 14 15 16 25 28 29];
Expt_nu = [1 2 17:20 23 24];
% Expt_nu = [13];4
% Expt_nu = [7]; %expt3
n_allunits = 96;

    all_mtrig_ampgram = zeros(length(lags),length(wfreqs),length(use_lfps));
    all_mtrig_V = zeros(length(lags),length(use_lfps));
    all_mn_cnts = zeros(length(lags),1);
    all_motrig_ampgram = zeros(length(lags),length(wfreqs),length(use_lfps));
    all_motrig_V = zeros(length(lags),length(use_lfps));
    all_mon_cnts = zeros(length(lags),1);

full_expt_vec = [];
full_trial_vec = [];
full_binned_spks = [];
full_trial_durs = [];
full_t = [];
full_im_flips = [];
full_alpha_phase = [];
for ee = 1:length(Expt_nu)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    n_mus = 96;
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    
    Trial_starts = [Expts{Expt_nu(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{Expt_nu(ee)}.Trials(:).End]/1e4;
            
    endEvents = [Expts{Expt_nu(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
%     used_trials = find(Trial_durs > min_trial_dur);
used_trials = 1:length(Trial_durs);    
    tt = nan(length(used_trials),1);
    for i = 1:length(used_trials)
        cur_n_tbins = floor(Trial_durs(used_trials(i))/dt);
        cur_n_shifts = ceil(cur_n_tbins/frames_per_jump);
        tt(i) = cur_n_shifts;
        
        cur_t_edges = Trial_starts(used_trials(i)):dt:(Trial_starts(used_trials(i)) + dt*cur_n_tbins);
        cur_t_cents = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_im_flips = Trial_starts(used_trials(i)):dt*frames_per_jump:Trial_ends(used_trials(i));
        
        cur_binned_spks = nan(n_mus,length(cur_t_cents));
        for j = 1:n_mus
            temp = histc(Clusters{j}.times,cur_t_edges);
            cur_binned_spks(j,:) = temp(1:end-1);
        end
        
        full_expt_vec = [full_expt_vec; ones(length(cur_t_cents),1)*Expt_nu(ee)];
        full_trial_vec = [full_trial_vec; ones(length(cur_t_cents),1)*used_trials(i)];
        full_binned_spks = [full_binned_spks; cur_binned_spks'];
        full_trial_durs = [full_trial_durs; ones(length(cur_t_cents),1)*Trial_durs(used_trials(i))];
        full_t = [full_t cur_t_cents];
        full_im_flips = [full_im_flips cur_im_flips];
    end
    
    ampgram = [];
    phasegram = [];
    all_V = [];
    all_alpha_phase = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V*FullV.intscale(1)/FullV.intscale(2);
        V = V(:);
%         cur_cwt = cwt(V,scales,'cmor1-1')';
%         ampgram(:,:,ll) = abs(cur_cwt);
%         phasegram(:,:,ll) = angle(cur_cwt);
%         V = filtfilt(b,a,V);
        all_V(:,ll) = V;
%         cur_alpha = filtfilt(b,a,V);
%         all_alpha_phase(:,ll) = angle(hilbert(cur_alpha));
    end
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
%     ampgram = nanzscore(ampgram);
    all_V = nanzscore(all_V);
    cur_set = find(full_expt_vec == Expt_nu(ee));
%     interp_alpha_phase = interp1(t_ax,all_alpha_phase,full_t(cur_set));
%     full_alpha_phase = [full_alpha_phase; interp_alpha_phase];
% end


%     binned_spks = nan(length(t_ax),length(use_lfps));
%     smoothed_spks = nan(length(t_ax),length(use_lfps));
%     for cc = 1:length(use_lfps)
%         binned_spks(:,cc) = hist(Clusters{use_lfps(cc)}.times,t_ax);
%         smoothed_spks(:,cc) = smooth(binned_spks(:,cc),spk_smwin);
%     end
%     smoothed_spks = zscore(smoothed_spks);

    %%
% %     sm_fac = 5;
%     clear cur_binned_spks
%     for j = 1:n_mus
%         temp = histc(Clusters{j}.times,t_ax);
% %         cur_binned_spks(j,:) = smooth(temp,sm_fac);
%         cur_binned_spks(j,:) = temp;
%     end
% %     cur_binned_spks = nanzscore(cur_binned_spks');
    
    Trial_start_inds = round(interp1(t_ax,1:length(t_ax),Trial_starts));
    Trial_end_inds = round(interp1(t_ax,1:length(t_ax),Trial_ends));
    bad_trials = find(Trial_start_inds < backlag | Trial_end_inds > length(t_ax) - forwardlag | ...
        isnan(Trial_start_inds) | isnan(Trial_end_inds));
    Trial_start_inds(bad_trials) = [];
    Trial_end_inds(bad_trials) = [];
    Trial_starts(bad_trials) = [];
    Trial_ends(bad_trials) = [];
    
    intrial_inds = [];
    for i =1 :length(Trial_start_inds)
        intrial_inds = [intrial_inds Trial_start_inds(i):Trial_end_inds(i)];
    end
    outtrial_inds = setdiff(1:length(t_ax),intrial_inds);
    
    
    cur_msac_times = [];
    cur_msac_amps = [];
    for i = 1:length(Trial_starts)
        cur_set = find(all_sac_start_times >= Trial_starts(i) & all_sac_start_times <= Trial_ends(i));
        cur_msac_times = [cur_msac_times all_sac_stop_times(cur_set)];
        cur_msac_amps = [cur_msac_amps all_sac_amps(cur_set)];
    end
    cur_mosac_times = [];
    cur_mosac_amps = [];
    for i = 1:length(Trial_starts)-1
        cur_set = find(all_sac_start_times >= Trial_ends(i) & all_sac_start_times <= Trial_starts(i+1));
        cur_mosac_times = [cur_mosac_times all_sac_stop_times(cur_set)];
        cur_mosac_amps = [cur_mosac_amps all_sac_amps(cur_set)];
    end
    cur_msac_inds = round(interp1(t_ax,1:length(t_ax),cur_msac_times));
    cur_msac_inds(cur_msac_inds < backlag | cur_msac_inds > length(t_ax) - forwardlag) = [];
    cur_mosac_inds = round(interp1(t_ax,1:length(t_ax),cur_mosac_times));
    cur_mosac_inds(cur_mosac_inds < backlag | cur_mosac_inds > length(t_ax) - forwardlag) = [];
    
%     trig_ampgram = zeros(length(lags),length(wfreqs),length(use_lfps));
%     trig_V = zeros(length(lags),length(use_lfps));
%     n_cnts = zeros(length(lags),1);
%     for i = 1:length(Trial_start_inds)
%         cur_inds = (Trial_start_inds(i)-backlag):(Trial_start_inds(i)+forwardlag);
%         cur_inds(cur_inds > Trial_end_inds(i)) = [];
%         cl = length(cur_inds);
%         trig_ampgram(1:cl,:,:) = ampgram(cur_inds,:,:);
%         trig_V(1:cl,:) = all_V(cur_inds,:);
%         n_cnts(1:cl) = n_cnts(1:cl) + 1;
%     end
%     trig_ampgram = bsxfun(@rdivide,trig_ampgram,n_cnts);
%     trig_V = bsxfun(@rdivide,trig_V,n_cnts);

    mtrig_ampgram = zeros(length(lags),length(wfreqs),length(use_lfps));
    mtrig_V = zeros(length(lags),length(use_lfps));
    mn_cnts = zeros(length(lags),1);
    for i = 1:length(cur_msac_inds)
        cur_inds = (cur_msac_inds(i)-backlag):(cur_msac_inds(i)+forwardlag);
        next_tend = find(Trial_end_inds > cur_msac_inds(i),1,'first');
        if ~isempty(next_tend)
        cur_inds(cur_inds > Trial_end_inds(next_tend)) = [];
        end
        if i < length(cur_msac_inds)
            cur_inds(cur_inds > cur_msac_inds(i+1)) = [];
        end
        cl = length(cur_inds);
%         mtrig_ampgram(1:cl,:,:) = ampgram(cur_inds,:,:);
        mtrig_V(1:cl,:) = mtrig_V(1:cl,:)+all_V(cur_inds,:);
        mn_cnts(1:cl) = mn_cnts(1:cl) + 1;
    end
%     mtrig_ampgram = bsxfun(@rdivide,mtrig_ampgram,mn_cnts);
%     mtrig_V = bsxfun(@rdivide,mtrig_V,mn_cnts);

%         motrig_ampgram = zeros(length(lags),length(wfreqs),length(use_lfps));
    motrig_V = zeros(length(lags),length(use_lfps));
    mon_cnts = zeros(length(lags),1);
    for i = 1:length(cur_mosac_inds)
        cur_inds = (cur_mosac_inds(i)-backlag):(cur_mosac_inds(i)+forwardlag);
        next_tstart = find(Trial_start_inds > cur_mosac_inds(i),1,'first');
        if ~isempty(next_tstart)
        cur_inds(cur_inds > Trial_start_inds(next_tstart)) = [];
        end
        if i < length(cur_mosac_inds)
            cur_inds(cur_inds > cur_mosac_inds(i+1)) = [];
        end
        cl = length(cur_inds);
%         motrig_ampgram(1:cl,:,:) = ampgram(cur_inds,:,:);
        motrig_V(1:cl,:) = motrig_V(1:cl,:)+all_V(cur_inds,:);
        mon_cnts(1:cl) = mon_cnts(1:cl) + 1;
    end
%     motrig_ampgram = bsxfun(@rdivide,motrig_ampgram,mon_cnts);
%     motrig_V = bsxfun(@rdivide,motrig_V,mon_cnts);

%     all_trig_V(ee,:,:) = trig_V;
%     all_trig_ampgram(ee,:,:,:) = trig_ampgram;
    all_mtrig_V = all_mtrig_V + mtrig_V;
    all_mn_cnts = all_mn_cnts + mn_cnts;
%     all_mtrig_ampgram(ee,:,:,:) = mtrig_ampgram;

    %     all_trig_ampgram(ee,:,:,:) = trig_ampgram;
    all_motrig_V = all_motrig_V + motrig_V;
    all_mon_cnts = all_mon_cnts + mon_cnts;
%     all_mtrig_ampgram(ee,:,:,:) = mtrig_ampgram;

%      trial_durs_on = (Trial_end_inds-Trial_start_inds)/Fsd;
%      uset = find(trial_durs_on > 1);
%    sMarkers_on = [Trial_start_inds(uset)' Trial_end_inds(uset)'];
%     
%       trial_durs_off = (Trial_start_inds(2:end)-Trial_end_inds(1:end-1))/Fsd;
%      uset = find(trial_durs_off > 1);
%      temp_end = Trial_end_inds(1:end-1); temp_start = Trial_start_inds(2:end);
%    sMarkers_off = [temp_end(uset)' temp_start(uset)'];
%     params.Fs = Fsd;
%     params.tapers = [2 3];
%     movingwin = [1 1];
%     clear Son Soff
%     for ll =1 :length(use_lfps)
%         [ Son(:,ll), f ]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers_on );
%         [ Soff(:,ll), f ]= mtspectrumc_unequal_length_trials( all_V(:,ll), movingwin, params, sMarkers_off );
%     end
%     
%     params.tapers = [3 5];
%     [Cmn_on,Phimn,Smn,Smm,f] = coherencyc_unequal_length_trials( all_V, movingwin, params, sMarkers_on );
%     [Cmn_off,Phimn,Smn,Smm,f] = coherencyc_unequal_length_trials( all_V, movingwin, params, sMarkers_off);
%     
%     all_Son(ee,:,:) = Son;
%     all_Soff(ee,:,:) = Soff;
%     all_Cmn_on(ee,:,:) = Cmn_on;
%     all_cmn_off(ee,:,:) = Cmn_off;
%     
%     phase_lock = nan(96,length(wfreqs));
%     pref_phase = nan(96,length(wfreqs));
%     for ll = 1:96
%         cur_spikebins = convert_to_spikebins(cur_binned_spks(ll,:)');
%        for ww = 1:length(wfreqs)
%             cur_phases = squeeze(phasegram(cur_spikebins,ww,nearest_lfps(ll)));
%             phase_lock(ll,ww) = circ_kappa(cur_phases);
%             pref_phase(ll,ww) = circ_mean(cur_phases);
%        end
%     end
%     
%     all_phase_lock(ee,:,:) = phase_lock';
%     all_pref_phase(ee,:,:) = pref_phase';
end

all_motrig_V = bsxfun(@rdivide,all_motrig_V,all_mon_cnts);
all_mtrig_V = bsxfun(@rdivide,all_mtrig_V,all_mn_cnts);

%%
% for i = 1:48
%     shadedErrorBar(lags/Fsd,squeeze(nanmean(all_mtrig_V(:,:,i))),squeeze(nanstd(all_mtrig_V(:,:,i)))/sqrt(8))
%     pause
%     clf
% end
%%
close all
figure
shadedErrorBar(lags/Fsd,nanmean(all_mtrig_V,2),nanstd(all_mtrig_V,[],2),{'color','g'})
% hold on
% shadedErrorBar(lags/Fsd,nanmean(all_motrig_V,2),nanstd(all_motrig_V,[],2),{'color','r'})

%%
expt1_lags = lags/Fsd;
save expt1_trig_avgs expt1_lags all_mtrig_V all_mn_cnts
%%
    load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
    X_pos = ArrayConfig.X;
    Y_pos = ArrayConfig.Y;

    trial_start_inds = 1+find(diff(full_trial_vec)~=0);
    n_trials = length(trial_start_inds);
    forwardlag = round(0.45/dt);
phasemat = nan(n_trials,forwardlag+1,length(use_lfps));
for i = 1:n_trials
    cur_inds = trial_start_inds(i):(trial_start_inds(i)+forwardlag);
    cl = length(cur_inds);
    phasemat(i,1:cl,:) = full_alpha_phase(cur_inds,:);
end

elec_dist = pdist2([X_pos(use_lfps)' Y_pos(use_lfps)'],[X_pos(use_lfps)' Y_pos(use_lfps)']);
n_pairs = length(use_lfps)*(length(use_lfps)-1)/2;
pair_relphase = nan(n_trials,forwardlag+1,n_pairs);
cnt = 1;
for i = 1:length(use_lfps)-1
    for j = i+1:length(use_lfps)
       pair_relphase(:,:,cnt) = phasemat(:,:,i)-phasemat(:,:,j);
       for ll = 1:forwardlag
           uset = find(~isnan(pair_relphase(:,ll,cnt)));
          rel_kappa(cnt,ll) = circ_r(squeeze(pair_relphase(uset,ll,cnt))); 
          cur_elec_dist(cnt) = elec_dist(i,j);
       end
       cnt = cnt + 1;
    end
end
[~,dist_ord] = sort(cur_elec_dist);

%%
    % save expt1_lfp_analysis all_* f wfreqs 

%%
close all
pow_diff = log(all_Son)-log(all_Soff);
for i = 1:length(use_lfps) 
shadedErrorBar(f,squeeze(mean(pow_diff(:,:,i))),squeeze(std(pow_diff(:,:,i)))/sqrt(8))
set(gca,'xscale','log')
xlim(f([2 end]))
xl = f([2 end])
line(xl,[0 0],'color','r')
pause
clf
end

%%
C_diff = log(all_Cmn_on)-log(all_cmn_off);
for i = 1:size(C_diff,3)
shadedErrorBar(f,squeeze(mean(C_diff(:,:,i))),squeeze(std(C_diff(:,:,i)))/sqrt(8))
set(gca,'xscale','log')
xlim(f([2 end]))
xl = f([2 end])
line(xl,[0 0],'color','r')
pause
clf
end

%%
for i = 1:96
    subplot(2,1,1)
    shadedErrorBar(wfreqs,mean(all_phase_lock(:,:,i)),std(all_phase_lock(:,:,i))/sqrt(8))
    set(gca,'xscale','log')
    xlim(wfreqs([end 1]))
    xl = wfreqs([end 1])
    % line(xl,[0 0],'color','r')
    subplot(2,1,2)
    cur_mean = zeros(length(wfreqs),1);
    for t = 1:length(wfreqs)
        cur_mean(t) = circ_mean(squeeze(all_pref_phase(:,t,i)));
    end
    plot(wfreqs,cur_mean)
     set(gca,'xscale','log')
    xlim(wfreqs([end 1]))
    xl = wfreqs([end 1])
   
    pause
    clf
end

%%
use_el = 1:24;
use_el([3 9]) = [];

load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
el_pos = [X_pos(:) Y_pos(:)];
el_dist = squareform(pdist(el_pos(use_lfps(use_el),:)));
el_dist = el_dist(:);
beta_0 = [0.5 2];
LB = [0 1];
UB = [1 100];

betas_coh = nan(length(f),2);
for ww = 1:length(f)
    id = 1;
    cur_corrmat = nan(24,24);
    cur_corrmat(1,1) = 1;
    for j = 2:24
        for k = 1:j-1
            cur_corrmat(j,k) = avg_Cmn_on(ww,id);
            cur_corrmat(k,j) = avg_Cmn_on(ww,id);
            id = id + 1;
        end
        cur_corrmat(j,j) = 1;
    end
    cur_corrmat = cur_corrmat(use_el,use_el);
    betas_coh(ww,:) = lsqnonlin(@(X) cur_corrmat(:)-expon_decay_fun(X,el_dist),beta_0,LB,UB);
end
