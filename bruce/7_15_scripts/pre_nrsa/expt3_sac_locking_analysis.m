%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

% load ./Expt2_compiled_windata_d1p5.mat
% resh_all_stims = resh_all_stims/std(resh_all_stims(:));
load ./Expt3_fixbased_data.mat

fix_win_dur = 0.15;

Pix2Deg = 0.018837;
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

nat_set = 1:685;
white_set = 686:913;
sim_set = 914:1141;
obj_set = 1142:1369;
noise_set = [white_set sim_set];

Fs = 3e4;
dsf = 150;Fsd = Fs/dsf;
[b_alpha,a_alpha] = butter(2,[6 12]/(Fsd/2));
scales = logspace(log10(20),log10(250),25);
% scales = [scales 70 80 90 100 110 120 135 150 165 180 195 210];
scales = scales*30/dsf;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);
use_lfps = [1:2:96];

forwardlag = round(Fsd*0.5);
backlag = round(Fsd*0.5);
lags = -backlag:forwardlag;

%% eliminate image flip fixations
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
% Expt_nu = [13 14 15 16 25 28 29];
Expt_nu = [7:12]; %expt 3 34
n_expts = length(Expt_nu);

full_ampgrams = [];
full_tax = [];
full_V = [];
for ee = 1:n_expts
    fprintf('Expt %d of %d\n',ee,n_expts);
    clear ampgram all_V* alpha_phase
    filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(1));
    load(filename);
    V = double(FullV.V);
    V = decimate(V,dsf);
    V = V(:);
    vlen = length(V);
    V_all = nan(vlen,length(use_lfps));
    %     alpha_phase = nan(vlen,length(use_lfps));
    ampgram = nan(vlen,length(wfreqs),length(use_lfps));
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        V_all(:,ll) = V;
        ampgram(:,:,ll) = cwt(V,scales,'cmor1-1')';
    end
    
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    full_ampgrams = [full_ampgrams; ampgram];
    full_tax = [full_tax; t_ax(:)];
    full_V = [full_V; V_all];
end
full_abs_ampgrams = abs(full_ampgrams);
full_abs_ampgrams = nanzscore(full_abs_ampgrams);

%%
% cd /home/james/Data/bruce/7_15_12/G034
% save expt3_sac_locking_analysis full_ampgrams full_tax full_V full_abs_ampgrams wfreqs lags Fsd

%%
cNT = length(full_tax);
rel_fix_start_inds = round(interp1(full_tax,1:length(full_tax),trial_start_times));
rel_fix_end_inds = round(interp1(full_tax,1:length(full_tax),trial_stop_times));

use = find(diff(trial_imnum) ~= 0)+1;
use2 = find(trial_imnum(2:end)-trial_imnum(1:end-1) == 0)+1;
use2(~ismember(use2,use+1)) = [];
use3 = find(trial_imnum(2:end)-trial_imnum(1:end-1) == 0)+1;
use3(~ismember(use3,use+2)) = [];
use4 = find(trial_imnum(2:end)-trial_imnum(1:end-1) == 0)+1;
use4(~ismember(use4,use+3)) = [];
usea = [use2; use3; use4];

noise_trials = find(ismember(trial_imnum,noise_set));
nat_trials = find(ismember(trial_imnum,nat_set));
obj_trials = find(ismember(trial_imnum,obj_set));

rel_fix_start_inds = rel_fix_start_inds(usea);
rel_fix_end_inds = rel_fix_end_inds(usea);

bad_fixs = find(rel_fix_start_inds <= backlag | rel_fix_end_inds >= cNT-forwardlag);
rel_fix_start_inds(bad_fixs) = [];
rel_fix_end_inds(bad_fixs) = [];
%%
n_fixs = length(rel_fix_end_inds);
fix_trig_anggram = zeros(length(lags),length(wfreqs),length(use_lfps)); 
% fix_trig_anggram = nan(n_fixs,length(lags),length(wfreqs),length(use_lfps)); 
fix_trig_ampgram = nan(n_fixs,length(lags),length(wfreqs),length(use_lfps)); 
% fix_trig_ampgram = zeros(length(lags),length(wfreqs),length(use_lfps)); 
fix_trig_arrayC = zeros(length(lags),length(wfreqs)); 
n_cnts = zeros(length(lags),1);
for ii = 1:n_fixs
    cur_inds = (rel_fix_start_inds(ii)-backlag):(rel_fix_start_inds(ii)+forwardlag);
%     cur_inds(cur_inds > rel_fix_end_inds(ii)) = [];
    cl = length(cur_inds);
%     fix_trig_ampgram(1:cl,:,:) = fix_trig_ampgram(1:cl,:,:) + full_abs_ampgrams(cur_inds,:,:);
    fix_trig_ampgram(ii,1:cl,:,:) = full_abs_ampgrams(cur_inds,:,:);
    
    cur_anggram = full_ampgrams(cur_inds,:,:);
    cur_anggram = bsxfun(@rdivide,cur_anggram,abs(cur_anggram));
    fix_trig_anggram(1:cl,:,:) = fix_trig_anggram(1:cl,:,:) + cur_anggram;    
%     fix_trig_anggram(ii,1:cl,:,:) = cur_anggram;    
    
    cur_array_coh = abs(squeeze(mean(cur_anggram,3)));
    fix_trig_arrayC(1:cl,:) = fix_trig_arrayC(1:cl,:) + cur_array_coh;    

    n_cnts(1:cl) = n_cnts(1:cl) + 1;
end
% fix_trig_ampgram = bsxfun(@rdivide,fix_trig_ampgram,n_cnts);
fix_trig_arrayC = bsxfun(@rdivide,fix_trig_arrayC,n_cnts);
% fix_trig_anggram = bsxfun(@rdivide,fix_trig_anggram,n_cnts);
fix_trig_anggram = angle(fix_trig_anggram);
% fix_trig_phaselock = abs(fix_trig_anggram);

%% COMPUTE WAVELET COHERENCE SPECTRUM
n_fixs = length(rel_fix_end_inds);
fix_trig_mat = zeros(length(lags),length(wfreqs),length(use_lfps),length(use_lfps)); 
n_cnts = zeros(length(lags),1);
for ii = 1:n_fixs
    cur_inds = (rel_fix_start_inds(ii)-backlag):(rel_fix_start_inds(ii)+forwardlag);
    cur_inds(cur_inds > rel_fix_end_inds(ii)) = [];
    cl = length(cur_inds);
    temp =  full_ampgrams(cur_inds,:,:);
    for cc = 1:length(use_lfps)
        cur_mat = bsxfun(@times,temp,conj(temp(:,:,cc)));
        fix_trig_mat(1:cl,:,:,cc) = fix_trig_mat(1:cl,:,:,cc) + cur_mat;
    end
    n_cnts(1:cl) = n_cnts(1:cl) + 1;
end
fix_trig_mat = bsxfun(@rdivide,fix_trig_mat,n_cnts);

coh_mat = fix_trig_mat;
avg_cohmat = zeros(length(lags),length(wfreqs));
cnt = 0;
for ii = 1:length(use_lfps)
    for jj = 1:length(use_lfps)
        if ii ~= jj
            coh_mat(:,:,ii,jj) = abs(fix_trig_mat(:,:,ii,jj))./sqrt(fix_trig_mat(:,:,ii,ii).*fix_trig_mat(:,:,jj,jj));
            avg_cohmat = avg_cohmat + squeeze(coh_mat(:,:,ii,jj));
            cnt = cnt + 1;
        else
            coh_mat(:,:,ii,jj) = 1;
        end
    end
end
avg_cohmat = avg_cohmat/cnt;

%%
use_el = 1:48;
use_el([5 6 18]) = [];

load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
el_pos = [X_pos(:) Y_pos(:)];
el_dist = squareform(pdist(el_pos(use_lfps(use_el),:)));
el_dist = el_dist(:);
beta_0 = [0.5 2];
LB = [0 1];
UB = [1 100];

% circ_dfun = @(xi,xj)(abs(mean(bsxfun(@minus,xi,xj),2)));
% circ_dfun = @(xi,xj)(mean(mod(bsxfun(@minus,xi,xj),pi),2));
circ_dfun = @(xi,xj)(mean(angle(bsxfun(@rdivide,1i*xi,1i*xj)),2));
% angle(exp(1i*x)./exp(1i*y))
% betas = nan(length(wfreqs),length(lags),2);
betas_coh = nan(length(wfreqs),length(lags),2);
% betas_phase = nan(length(wfreqs),length(lags),2);
for ww = 1:length(wfreqs)
    for tt = 1:length(lags);
        cur_corrmat = squeeze(coh_mat(tt,ww,use_el,use_el));

        %         cur_set = squeeze(fix_trig_ampgram(:,tt,ww,use_el));
%         cur_set = squeeze(fix_trig_anggram(:,tt,ww,use_el));

%         cur_set(isnan(cur_set(:,1)),:) = [];
%         cur_corrmat = corr(cur_set);
%         cur_corrmat = squareform(pdist(cur_set',@(Xi,Xj) circ_dfun(Xi,Xj)));
%         cur_corrmat = 1-cur_corrmat/pi;
%         betas(ww,tt,:) = lsqnonlin(@(X) cur_corrmat(:)-expon_decay_fun(X,el_dist),beta_0,LB,UB);
%         betas_phase(ww,tt,:) = lsqnonlin(@(X) cur_corrmat(:)-expon_decay_fun(X,el_dist),beta_0,LB,UB);
        betas_coh(ww,tt,:) = lsqnonlin(@(X) cur_corrmat(:)-expon_decay_fun(X,el_dist),beta_0,LB,UB);
        end
    ww
end
%%
pvals = nan(length(use_lfps),length(lags),length(wfreqs));
zvals = pvals;
for ch = 1:length(use_lfps);
    for ll = 1:length(lags)
        for ww = 1:length(wfreqs)
            cur_set = squeeze(fix_trig_anggram(:,ll,ww,ch));
            cur_set = cur_set(~isnan(cur_set));
            [pvals(ch,ll,ww),zvals(ch,ll,ww)] = circ_rtest(cur_set);
        end
    end
end
