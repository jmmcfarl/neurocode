%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt2_compiled_windata_d1p5.mat
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
el_pos = [X_pos(:) Y_pos(:)];
grid_ids = nan(10,10);
for j = 1:10
    for i = 1:10
        cur = find(X_pos==i & Y_pos==j);
        if ~isempty(cur)
            grid_ids(j,i) = cur;
        end
    end
end
use_lfps = grid_ids(:);
use_lfps = downsample(use_lfps,4);
use_lfps(isnan(use_lfps)) = [];

fix_win_dur = 0.15;

Pix2Deg = 0.018837;
NT = length(full_fix_ids);
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

Fs = 3e4;
% dsf = 150;Fsd = Fs/dsf;
dsf = 120;Fsd = Fs/dsf;
[b_alpha,a_alpha] = butter(2,[6 12]/(Fsd/2));
scales = logspace(log10(15),log10(250),30); 
% scales = [scales 300 400];
% scales = [scales 70 80 90 100 110 120 135 150 165 180 195 210];
scales = scales*30/dsf;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);
% use_lfps = [1:2:96];
% use_lfps = [1];

forwardlag = round(Fsd*0.6);
backlag = round(Fsd*0.2);
lags = -backlag:forwardlag;

% load ./expt2_parsed_g034 all_sac_times
%% eliminate image flip fixations
%used_fixs = 1:length(fix_is_sac);
used_fixs = find(fix_is_sac==1);
flip_fixs = find(fix_is_sac==0);
full_fix_ids(ismember(full_fix_ids,flip_fixs)) = nan;
full_fix_starts = full_fix_starts(used_fixs);
full_fix_wends = full_fix_wends(used_fixs);
full_fix_ends = full_fix_ends(used_fixs);
fix_start_inds = fix_start_inds(used_fixs);
fix_start_times = fix_start_times(used_fixs);
fix_end_inds = fix_end_inds(used_fixs);
fix_end_times = fix_end_times(used_fixs);
% 
% all_sac_start_times = all_sac_times(:,1);
% 
% fix_sac_times = nan(size(fix_start_times));
% for i = 1:length(fix_start_times)
%     prev_sac = find(all_sac_start_times < fix_start_times(i),1,'last');
%     fix_sac_times(i) = all_sac_start_times(prev_sac);
% end
% fix_start_times = fix_sac_times;
%%
%%
% Expt_nu = [13 14 15];
Expt_nu = [13 14 15 16 25 28 29];
% Expt_nu = [13];
n_expts = length(Expt_nu);

n_fixs = length(full_fix_wends);
fix_expt_num = nan(n_fixs,1);
for i = 1:n_fixs
    cur_inds = full_fix_starts(i):full_fix_wends(i);
    fix_expt_num(i) = unique(full_expt_vec(cur_inds));
end
use_fixs = find(fix_expt_num <= 3);
full_fix_starts = full_fix_starts(use_fixs);
full_fix_wends = full_fix_wends(use_fixs);
fix_start_times = fix_start_times(use_fixs);
fix_end_times = fix_end_times(use_fixs);
n_fixs = length(use_fixs);

% fix_start_log = zeros(size(full_t));
% fix_end_log = fix_start_log;
% fix_start_log(full_fix_starts) = 1;
% fix_end_log(full_fix_ends) = 1;

full_ampgrams = [];
full_abs_ampgrams = [];
% full_alpha_phase = [];
full_tax = [];
for ee = 1:n_expts
    fprintf('Expt %d of %d\n',ee,n_expts);
    clear ampgram all_V* alpha_phase
    filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(1));
    load(filename);
    V = double(FullV.V);
    V = decimate(V,dsf);
    V = V(:);
    vlen = length(V);
    V_alpha_all = nan(vlen,length(use_lfps));
    %     alpha_phase = nan(vlen,length(use_lfps));
    ampgram = nan(vlen,length(wfreqs),length(use_lfps));
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        
%         ord = 4;
%         acorr_seq = xcov(V,ord);
%         acorr_seq = acorr_seq(ord+1:end);
%         A = levinson(acorr_seq(1:ord),ord);
%         V = filter(-A,1,V);
        
        ampgram(:,:,ll) = cwt(V,scales,'cmor1-1')';
    end
    
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);

    full_ampgrams = [full_ampgrams; ampgram];
    full_abs_ampgrams = [full_abs_ampgrams; nanzscore(abs(ampgram))];
%     full_abs_ampgrams = [full_abs_ampgrams; abs(ampgram)];
    full_tax = [full_tax; t_ax(:)];
    clear ampgram
end
% full_abs_ampgrams = abs(full_ampgrams);
% full_abs_ampgrams = nanzscore(full_abs_ampgrams);

%%
cNT = size(full_ampgrams,1);
rel_fix_start_inds = round(interp1(full_tax,1:length(full_tax),fix_start_times));
rel_fix_end_inds = round(interp1(full_tax,1:length(full_tax),fix_end_times));

% big_sacs = find(fix_sac_amps(used_fixs) > 1);
% small_sacs = find(fix_sac_amps(used_fixs) < 1);

% rel_fix_start_inds = rel_fix_start_inds(big_sacs);
% rel_fix_end_inds = rel_fix_end_inds(big_sacs);
% rel_fix_start_inds = rel_fix_start_inds(small_sacs);
% rel_fix_end_inds = rel_fix_end_inds(small_sacs);

% use_inds = find(~isnan(full_fix_ids));
% rel_fix_start_inds = find(fix_start_log(use_inds) == 1);
% rel_fix_end_inds = find(fix_end_log(use_inds) == 1);
bad_fixs = find(rel_fix_start_inds <= backlag | rel_fix_end_inds >= cNT-forwardlag);
rel_fix_start_inds(bad_fixs) = [];
rel_fix_end_inds(bad_fixs) = [];

%%
% fix_trig_real = nan(n_fixs,length(lags),length(wfreqs),length(use_lfps)); 
fix_trig_anggram = zeros(length(lags),length(wfreqs),length(use_lfps));
% fix_trig_anggram = nan(n_fixs,length(lags),length(wfreqs),length(use_lfps));
fix_trig_ampgram = zeros(length(lags),length(wfreqs),length(use_lfps));
% fix_trig_ampgram = nan(n_fixs,length(lags),length(wfreqs),length(use_lfps));
fix_trig_arrayC = zeros(length(lags),length(wfreqs));
n_cnts = zeros(length(lags),1);
for ii = 1:n_fixs
    cur_inds = (rel_fix_start_inds(ii)-backlag):(rel_fix_start_inds(ii)+forwardlag);
    cur_inds(cur_inds > rel_fix_end_inds(ii)) = [];
    cl = length(cur_inds);
    fix_trig_ampgram(1:cl,:,:) = fix_trig_ampgram(1:cl,:,:) + full_abs_ampgrams(cur_inds,:,:);
    %     fix_trig_ampgram(ii,1:cl,:,:) = full_abs_ampgrams(cur_inds,:,:);
    
    %     fix_trig_real(ii,1:cl,:,:) = real(full_ampgrams(cur_inds,:,:));
    
    cur_anggram = full_ampgrams(cur_inds,:,:);
    cur_anggram = bsxfun(@rdivide,cur_anggram,abs(cur_anggram));
    fix_trig_anggram(1:cl,:,:) = fix_trig_anggram(1:cl,:,:) + cur_anggram;
    %     fix_trig_anggram(ii,1:cl,:,:) = cur_anggram;
    
    %     cur_array_coh = abs(squeeze(mean(cur_anggram,3)));
    %     fix_trig_arrayC(1:cl,:) = fix_trig_arrayC(1:cl,:) + cur_array_coh;
    
    n_cnts(1:cl) = n_cnts(1:cl) + 1;
end
fix_trig_ampgram = bsxfun(@rdivide,fix_trig_ampgram,n_cnts);
fix_trig_arrayC = bsxfun(@rdivide,fix_trig_arrayC,n_cnts);
fix_trig_anggram = bsxfun(@rdivide,fix_trig_anggram,n_cnts);
% fix_trig_anggram = angle(fix_trig_anggram);
fix_trig_phaselock = abs(fix_trig_anggram);

% n = repmat(n_cnts,[1 length(wfreqs) length(use_lfps)]);

d_p = 0.01;
thresh_r = (log(d_p) + 1 - 2*n_cnts).^2 - 1 - 4*n_cnts - 4*n_cnts.^2;
thresh_r = sqrt(thresh_r/4)./n_cnts;
thresh_r = repmat(thresh_r,[1 length(wfreqs)]);

avg_phaselock = squeeze(mean(fix_trig_phaselock,3));
avg_phaselock(avg_phaselock < thresh_r) = nan;

cd ~/Data/bruce/7_15_12/G034
save fv_fix_trig_phaselock wfreqs lags Fsd thresh_r avg_phaselock fix_trig_phaselock fix_trig_ampgram

% time_means = squeeze(nanmean(fix_trig_real));
% time_stds = squeeze(nanstd(fix_trig_real));
% fix_trig_real = reshape(fix_trig_real,[n_fixs*length(lags) length(wfreqs) length(use_lfps)]);
% ov_means = squeeze(nanmean(fix_trig_real));
% ov_stds = squeeze(nanstd(fix_trig_real));
% ov_means = repmat(shiftdim(ov_means,-1),[length(lags) 1 1]);
% ov_stds = repmat(shiftdim(ov_stds,-1),[length(lags) 1 1]);
% stim_gkl_div = 1./(2*ov_stds.^2).*((time_means - ov_means).^2 + ...
%         time_stds.^2 - ov_stds.^2) + log(ov_stds./time_stds);


%% COMPUTE WAVELET COHERENCE SPECTRUM
bad_lfps = [9 10 11 16 35 92];
cur_use_lfps = find(~ismember(use_lfps,bad_lfps));

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
use_lags = find(lags/Fsd > 0 & lags/Fsd < 0.5);
coh_mat2 = coh_mat(use_lags,:,cur_use_lfps,cur_use_lfps);
tavg_coh_mat = squeeze(nansum(bsxfun(@times,coh_mat2,n_cnts(use_lags))))/sum(n_cnts(use_lags));
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
xpos = X_pos(use_lfps(cur_use_lfps))*0.4;
ypos = Y_pos(use_lfps(cur_use_lfps))*0.4;
distmat = squareform(pdist([xpos(:) ypos(:)]));
distvec = reshape(distmat,length(cur_use_lfps)^2,1);
tavg_coh_r = reshape(tavg_coh_mat,[length(wfreqs) length(cur_use_lfps)^2]);
  
close all
n_bins = 15;
% dist_bin_edges = linspace(0.2,4.3,n_bins+1);
uset = find(distvec > 0);
dist_bin_edges = prctile(distvec(uset),linspace(0,100,n_bins+1));
dist_bin_edges(find(diff(dist_bin_edges) < 0.01)+1) = [];
dist_bins = 0.5*dist_bin_edges(1:end-1) + 0.5*dist_bin_edges(2:end);
n_bins = length(dist_bins);
bin_avg_coh = nan(n_bins,length(wfreqs));
for i = 1:n_bins
    cur_set = find(distvec >= dist_bin_edges(i) & distvec < dist_bin_edges(i+1));
    bin_avg_coh(i,:) = nanmean(tavg_coh_r(:,cur_set),2);
    nset(i) = length(cur_set);
end

% close all
% uset = find(distvec > 0);
% dist_bin_cents = unique(distvec(uset));
% nbins = length(dist_bin_cents);
% bin_avg_coh = nan(n_bins,length(wfreqs));
% for i = 1:n_bins
%     cur_set = find(distvec == dist_bin_cents(i));
%     bin_avg_coh(i,:) = nanmean(tavg_coh_r(:,cur_set),2);
%     nset(i) = length(cur_set);
% end

pcolor(dist_bins,wfreqs,bin_avg_coh');shading flat

ov_avg_coh = mean(tavg_coh_r(:,uset),2);
bin_avg_coh_norm = bsxfun(@rdivide,bin_avg_coh,ov_avg_coh');
% figure
% pcolor(dist_bins,wfreqs,bin_avg_coh_norm');shading flat

%%
% use_el = 1:48;
% use_el([5 6 18]) = [];
% use_el = 1:24;
% use_el([3 9]) = [];

el_dist = squareform(pdist(el_pos(use_lfps(use_el),:)));
el_dist = el_dist(:);
beta_0 = [0.5 2];
LB = [0 1];
UB = [1 100];

circ_dfun = @(xi,xj)(abs(mean(bsxfun(@minus,xi,xj),2)));

% circ_dfun = @(xi,xj)(mean(mod(bsxfun(@minus,xi,xj),pi),2));
% circ_dfun = @(xi,xj)(mean(angle(bsxfun(@rdivide,1i*xi,1i*xj)),2));

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

save coh_dist_data betas_coh wfreqs lags Fsd 

%%
% fix_trig_phaselock = abs(fix_trig_anggram);
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
