%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt2_compiled_windata_d1p5.mat
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

fix_win_dur = 0.15;

Pix2Deg = 0.018837;
NT = length(full_fix_ids);
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

Fs = 3e4;
dsf1 = 50;
dsf2 = 5;
dsf = dsf1*dsf2;
Fsd = Fs/dsf;
Fsd1 = Fs/dsf1;
[b_alpha,a_alpha] = butter(2,[6 12]/(Fsd/2));
scales = logspace(log10(5),log10(250),40);
% scales = [scales 70 80 90 100 110 120 135 150 165 180 195 210];
scales = scales*30/dsf1;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd1);
nwfreqs = length(wfreqs);
use_lfps = [1:4:96];

forwardlag = round(Fsd*0.6);
backlag = round(Fsd*0.2);
lags = -backlag:forwardlag;

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
%%
n_fixs = length(full_fix_wends);
fix_expt_num = nan(n_fixs,1);
for i = 1:n_fixs
    cur_inds = full_fix_starts(i):full_fix_wends(i);
    fix_expt_num(i) = unique(full_expt_vec(cur_inds));
end
%%
Expt_nu = [13 14 15 16 25 28 29];
n_expts = length(Expt_nu);

% fix_start_log = zeros(size(full_t));
% fix_end_log = fix_start_log;
% fix_start_log(full_fix_starts) = 1;
% fix_end_log(full_fix_ends) = 1;

full_ampgrams = [];
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
        V = decimate(V,dsf1);
        V = V(:);
        cur_ampgram = abs(cwt(V,scales,'cmor1-1')');
        ampgram(:,:,ll) = downsample(cur_ampgram,dsf2);
    end
    
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);

    full_ampgrams = [full_ampgrams; nanzscore(ampgram)];
    full_tax = [full_tax; t_ax(:)];
end

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
n_fixs = length(rel_fix_end_inds);
fix_trig_ampgram = zeros(length(lags),length(wfreqs),length(use_lfps)); 
n_cnts = zeros(length(lags),1);
for ii = 1:n_fixs
    cur_inds = (rel_fix_start_inds(ii)-backlag):(rel_fix_start_inds(ii)+forwardlag);
    cur_inds(cur_inds > rel_fix_end_inds(ii)) = [];
    cl = length(cur_inds);
    fix_trig_ampgram(1:cl,:,:) = fix_trig_ampgram(1:cl,:,:) + full_ampgrams(cur_inds,:,:);
%     fix_trig_ampgram(ii,1:cl,:,:) = full_abs_ampgrams(cur_inds,:,:);
    
    n_cnts(1:cl) = n_cnts(1:cl) + 1;
end
fix_trig_ampgram = bsxfun(@rdivide,fix_trig_ampgram,n_cnts);

