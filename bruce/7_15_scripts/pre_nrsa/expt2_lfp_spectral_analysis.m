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
use_lfps = [1:4:96];

Fs = 3e4;
dsf = 60;Fsd = Fs/dsf;
scales = logspace(log10(2.5),log10(85),40);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);
use_lfps = [1:4:96];
[b,a] = butter(2,[5 12]/(Fsd/2));

alpha_dp = 0.15;

forwardlag = round(Fsd*0.6);
backlag = round(Fsd*0.4);
lags = -backlag:forwardlag;

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
 cd ~/Data/bruce/7_15_12/G034/
% Expt_nu = [3 8 15 19 26 29];
Expt_nu = [13 14 15 16 25 28 29];
n_allunits = 96;
for ee = 1:length(Expt_nu)
   
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
        
    filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(1));
    load(filename);
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);

    %%
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
        V = V(:);
        cur_cwt = cwt(V,scales,'cmor1-1')';
        ampgram(:,:,ll) = abs(cur_cwt);
%         phasegram(:,:,ll) = angle(cur_cwt);
        all_V(:,ll) = V;
%         cur_alpha = filtfilt(b,a,V);
%         all_alpha_phase(:,ll) = angle(hilbert(cur_alpha));
    end
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
%         ampgram = nanzscore(ampgram);

    %%
    cur_fixs = find(fix_expt_num == ee);
    cur_n_fixs = length(cur_fixs);
    
    fix_start_inds = round(interp1(t_ax,1:length(t_ax),full_t(full_fix_starts(cur_fixs))));
    fix_stop_inds = round(interp1(t_ax,1:length(t_ax),full_t(full_fix_ends(cur_fixs))));

    %%
    trig_ampgram = zeros(length(lags),length(wfreqs),length(use_lfps));
    trig_V = zeros(length(lags),length(use_lfps));
    n_cnts = zeros(length(lags),1);
    for i = 1:length(fix_start_inds)
        cur_inds = (fix_start_inds(i)-backlag):(fix_start_inds(i)+forwardlag);
        cur_inds(cur_inds > fix_stop_inds(i)) = [];
        cl = length(cur_inds);
        trig_ampgram(1:cl,:,:) = trig_ampgram(1:cl,:,:) + ampgram(cur_inds,:,:);
        trig_V(1:cl,:) = trig_V(1:cl,:) + all_V(cur_inds,:);
        n_cnts(1:cl) = n_cnts(1:cl) + 1;
    end
    trig_ampgram = bsxfun(@rdivide,trig_ampgram,n_cnts);
    trig_V = bsxfun(@rdivide,trig_V,n_cnts);
        
    
end

%%
