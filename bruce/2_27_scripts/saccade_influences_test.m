clear all
% close all

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-5 5;-5 5];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;
% nlags = 4;


Pix2Deg = 0.018837;

% down-sampling fraction for image
dsfrac = 2;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% [XX,YY] = meshgrid(xax,yax);

% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.75 6.25; -3.5 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

cd /Users/James/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_d2
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

SDIM = length(xpatch_inds);
cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_image_patches_corrected

%%
NT = length(used_fixs);
X_resh = reshape(X,NT,SDIM^2);
dt = .005;
spikes_binned = [];
time_since_fix = [];
stim_outs = [];

%%
n_used_cells = size(spk_cnts,2);

%%
fix_stimout = zeros(n_used_cells,length(used_fixs));
all_model_time_axis = [];
all_model_blockids = [];
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    cur_set = find(blockids(used_fixs) == blockid);
    
    for i = 1:length(cur_set)
        start_T = all_fix_start_times(used_fixs(cur_set(i)));
        end_T = all_fix_stop_times(used_fixs(cur_set(i)));
        cur_tedges = start_T:dt:end_T;
        cur_tcents = 0.5*cur_tedges(1:end-1)+0.5*cur_tedges(2:end);
        temp_binned = zeros(10,length(cur_tcents));
        temp_modouts = zeros(10,length(cur_tcents));
        for c = 1:n_used_cells
            temp = histc(Blocks{blockid}.spktimes{c},cur_tedges);
            temp_binned(c,:) = temp(1:end-1);
        end
        spikes_binned = [spikes_binned; temp_binned'];
        
        time_since_fix = [time_since_fix cur_tcents-start_T];
        all_model_time_axis = [all_model_time_axis cur_tcents];
        all_model_blockids = [all_model_blockids blockid*ones(1,length(cur_tcents))];
        
    end
end

%%
max_tsf = 0.6; nbins = 70;
used_tsf = time_since_fix;
used_tsf(used_tsf > max_tsf) = max_tsf;
tax = prctile(used_tsf,linspace(100/nbins,100-100/nbins,nbins));

sac_trg_rate = zeros(n_used_cells,length(tax));
taxb = [0 tax];
for i = 1:length(tax)
    curset = find(time_since_fix >= taxb(i) & time_since_fix < taxb(i+1));
    for c = 1:n_used_cells
        sac_trg_rate(c,i) = mean(spikes_binned(curset,c))/dt;
    end
    n_occ(i) = length(curset);
end
avg_rates = mean(spikes_binned)/dt;
all_sac_dir = atan2(all_sac_dy,all_sac_dx);
%%
post_sac_win = 0.1;
post_sac_win2 = 0.2;
post_sac_spkcnts = zeros(10,length(used_fixs));
post_sac_spkcnts2 = zeros(10,length(used_fixs));
for blockid = 1:4
    fprintf('Block %d of %d\n',blockid,4);
    cur_set = find(blockids(used_fixs) == blockid);
    for i = 1:length(cur_set)
        start_T = all_fix_start_times(used_fixs(cur_set(i)));
        end_T = start_T + post_sac_win;        
        end_T2 = start_T + post_sac_win2;        
        for c = 1:10
            post_sac_spkcnts(c,cur_set(i)) = ...
                sum(Blocks{blockid}.spktimes{c} >= start_T & Blocks{blockid}.spktimes{c} <= end_T);
            post_sac_spkcnts2(c,cur_set(i)) = ...
                sum(Blocks{blockid}.spktimes{c} >= end_T & Blocks{blockid}.spktimes{c} <= end_T2);
        end
    end
end


%%
for c = 1:10
    [modsacdir_corr(c),modsacdir_p(c)] = circ_corrcl(mod(all_sac_dir(used_fixs),pi),post_sac_spkcnts(c,:));
    [sacdir_corr(c),sacdir_p(c)] = circ_corrcl(all_sac_dir(used_fixs),post_sac_spkcnts(c,:));   
    [sacamp_corr(c),sacdir_p(c)] = corr(all_sac_amps(used_fixs),post_sac_spkcnts(c,:)','type','spearman');   
    [sacdur_corr(c),sacdir_p(c)] = corr(all_sac_durs(used_fixs),post_sac_spkcnts(c,:)','type','spearman');   
 
    [modsacdir_corr2(c),modsacdir_p2(c)] = circ_corrcl(mod(all_sac_dir(used_fixs),pi),post_sac_spkcnts2(c,:));
    [sacdir_corr2(c),sacdir_p2(c)] = circ_corrcl(all_sac_dir(used_fixs),post_sac_spkcnts2(c,:));   
    [sacamp_corr2(c),sacdir_p2(c)] = corr(all_sac_amps(used_fixs),post_sac_spkcnts2(c,:)','type','spearman');   
    [sacdur_corr2(c),sacdir_p2(c)] = corr(all_sac_durs(used_fixs),post_sac_spkcnts2(c,:)','type','spearman');   
end

%%
clear sacdiravg
nbins = 20;
phase_bin_edges = prctile(all_sac_dir,linspace(0,100,nbins+1));
for i = 1:nbins
    curset = find(all_sac_dir(used_fixs) >= phase_bin_edges(i) & all_sac_dir(used_fixs) < phase_bin_edges(i+1));
    sacdiravg(:,i) = mean(post_sac_spkcnts(:,curset),2);
end

