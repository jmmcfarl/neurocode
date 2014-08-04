clear all
close all
addpath('~/James_scripts/autocluster/');

global data_dir base_save_dir init_save_dir Expt_name Vloaded n_probes loadedData
Expt_name = 'M297';

Expt_num = str2num(Expt_name(2:end));
if Expt_name(1) == 'M'
    rec_type = 'LP';
elseif Expt_name(1) == 'G'
    rec_type = 'UA';
end

if Expt_num >= 280
    data_loc = '/media/NTlab_data3/Data/bruce/';
else
    data_loc = '/home/james/Data/bruce/';
end

if Expt_num >= 281
    data_dir2 = ['/media/NTlab_data3/Data/bruce/' Expt_name];
else
    data_dir2 = ['/home/james/Data/bruce/' Expt_name];
end

base_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
init_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering/init'];
spkdata_dir = [data_loc Expt_name '/spikes/'];

if ~exist(base_save_dir,'dir');
    mkdir(base_save_dir);
end
if ~exist(init_save_dir,'dir');
    mkdir(init_save_dir);
end

%location of FullV files
data_dir = [data_loc Expt_name];


Vloaded = nan;
cd(data_dir2);
if Expt_name(1) == 'G';
    if strcmp(Expt_name,'G029')
        load('G029Expts.mat');
    else
        load(sprintf('jbe%sExpts.mat',Expt_name));
    end
    n_probes = 96;
elseif Expt_name(1) == 'M'
    load(sprintf('lem%sExpts.mat',Expt_name));
    n_probes = 24;
end

target_probes = 1:n_probes;
force_new_clusters = false; %if you want to

elen = cellfun(@(x) length(x),Expts);
target_blocks = find(elen > 0);

if strcmp(Expt_name,'M232')
    target_blocks = [37    38    39    43    46    47];
elseif strcmp(Expt_name,'M235')
    target_blocks = [52    53    54    55    56    57    58    59    60];
elseif strcmp(Expt_name,'M239')
    target_blocks = [37    38    39    43    44    49];
elseif strcmp(Expt_name,'M281')
    target_blocks(target_blocks == 13) = [];
elseif strcmp(Expt_name,'M289')
    target_blocks(target_blocks == 14) = [];
elseif strcmp(Expt_name,'G087')
    target_blocks(target_blocks == 15) = [];
elseif strcmp(Expt_name,'G093')
    target_blocks(target_blocks == 28 | target_blocks == 52) = [];
end


%don't apply to blocks where we dont have the FullV data
missing_Vdata = [];
for bb = target_blocks
    if Expt_name(1) == 'M'
        check_name = [data_dir sprintf('/Expt%dFullV.mat',bb)];
    else
        check_name = [data_dir sprintf('/Expt%d.p1FullV.mat',bb)];
    end
    if ~exist(check_name,'file')
        missing_Vdata = [missing_Vdata bb];
    end
end
target_blocks(ismember(target_blocks,missing_Vdata)) = [];

%% load in basic cluster stats for all blocks
rclust_data = [base_save_dir '/Ref_Clusters.mat'];
load(rclust_data,'RefClusters');

for bb = target_blocks
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',bb)];
    load(cur_dat_name,'Clusters');
    all_block_dprimes(bb,:) = cellfun(@(x) x.dprime,Clusters);
    all_block_LL(bb,:) = cellfun(@(x) x.LL,Clusters);
    for pp = 1:n_probes
        temp_Lratios = Clusters{pp}.Lratios;
        temp_isodists = Clusters{pp}.iso_dists;
        temp_nspks = Clusters{pp}.n_spks;
        all_block_best_Lratio(bb,pp) = nanmin(temp_Lratios);
        all_block_best_isodist(bb,pp) = nanmax(temp_isodists);
    end
    clear Clusters
end
all_block_LL(all_block_LL == 0) = nan;

global full_save_dir

full_save_dir = [base_save_dir '/full'];
if ~exist(full_save_dir,'dir');
    mkdir(full_save_dir);
end

% figure
% imagesc(all_block_dprimes);
% figure
% imagesc(abs(all_block_best_isodist));
%% load stats for each unique cluster

clear SU_clust_data
su_cnt = 1;
for pp = 1:n_probes
    temp_Lratios = RefClusters{pp}.Lratios;
    temp_isodists = RefClusters{pp}.iso_dists;
    temp_nspks = RefClusters{pp}.n_spks;
    cur_n_sus = length(temp_Lratios);
    for jj = 1:cur_n_sus
        SU_clust_data(su_cnt).Lratio = temp_Lratios(jj);
        SU_clust_data(su_cnt).iso_dist = temp_isodists(jj);
        SU_clust_data(su_cnt).spk_rate = temp_nspks(jj+1)/RefClusters{pp}.recDur;
        SU_clust_data(su_cnt).n_spks = temp_nspks(jj+1);
        SU_clust_data(su_cnt).base_block_num = RefClusters{pp}.base_block;
        SU_clust_data(su_cnt).probe_num = pp;
        SU_clust_data(su_cnt).cluster_label = jj+1;
        su_cnt = su_cnt + 1;
    end
end

lrat_thresh = 1e3;
iso_thresh = 2;
su_pnums = [SU_clust_data(:).probe_num];
su_cnums = [SU_clust_data(:).cluster_label];
% good_SUs = find([SU_clust_data(:).Lratio] <= lrat_thresh | [SU_clust_data(:).iso_dist] >= iso_thresh);
good_SUs = find([SU_clust_data(:).iso_dist] >= iso_thresh);
good_SU_pnums = su_pnums(good_SUs);
good_SU_cnums = su_cnums(good_SUs);

%% LOOK at best clustering block projections for all potential sus
n_good_units = length(good_SUs);
n_cols = ceil(sqrt(n_good_units)); n_rows = ceil(n_good_units/n_cols);
f1 = figure();
f2 = figure();
for ii = 1:n_good_units
    su_inds = RefClusters{good_SU_pnums(ii)}.spike_clusts == good_SU_cnums(ii);
    all_inds = RefClusters{good_SU_pnums(ii)}.spike_clusts > 0;
    cur_block = RefClusters{good_SU_pnums(ii)}.base_block;
    spike_xy = RefClusters{good_SU_pnums(ii)}.spike_xy;
    
    set(0,'CurrentFigure',f1);
    subplot(n_cols,n_rows,ii);
    hold on
    plot(spike_xy(all_inds,1),spike_xy(all_inds,2),'k.');
    plot(spike_xy(su_inds,1),spike_xy(su_inds,2),'r.');
    title(sprintf('U%d P%d B%d\n',good_SUs(ii),good_SU_pnums(ii),cur_block));
    axis tight
    
    set(0,'CurrentFigure',f2);
    subplot(n_cols,n_rows,ii);hold on
    [handles, details] = DensityPlot_jmm(spike_xy(:,1),spike_xy(:,2),'sqrtsc','ynormal','sd',[1 1]);
    title(sprintf('Unit %d\n',good_SUs(ii)));
    axis tight
end

%% compute iso dist and Lratio for each cluster in each block
% target_blocks = [3 5 11];
pSU_Lratio_mat = nan(max(target_blocks),length(SU_clust_data));
pSU_isodist_mat = nan(max(target_blocks),length(SU_clust_data));
pSU_refract_mat = nan(max(target_blocks),length(SU_clust_data));
for bb = target_blocks
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',bb)];
    load(cur_dat_name);
    for ii = 1:length(SU_clust_data)
        if length(Clusters{su_pnums(ii)}.Lratios) >= su_cnums(ii)-1
            if ~isnan(Clusters{su_pnums(ii)}.Lratios)
                pSU_Lratio_mat(bb,ii) = Clusters{su_pnums(ii)}.Lratios(su_cnums(ii)-1)/Clusters{su_pnums(ii)}.n_spks(su_cnums(ii));
                pSU_isodist_mat(bb,ii) = abs(Clusters{su_pnums(ii)}.iso_dists(su_cnums(ii)-1));
                otherInds = setdiff(1:length(Clusters{su_pnums(ii)}.n_spks),su_cnums(ii));
                N_otherSpikes = sum(Clusters{su_pnums(ii)}.n_spks(otherInds));
                N_clustSpikes = Clusters{su_pnums(ii)}.n_spks(su_cnums(ii));
                
                %Isodist not defined when there are more spikes in the
                %cluster than outside it
                if N_otherSpikes<N_clustSpikes
                    pSU_isodist_mat(bb,ii) = nan;
                end
            end
        end
    end
end

% figure;
% imagescnan(pSU_Lratio_mat)
% set(gca,'xtick',1:length(SU_clust_data),'xticklabel',1:length(SU_clust_data));
% xlabel('SU number','fontsize',12);
% caxis([0 2]);
% 
% figure;
% imagescnan(pSU_isodist_mat)
% set(gca,'xtick',1:length(SU_clust_data),'xticklabel',1:length(SU_clust_data));
% xlabel('SU number','fontsize',12);
% ca = caxis();
% caxis([2 ca(2)]);

%% CHECK SPIKE CORRELATIONS
block_num = 15;
cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',block_num)];
load(cur_dat_name,'Clusters');
if Expt_name(1) == 'G'
    sfile = [data_dir sprintf('/Expt%d.p%dFullV.mat',block_num,1)];
    [loadedData.V,loadedData.Vtime,loadedData.Fs] = Load_FullV(sfile, false, [100 nan],1);
else
    sfile_name = [data_dir sprintf('/Expt%dFullV.mat',block_num)];
    if Vloaded ~= block_num
        fprintf('Loading data file %s\n',sfile_name);
        [loadedData.V,loadedData.Vtime,loadedData.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
        Vloaded = block_num;
    end
end

bin_width = 0.0015;
t_axis = loadedData.Vtime(1):bin_width:loadedData.Vtime(end);
binned_spikes = [];
for ii = 1:length(good_SUs)
    cur_su_ind = good_SUs(ii);
    cur_probe = SU_clust_data(cur_su_ind).probe_num;
    cur_spk_set = Clusters{cur_probe}.spike_clusts == SU_clust_data(cur_su_ind).cluster_label;
    cur_spk_times = Clusters{cur_probe}.times(cur_spk_set);
    cur_binned_spikes = hist(cur_spk_times,t_axis);
    binned_spikes = [binned_spikes cur_binned_spikes'];
end
[II,JJ] = meshgrid(1:length(good_SUs));
full_bs_corrmat = corr(binned_spikes);
bs_corrmat = full_bs_corrmat;
bs_corrmat(II <= JJ) = nan;

cond_p1 = nan(length(good_SUs),length(good_SUs));
cond_p2 = cond_p1;
for ii = 1:length(good_SUs)-1
    for jj = (ii+1):length(good_SUs)
        f_is_spike = binned_spikes(:,ii) > 0;
        s_is_spike = binned_spikes(:,jj) > 0;
        cond_p1(ii,jj) = sum(binned_spikes(f_is_spike,jj) > 0)/sum(f_is_spike);
        cond_p2(ii,jj) = sum(binned_spikes(s_is_spike,ii) > 0)/sum(s_is_spike);
    end
end
max_cond_prob = max(cond_p1,cond_p2);

figure
imagesc(max_cond_prob); set(gca,'ydir','normal');
set(gca,'xtick',1:length(good_SUs),'xticklabel',good_SUs);
set(gca,'ytick',1:length(good_SUs),'yticklabel',good_SUs);
caxis([0 1]);
colorbar;

clear binned_spikes
%% COMPARE spike waveforms for pair of clusters on a given pair of adjacent probes
pair = [29 31];
spk_pts = [-12:27];

probe1 = SU_clust_data(pair(1)).probe_num;
probe2 = SU_clust_data(pair(2)).probe_num;

spk_inds1 = Clusters{probe1}.spk_inds;
spk_inds1(spk_inds1 < -spk_pts(1)) = [];
cur_spk_set = Clusters{probe1}.spike_clusts == SU_clust_data(pair(1)).cluster_label;
Spikes = getSpikeSnippets(loadedData.V,loadedData.Vtime,spk_inds1(cur_spk_set)',spk_pts,probe1);
trig_avg_1 = squeeze(mean(Spikes.V));
trig_avg_1 = trig_avg_1(:,[probe1 probe2]);

spk_inds2 = Clusters{probe2}.spk_inds;
spk_inds2(spk_inds2 < -spk_pts(1)) = [];
cur_spk_set = Clusters{probe2}.spike_clusts == SU_clust_data(pair(2)).cluster_label;
Spikes = getSpikeSnippets(loadedData.V,loadedData.Vtime,spk_inds2(cur_spk_set)',spk_pts,probe2);
trig_avg_2 = squeeze(mean(Spikes.V));
trig_avg_2 = trig_avg_2(:,[probe1 probe2]);

pair_wvfrm_corr = corr(trig_avg_1(:),trig_avg_2(:));
figure;
subplot(2,1,1);hold on
plot(trig_avg_1(:,1),'k');
plot(trig_avg_1(:,2),'r');
title(sprintf('SU %d',pair(1)));
subplot(2,1,2); hold on
plot(trig_avg_2(:,1),'k');
plot(trig_avg_2(:,2),'r');
title(sprintf('SU %d',pair(2)));
fprintf('Waveform correlation %.3f\n',pair_wvfrm_corr);
ginds = find(ismember(good_SUs,pair));
fprintf('Spike time correlation %.3f\n',max_cond_prob(ginds(1),ginds(2)));
fprintf('SU %d, N%d, detected on Ch %d\n',pair(1),SU_clust_data(pair(1)).cluster_label,probe1);
fprintf('SU %d, N%d detected on Ch %d\n',pair(2),SU_clust_data(pair(2)).cluster_label,probe2);

%%
init_use_SUs = [];
switch Expt_name
    case 'M232'
        init_use_SUs = [2 4 5 8 14 17 19 21 28]; %CHECKED
    case 'M235'
        init_use_SUs = [2 3 7 14 18 22]; %CHECKED
    case 'M239'
        init_use_SUs = [3 4 6 9 10 13 16 17 20 26 28 30]; %CHECKED
    case 'M266'
        init_use_SUs = [13 14 24]; %M266 %CHECKED
    case 'M270'
        init_use_SUs = [10 12 17 18 24 27]; %for M270 CHECKED
    case 'M275'
        init_use_SUs = [3 6 17 19]; %M275 CHECKED
    case 'M277'
            init_use_SUs = [4 5 7 12 16 18]; %M277 CHECKED
    case 'M281'
        init_use_SUs = [4 5 6 7 8 9 10 13 17 18 20 22 24 26 27 28]; %M281
    case 'M287'
        init_use_SUs = [2 6 7 8 9 11 12 13 16 18 19 21 24]; %M287
    case 'M289'
        init_use_SUs = [3 8 9 17 19 20]; %M289
    case 'M294'
        init_use_SUs = [3 6 8 9 10 11 13 15 18 22]; %M294
    case 'M296'
        init_use_SUs = [11 13 15 16 17 18 21 23 26];
    case 'M297'
        init_use_SUs = [1 2 3 4 5 6 10 13 14 18 23 24 25 31 32];

    case 'G029'
        init_use_SUs = [2 4 5 9 14 23 24 31 39 47 49 55 63 66 70 71 80 81]; %G029 %CHECKED
    case 'G081'
        init_use_SUs = [9 18 37 42 50 56 72 80 86 94]; %for G081 %CHECKED
    case 'G085'
        init_use_SUs = [18 19 22 29 48 55 89]; %for G085 %CHECKED
    case 'G086'
         init_use_SUs = [12 42 43 48 49 91 94]; %for G086 CHECKED
   case 'G087'
        init_use_SUs = [19 36 48 64 94]; %for G087 CHECKED
    case 'G088'
        init_use_SUs = [17 44 55 64 85]; %for G088 CHECKED
    case 'G089'
        init_use_SUs = [25 48]; %for G089 CHECKED
    case 'G091'
        init_use_SUs = [18 42 43]; %for G091 CHECKED
    case 'G093'
        init_use_SUs = [22 34 48 74 76 86 94]; %for G093 CHECKED
    case 'G095'
        init_use_SUs = [7 9 12 18 44 45 56 60 75 77 90]; %for G095 CHECKED
    case 'G103'
        init_use_SUs = [38 94]; %G103 CHECKED
end

%% M232
% init_use_SUs = [2 4 5 8 14 17 19 21 28]; %CHECKED

%% M235
% init_use_SUs = [2 3 7 14 18 22]; %CHECKED

%% M239
% init_use_SUs = [3 4 6 9 10 13 16 17 20 26 28 30]; %CHECKED

%% newer LEM data
% init_use_SUs = [13 14 24]; %M266 %CHECKED
% init_use_SUs = [10 12 17 18 24 27]; %for M270 CHECKED
% init_use_SUs = [3 6 17 19]; %M275 CHECKED
% init_use_SUs = [4 5 7 12 16 18]; %M277 CHECKED
% init_use_SUs = [4 5 6 7 8 9 10 13 17 18 20 22 24 26 27 28]; %M281
% init_use_SUs = [2 6 7 8 9 11 12 13 16 18 19 21 24]; %M287
% init_use_SUs = [3 8 9 17 19 20]; %M289
% init_use_SUs = [3 6 8 9 10 11 13 15 18 22]; %M294
%% JBE data
% init_use_SUs = [2 4 5 9 14 23 24 31 39 47 49 55 63 66 70 71 80 81]; %G029 %CHECKED
% init_use_SUs = [9 18 37 42 50 56 72 80 86 94]; %for G081 %CHECKED
% init_use_SUs = [18 19 22 29 48 55 89]; %for G085 %CHECKED
% init_use_SUs = [12 42 43 48 49 91 94]; %for G086 CHECKED
% init_use_SUs = [19 36 48 64 94]; %for G087 CHECKED
% init_use_SUs = [17 44 55 64 85]; %for G088 CHECKED
% init_use_SUs = [25 48]; %for G089 CHECKED
% init_use_SUs = [18 42 43]; %for G091 CHECKED
% init_use_SUs = [22 34 48 74 76 86 94]; %for G093 CHECKED
% init_use_SUs = [7 9 12 18 44 45 56 60 75 77 90]; %for G095 CHECKED
% init_use_SUs = [38 94]; %G103 CHECKED

%% ASSIGN SU NUMBERS AND MANUALLY ELIMINATE BLOCKS WHERE CLUSTERING ISN"T CLEAR
iso_thresh = 2.25;
Lratio_thresh = 1e3;

SU_numbers = 1:length(init_use_SUs);
SU_ID_mat = nan(max(target_blocks),length(SU_clust_data));
SU_ID_mat(:,init_use_SUs) = repmat(SU_numbers,max(target_blocks),1);
SU_ID_mat(setdiff(1:max(target_blocks),target_blocks),:) = nan;

%iso dist must be bigger than thresh in each block
SU_ID_mat(pSU_isodist_mat < iso_thresh) = nan;
no_iso_dist = isnan(pSU_isodist_mat); %instances where iso distwas not defined, use Lratio
SU_ID_mat(no_iso_dist(pSU_Lratio_mat(no_iso_dist) > Lratio_thresh)) = nan;
bad = find(isnan(pSU_Lratio_mat) & isnan(pSU_isodist_mat)); %don't use cluster when neither was defined
SU_ID_mat(bad) = nan;

%secondary check based on visual inspection, exclude cluster/blocks that don't look
%good enough
switch Expt_name
    
    case 'M266'
        SU_ID_mat(1:4,14) = 1;
        SU_ID_mat(5:end,14) = nan;
        
    case 'M270'
        SU_ID_mat(21:25,12) = nan;
        SU_ID_mat([9 14 21 22 24],17) = nan;
        SU_ID_mat([1 6:24],18) = nan;
        SU_ID_mat([1:22],27) = nan;
        
    case 'M275'
        SU_ID_mat(1,6)  = nan;
        
    case 'M277'
        SU_ID_mat(27:end,4) = nan;
        SU_ID_mat(1:26,5) = nan;
        SU_ID_mat(27:end,5) = 1;
        SU_ID_mat(1:27,7) = nan;
        SU_ID_mat(27:end,12) = nan;
        SU_ID_mat(27:end,16) = nan;
        
    case 'M281'
        SU_ID_mat([11],5) = nan;
        SU_ID_mat([16 17 27 28],6) = nan;
        SU_ID_mat([5 7 8 10 25 27 28],7) = nan;
        SU_ID_mat([15 16 17 18],9) = nan;
        SU_ID_mat([15],10) = nan;
        SU_ID_mat([6 7 14],17) = nan;
        SU_ID_mat([19:end],18) = nan;
        SU_ID_mat([7 12],19) = nan;
        SU_ID_mat([16:end],22) = nan;
        SU_ID_mat([1:12],24) = nan;
        SU_ID_mat([1:18],26) = nan;
        SU_ID_mat([1:12 19:end],27) = nan;
        SU_ID_mat([1:17],28) = nan;
        
    case 'M287'
        SU_ID_mat([1:10],2) = nan;
        SU_ID_mat([12:end],6) = nan;
        SU_ID_mat([1:10],7) = nan;
        SU_ID_mat([1:10],8) = nan;
        SU_ID_mat([1:31],9) = nan;
        SU_ID_mat([10:end],11) = nan;
        SU_ID_mat([9:end],13) = nan;
        SU_ID_mat([1:23 28],18) = nan;
        SU_ID_mat([11:end],19) = nan;
        SU_ID_mat([1:10],21) = nan;
        SU_ID_mat([1:26],24) = nan;
        
    case 'M289'
        SU_ID_mat([1:18],3) = nan;
        SU_ID_mat([1:5],8) = nan;
        SU_ID_mat([1:5],20) = nan;
        
    case 'M294'
        SU_ID_mat([24],3) = nan;
        SU_ID_mat([24 30:end],6) = nan;
        SU_ID_mat([30:end],8) = nan;
        SU_ID_mat([30:end],9) = nan;
        SU_ID_mat([1:22],10) = nan;
        SU_ID_mat([30:end],13) = nan;
        SU_ID_mat([1:29],17) = nan;
        SU_ID_mat([1:12],22) = nan;
        
    case 'M296'
        SU_ID_mat([1],13) = nan;
        SU_ID_mat([12 13],15) = nan;
        SU_ID_mat([25 30 31],17) = nan;
        SU_ID_mat([4 5 6],21) = nan;
        SU_ID_mat([1:4 12 13],23) = nan;
        SU_ID_mat([14:end],26) = nan;
    case 'M297'
        SU_ID_mat([1:2 26:29],1) = nan;
        SU_ID_mat([1:2 6 30],2) = nan;
        SU_ID_mat(~isnan(SU_ID_mat(:,2)),2) = 1; %probe 2 and 1 are picking up the same unit
        SU_ID_mat([2],3) = nan;
        SU_ID_mat([1:16],4) = nan;
        SU_ID_mat([2],5) = nan;
        SU_ID_mat([1 2],6) = nan;
        SU_ID_mat([2 7],10) = nan;
        SU_ID_mat([1 2],13) = nan;
        SU_ID_mat([2],14) = nan;
        SU_ID_mat([2 5 6 7 16:29],18) = nan;
        SU_ID_mat([1:9],23) = nan;
        SU_ID_mat([7],24) = nan;
        SU_ID_mat([1 3 18:30],25) = nan;
        SU_ID_mat([1:5],31) = nan;
        SU_ID_mat([1:11],32) = nan;
end

% figure;
% imagescnan(SU_ID_mat);

final_SU_set = unique(SU_ID_mat(~isnan(SU_ID_mat)));
fprintf('Finalized %d unique SUs\n',length(final_SU_set));

%%
[grid_CLUSTnums,grid_BLOCKS] = meshgrid(1:length(SU_clust_data),1:max(target_blocks));

clust_avg_rates = nan(length(final_SU_set),max(target_blocks),1);
clust_X_prcs = nan(length(final_SU_set),max(target_blocks),3);
back_X_prcs = nan(length(final_SU_set),max(target_blocks),3);
clust_iso_dists = nan(length(final_SU_set),max(target_blocks));
clust_Lratios = nan(length(final_SU_set),max(target_blocks));
clust_dprimes = nan(length(final_SU_set),max(target_blocks));
clust_refracts = nan(length(final_SU_set),max(target_blocks));
clust_iso_reliable = nan(length(final_SU_set),max(target_blocks));
for bb = 1:length(target_blocks)
    fprintf('Processing block %d of %d\n',bb,length(target_blocks));
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',target_blocks(bb))];
    load(cur_dat_name,'Clusters');
    
    for cc = 1:length(final_SU_set);
        cur_SU_num = final_SU_set(cc);
        uu = find(SU_ID_mat == cur_SU_num);
        used_clusts = grid_CLUSTnums(uu);
        cur_SU_probes = su_pnums(used_clusts);
        cur_SU_clusts = su_cnums(used_clusts);
        cur_SU_blocks = grid_BLOCKS(uu);
        best_probe = mode(cur_SU_probes);
        
        temp = find(cur_SU_blocks == target_blocks(bb));
        if isempty(temp)
            spike_xy = Clusters{best_probe}.spike_xy;
            clust_spikes = [];
            back_X_prcs(cc,target_blocks(bb),:) = prctile(spike_xy(:,1),[10 50 90]);
        else
            spike_xy = Clusters{cur_SU_probes(temp)}.spike_xy;
            clust_spikes = find(Clusters{cur_SU_probes(temp)}.spike_clusts == cur_SU_clusts(temp));
            non_clust_spikes = find(Clusters{cur_SU_probes(temp)}.spike_clusts > 0);
            non_clust_spikes(ismember(non_clust_spikes,clust_spikes)) = [];
            back_spikes = find(Clusters{cur_SU_probes(temp)}.spike_clusts == 1);
            
            clust_avg_rates(cc,target_blocks(bb)) = length(clust_spikes)/Clusters{cur_SU_probes(temp)}.recDur;
            clust_X_prcs(cc,target_blocks(bb),:) = prctile(spike_xy(clust_spikes,1),[10 50 90]);
            back_X_prcs(cc,target_blocks(bb),:) = prctile(spike_xy(:,1),[10 50 90]);
        
            clust_spike_times = Clusters{cur_SU_probes(temp)}.times(clust_spikes);
            isis = diff(clust_spike_times);
            refract_ratio = sum(isis < 1e-3)/length(isis);
            clust_refracts(cc,target_blocks(bb)) = refract_ratio;
            
            clust_mean = mean(spike_xy(clust_spikes,1));
            back_mean = mean(spike_xy(back_spikes,1));
            clust_var = var(spike_xy(clust_spikes,1));
            back_var = var(spike_xy(back_spikes,1));
            dprime = abs((clust_mean - back_mean))./sqrt(0.5*(clust_var + back_var));
            
            clust_iso_dists(cc,target_blocks(bb)) = Clusters{cur_SU_probes(temp)}.iso_dists(cur_SU_clusts(temp)-1);
            clust_Lratios(cc,target_blocks(bb)) = Clusters{cur_SU_probes(temp)}.Lratios(cur_SU_clusts(temp)-1);
            clust_dprimes(cc,target_blocks(bb)) = dprime;
            clust_iso_reliable(cc,target_blocks(bb)) = length(clust_spikes) < length(non_clust_spikes);
        end
    end
end
tempr = ones(length(final_SU_set),1);
tempc = ones(max(target_blocks),1);
SU_allBlock_Data = struct('Lratios',mat2cell(clust_Lratios,tempr,tempc),'isoDists',mat2cell(clust_iso_dists,tempr,tempc),...
    'isoReliable',mat2cell(clust_iso_reliable,tempr,tempc),'dprimes',mat2cell(clust_dprimes,tempr,tempc),...
    'refract',mat2cell(clust_refracts,tempr,tempc));
%%
fname = [base_save_dir '/final_cluster.mat'];
fprintf('Saving final clustering to %s\n',fname);
SU_target_blocks = target_blocks;
save(fname,'SU_clust_data','SU_ID_mat','SU_target_blocks','SU_allBlock_Data');


%% CREATE SCATTERPLOTS FOR EACH BLOCK [STILL DOESN"T WORK FOR SUS THAT ARE CLUSTERED ON MULTIPLE PROBES ACROSS BLOCKS]
fin_save_dir = [base_save_dir '/final'];
if ~exist(fin_save_dir,'dir')
    mkdir(fin_save_dir);
end
n_blocks = length(target_blocks);
n_cols = ceil(sqrt(n_blocks)); n_rows = ceil(n_blocks/n_cols);
close all

for cc = 1:length(final_SU_set)
    
    cur_SU_num = final_SU_set(cc);
    uu = find(SU_ID_mat == cur_SU_num);
    used_clusts = grid_CLUSTnums(uu);
    if ~isempty(used_clusts)
        cur_SU_probes = su_pnums(used_clusts);
        cur_SU_clusts = su_cnums(used_clusts);
        cur_SU_blocks = grid_BLOCKS(uu);
        best_probe = mode(cur_SU_probes);
        best_clust = mode(cur_SU_clusts);
        %     probe_num = su_pnums(ii);
        %     clust_num = su_cnums(ii);
        fprintf('Printing final cluster figures for Probe %d Unit %d\n',best_probe,cur_SU_num);
        
        ref_block = RefClusters{best_probe}.base_block;
        probe_nums = ones(max(target_blocks),1)*best_probe;
        probe_nums(cur_SU_blocks) = cur_SU_probes;
        
        clust_nums = ones(max(target_blocks),1)*best_clust;
        clust_nums(cur_SU_blocks) = cur_SU_clusts;
        
        clust_ids = ones(max(target_blocks),1)*mode(used_clusts);
        clust_ids(cur_SU_blocks) = used_clusts;
        %     [full_scatter_fig] = regenerate_allblock_xyscatters_v2(best_probe,target_blocks,best_clust,0);
        [full_scatter_fig] = regenerate_allblock_xyscatters_v2(ref_block,probe_nums,target_blocks,clust_nums);
        set(0,'CurrentFigure',full_scatter_fig);
        for bb = 1:length(target_blocks)
            subplot(n_cols,n_rows,bb)
            if ~isnan(SU_ID_mat(target_blocks(bb),clust_ids(bb)))
                title(sprintf('Block %d',target_blocks(bb)),'Color','r');
            else
                title(sprintf('Block %d',target_blocks(bb)),'Color','k');
            end
        end
        
        fillPage(full_scatter_fig,'papersize',[14 14]);
        if length(unique(probe_nums)) > 1
            pfname_sc = [fin_save_dir sprintf('/Unit%d_Probe',cur_SU_num)];
            un_pnums = unique(probe_nums);
            for jj = 1:length(un_pnums)
               pfname_sc = strcat(pfname_sc,sprintf('%d_',un_pnums(jj))); 
            end
            pfname_sc = strcat(pfname_sc,'scatter');
        else
            pfname_sc = [fin_save_dir sprintf('/Unit%d_Probe%d_scatter',cur_SU_num,best_probe)];
        end
        print(full_scatter_fig,pfname_sc,'-dpng');
        close(full_scatter_fig);
    end
end


%%
fin_save_dir = [base_save_dir '/final'];
if ~exist(fin_save_dir,'dir')
    mkdir(fin_save_dir);
end
xl = [0 max(target_blocks)];

for cc = 1:length(final_SU_set)
    cur_SU_num = final_SU_set(cc);
    uu = find(SU_ID_mat == cur_SU_num);
    used_clusts = grid_CLUSTnums(uu);
    cur_SU_probes = su_pnums(used_clusts);
    cur_SU_clusts = su_cnums(used_clusts);
    cur_SU_blocks = grid_BLOCKS(uu);
    best_probe = mode(cur_SU_probes);
    bax = 1:max(target_blocks);

    cur_fig = figure('visible','off');
    subplot(2,2,1);hold on
    errorbar(bax,back_X_prcs(cc,:,2),back_X_prcs(cc,:,2)-back_X_prcs(cc,:,1),...
        back_X_prcs(cc,:,3)-back_X_prcs(cc,:,2),'k','linewidth',0.5)
    errorbar(bax+0.25,clust_X_prcs(cc,:,2),clust_X_prcs(cc,:,2)-clust_X_prcs(cc,:,1),...
        clust_X_prcs(cc,:,3)-clust_X_prcs(cc,:,2),'r','linewidth',2)
    ylabel('Spike Xdists');
    xlabel('Block number');
    xlim(xl);

    subplot(2,2,3);hold on
    plot(bax,clust_dprimes(cc,:),'ko-');
    xlim(xl);
    ylabel('DPrime');
    xlabel('Block number');
    
    subplot(2,2,2);hold on
    [AX,H1,H2] = plotyy(bax,clust_iso_dists(cc,:),1:max(target_blocks),-log10(clust_Lratios(cc,:)));
    set(get(AX(1),'Ylabel'),'String','Iso distance');
    set(get(AX(2),'Ylabel'),'String','Negative Log Lratio');
    set([H1 H2],'LineWidth',2)
    xlabel('Block number');
    hold(AX(1),'on'); hold(AX(2),'on');
    plot(AX(1),xl,[iso_thresh iso_thresh],'b--');
    xlim(xl);
    plot(AX(2),xl,-log10([Lratio_thresh Lratio_thresh]),'g--');
    uset = find(clust_iso_reliable(cc,:));
    plot(AX(1),bax(uset),clust_iso_dists(cc,uset),'o');
    plot(AX(2),bax,-log10(clust_Lratios(cc,:)),'o');
    xlim(AX(2),xl);
     
    subplot(2,2,4);hold on
    plot(bax,clust_refracts(cc,:)*100,'ro-','linewidth',2);
    ylabel('Refractory percentage');
    xlabel('Block number');
    xlim(xl);
    
    fillPage(cur_fig,'papersize',[11 7]);
    pfname = [fin_save_dir sprintf('/Unit%d_Probe%d_sep',cur_SU_num,best_probe)];
    print(cur_fig,pfname,'-dpng');
    close(cur_fig);    

end
%%
spkdata_dir = [data_loc Expt_name '/spikes/'];
blocks_with_clusters = find(any(~isnan(SU_ID_mat),2));
%     spk_data_name = [spkdata_dir Expt_name sprintf('_p%d_blk%d.mat',cur_probe_num,bb)];

clear all_avg_waveforms
for bb = blocks_with_clusters'
    cur_clust_set = find(~isnan(SU_ID_mat(bb,:)));
    cur_probe_set = su_pnums(cur_clust_set);
    
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',bb)];
    load(cur_dat_name,'Clusters');
    
    for cc = 1:length(cur_probe_set)
        cur_probe_num = cur_probe_set(cc);
        cur_su_num = SU_ID_mat(bb,cur_clust_set(cc));
        cur_su_ind = find(final_SU_set == cur_su_num);
        cur_Cluster = Clusters{cur_probe_num};
        
        
        if Expt_name(1) == 'G'
            loadedData = [data_dir sprintf('/Expt%d.p%dFullV.mat',bb,probe_num)];
        else
            sfile_name = [data_dir sprintf('/Expt%dFullV.mat',bb)];
            if Vloaded ~= bb
                fprintf('Loading data file %s\n',sfile_name);
                [loadedData.V,loadedData.Vtime,loadedData.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
                Vloaded = bb;
            end
        end
    
        V = loadedData.V(:,cur_Cluster.use_chs);
        Vtime = loadedData.Vtime;
        [spk_id, trig_thresh] = triggerSpikes(V(:,cur_Cluster.trig_ch),cur_Cluster.params.thresh_sign,[],cur_Cluster.trig_thresh);
        spk_id(spk_id <= abs(cur_Cluster.params.spk_pts(1)) | spk_id >= length(V)-cur_Cluster.params.spk_pts(end)) = []; %get rid of spikes at the edges
       
        %extract spike snippets
        Spikes = getSpikeSnippets(V,Vtime,spk_id,cur_Cluster.params.spk_pts,cur_Cluster.trig_ch);
        clear V Vtime
        
        % artifact detection
        artifact_ids = find_spike_artifacts(Spikes,cur_Cluster.params);
        Spikes.V(artifact_ids,:,:) = [];
        Spikes.times(artifact_ids) = [];
        Spikes.trig_vals(artifact_ids) = [];
        spk_id(artifact_ids) = [];
        Spikes.spk_inds = spk_id(:);
        
        if length(cur_Cluster.times) ~= length(Spikes.times)
            error('Spike alignment problem!');
        end
        
        cluster_stats = get_cluster_stats(Spikes.V(:,:,cur_Cluster.trig_ch),cur_Cluster.spike_clusts);
        
        all_avg_waveforms(bb,cur_su_ind).mean_spike = cluster_stats.mean_spike;
        all_avg_waveforms(bb,cur_su_ind).std_spike = cluster_stats.std_spike;
    end
end

%%
fin_save_dir = [base_save_dir '/final'];
if ~exist(fin_save_dir,'dir')
    mkdir(fin_save_dir);
end
n_blocks = length(target_blocks);
n_cols = ceil(sqrt(n_blocks)); n_rows = ceil(n_blocks/n_cols);
close all

for cc = 1:length(final_SU_set);
    
    cur_SU_num = final_SU_set(cc);
    uu = find(SU_ID_mat == cur_SU_num);
    used_clusts = grid_CLUSTnums(uu);
    cur_SU_probes = su_pnums(used_clusts);
    cur_SU_clusts = su_cnums(used_clusts);
    cur_SU_blocks = grid_BLOCKS(uu);
    best_probe = mode(cur_SU_probes);
    
    full_wvfrm_fig = figure('visible','off');
    for bb = 1:length(target_blocks)
        cur_probe = find(SU_ID_mat(target_blocks(bb),:) == final_SU_set(cc));
        if ~isempty(cur_probe)
            cur_cnum = su_cnums(cur_probe);
            n_pts = size(all_avg_waveforms(target_blocks(bb),cc).mean_spike,1);
            subplot(n_cols,n_rows,bb);hold on
            plot(1:n_pts,all_avg_waveforms(target_blocks(bb),cc).mean_spike(:,1),'k','linewidth',2);
            if size(all_avg_waveforms(target_blocks(bb),cc).mean_spike,2) > 1
                plot(1:n_pts,all_avg_waveforms(target_blocks(bb),cc).mean_spike(:,2:end),'b','linewidth',1);
                errorbar(1:n_pts,all_avg_waveforms(target_blocks(bb),cc).mean_spike(:,cur_cnum),...
                    all_avg_waveforms(target_blocks(bb),cc).std_spike(:,cur_cnum),'r','linewidth',1);
            else
                fprintf('Warning: no waveform available but marked as clustered!\n');
            end
            title(sprintf('Block %d',target_blocks(bb)));
            axis off
        end
    end
    
    pfname = [fin_save_dir sprintf('/Unit%d_Probe%d_wvfrm',cur_SU_num,best_probe)];
    print(full_wvfrm_fig,pfname,'-dpng');
    close(full_wvfrm_fig);
end
