close all

elen = cellfun(@(x) length(x),Expts);
target_blocks = find(elen > 0);

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
figure
imagesc(abs(all_block_best_isodist));
%%

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

%% CHECK SPIKE CORRELATIONS
block_num = 17;
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

bin_width = 0.0005;
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
full_bs_corrmat2 = corr(binned_spikes);
bs_corrmat2 = full_bs_corrmat2;
bs_corrmat2(II <= JJ) = nan;

figure
subplot(2,1,1);
imagesc(bs_corrmat); set(gca,'ydir','normal');
set(gca,'xtick',1:length(good_SUs),'xticklabel',good_SUs);
set(gca,'ytick',1:length(good_SUs),'yticklabel',good_SUs);
caxis([0 0.1]);
colorbar;
subplot(2,1,2);
imagesc(bs_corrmat2); set(gca,'ydir','normal');
set(gca,'xtick',1:length(good_SUs),'xticklabel',good_SUs);
set(gca,'ytick',1:length(good_SUs),'yticklabel',good_SUs);
caxis([0 0.1]);
colorbar;

clear binned_spikes
%% COMPARE 
pair = [21 24];
spk_pts = [-12:27];

probe1 = SU_clust_data(pair(1)).probe_num;
probe2 = SU_clust_data(pair(2)).probe_num;

spk_inds1 = Clusters{probe1}.spk_inds;
cur_spk_set = Clusters{probe1}.spike_clusts == SU_clust_data(pair(1)).cluster_label;
Spikes = getSpikeSnippets(loadedData.V,loadedData.Vtime,spk_inds1(cur_spk_set)',spk_pts,probe1);
trig_avg_1 = squeeze(mean(Spikes.V));
trig_avg_1 = trig_avg_1(:,[probe1 probe2]);

spk_inds2 = Clusters{probe2}.spk_inds;
cur_spk_set = Clusters{probe2}.spike_clusts == SU_clust_data(pair(2)).cluster_label;
Spikes = getSpikeSnippets(loadedData.V,loadedData.Vtime,spk_inds2(cur_spk_set)',spk_pts,probe2);
trig_avg_2 = squeeze(mean(Spikes.V));
trig_avg_2 = trig_avg_2(:,[probe1 probe2]);

pair_wvfrm_corr = corr(trig_avg_1(:),trig_avg_2(:));
figure;
subplot(2,1,1)
plot(trig_avg_1);
title(sprintf('SU %d',pair(1)));
subplot(2,1,2)
plot(trig_avg_2);
title(sprintf('SU %d',pair(2)));
fprintf('Waveform correlation %.3f\n',pair_wvfrm_corr);

%%
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
    subplot(n_cols,n_rows,ii);hold on
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


%%
% target_blocks = [3 5 11];
pSU_Lratio_mat = nan(max(target_blocks),length(SU_clust_data));
pSU_isodist_mat = nan(max(target_blocks),length(SU_clust_data));
for bb = target_blocks
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',bb)];
    load(cur_dat_name);
    for ii = 1:length(SU_clust_data)
        if length(Clusters{su_pnums(ii)}.Lratios) >= su_cnums(ii)-1
        pSU_Lratio_mat(bb,ii) = Clusters{su_pnums(ii)}.Lratios(su_cnums(ii)-1)/Clusters{su_pnums(ii)}.n_spks(su_cnums(ii));
        pSU_isodist_mat(bb,ii) = abs(Clusters{su_pnums(ii)}.iso_dists(su_cnums(ii)-1));
        end
    end
end

figure;
imagescnan(pSU_Lratio_mat)
set(gca,'xtick',1:length(SU_clust_data),'xticklabel',1:length(SU_clust_data));
xlabel('SU number','fontsize',12);
caxis([0 2]);

figure;
imagescnan(pSU_isodist_mat)
set(gca,'xtick',1:length(SU_clust_data),'xticklabel',1:length(SU_clust_data));
xlabel('SU number','fontsize',12);
ca = caxis();
caxis([2 ca(2)]);

%%
% init_use_SUs = [10 12 17 18 21 24 27]; %for M270
% init_use_SUs = [9 18 37 42 50 56 72 80 86 94]; %for G081
% init_use_SUs = [18 19 22 29 48 55 89]; %for G085
% init_use_SUs = [7 9 12 18 44 45 56 60 75 77 90]; %for G095
init_use_SUs = [12 42 43 48 49 91 94]; %for G086
SU_numbers = 1:length(init_use_SUs);
SU_ID_mat = nan(max(target_blocks),length(SU_clust_data));
SU_ID_mat(:,init_use_SUs) = repmat(SU_numbers,max(target_blocks),1);

% SU_ID_mat(pSU_Lratio_mat > 1e3) = nan;
SU_ID_mat(pSU_isodist_mat < 2.25) = nan;
no_iso_dist = isnan(pSU_isodist_mat);
SU_ID_mat(no_iso_dist(pSU_Lratio_mat(no_iso_dist) > 1e3)) = nan;
% SU_ID_mat(~(pSU_Lratio_mat < 1e3) & ~(pSU_isodist_mat > 2)) = nan;

figure;
imagescnan(SU_ID_mat);

%%
% NOTE THIS ONLY WORKS WHEN SUS ARE ON THE SAME PROBE ACROSS ALL BLOCKS
n_cols = ceil(sqrt(n_blocks)); n_rows = ceil(n_blocks/n_cols);
close all
for ii = init_use_SUs(1:end)
    probe_num = su_pnums(ii);
    fprintf('Probe %d\n',probe_num);
    %look for existing full scatter figure and open if it exists
    pfname_sc = [full_save_dir sprintf('/Probe%d_fullclust_scatter.fig',probe_num)];
    if exist(pfname_sc,'file') %&& ~ishandle(full_scatter_fig)
        full_scatter_fig = open(pfname_sc); set(full_scatter_fig,'visible','on');
    end
    set(0,'CurrentFigure',full_scatter_fig);
    for bb = target_blocks
    subplot(n_cols,n_rows,bb)
        if ~isnan(SU_ID_mat(bb,ii))
           title(sprintf('Block %d',bb),'Color','r'); 
        else
           title(sprintf('Block %d',bb),'Color','k');             
        end
    end

    resp = input('Enter vector of block numbers to flip\n');
    if ~isempty(resp)
        to_good = resp(isnan(SU_ID_mat(resp,ii)));
        to_bad = resp(~isnan(SU_ID_mat(resp,ii)));
        SU_ID_mat(to_good,ii) = find(init_use_SUs==ii);
        SU_ID_mat(to_bad,ii) = nan;
        
        set(0,'CurrentFigure',full_scatter_fig);
        for bb = target_blocks
            subplot(n_cols,n_rows,bb)
            if ~isnan(SU_ID_mat(bb,ii))
                title(sprintf('Block %d',bb),'Color','r');
            else
                title(sprintf('Block %d',bb),'Color','k');
            end
        end
       
    end
    
    pause
    close all
end

%%
fname = [base_save_dir '/final_cluster.mat'];
fprintf('Saving final clustering to %s\n',fname);
SU_target_blocks = target_blocks;
save(fname,'SU_clust_data','SU_ID_mat','SU_target_blocks');

%%
fin_save_dir = [base_save_dir '/final'];
if ~exist(fin_save_dir,'dir')
    mkdir(fin_save_dir);
end

n_cols = ceil(sqrt(n_blocks)); n_rows = ceil(n_blocks/n_cols);
close all
for ii = init_use_SUs
    
    probe_num = su_pnums(ii);
    clust_num = su_cnums(ii);
    fprintf('Printing final cluster figures for Probe %d Clust %d\n',probe_num,clust_num);
    
    [full_scatter_fig] = regenerate_allblock_xyscatters(probe_num,target_blocks,clust_num,0);    
    set(0,'CurrentFigure',full_scatter_fig);
    for bb = target_blocks
        subplot(n_cols,n_rows,bb)
        if ~isnan(SU_ID_mat(bb,ii))
            title(sprintf('Block %d',bb),'Color','r');
        else
            title(sprintf('Block %d',bb),'Color','k');
        end
    end

    pfname_sc = [fin_save_dir sprintf('/Unit%d_Probe%d_scatter',ii,probe_num)];
    print(full_scatter_fig,pfname_sc,'-dpng');
    close(full_scatter_fig);
end

