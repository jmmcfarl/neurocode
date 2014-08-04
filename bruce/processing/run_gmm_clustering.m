clear all
close all

% Expt_name = 'G081';
% target_probes = [9 12 13 18 37 42 50 72 80 86 94];
% Expt_name = 'G085';
% target_probes = [9 18 19 22 29 33 42 48 55 58 89];
Expt_name = 'G086';
% target_probes = [12 36 42 43 48 49 51 91 94];
target_probes = [36];
% Expt_name = 'G087';
% target_probes = [9 19 36 48 64 94];
% Expt_name = 'G088';
% target_probes = [9 19 44 48 55 60 64 65 70 85];
% Expt_name = 'G089';
% target_probes = [9 17 25 48 60 70 85 94 96];
% Expt_name = 'G091';
% target_probes = [9 18 42 43 47 64 65 94];
% Expt_name = 'G092';
% target_probes = [22 43 49 74 85 94];
% Expt_name = 'G093';
% target_probes = [9 22 33 44 48 74 76 86 94];
% Expt_name = 'G095';
% target_probes = [7 9 18 33 44 49 55 58 59 74 76 89];
% Expt_name = 'M266';
% target_probes = [13 14 24];
% Expt_name = 'M270';
% target_probes = [10 11 12 15 16 17 18 21:24];

dir_prefix = '~';
% dir_prefix = '/Volumes/james';
% dir_prefix = '/media/NTlab_data1';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
cd(data_dir);

if Expt_name(1) == 'G'
    load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
n_probes = 96;
elseif Expt_name(1) == 'M'
    load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
n_probes = 24;
end
n_expts = length(Expts);
bad_expt = [];

%%
init_mahal_mat = nan(n_expts,n_probes,4);
init_marked = nan(n_expts,n_probes);
for ee = 1:n_expts
    fprintf('Loading Expt%dClusterTimes\n',ee);
    fname = sprintf('Expt%dClusterTimes.mat',ee);
    if ~exist(fname,'file')
        bad_expt = [bad_expt ee];
    else
        load(fname);
        for cc = 1:n_probes
            init_mahal_mat(ee,cc,:) = Clusters{cc}.mahal;
            if isfield(Clusters{cc},'marked')
            init_marked(ee,cc) = Clusters{cc}.marked;
            end
        end
    end
end
imagesc(squeeze(init_mahal_mat(:,:,1)));
caxis([0 5]);
%% PERFORM INITIAL 2-COMPONENT GMM FITS

n_comps = 2;
n_reps = 1;
avg_wvfrm = zeros(40,n_comps);
std_wvfrm = zeros(40,n_comps);
sigmas = nan(1,n_comps);
if ~exist('stored_spkxy','var')
    stored_spkxy = cell(n_expts,length(target_probes));
    stored_spkt = cell(n_expts,length(target_probes));
end
target_expts = 1:n_expts;
target_expts(ismember(target_expts,bad_expt)) = [];
for ee = target_expts
    
    if isempty(stored_spkxy{ee,1})
        fprintf('Loading Expt%dClusterTimesDetails\n',ee);
        fname = sprintf('Expt%dClusterTimesDetails.mat',ee);
        load(fname);
    end
    for pp = 1:length(target_probes)
        
        fprintf('Loading Expt%d P%d Spikes\n',ee,target_probes(pp));
        fname = sprintf('Spikes/nby%s.p%dt%d.mat',Expt_name,target_probes(pp),ee);
        load(fname);
        
        if isempty(stored_spkxy{ee,pp})
            n_clst_spks = length(ClusterDetails{target_probes(pp)}.t);
            if n_clst_spks ~= length(Spikes.times)
                use_set = find(ismember(Spikes.times,ClusterDetails{target_probes(pp)}.t));
                use_set2 = find(ismember(ClusterDetails{target_probes(pp)}.t,Spikes.times));
                fprintf('Warning, no match for %d of %d spikes\n',length(setdiff(1:n_clst_spks,use_set)),n_clst_spks);
            else
                use_set = 1:n_clst_spks;
                use_set2 = 1:n_clst_spks;
            end
            stored_spkxy{ee,pp} = ClusterDetails{target_probes(pp)}.xy;
            stored_spkt{ee,pp} = ClusterDetails{target_probes(pp)}.t;
        end
        if size(Spikes.values,1) ~= length(Spikes.times)
            Spikes.values = Spikes.values';
        end
        
        spk_values = double(Spikes.values(use_set,:));
       
        clust_gmm = gmdistribution.fit(stored_spkxy{ee,pp},n_comps,'Replicates',n_reps,'Regularize',1e-15,'Options',statset('MaxIter',1000));
        [clust_idx,nlogl] = cluster(clust_gmm,stored_spkxy{ee,pp});
        D = mahal(clust_gmm,clust_gmm.mu);
        mahal_d = sqrt(2./((1./D(1,2))+(1./D(2,1))));
        for j = 1:n_comps
            sigmas(j) = sqrt(sum(diag(clust_gmm.Sigma(:,:,j)))); %var for this dimemsion
        end
        nsd = diff(clust_gmm.mu)./sigmas;
        dprime = sqrt(sum(nsd.^2));
        
        mv1 = mean(stored_spkxy{ee,pp}(clust_idx == 1,:));
        mv2 = mean(stored_spkxy{ee,pp}(clust_idx == 2,:));
        cov1 = cov(stored_spkxy{ee,pp}(clust_idx == 1,:));
        cov2 = cov(stored_spkxy{ee,pp}(clust_idx == 2,:));
        D1 = (mv1-mv2)*inv(cov2)*(mv1-mv2)';
        
        D2 = (mv1-mv2)*inv(cov1)*(mv1-mv2)';
        tmahal_d = sqrt(2./((1./D1)+(1./D2)));
        
        clear avg_wvfrm std_wvfrm
        
        avg_wvfrm(:,1) = mean(spk_values(clust_idx(use_set2)==1,:),1);
        avg_wvfrm(:,2) = mean(spk_values(clust_idx(use_set2)==2,:),1);
        
        smallest_vals = min(avg_wvfrm);
        [~,clust_ord] = sort(smallest_vals,'ascend');
        avg_wvfrm = avg_wvfrm(:,clust_ord);
        clust_idx = clust_ord(clust_idx);
        
        std_wvfrm(:,1) = std(spk_values(clust_idx(use_set2)==1,:),[],1);
        std_wvfrm(:,2) = std(spk_values(clust_idx(use_set2)==2,:),[],1);
        avg_wvfrm = avg_wvfrm*Spikes.maxv/Spikes.maxint;
        std_wvfrm = std_wvfrm*Spikes.maxv/Spikes.maxint;
        
        autoclust(ee,pp).probenum = target_probes(pp);
        autoclust(ee,pp).gmm = clust_gmm;
        autoclust(ee,pp).n_comps = n_comps;
        autoclust(ee,pp).idx = clust_idx;
        autoclust(ee,pp).times = ClusterDetails{target_probes(pp)}.t;
        autoclust(ee,pp).avg_wvfrm = avg_wvfrm;
        autoclust(ee,pp).std_wvfrm = std_wvfrm;
        autoclust(ee,pp).nlogl = nlogl;
        autoclust(ee,pp).mahal_d = mahal_d;
        autoclust(ee,pp).tmahal_d = tmahal_d;
        autoclust(ee,pp).dprime = dprime;
        autoclust(ee,pp).man_code = nan;
        autoclust(ee,pp).cluster_id = pp;
    end
end
for ee = 1:length(bad_expt)
    for pp = 1:length(target_probes)
    autoclust(bad_expt(ee),pp).probenum = nan;
    autoclust(bad_expt(ee),pp).gmm = nan;
    autoclust(bad_expt(ee),pp).n_comps = nan;
    autoclust(bad_expt(ee),pp).idx = nan;
    autoclust(bad_expt(ee),pp).times = nan;
    autoclust(bad_expt(ee),pp).avg_wvfrm = nan;
    autoclust(bad_expt(ee),pp).std_wvfrm = nan;
    autoclust(bad_expt(ee),pp).nlogl = nan;
    autoclust(bad_expt(ee),pp).mahal_d = nan;
    autoclust(bad_expt(ee),pp).tmahal_d = nan;
    autoclust(bad_expt(ee),pp).dprime = nan;
    autoclust(bad_expt(ee),pp).man_code = nan;
    autoclust(bad_expt(ee),pp).cluster_id = nan;
end
end
sname = 'cur_aclust_data';
save(sname,'autoclust','target_probes');

%% VISUALIZE SPECIFIED PROBE
close all
probe_num = 1;
expt_nums = 1:n_expts;
% expt_nums = 20;
use_plots = [1 0 1];
dsize = 5;
expt_nums(ismember(expt_nums,bad_expt)) = [];
% use_plots = [plot_dens plot_wvfrm plot_scat];
fprintf('Plotting probe %d of %d\n',probe_num,length(target_probes));
visualize_probe_clustering(probe_num,expt_nums,stored_spkxy,autoclust,use_plots,dsize);

%% REFINE TARGET GMMS
probe_num = 11;
poss_n_comps = 2:3;
n_reps = 1;
target_expts = 1:n_expts;
% target_expts = 29:n_expts;
target_expts(ismember(target_expts,bad_expt)) = [];
dsize = 5;
reload = 0;
new_autoclust = refine_gmms(probe_num,stored_spkxy,stored_spkt,autoclust,poss_n_comps,n_reps,Expt_name,target_probes,target_expts,dsize,reload);
have_char = 0;
while have_char == 0
    ui = input('Use new aclust (y/n)?\n','s');
    if strcmpi(ui,'y')
        autoclust = new_autoclust;
        clear new_autoclust
        have_char = 1;
    end
    if have_char == 0
        fprintf('Invalid input\n');
    end
end
close all

%% FLIP TARGET CIDX
probe_num = 9;
poss_expts = 1:n_expts;poss_expts(ismember(poss_expts,bad_expt)) = [];
target_expts = [60 68 71];
% target_inds = find(ismember(poss_expts,target_expts));
autoclust = flip_cluster_assignment(probe_num,autoclust,target_expts,stored_spkxy);

%% SUMMARY PLOTS
close all
mah_thresh = 2;
min_used_blocks = 5;

used_expts = 1:n_expts;
for ee = 1:n_expts
if isempty(autoclust(ee,1).probenum)
    used_expts(ee) = [];
end
end

aclust_mahals = nan(n_expts,length(target_probes));
aclust_dprimes = nan(n_expts,length(target_probes));
for pp = 1:length(target_probes)    
    aclust_dprimes(used_expts,pp) = [autoclust(used_expts,pp).dprime];
    aclust_mahals(used_expts,pp) = [autoclust(used_expts,pp).mahal_d]';
end

aclust_CellList = zeros(n_expts,n_probes);
for pp = 1:length(target_probes)
    for ee = 1:n_expts
        if ~isempty(autoclust(ee,pp).cluster_id)
       aclust_CellList(ee,target_probes(pp)) = autoclust(ee,pp).cluster_id;
        end
    end
end

cmap = jet(length(target_probes));
plot(aclust_mahals,'k');
hold on
for pp = 1:length(target_probes)
    good_set = find(aclust_mahals(:,pp) > mah_thresh);
    if length(good_set) < min_used_blocks
        good_set = [];
    end
%     aclust_CellList(good_set,target_probes(pp)) = pp;
bad_set = setdiff(1:n_expts,good_set);
    aclust_CellList(bad_set,target_probes(pp)) = 0;
    plot(good_set,aclust_mahals(good_set,pp),'.','color',cmap(pp,:),'markersize',16);
end

figure
imagesc(aclust_CellList);

%% MERGE CLUSTER IDS
merge_cluster_ids = [1 2 3];
for ee = 1:n_expts
    cur_probes = find(ismember(aclust_CellList(ee,:),merge_cluster_ids));
    aclust_CellList(ee,cur_probes) = 0;
    [~,best_probe] = max(aclust_mahals(ee,ismember(target_probes,cur_probes)),[],2);
    best_probe = cur_probes(best_probe);
    aclust_CellList(ee,best_probe) = merge_cluster_ids(1);
end
figure
imagesc(aclust_CellList);

%% VISUALIZE SPECIFIED CLUSTER ID
clust_id = 8;
expt_nums = 1:n_expts;
expt_nums = 20;
use_plots = [1 1 1];
dsize = 5;
expt_nums(ismember(expt_nums,bad_expt)) = [];

expt_nums = expt_nums(any(aclust_CellList(expt_nums,:) == clust_id,2));
probe_nums = nan(size(expt_nums));
for ee = 1:length(expt_nums)
    probe_nums(ee) = find(aclust_CellList(expt_nums(ee),target_probes) == clust_id);
end

% use_plots = [plot_dens plot_wvfrm plot_scat];
fprintf('Plotting probe %d of %d\n',probe_num,length(target_probes));
visualize_probe_clustering(probe_nums,expt_nums,stored_spkxy,autoclust,use_plots,dsize);

%% SAVE FINAL ACLUST
save_dir = ['~/Data/bruce/' Expt_name];
cd(save_dir);
sname = 'fin_aclust_data';
save(sname,'autoclust','aclust_*','target_probes');

