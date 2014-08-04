cd ~/Data/bruce/M232/
load CellList
n_probes = 24;
bar_expts = [37 38 39 43 46 47];
all_sus = squeeze(CellList(bar_expts,:,:));
all_sus = all_sus(:);
unique_su_nums = unique(all_sus(all_sus > 0));

clear su_*
probe_has_su = zeros(length(bar_expts),n_probes);
for ii = 1:length(unique_su_nums)
    curset = find(CellList(bar_expts,:,:) == unique_su_nums(ii));
    [su_block{ii},su_probe{ii},su_extra{ii}] = ind2sub(size(CellList(bar_expts,:,:)),curset);
    probe_has_su(su_block{ii},su_probe{ii}) = 1;
end

%%
global data_dir base_save_dir init_save_dir Expt_name Vloaded n_probes loadedData
Expt_name = 'M232';

% data_loc = '/media/NTlab_data1/Data/bruce/';
data_loc = '/home/james/Data/bruce/';

base_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
init_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering/init'];

if ~exist(base_save_dir,'dir');
    mkdir(base_save_dir);
end
if ~exist(init_save_dir,'dir');
    mkdir(init_save_dir);
end

%location of FullV files
data_dir = [data_loc Expt_name];
Vloaded = nan;
%%
block_num = 47;
probe_num = 8;

fprintf('Loading block %d Clusters\n',block_num);
cur_clust_data = [base_save_dir sprintf('/Block%d_Clusters.mat',block_num)];
load(cur_clust_data,'Clusters');

me_clust = Clusters{probe_num};

bclust_data = sprintf('Expt%dClusterTimes.mat',block_num);
bdat = load(bclust_data);
br_clust = bdat.Clusters{probe_num};

sfile_name = [data_dir sprintf('/Expt%dFullV.mat',block_num)];
if Vloaded ~= block_num
    fprintf('Loading data file %s\n',sfile_name);
    [loadedData.V,loadedData.Vtime,loadedData.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
    Vloaded = block_num;
end

% use_chs = me_clust.use_chs;
use_chs = probe_num;
% trig_ch = me_clust.trig_ch;
trig_ch = 1;
me_Nspks = length(me_clust.spike_clusts);
target_Nspks = 2*me_Nspks;
V = loadedData.V(:,use_chs);
Vtime = loadedData.Vtime;
Fs = loadedData.Fs;
[spk_id, trig_thresh,noise_sigma] = triggerSpikes(V(:,trig_ch),me_clust.params.thresh_sign,target_Nspks);
spk_id(spk_id <= abs(me_clust.params.spk_pts(1)) | spk_id >= length(V)-me_clust.params.spk_pts(end)) = []; %get rid of spikes at the edges
Spikes = getSpikeSnippets(V,Vtime,spk_id,me_clust.params.spk_pts,trig_ch);

br_spike_inds = round(interp1(loadedData.Vtime,1:length(loadedData.Vtime),br_clust.times));
br_set = find(ismember(spk_id,br_spike_inds));

% Spikes2 = getSpikeSnippets(V,Vtime,br_spike_inds,me_clust.params.spk_pts,trig_ch);

[N_spks,D,N_chs] = size(Spikes.V);
if N_chs > 1
    AllV = reshape(Spikes.V,N_spks,D*N_chs);
else
    AllV = Spikes.V;
end
C = cov(AllV);
[pc_coeffs, E] = eig(C);
[a,b] = max(diag(E));
Eval = diag(E);
if b == 1
    fprintf('Max Eigernvector First\n');
else
    pc_coeffs = fliplr(pc_coeffs); %put largest first;
    Eval = Eval(end:-1:1);
end

%arbitrary sign convention just so that its consistent when re-applying
for j = 1:size(pc_coeffs,1)
    if sum(pc_coeffs(:,j)) < 0
        pc_coeffs(:,j) = -pc_coeffs(:,j);
    end
end
pc_scores = AllV*pc_coeffs;

n_comps = 2;
pc_ks = nan(D,1);
for ii = 1:D
    pc_ks(ii) = lillie_KSstat(pc_scores(:,ii));
end
[~,use_pcs] = sort(pc_ks,'descend');
use_pcs = use_pcs(1:4);


tdim_ks = nan(D*N_chs,1);
for ii = 1:D*N_chs
    tdim_ks(ii) = lillie_KSstat(AllV(:,ii));
end
peak_locs = find(diff(sign(diff(tdim_ks))) < 0);
if ~isempty(peak_locs)
    peak_locs = peak_locs + 1;
    peak_amps = tdim_ks(peak_locs);
    [~,use_peaks] = sort(peak_amps,'descend');
    peak_locs = peak_locs(use_peaks);
    if length(peak_locs) > 4
        peak_locs = peak_locs(1:4);
    end
else
    peak_locs = [first_peak second_peak round((second_peak+first_peak)/2) second_peak + round((second_peak-first_peak)/2)];
end
use_tdims = peak_locs;
use_tdims(use_tdims > D*N_chs) = [];







