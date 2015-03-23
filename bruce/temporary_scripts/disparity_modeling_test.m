clear all
close all

monName = 'jbe';
Expt_name = 'M008';
rec_type = 'LP';
data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];


%%

min_trial_dur = 2;
trial_dur = 4;
use_nPix = 20;
dt = 0.01;
flen = 12;
Fr = 1;
n_probes = 24;
spatial_usfac = 1;
%%
cd(data_dir);
ExptFileName = strcat(monName,Expt_name,'Expts.mat');
fprintf('Loading existing file %s\n',ExptFileName);
load(ExptFileName);
cd stims/
load('./stim_data.mat');
load('./expt_data.mat');

%%
include_expts = {'rls.dpXceXFrRC'};
usable_expts = find(cellfun(@(x) ~isempty(x),Expts));
expt_names = cellfun(@(x) x.Header.expname,Expts(usable_expts),'uniformoutput',0);
        
cur_block_set = usable_expts(cellfun(@(x) any(strcmp(x,include_expts)),expt_names));
n_blocks = length(cur_block_set);

full_nPix = unique(expt_npix(cur_block_set));
if length(full_nPix) > 1
    warning('multiple npix detected');
end

%%
load([data_dir '/CellList']);
good_sus = find(all(CellList(cur_block_set,:,1) > 0));

%% Spatial resolution
all_dws = cellfun(@(x) x.Stimvals.dw,Expts(cur_block_set));
base_sp_dx = mode(all_dws);
if length(unique(all_dws)) > 1
    fprintf('Warning, multiple different dws detected, using %.3f\n',base_sp_dx);
end
sp_dx = base_sp_dx/spatial_usfac; %model dx in deg
ov_RF_pos = Expts{cur_block_set(1)}.Stimvals.rf(1:2);

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
all_Lstim_mat = [];
all_Rstim_mat = [];
all_frame_ce = [];
all_frame_dp = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_blk = [];
all_trial_Fr = [];
all_nframes = [];
all_nstims = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_trial_rptframes = [];

trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
cur_spkind_offset = 0;
for ee = 1:n_blocks;
    fprintf('Block %d of %d;  UNMATCHED EXPT TYPE\n',ee,n_blocks);
    cur_block = cur_block_set(ee);
    
    fname = sprintf('Expt%dClusterTimes.mat',cur_block_set(ee));
    load(fname);
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times' + cur_toffset);
%         all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},ones(size(Clusters{cc}.times')));
    end
    
    trial_start_times = [Expts{cur_block}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_block}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_block}.Trials(:).dur]/1e4;
    trial_ids = [Expts{cur_block}.Trials(:).id];
    [un_ids,id_inds] = unique(trial_ids);
    rpt_trials = false;
    if length(un_ids) < length(trial_ids)
        rpt_trials = true;
        fprintf('Warning, repeat trial inds detected!\n');
        use_trials = [];
    else
        use_trials = find(trial_durs >= min_trial_dur);
    end
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    trial_Se = [Expts{cur_block}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    trial_Fr = [Expts{cur_block}.Trials(:).Fr];
    trial_Fr = trial_Fr(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    all_trial_Fr = cat(1,all_trial_Fr,trial_Fr(use_trials)');
    all_trial_blk = cat(1,all_trial_blk,ones(length(use_trials),1)*ee);
        
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    buffer_pix = floor((expt_npix(cur_block) - full_nPix)/2);
%     if buffer_pix == -1
%         for ii = 1:length(left_stim_mats)
%             left_stim_mats{ii} = [zeros(size(left_stim_mats{ii},1),1) left_stim_mats{ii} zeros(size(left_stim_mats{ii},1),1)];
%         end
%         buffer_pix = 0;
%     end
    cur_use_pix = (1:full_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    cur_nrpt_frames = zeros(n_trials,1);
    [trial_nframes,trial_nstims] = deal(nan(n_trials,1));
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
        cur_stim_ce = Expts{cur_block}.Trials(use_trials(tt)).ce;
        cur_stim_dp = Expts{cur_block}.Trials(use_trials(tt)).dp;
        
        if trial_Fr(use_trials(tt)) == 3
           cur_stim_times = repmat(cur_stim_times,[3 1]);
           cur_stim_times = cur_stim_times(:); cur_stim_times = cur_stim_times(1:301)';
           cur_stim_ce = repmat(cur_stim_ce',[3 1]);
           cur_stim_ce = cur_stim_ce(:); cur_stim_ce = cur_stim_ce(1:301);
           cur_stim_dp = repmat(cur_stim_dp',[3 1]);
           cur_stim_dp = cur_stim_dp(:); cur_stim_dp = cur_stim_dp(1:301); 
        end
        
        
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        trial_nframes(tt) = n_frames;
        trial_nstims(tt) = length(cur_stim_times);
        if isfield(Expts{cur_block}.Trials(use_trials(tt)),'rptframes')
            cur_nrpt_frames(tt) = length(Expts{cur_block}.Trials(use_trials(tt)).rptframes);
        end
        if n_frames > 0
            if length(cur_stim_times) == 1
                cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
            end
            cur_t_edges = [cur_stim_times'; cur_stim_times(end) + dt*Fr];
        end
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt & trial_nframes(tt) == trial_nstims(tt)
            use_frames = min(length(cur_stim_times),n_frames);
            cur_Lstim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            cur_Rstim_mat = double(right_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times+cur_toffset < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_nframes = [all_nframes; trial_nframes(:)];
            all_nstims = [all_nstims; trial_nstims(:)];
            all_stim_times = [all_stim_times; cur_stim_times' + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_Lstim_mat = [all_Lstim_mat; cur_Lstim_mat];
            all_Rstim_mat = [all_Rstim_mat; cur_Rstim_mat];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_frame_ce = [all_frame_ce; cur_stim_ce];
            all_frame_dp = [all_frame_dp; cur_stim_dp];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
    all_trial_rptframes = [all_trial_rptframes; cur_nrpt_frames];
    
    %need to keep track of block time offsets for LP recordings
    if strcmp(rec_type,'LP')
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
        cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
    end
end

%% BIN SPIKES FOR MU AND SU
clust_params.n_probes = n_probes;
all_binned_mua = get_quickbinned_mu(all_spk_times,all_clust_ids,all_t_axis,all_t_bin_edges,all_bin_edge_pts,clust_params);

%%
full_nPix_us = spatial_usfac*full_nPix;
% if spatial_usfac > 1
%     all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
%     for ii = 1:size(all_stim_mat,2)
%         for jj = 1:spatial_usfac
%             all_stimmat_up(:,spatial_usfac*(ii-1)+jj) = all_stim_mat(:,ii);
%         end
%     end
% elseif spatial_usfac == 1
%     all_stimmat_up = all_stim_mat;
% end


stim_params_us = NMMcreate_stim_params([flen full_nPix_us*2],dt);

%%
flip_frames = all_frame_ce == -1;
all_Rstim_mat(flip_frames,:) = -all_Rstim_mat(flip_frames,:);
% all_Lstim_mat(flip_frames,:) = -all_Lstim_mat(flip_frames,:);

%%
blank_frames = all_frame_ce == -1009;
all_Rstim_mat(blank_frames,:) = 0;
all_Lstim_mat(blank_frames,:) = 0;
%% CREATE STIMULUS MATRIX
stim_params = NIMcreate_stim_params([flen full_nPix_us*2],dt);
all_Xmat = create_time_embedding([all_Lstim_mat all_Rstim_mat],stim_params);

%%
all_Xmat(all_Xmat == -128) = nan;
bad_trials = [];
good_inds = [];
un_trials = unique(all_trialvec);
for tt = 1:length(un_trials)
    cur_set = find(all_trialvec == un_trials(tt));
    if any(isnan(reshape(all_Xmat(cur_set,:),[],1)))
        bad_trials = [bad_trials un_trials(tt)];
    else
        good_inds = [good_inds; cur_set(:)];
    end
end
%%


% DEFINE DATA USED FOR ANALYSIS
beg_buffer = 0.15;
end_buffer = 0.05;
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
NT = length(used_inds);

used_inds(ismember(all_trialvec(used_inds),bad_trials)) = [];

%%
silent = 0;
base_lambda_d2XT = 100;
base_lambda_L1 = 1;

n_stim_filts = 3;
mod_signs = ones(1,n_stim_filts);
NL_types = [{'lin' 'quad','quad'}];
init_d2XT = [ones(n_stim_filts,1)];
init_L2 = [zeros(n_stim_filts,1);];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',base_lambda_d2XT,'lambda_L1',base_lambda_L1,'boundary_conds',[0 0 0]);

for ss = 22
    fprintf('Fitting model for MU %d of %d\n',ss,n_probes);
    Robs = all_binned_mua(used_inds,ss);
    
    if sum(Robs) > 0
        gqm1 = NMMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params);
%         for ii = 1:length(gqm1.mods)
%             gqm1.mods(ii).filtK(:) = 0;
%         end
        gqm1 = NMMfit_filters(gqm1,Robs,all_Xmat(used_inds,:),[],[],silent);
        
        all_mod_fits(ss) = gqm1;
        
        [LL, penLL, pred_rate, G, gint, fgint, nullLL] = NMMmodel_eval(all_mod_fits(ss),Robs,all_Xmat(used_inds,:));
        all_mod_LLimp(ss) = (LL-nullLL)/log(2);
    end
end

%%
for ss = 1:n_probes
    all_filts(ss,:,:) = reshape(all_mod_fits(ss).mods(1).filtK,flen,use_nPix);
end
temp_profiles = squeeze(std(all_filts,[],3));
space_profiles = squeeze(std(all_filts,[],2));
[~,best_lags] = max(temp_profiles,[],2);
[~,best_pos] = max(space_profiles,[],2);
for ss = 1:n_probes
    best_temp_profiles(ss,:) = squeeze(all_filts(ss,:,best_pos(ss)));
    best_space_profiles(ss,:) = squeeze(all_filts(ss,best_lags(ss),:));
end
figure
imagesc(0:dt*1e3:(flen-1)*dt*1e3,1:n_probes,best_temp_profiles);
xlabel('Time (ms)','fontsize',12);
ylabel('Probe','fontsize',12);

figure
imagesc(1:use_nPix,1:n_probes,best_space_profiles);
xlabel('Bar position','fontsize',12);
ylabel('Probe','fontsize',12);