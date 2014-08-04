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
nlags = 1;
fix_win = 0.2;

cd /Users/James/James_scripts/bruce
load modfits_dualeye_25_ds4
all_fix_start_times(length(all_fix_stop_times)+1:end) = [];
blockids(length(all_fix_stop_times)+1:end) = [];
%%
X = [];
all_spike_bin_vecs = [];
all_time_since_fix = [];
all_eyepos = [];
all_fix_range = [];
all_time_inds = [];
all_block_inds = [];
all_sac_inds = [];
for blockid = 1:4
    
    fprintf('Processing block %d of %d...\n',blockid,4);
    
    block_times = Blocks{blockid}.blocktimes;
    stim_times = Blocks{blockid}.stimtime;
    stimID = Blocks{blockid}.stimids;
    
    n_stims = length(stim_times);
    
    cd ~/Data/bruce/2_27_12/saccades/
    load(sprintf('lemM232.5%d.em.sac.mat',blockid))
    
    % identify saccade start and stop times
    EyeStartT = Expt.Trials.Start/10000; % time of first eye sample
    EyeEndT = Expt.Trials.End/10000; % time of last eye sample
    Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
    eyets = EyeStartT:Eyedt:EyeEndT; %eye tracking time axis (sec)
    sac_buffer_inds = round(sac_buffer/Eyedt);
    
    reye_pos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
    leye_pos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];
    
    % load patch video data
    cd ~/Data/bruce/2_27_12/stimrecon/
    sname = sprintf('image_patch_block%d_leye_25_dsf4',blockid);
    load(sname)
    
    % interpolate eye signal onto stimulus time axis
%     eye_interp = interp1(eyets(1:end-1),reye_pos,recon_t);
        eye_interp = interp1(eyets(1:end-1),leye_pos,recon_t);
    %     eye_interp = interp1(eyets(1:end-1),leye_pos+reye_pos,recon_t)/2;
    
    
    %compute binned spike count vectors at stimulus resolution
    time_bin_edges = [(recon_t(1)-stimres/2) (recon_t+stimres/2)];
    spike_bin_vecs = zeros(10,length(recon_t));
    for cellid = 1:10
        spike_times = Blocks{blockid}.spktimes{cellid};
        
        spikes_binned = histc(spike_times,time_bin_edges);
        spikes_binned(end) = [];
        spike_bin_vecs(cellid,:) = spikes_binned;
    end
    
    cur_fix_set = find(blockids==blockid);
    cur_fix_start_inds = round(interp1(recon_t,1:length(recon_t),all_fix_start_times(cur_fix_set)));
    cur_fix_stop_inds = round(interp1(recon_t,1:length(recon_t),all_fix_stop_times(cur_fix_set)));
    bad_fixs = find(isnan(cur_fix_start_inds) | isnan(cur_fix_stop_inds));
    cur_fix_start_inds(bad_fixs) = [];
    cur_fix_stop_inds(bad_fixs) = [];
    
    used_samps = [];
    for i = 1:length(cur_fix_start_inds)
        fprintf('Fixation %d of %d\n',i,length(cur_fix_start_inds));
        cur_set = cur_fix_start_inds(i):cur_fix_stop_inds(i);
        %         cur_set = cur_fix_start_inds(i):(cur_fix_start_inds(i)+round(fix_win/stimres));
        used_samps = [used_samps; cur_set(:)];
        cur_ims = ov_im_patch(cur_set,:,:);
        cur_ims = cur_ims - mean(cur_ims(:));
%         cur_ims = dist_shift3d(cur_ims,-round(right_dY_seqs(end,cur_fix_set(i))*Fsd),2);
%         cur_ims = dist_shift3d(cur_ims,-round(right_dX_seqs(end,cur_fix_set(i))*Fsd),3);
        cur_ims = dist_shift3d(cur_ims,-round(left_dY_seqs(end,cur_fix_set(i))*Fsd),2);
        cur_ims = dist_shift3d(cur_ims,-round(left_dX_seqs(end,cur_fix_set(i))*Fsd),3);
        
        %         S = makeStimRows(cur_ims,nlags,0);
        %         X = [X; S];
        X = [X; reshape(cur_ims,length(cur_set),size(cur_ims,2)*size(cur_ims,3))];
        all_time_since_fix = [all_time_since_fix; (1:length(cur_set))'];
        all_fix_range = [all_fix_range; range(eye_interp(cur_set,:))];
    end
    
    all_spike_bin_vecs = [all_spike_bin_vecs; spike_bin_vecs(:,used_samps)'];
    all_eyepos = [all_eyepos; eye_interp(used_samps,:)];
    all_time_inds = [all_time_inds; recon_t(used_samps)'];
    all_block_inds = [all_block_inds; ones(size(used_samps))*blockid];
    all_sac_inds = [all_sac_inds; cur_fix_set];
end

%% convert spike data to bin indices
spikebins = cell(10,1);
for cellid = 1:10
    un_spk_cnts = unique(all_spike_bin_vecs(:,cellid));
    cur_spikebins = [];
    for i = 1:length(un_spk_cnts)
        cur_set = find(all_spike_bin_vecs(:,cellid) == un_spk_cnts(i));
        cur_spikebins = [cur_spikebins; repmat(cur_set(:),un_spk_cnts(i),1)];
    end
    cur_spikebins = sort(cur_spikebins);
    spikebins{cellid} = [spikebins{cellid}; cur_spikebins ];
end

%% shift all spike bins forward to account for desired temporal lag in response
desired_lag = 2;
for cellid = 1:10
    spikebins{cellid} = spikebins{cellid} + desired_lag;
    spikebins{cellid}(spikebins{cellid} > length(all_time_since_fix)) = [];
end

%%
backwin = 0;
forwardwin = 0.75;
T = backwin+forwardwin+stimres;
tax = -backwin:stimres:forwardwin;
temp = all_time_since_fix;
temp(temp > length(tax)) = nan;
Tmat = tbrep(temp,1:length(tax));

%%
SDIM = size(ov_im_patch,2);
[TT,XX,YY] = meshgrid(1:nlags,1:SDIM,1:SDIM);
XX = permute(XX,[2 1 3]);
YY = permute(YY,[2 1 3]);
XX = XX(:);
YY = YY(:);
uset = find(XX >= 1 & XX <= 33 & YY >= 1 & YY <= 33);
usdim = sqrt(length(uset)/nlags);

%%
% load modfits_avgeye_25_ds4_v2

% cellid = 5;
for cellid = 1:5;
    
    %     SDIM = size(ov_im_patch,2); %number of pixels (if 2d image, must be a square array at this point)
    SDIM = usdim; %number of pixels (if 2d image, must be a square array at this point)
    flen = nlags; %number of time lags
    kern_l = flen*SDIM^2; %filter dimensionality.  Note, this is flen*SDIM for 1 spatial dim.  Also note this must be equal to size(X,2)
    nmods = 1; %number of filters
    
    %basic filter initialization parameters (much of which is vestigial and not
    %used in this code, such as the internal NL tent-basis rep)
    n_nlx_samps = 11; %number of samples of the internal NL x-axis
    defmod.nlx = norminv(linspace(.001,0.999,n_nlx_samps)); %internal NL x-axis
    defmod.nly = defmod.nlx; %internal NL y-axis
    defmod.w = 1; %model weight
    defmod.SDIM = SDIM; %spatial dimensions
    defmod.fsdim = SDIM^2; %spatial dimensions if 2d stim
    defmod.pids = 1:defmod.fsdim; %pixels that are used for 2d stim (can be used to define ROIs)
    defmod.h = 1; %setting PSC kernel to delta
    defmod.locLambda = 0; %localization penalty
    defmod.lambda_dX = 5e6; %kernel spatial smoothness penalty
    defmod.lambda_dT = 0; %kernel temporal smoothness penalty
    defmod.lambda_L1x = 2e3; %kernel sparseness penalty
    basis = 'pix'; %basis type, can be either 'white' or 'pix'.  Set to 'pix'.
    
    init_kerns = randn(kern_l,nmods); %initial estimate of filter coeffs, in this case initialized randomly
    init_signs = ones(nmods,1); %vector specifying the 'sign' of each filter. 1 for excitatory, -1 for inhibitory
    init_betas = 2*ones(nmods,1); %vector specifying the beta parameters for the logexp internal NLs
    
    glm0 = createGLM_tlin(init_kerns,init_signs,defmod,basis);
    glm0.image_type = '2d'; %MUST specify '1d' or '2d' here indicating spatial dimensionality of stimulus
    
    glm0 = normalizeRFs_full(glm0,X(:,uset)); %This normalizes filter coefs so that the filtered stimulus distribution is standard normal (not really necessary)
    %     glm0 = mod_fit(end,cellid);
    glm0.const = -5;
    glm0.sacmod = zeros(size(tax));
    
    lexp_mod_sachist(cellid) = fitGLM_lexp_tempmod(glm0,X(:,uset),spikebins{cellid},Tmat,'tots'); %fit the model
    lexp_mod(cellid) = fitGLM_lexp(lexp_mod_sachist(cellid),X(:,uset),spikebins{cellid},'tots'); %fit the model
    glm0 = lexp_mod_sachist(cellid);
    sachist_mod(cellid) = fitGLM_temponly(glm0,Tmat,spikebins{cellid},'tots');
    
    %     lexp_mod_space(cellid) = fitGLM_lexp(glm0,X(:,uset),spikebins{cellid},'tots',50); %fit the model
    %     lexp_mod_space(cellid) = fitGLM_lexp(lexp_mod_space(cellid),X(:,uset),spikebins{cellid},'tots',20); %fit the model
end


%%
cd /Users/James/James_scripts/bruce
save onefilt_models_25_50lag_lefteye *_mod all_* spikebins

%%
clear spatiotempmod_LL spacemod_LL tempmod_LL
n_range_bins = 10;
range_bin_edges = prctile(all_time_since_fix,linspace(0,100,n_range_bins+1));
range_bin_centers = (range_bin_edges(1:end-1)+range_bin_edges(2:end))/2;
for cellid = 6:10;
    for ii = 1:n_range_bins
        time_set = find(all_time_since_fix >= range_bin_edges(ii) & all_time_since_fix < range_bin_edges(ii+1));
        un_spk_cnts = unique(all_spike_bin_vecs(time_set,cellid));
        cur_spikebins = [];
        for i = 1:length(un_spk_cnts)
            cur_set = find(all_spike_bin_vecs(time_set,cellid) == un_spk_cnts(i));
            cur_spikebins = [cur_spikebins; repmat(cur_set(:),un_spk_cnts(i),1)];
        end
        cur_spikebins = sort(cur_spikebins);
        
        spatiotempmod_LL(cellid,ii) = getLLGLM_lexp_tempmod(lexp_mod_sachist(cellid),X(time_set,uset),cur_spikebins,Tmat(time_set,:),'none');
        spacemod_LL(cellid,ii) = getLLGLM_lexp(lexp_mod(cellid),X(time_set,uset),cur_spikebins,'none');
        tempmod_LL(cellid,ii) = getLLGLM_temponly(sachist_mod(cellid),Tmat(time_set,:),cur_spikebins,'none');
        
    end
end

for cellid = 6:10
    sac_mod(cellid,:) = sachist_mod(cellid).sacmod;
end
%%
for cellid = 1:10
    filtX(cellid,:) = X(:,uset)*get_k_mat(lexp_mod(cellid));
end

%%
cellid = 6;
k_mat = get_k_mat(lexp_mod_sachist(cellid));
out = X(:,uset)*k_mat;
out_filt = out;
out_filt(out_filt < 0) = 0;
out_temp = Tmat*lexp_mod_sachist(cellid).sacmod';

ov_g = out_temp + out_filt + lexp_mod_sachist(cellid).const;
pred_rate = log(1+exp(ov_g));

figure
plot(out_filt)
hold on
plot(out_temp,'k')
plot(pred_rate,'g')
plot(smooth(all_spike_bin_vecs(:,cellid),3),'r')
sac_inds = find(all_time_since_fix==1);
plot(sac_inds,ones(size(sac_inds)),'r*')

%xlim([1.051 1.066]*1e4)

%%
cellid =2;
plot(smooth(all_spike_bin_vecs(:,cellid),10),'r')
hold on
plot(all_block_inds,'k')

