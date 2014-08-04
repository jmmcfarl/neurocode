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

%%
X = [];
spk_cnts = [];
all_sac_amps = [];
all_sac_vels = [];
all_sac_durs = [];
all_fix_start_times = [];
all_fix_stop_times = [];
all_stim_filtered = [];
all_left_locs = [];
all_right_locs = [];
all_left_stds = [];
all_right_stds = [];
all_is_blinks = [];
blockids = [];
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
    
    %if we haven't already computed saccade times do so now
    sname = sprintf('sac_times_corf_block%d',blockid);
    if ~exist([sname '.mat'],'file')
        fprintf('Computing saccade times\n');
        avg_eyepos = (reye_pos + leye_pos)/2;
        clear sm_avg_eyepos eye_vel
        sm_avg_eyepos(:,1) = smooth(avg_eyepos(:,1),3);
        sm_avg_eyepos(:,2) = smooth(avg_eyepos(:,2),3);
        eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/Eyedt;
        eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/Eyedt;
        eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
        
        
        vergence = reye_pos - leye_pos;
        vervel = [0 0;diff(vergence)]/Eyedt;
        %         verspeed = abs(vervel);
        verspeed = sqrt(vervel(:,1).^2+vervel(:,2).^2);
        %         verspeed = max(vervel,[],2); %take the larger of the hor and ver vergence speeds (Read and Cummin 2003).
        blink_on = find(verspeed(2:end) > blink_thresh & verspeed(1:end-1) <= blink_thresh);
        blink_off = find(verspeed(2:end) <= blink_thresh & verspeed(1:end-1) > blink_thresh);
        blink_dur = eyets(blink_off) - eyets(blink_on);
        is_blink = find(blink_dur > min_blink_dur);
        blink_times = eyets(round((blink_on(is_blink)+blink_off(is_blink))/2));
        
        %find saccade start and stop indices
        sac_inds = find(eye_speed(1:end-1) < sac_eyespeed & eye_speed(2:end) > sac_eyespeed);
        
        sac_start_inds = nan(size(sac_inds));
        sac_stop_inds = nan(size(sac_inds));
        for i = 1:length(sac_inds)
            temp = find(eye_speed(1:sac_inds(i)) < thresh_eyespeed,1,'last');
            if ~isempty(temp)
                sac_start_inds(i) = temp;
            end
            temp = find(eye_speed(sac_inds(i)+1:end) < thresh_eyespeed,1,'first');
            if ~isempty(temp)
                sac_stop_inds(i) = sac_inds(i)+temp-1;
            end
        end
        
        %identify start and stop times of unique saccades
        sac_vec = zeros(size(reye_pos,1),1);
        for i = 1:length(sac_start_inds)
            if ~isnan(sac_start_inds(i)) & ~isnan(sac_stop_inds(i))
                sac_vec(sac_start_inds(i):sac_stop_inds(i)+sac_buffer_inds) = 1;
            end
        end
        sac_vec(length(eyets)+1:end) = [];
        sac_vec([1 end]) = 0;
        sac_start_indsn = find(sac_vec(1:end-1) == 0 & sac_vec(2:end) == 1);
        sac_stop_indsn = find(sac_vec(1:end-1) == 1 & sac_vec(2:end) == 0);
        if length(sac_start_indsn) ~= length(sac_stop_indsn)
            error('saccade mis-alignment');
        end
        
        sac_start_times = eyets(sac_start_indsn);
        sac_stop_times = eyets(sac_stop_indsn);
        
        %compute saccade amplitudes, velocities, and durations
        sac_dx = avg_eyepos(sac_stop_indsn,1) - avg_eyepos(sac_start_indsn,1);
        sac_dy = avg_eyepos(sac_stop_indsn,2) - avg_eyepos(sac_start_indsn,2);
        sac_amps = sqrt(sac_dx.^2+sac_dy.^2);
        sac_peakvel = zeros(size(sac_start_indsn));
        for i = 1:length(sac_start_indsn)
            sac_peakvel(i) = max(eye_speed(sac_start_indsn(i):sac_stop_indsn(i)));
        end
        sac_durs = sac_stop_times-sac_start_times;
        
        cd ~/Data/bruce/2_27_12/saccades/
        save(sname,'sac_start_times','sac_stop_times','sac_start_inds','sac_stop_inds','sac_vec',...
            'sac_amps','sac_peakvel','sac_durs','blink_times')
    else
        load(sname)
    end
    
    % load patch video data
    cd ~/Data/bruce/2_27_12/stimrecon/
    sname = sprintf('image_patch_block%d_avgeye_25_dsf4',blockid);
    load(sname)
    
    % interpolate eye signal onto stimulus time axis
    reye_interp = interp1(eyets(1:end-1),reye_pos,recon_t);
    leye_interp = interp1(eyets(1:end-1),leye_pos,recon_t);
    
    avg_interp = 0.5*reye_interp + 0.5*leye_interp;
    
    %create an interpolated 0-1 vector of saccade times
    sac_vec_interp = zeros(size(recon_t));
    for i = 1:length(sac_start_times)
        cur_set = find(recon_t >= sac_start_times(i) & recon_t <= sac_stop_times(i));
        sac_vec_interp(cur_set) = 1;
    end
    sac_vec_interp([1 end]) = 0;
    
    %indices of the start and stops of fixations
    fix_start_inds = 1+[0 find(sac_vec_interp(1:end-1) == 1 & sac_vec_interp(2:end) == 0)];
    fix_stop_inds = [find(sac_vec_interp(1:end-1) == 0 & sac_vec_interp(2:end) == 1) length(recon_t)-1];
    if length(fix_start_inds) ~= length(fix_stop_inds)
        error('Fixation mis-alignment');
    end
    fix_durs = (fix_stop_inds - fix_start_inds)*stimres;
    too_short = find(fix_durs < min_fix_dur); %only keep fixations that are minimum duration
    
    fprintf('%d of %d fixations too short\n',length(too_short),length(fix_durs));
    fix_start_inds(too_short) = []; fix_stop_inds(too_short) = []; fix_durs(too_short) = [];
    fix_start_times = recon_t(fix_start_inds);
    
    fix_dx = reye_interp(fix_stop_inds,1) - reye_interp(fix_start_inds,1);
    fix_dy = reye_interp(fix_stop_inds,2) - reye_interp(fix_start_inds,2);
    right_drift_amps = sqrt(fix_dx.^2+fix_dy.^2);
    fix_dx = leye_interp(fix_stop_inds,1) - leye_interp(fix_start_inds,1);
    fix_dy = leye_interp(fix_stop_inds,2) - leye_interp(fix_start_inds,2);
    left_drift_amps = sqrt(fix_dx.^2+fix_dy.^2);
    
    %determine correspondence between fixation starts and previous saccades
    fix_prev_sac_inds = zeros(size(fix_start_inds));
    for i = 1:length(fix_start_inds)
        [~,fix_prev_sac_inds(i)] = min(abs(sac_stop_times-fix_start_times(i)));
    end
    
    % create vector determining whether current stimulus was high-pass filtered
    filt_stims = mod(0:500,4) + 1 >2;
    is_stim_filtered = filt_stims(stimID(stim_num));
    
    %find times where we skip stimulus presentation data
    %also find blinks
    bad_times = [0 diff(recon_t)] > 0.05;
    fix_skips = zeros(size(fix_start_inds));
    blinks = zeros(size(fix_start_inds));
    for i = 1:length(fix_start_inds)
        if any(bad_times(fix_start_inds(i):fix_stop_inds(i)))
            fix_skips(i) = 1;
        end
        if min(abs(blink_times - fix_start_times(i))) < 0.1
            blinks(i) = 1;
        end
    end
    
    % find times where eye-position is within central window
    in_window = (avg_interp(:,1) >= accept_window(1,1) & avg_interp(:,1) <= accept_window(1,2) & ...
        avg_interp(:,2) >= accept_window(2,1) & avg_interp(:,2) <= accept_window(2,2));
    
    off_screen = max(any(isnan(ov_im_patch),3),[],2);
    
    
    %     used_fix_start_inds = find(in_window(fix_start_inds)==1 & off_screen(fix_start_inds) == 0 & fix_skips'==0 & blinks'==0);
%     used_fix_start_inds = find(in_window(fix_start_inds)==1 & off_screen(fix_start_inds) == 0 & fix_skips'==0);
    used_fix_start_inds = find(in_window(fix_start_inds)==1 & off_screen(fix_start_inds) == 0 & fix_skips'==0 & is_stim_filtered(fix_start_inds)'==0);
    %     used_fix_start_inds = find(in_window(fix_start_inds)==1 & is_stim_filtered(fix_start_inds)' == 0);
    n_used_fixs = length(used_fix_start_inds);
    fprintf('Out of window: %d,  off-screen:  %d, gray-stim: %d,  blinks: %d\n',sum(in_window(fix_start_inds)==0),sum(off_screen(fix_start_inds)==1),sum(fix_skips==1),sum(blinks==1));
    fprintf('analyzing %d of %d fixations\n',n_used_fixs,length(fix_start_inds));
    
    %compute binned spike vectors
    spike_bin_vecs = [];
    Nimage = length(stim_times);
    for ii = 1:Nimage
        onsetT = stim_times(ii); %stimulus onset time (s)
        if ii < Nimage %for the last image, take the end time as the block end time
            endT = stim_times(ii+1)-stimres; %set end time as one time bin before the following stimulus presentation
            endT = endT - 0.2; %last 200ms doesn't have any image presentation!
        else
            endT = Blocks{blockid}.blocktimes(2,end);
        end
        cur_taxis = onsetT:stimres:endT;
        %compute binned spike count vectors at stimulus resolution
        time_bin_edges = [(cur_taxis(1)-stimres/2) (cur_taxis+stimres/2)];
        cur_spike_bin_vecs = zeros(10,length(cur_taxis));
        for cellid = 1:10            
            spikes_binned = histc(Blocks{blockid}.spktimes{cellid},time_bin_edges);
            spikes_binned(end) = [];
            cur_spike_bin_vecs(cellid,:) = spikes_binned;
        end
        spike_bin_vecs = [spike_bin_vecs; cur_spike_bin_vecs'];
    end
    spike_bin_vecs = spike_bin_vecs';
    
%     %compute binned spike count vectors at stimulus resolution
%     time_bin_edges = (recon_t(1)-stimres/2):stimres:(recon_t(end)+stimres/2);
%     recon_t2 = recon_t(1):stimres:recon_t(end);
%     cur_spike_bin_vecs = zeros(10,length(recon_t2));
%     for cellid = 1:10
%         spikes_binned = histc(Blocks{blockid}.spktimes{cellid},time_bin_edges);
%         spikes_binned(end) = [];
%         cur_spike_bin_vecs(cellid,:) = spikes_binned;
%     end
    
    
    cur_fix_image_set = zeros(n_used_fixs,size(ov_im_patch,2),size(ov_im_patch,3));
    cur_fix_spkcnt_set = zeros(n_used_fixs,10);
    right_fix_locs = zeros(n_used_fixs,2);
    left_fix_locs = zeros(n_used_fixs,2);
    right_fix_std = zeros(n_used_fixs,2);
    left_fix_std = zeros(n_used_fixs,2);
    for i = 1:n_used_fixs
        cur_set = fix_start_inds(used_fix_start_inds(i)) : (fix_start_inds(used_fix_start_inds(i))+round(min_fix_dur/stimres)); %use only window at begenning of fixation
        cur_avg_image = nanmean(ov_im_patch(cur_set,:,:));
        right_fix_locs(i,:) = mean(reye_interp(cur_set,:));
        left_fix_locs(i,:) = mean(leye_interp(cur_set,:));
        left_fix_std(i,:) = std(leye_interp(cur_set,:));
        right_fix_std(i,:) = std(reye_interp(cur_set,:));
        cur_fix_image_set(i,:,:) = cur_avg_image - nanmean(cur_avg_image(:)); %within fixation mean subtraction
        cur_fix_spkcnt_set(i,:) = sum(spike_bin_vecs(:,cur_set),2); %total spike counts within fixation window
    end
        
    cur_fix_start_times = sac_stop_times(fix_prev_sac_inds(used_fix_start_inds));
    cur_fix_stop_times = [sac_start_times(fix_prev_sac_inds(used_fix_start_inds(1:end-1))+1) recon_t(end)];
    
    X = [X; cur_fix_image_set];
    spk_cnts = [spk_cnts; cur_fix_spkcnt_set];
    all_sac_amps = [all_sac_amps; sac_amps(fix_prev_sac_inds(used_fix_start_inds))];
    all_sac_durs = [all_sac_durs; sac_durs(fix_prev_sac_inds(used_fix_start_inds))'];
    all_sac_vels = [all_sac_vels; sac_peakvel(fix_prev_sac_inds(used_fix_start_inds))];
    all_left_locs = [all_left_locs; left_fix_locs];
    all_right_locs = [all_right_locs; right_fix_locs];
    all_left_stds = [all_left_stds; left_fix_std];
    all_right_stds = [all_right_stds; right_fix_std];
    all_fix_start_times = [all_fix_start_times; cur_fix_start_times'];
    all_fix_stop_times = [all_fix_stop_times; cur_fix_stop_times'];
%     all_fix_stop_times = [all_fix_stop_times; recon_t(end)];
    blockids = [blockids; ones(length(used_fix_start_inds),1)*blockid];
    all_stim_filtered = [all_stim_filtered; is_stim_filtered(fix_start_inds(used_fix_start_inds))'];
    all_is_blinks = [all_is_blinks; blinks(used_fix_start_inds)'];
end

%%
disparity = max(all_left_locs - all_right_locs,[],2);
max_allowed_disparity = 1;
bad_fix = find(disparity > max_allowed_disparity);
X(bad_fix,:) = [];

spk_cnts(bad_fix,:) = [];
all_left_locs(bad_fix,:) = [];
all_right_locs(bad_fix,:) = [];
all_left_stds(bad_fix,:) = [];
all_right_stds(bad_fix,:) = [];
all_sac_amps(bad_fix) = [];
all_sac_durs(bad_fix) = [];
all_sac_vels(bad_fix) = [];
all_fix_start_times(bad_fix) = [];
all_fix_stop_times(bad_fix) = [];
blockids(bad_fix) = [];
all_stim_filtered(bad_fix) = [];
all_is_blinks(bad_fix) = [];

%% convert spike data to bin indices
spikebins = cell(10,1);
for cellid = 1:10
    un_spk_cnts = unique(spk_cnts(:,cellid));
    cur_spikebins = [];
    for i = 1:length(un_spk_cnts)
        cur_set = find(spk_cnts(:,cellid) == un_spk_cnts(i));
        cur_spikebins = [cur_spikebins; repmat(cur_set(:),un_spk_cnts(i),1)];
    end
    cur_spikebins = sort(cur_spikebins);
    spikebins{cellid} = [spikebins{cellid}; cur_spikebins ];
end

%%
Xmat = X;
SDIM = size(ov_im_patch,2);
kern_l = SDIM^2;
nmods = 1;
init_kerns = randn(kern_l,nmods);
init_signs = ones(nmods,1);
init_betas = 2*ones(nmods,1);

cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat
defmod.fsdim = SDIM^2;
defmod.pids = 1:defmod.fsdim;
defmod.h = 1;
defmod.SDIM = SDIM;
defmod.locLambda = 0;
defmod.lambda_dX = 1e6;
defmod.lambda_dT = 0;
defmod.lambda_L1x = 1e3;
basis = 'pix';

% cellid = 7;
for cellid = 1:10
    fprintf('Fitting cell %d of %d\n',cellid,10);
%     glm0 = createGLM_lexp(init_kerns,init_signs,init_betas,defmod,basis,[],[]);
    glm0 = createGLM_tlin(init_kerns,init_signs,defmod,basis,[],[]);
%     glm0 = createGLM_quad(init_kerns,init_signs,defmod,basis,[],[]);
    glm0.image_type = '2d';
    
    glm0 = normalizeRFs_full(glm0,Xmat);
    ml_sta_orig(cellid) = fitGLM_lexp(glm0,Xmat,spikebins{cellid},'tots');
end

%%
mod_fit(1,:) = ml_sta_orig;
used_cells = [1:8];
% used_cells = [7];
n_iter = 4;

Pix2Deg = 0.018837;
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;

max_pert = 10;
deltax_range = max_pert;
deltay_range = max_pert;
delta_x = (-deltax_range:deltax_range);
delta_y = (-deltay_range:deltay_range);
[dX,dY] = meshgrid(delta_x,delta_y);
dX = dX(:);
dY = dY(:);
xax = linspace(-SDIM/2,SDIM/2,SDIM)/Fsd;
yax = linspace(-SDIM/2,SDIM/2,SDIM)/Fsd;
[tempx,tempy] = meshgrid(xax,yax);
K = length(dX);
NT = size(Xmat,1);

ROI = find(tempx >= -1 & tempx <= 1 & tempy >= -1 & tempy <= 1);

all_avg_locs = 0.5*all_left_locs + 0.5*all_right_locs;
lerr = all_left_locs - all_avg_locs;
rerr = all_right_locs - all_avg_locs;

dXl = bsxfun(@minus,dX/Fsd,lerr(:,1)');
dYl = bsxfun(@minus,dY/Fsd,lerr(:,2)');
dXr = bsxfun(@minus,dX/Fsd,rerr(:,1)');
dYr = bsxfun(@minus,dY/Fsd,rerr(:,2)');

error_sigma = 0.4; %WAS 0.5
% ov_pen = exp(-(dXl.^2+dYl.^2+dXr.^2+dYr.^2)/(2*error_sigma^2));
ov_pen = exp(-((dX/Fsd).^2+(dY/Fsd).^2)/(2*error_sigma^2));

delta_sigma = 0.05; %WAS 0.2*Fsd
cur_distmat = pdist2([dX dY]/Fsd,[dX dY]/Fsd);
cur_Amat = exp(-cur_distmat.^2/(2*delta_sigma^2));
cur_Amat = bsxfun(@times,cur_Amat,ov_pen');
cur_Amat = bsxfun(@rdivide,cur_Amat,sum(cur_Amat,2));

dX_seqs = zeros(n_iter,NT);
dY_seqs = zeros(n_iter,NT);
left_dX_seqs = zeros(n_iter,NT);
left_dY_seqs = zeros(n_iter,NT);
right_dX_seqs = zeros(n_iter,NT);
right_dY_seqs = zeros(n_iter,NT);

%%
for n = 2:n_iter+1
    fprintf('Iteration %d of %d\n',n-1,n_iter);    
  
    %%
    clear k_pert poss_LLs
    for cellid = used_cells
        fprintf('Computing shifts: cell %d of %d\n',cellid,10);
        k_vec = get_k_mat(mod_fit(n-1,cellid));
        k_mat = squeeze(reshape(k_vec,[SDIM,SDIM,nmods]));
        
        k_pert = nan(nmods,length(delta_y),length(delta_x),length(k_vec));
        for xx = 1:length(delta_x)
            for yy = 1:length(delta_y)
                for i = 1:nmods
                    k_shifted = dist_shift2d(squeeze(k_mat(:,:,i)), delta_x(xx),2,0);
                    k_shifted = dist_shift2d(k_shifted,delta_y(yy),1,0);
                    k_pert(i,yy,xx,:) = k_shifted(:);
                end
            end
        end
        k_pert = reshape(k_pert,[nmods length(delta_x)*length(delta_y) kern_l]);
        n_perts = size(k_pert,2);
        k_pert = permute(k_pert,[3 2 1]);
        
        min_x = min(mod_fit(n-1,cellid).mods(1).nlx);
        poss_gen_funs = ones(NT,K)*mod_fit(n-1,cellid).const;
        for i = 1:nmods
            poss_spatial_outputs = Xmat(:,ROI)*squeeze(k_pert(ROI,:,i));
            
            if strcmp(mod_fit(n-1,cellid).mods(i).nltype,'lexp')
            %take logexp of internal genfun
            poss_spatial_outputs = 1/mod_fit(n-1,cellid).mods(1).beta*log(1+exp(mod_fit(n-1,cellid).mods(1).beta*poss_spatial_outputs)) - ...
                1/mod_fit(n-1,cellid).mods(1).beta*log(1+exp(mod_fit(n-1,cellid).mods(1).beta*min_x));
            elseif strcmp(mod_fit(n-1,cellid).mods(i).nltype,'lin')
                
            elseif strcmp(mod_fit(n-1,cellid).mods(i).nltype,'quad')
                poss_spatial_outputs = poss_spatial_outputs.^2;               
            end

            poss_gen_funs = poss_gen_funs + poss_spatial_outputs*mod_fit(n-1,cellid).mods(i).w;
        end
        poss_rates = log(1+exp(poss_gen_funs));
        
        poss_LLs(cellid,:,:) = repmat(spk_cnts(:,cellid),1,n_perts).*log(poss_rates) - poss_rates - ...
            repmat(log(factorial(spk_cnts(:,cellid))),1,n_perts);
    end
    %%
    B = squeeze(exp(sum(poss_LLs(used_cells,:,:),1))); %convert to overall likes
    Pi = ov_pen(:,1)'/sum(ov_pen(:,1));
    
    %NOW FOR VITERBI SEQUENCE ESTIMATION
    fprintf('Computing viterbi sequence\n');
    %initialize variables
    delta=zeros(NT,K); %delta is the maximized log probability as a function of time
    psi=zeros(NT,K); %Psi stores the most likely preceeding state
    %initialization
    delta(1,:) = log(Pi)+log(B(1,:));    % Eq. 105(a) Rabiner
    for t=2:NT
%         cur_distmat = pdist2([dXl(:,t-1) dYl(:,t-1) dXr(:,t-1) dYr(:,t-1)],[dXl(:,t) dYl(:,t) dXr(:,t) dYr(:,t)]);
%         cur_Amat = exp(-cur_distmat.^2/(2*delta_sigma^2));
%         cur_Amat = bsxfun(@times,cur_Amat,ov_pen(:,t)');
%         cur_Amat = bsxfun(@rdivide,cur_Amat,sum(cur_Amat,2));
        temp = repmat(delta(t-1,:),K,1) + log(cur_Amat);
        [delta(t,:),psi(t,:)] = max(temp,[],2);
        delta(t,:) = delta(t,:) + log(B(t,:));
    end
    
    % Backtracking Viterbi
    state_seq = zeros(NT,1);
    [llik_best,state_seq(NT)] = max(delta(NT,:));
    for t=NT-1:-1:1,
        state_seq(t) = psi(t+1,state_seq(t+1));
    end
    
    dX_seqs(n,:) = dX(state_seq)';
    dY_seqs(n,:) = dY(state_seq)';

    left_dX_seqs(n,:) = lerr(:,1) - dX_seqs(n,:)'/Fsd;
    left_dY_seqs(n,:) = lerr(:,2) - dY_seqs(n,:)'/Fsd;
    right_dX_seqs(n,:) = rerr(:,1) - dX_seqs(n,:)'/Fsd;
    right_dY_seqs(n,:) = rerr(:,2) - dY_seqs(n,:)'/Fsd;
    
    disp('DONE')
    
    %%
    Xmat_resh = reshape(Xmat,NT,SDIM,SDIM);
    for i = 1:NT
        Xmat_resh(i,:,:) = dist_shift2d(squeeze(Xmat_resh(i,:,:)),-dX(state_seq(i)),2,0);
        Xmat_resh(i,:,:) = dist_shift2d(squeeze(Xmat_resh(i,:,:)),-dY(state_seq(i)),1,0);
    end    
    Xmat_new = reshape(Xmat_resh,NT,SDIM^2);
    
    %%
    for cellid = 1:10
        fprintf('Fitting cell %d of %d\n',cellid,10);
        mod_fit(n,cellid) = fitGLM_lexp(mod_fit(n-1,cellid),Xmat_new,spikebins{cellid},'tots');
    end
    
end

%%
cd /Users/James/James_scripts/bruce
save modfits_dualeye_25_ds4_v2 mod_fit dsfrac Fsd state_seq *_seqs all_* blockids

%%

figure
it_num = n_iter;

xax = linspace(-SDIM/2,SDIM/2,SDIM)/Fsd;
yax = linspace(-SDIM/2,SDIM/2,SDIM)/Fsd;

% xl = [-2 2];
% yl = [-2 2];
xl = [-1. 1.];
yl = [-1. 1.];

% for i = 1:10
%     k_vec = get_k_mat(mod_fit(it_num,i));
%     minz(i) = min(min(k_vec(:)),min(k_vec(:)));
%     maxz(i) = max(max(k_vec(:)),max(k_vec(:)));
% end
% minz = min(minz)/2;
% maxz = max(maxz)/2;

for cellid = 1:5
    k_vec = get_k_mat(mod_fit(1,cellid));
    k_mat = reshape(k_vec,SDIM,SDIM);
    
    k_vec = get_k_mat(mod_fit(it_num,cellid));
    k_matr = reshape(k_vec,SDIM,SDIM);
    
    minz = min(min(k_mat(:)),min(k_matr(:)));
    maxz = max(max(k_mat(:)),max(k_matr(:)));
    
    subplot(5,4,(cellid-1)*4+1)
    imagesc(xax,yax,k_mat); set(gca,'ydir','normal');
    caxis([minz maxz]);
    colormap(gray);
    xlim(xl); ylim(yl);
    
    subplot(5,4,(cellid-1)*4+2)
    imagesc(xax,yax,k_matr); set(gca,'ydir','normal');
    caxis([minz maxz]);
    colormap(gray);
    xlim(xl); ylim(yl);
    
end

for cellid = 6:10
    k_vec = get_k_mat(mod_fit(1,cellid));
    k_mat = reshape(k_vec,SDIM,SDIM);
    
    k_vec = get_k_mat(mod_fit(it_num,cellid));
    k_matr = reshape(k_vec,SDIM,SDIM);
    
    minz = min(min(k_mat(:)),min(k_matr(:)));
    maxz = max(max(k_mat(:)),max(k_matr(:)));
    
    subplot(5,4,(cellid-6)*4+3)
    imagesc(xax,yax,k_mat); set(gca,'ydir','normal');
    caxis([minz maxz]);
    colormap(gray);
    xlim(xl); ylim(yl);
    
    subplot(5,4,(cellid-6)*4+4)
    imagesc(xax,yax,k_matr); set(gca,'ydir','normal');
    caxis([minz maxz]);
    colormap(gray);
    xlim(xl); ylim(yl);
    
end

%%
nmods = 2;
init_kerns = randn(kern_l,nmods);
init_signs = ones(nmods,1);
% init_signs(end) = -1;
init_betas = 2*ones(nmods,1);

defmod.fsdim = SDIM^2;
defmod.pids = 1:defmod.fsdim;
defmod.h = 1;
defmod.SDIM = SDIM;
defmod.locLambda = 0;
defmod.lambda_dX = 1e6;
defmod.lambda_dT = 0;
% defmod.lambda_L1x = 1e3;
defmod.lambda_L1x = 500;
basis = 'pix';

cellids = 1:10;
for cellid = cellids
    fprintf('Fitting cell %d of %d\n',cellid,10);
%     glm0 = createGLM_lexp(init_kerns,init_signs,init_betas,defmod,basis,[],[]);
%     glm0.image_type = '2d';
%     glm0 = normalizeRFs_full(glm0,Xmat_new);
%     mod_fit3(cellid) = fitGLM_lexp(glm0,Xmat_new,spikebins{cellid},'tots');
 
        glm0 = createGLM_tlin(init_kerns,init_signs,defmod,basis,[],[]);
    glm0.image_type = '2d';
    glm0 = normalizeRFs_full(glm0,Xmat_new);
    mod_tlin(cellid) = fitGLM_lexp(glm0,Xmat_new,spikebins{cellid},'tots');
%     mod_tlin(cellid) = fitWeights_lexp_2(mod_tlin(cellid),Xmat_new*get_k_mat(mod_tlin(cellid)),spikebins{cellid});
%     mod_tlin(cellid) = fitGLM_lexp(mod_tlin(cellid),Xmat_new,spikebins{cellid},'tots');

%     glm0 = createGLM_quad(init_kerns,init_signs,defmod,basis,[],[]);
%     glm0.image_type = '2d';
%     glm0 = normalizeRFs_full(glm0,Xmat_new);
%     mod_quad3(cellid) = fitGLM_lexp(glm0,Xmat_new,spikebins{cellid},'tots');
end


%%
disparity = all_right_locs - all_left_locs;
plot(disparity)
axis tight
shg
set(gca,'fontsize',14)
box off
xlabel('Fixation number','fontsize',16)
ylabel('Disparity (degrees)','fontsize',16)
legend('Horizontal','Vertical')

%%
disparity = all_right_locs - all_left_locs;
plot(right_dX_seqs(n_iter,:),'b.-')
hold on
plot(left_dX_seqs(n_iter,:),'r.-')
plot(dX_seqs(n_iter,:)/Fsd,'g.-')
% plot(right_dX_seqs(end,:)-left_dX_seqs(end,:),'g')
plot(disparity(:,1),'k.-')
set(gca,'fontsize',14)
xlabel('Fixation number','fontsize',16)
ylabel('Eye error (degrees)','fontsize',16)
axis tight
% xlim([2060 2100])

%%
disparity = all_right_locs - all_left_locs;
plot(right_dY_seqs(n_iter,:),'b.-')
hold on
plot(left_dY_seqs(n_iter,:),'r.-')
plot(dY_seqs(n_iter,:)/Fsd,'g.-')
% plot(right_dX_seqs(end,:)-left_dX_seqs(end,:),'g')
plot(disparity(:,2),'k.-')
set(gca,'fontsize',14)
xlabel('Fixation number','fontsize',16)
ylabel('Eye error (degrees)','fontsize',16)
axis tight
% xlim([2060 2100])