clear all
% close all

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-7 7;-7 7];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;
% nlags = 4;

X = [];
spk_cnts = [];
all_sac_amps = [];
all_sac_vels = [];
all_sac_durs = [];
all_fix_start_times = [];
all_fix_stop_times = [];
all_stim_filtered = [];
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
%     if ~exist([sname '.mat'],'file')
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
%     else
%         load(sname)
%     end
    
    % load patch video data
    cd ~/Data/bruce/2_27_12/stimrecon/
    %     sname = sprintf('image_patch_block%d_reye_50_dsf4',blockid);
    %     sname = sprintf('image_patch_block%d_leye_50_dsf4',blockid);
    sname = sprintf('image_patch_block%d_avgeye_25_dsf4',blockid);
    load(sname)
    
    % interpolate eye signal onto stimulus time axis
%     eye_interp = interp1(eyets(1:end-1),reye_pos,recon_t);
%         eye_interp = interp1(eyets(1:end-1),leye_pos,recon_t);
        eye_interp = interp1(eyets(1:end-1),leye_pos+reye_pos,recon_t)/2;
    
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
    
    fix_dx = eye_interp(fix_stop_inds,1) - eye_interp(fix_start_inds,1);
    fix_dy = eye_interp(fix_stop_inds,2) - eye_interp(fix_start_inds,2);
    drift_amps = sqrt(fix_dx.^2+fix_dy.^2);
    
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
    in_window = (eye_interp(:,1) >= accept_window(1,1) & eye_interp(:,1) <= accept_window(1,2) & ...
        eye_interp(:,2) >= accept_window(2,1) & eye_interp(:,2) <= accept_window(2,2));
    
    off_screen = max(any(isnan(ov_im_patch),3),[],2);
    
    
    used_fix_start_inds = find(in_window(fix_start_inds)==1 & off_screen(fix_start_inds) == 0 & fix_skips'==0 & blinks'==0);
    %     used_fix_start_inds = find(in_window(fix_start_inds)==1 & is_stim_filtered(fix_start_inds)' == 0);
    n_used_fixs = length(used_fix_start_inds);
    fprintf('Out of window: %d,  off-screen:  %d, gray-stim: %d,  blinks: %d\n',sum(in_window(fix_start_inds)==0),sum(off_screen(fix_start_inds)==1),sum(fix_skips==1),sum(blinks==1));
    fprintf('analyzing %d of %d fixations\n',n_used_fixs,length(fix_start_inds));
    
    %compute binned spike count vectors at stimulus resolution
    time_bin_edges = [(recon_t(1)-stimres/2) (recon_t+stimres/2)];
    spike_bin_vecs = zeros(10,length(recon_t));
    for cellid = 1:10
        spike_times = Blocks{blockid}.spktimes{cellid};
        
        spikes_binned = histc(spike_times,time_bin_edges);
        spikes_binned(end) = [];
        spike_bin_vecs(cellid,:) = spikes_binned;
    end
    
    cur_fix_image_set = zeros(n_used_fixs,size(ov_im_patch,2),size(ov_im_patch,3));
    cur_fix_spkcnt_set = zeros(n_used_fixs,10);
    for i = 1:n_used_fixs
        %        cur_set = fix_start_inds(used_fix_start_inds(i)):fix_stop_inds(used_fix_start_inds(i)); %use entire fixation
        cur_set = fix_start_inds(used_fix_start_inds(i)) : (fix_start_inds(used_fix_start_inds(i))+round(min_fix_dur/stimres)); %use only window at begenning of fixation
        cur_avg_image = nanmean(ov_im_patch(cur_set,:,:));
        
        %        cur_fix_image_set(i,:,:) = cur_avg_image;
%         cur_fix_image_set(i,:,:) = cur_avg_image; %within fixation mean subtraction
        cur_fix_image_set(i,:,:) = cur_avg_image - nanmean(cur_avg_image(:)); %within fixation mean subtraction
        
        cur_fix_spkcnt_set(i,:) = sum(spike_bin_vecs(:,cur_set),2); %total spike counts within fixation window
    end
    
    X = [X; cur_fix_image_set];
    spk_cnts = [spk_cnts; cur_fix_spkcnt_set];
    all_sac_amps = [all_sac_amps; sac_amps(fix_prev_sac_inds(used_fix_start_inds))];
    all_sac_durs = [all_sac_durs; sac_durs(fix_prev_sac_inds(used_fix_start_inds))'];
    all_sac_vels = [all_sac_vels; sac_peakvel(fix_prev_sac_inds(used_fix_start_inds))];
    all_fix_start_times = [all_fix_start_times; sac_stop_times(fix_prev_sac_inds(used_fix_start_inds))'];
    all_fix_stop_times = [all_fix_stop_times; sac_start_times(fix_prev_sac_inds(used_fix_start_inds(1:end-1))+1)'];
    all_fix_stop_times = [all_fix_stop_times; recon_t(end)];
    blockids = [blockids; ones(length(used_fix_start_inds),1)*blockid];
    all_stim_filtered = [all_stim_filtered; is_stim_filtered(fix_start_inds(used_fix_start_inds))'];
end

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
% SDIM = 33;
%
% cellid = 3;
% temp = spk_cnts(:,cellid)'*reshape(X,[size(X,1) size(X,2)*size(X,3)]);
% temp = reshape(temp,size(X,2),size(X,3))/sum(spk_cnts(:,cellid));
%
% sta = temp - squeeze(mean(X));
%
% imagesc(sta);colormap(gray);
% % shg
%
%%
Xmat = reshape(X,[size(X,1) size(X,2)*size(X,3)]);
SDIM = 33;
kern_l = SDIM^2;
nmods = 1;
% init_kerns = sta(:);
% init_signs = 1;
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

for cellid = 1:10
    fprintf('Fitting cell %d of %d\n',cellid,10);
    glm0 = createGLM_lexp(init_kerns,init_signs,init_betas,defmod,basis,[],[]);
    glm0.image_type = '2d';
    
    glm0 = normalizeRFs_full(glm0,Xmat);
    ml_sta_orig(cellid) = fitGLM_lexp(glm0,Xmat,spikebins{cellid},'tots');
end

%%
n_iter = 5;

Pix2Deg = 0.018837;
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;

deltax_range = 10;
deltay_range = 10;
delta_x = -deltax_range:deltax_range;
delta_y = -deltay_range:deltay_range;
[dX,dY] = meshgrid(delta_x,delta_y);
dX = dX(:);
dY = dY(:);
xax = linspace(-SDIM/2,SDIM/2,SDIM)/Fsd;
yax = linspace(-SDIM/2,SDIM/2,SDIM)/Fsd;
[tempx,tempy] = meshgrid(xax,yax);
K = length(dX);
NT = size(Xmat,1);

ROI = find(tempx >= -1.0 & tempx <= 1.0 & tempy >= -1.0 & tempy <= 1.0);

dist_mat = squareform(pdist([dX dY]));
A = zeros(K);
% A_self = 0.6;
% A_other = 1-A_self;
% A(dist_mat == 0) = A_self;
% A(dist_mat == 1) = A_other;

sigma1 = 5;
sigma2 = 0.5;
w1 = 0.25;
w2 = 0.75;
A1 = 1/sqrt(2*pi*sigma1^2)*exp(-dist_mat.^2/(2*sigma1^2));
A2 = 1/sqrt(2*pi*sigma2^2)*exp(-dist_mat.^2/(2*sigma2^2));
A = w1*A1 + w2*A2;

% for i = 1:K
%     A(i,:) = A(i,:)/sum(A(i,:));
% end
A = A';

eps = 1e-200;

% compute posterior over hidden states

mod_fit(1,:) = ml_sta_orig;
cur_Xmat = Xmat;

delta_X_seqs = zeros(n_iter+1,size(Xmat,1));
delta_Y_seqs = zeros(n_iter+1,size(Xmat,1));
cur_state_seq = ones(size(Xmat,1),1)*(K+1)/2;

cur_Amat = shiftdim(repmat(A,[1 1 NT]),1);

for n = 2:n_iter+1
    fprintf('Iteration %d of %d\n',n-1,n_iter);
    %%
        
    for 1:NT
dist_mat = squareform(pdist([dX dY]));
       cur_pdist = squareform(pdist2( 
    end
    
    
    clear k_pert poss_LLs
    for cellid = 1:10
        fprintf('Computing shifts: cell %d of %d\n',cellid,10);
        k_vec = get_k_mat(mod_fit(n-1,cellid));
        k_mat = reshape(k_vec,SDIM,SDIM);
        %         spatial_output = cur_Xmat*k_vec;
        
        k_pert = nan(length(delta_y),length(delta_x),length(k_vec));
        for xx = 1:length(delta_x)
            for yy = 1:length(delta_y)
                k_shifted = dist_shift2d(k_mat, delta_x(xx), 2,0);
                k_shifted = dist_shift2d(k_shifted,delta_y(yy),1,0);
                k_pert(yy,xx,:) = k_shifted(:);
            end
        end
        k_pert = reshape(k_pert,length(delta_x)*length(delta_y),size(cur_Xmat,2));
        n_perts = size(k_pert,1);
        
        min_x = min(mod_fit(n-1,cellid).mods(1).nlx);
        poss_spatial_outputs = cur_Xmat(:,ROI)*k_pert(:,ROI)';
        %take logexp of internal genfun
        poss_spatial_outputs = 1/mod_fit(n-1,cellid).mods(1).beta*log(1+exp(mod_fit(n-1,cellid).mods(1).beta*poss_spatial_outputs)) - ...
            1/mod_fit(n-1,cellid).mods(1).beta*log(1+exp(mod_fit(n-1,cellid).mods(1).beta*min_x));
        
        poss_gen_funs = poss_spatial_outputs+mod_fit(n-1,cellid).const;
        poss_rates = log(1+exp(poss_gen_funs));
        
        poss_LLs(cellid,:,:) = repmat(spk_cnts(:,cellid),1,n_perts).*log(poss_rates) - poss_rates - ...
            repmat(log(factorial(spk_cnts(:,cellid))),1,n_perts);
    end
    %%
    B = squeeze(exp(sum(poss_LLs,1))); %convert to overall likes
        
    Pi = ones(1,K);
    
    %initialize forward and backward messages
    alpha=zeros(NT,K);
    beta=zeros(NT,K);
    gamma=zeros(NT,K);
    
    scale=zeros(NT,1); %initialize rescaling parameters
    
    %compute rescaled forward messages
    alpha(1,:)=Pi.*B(1,:);
    scale(1)=sum(alpha(1,:));
    alpha(1,:)=alpha(1,:)/scale(1);
    upA_padded = A;
    npads = K;
    upA_padded = cat(1,upA_padded,zeros(npads,size(upA_padded,2)));
    upA_padded = cat(1,zeros(npads,size(upA_padded,2)),upA_padded);
    upA_padded = cat(2,upA_padded,zeros(size(upA_padded,1),npads));
    upA_padded = cat(2,zeros(size(upA_padded,1),npads),upA_padded);
    x_used = (1:K)+npads;
    y_used = (1:K)+npads;
    for t=2:NT
%         alpha(t,:)=(alpha(t-1,:)*squeeze(cur_Amat(t,:,:))).*B(t,:);
        cur_yused = y_used + dX(cur_state_seq(t-1))*length(delta_y) + dY(cur_state_seq(t-1));
        cur_xused = x_used + dX(cur_state_seq(t))*length(delta_y) + dY(cur_state_seq(t));
        alpha(t,:)=(alpha(t-1,:)*upA_padded(cur_yused,cur_xused)).*B(t,:);
        scale(t) = sum(alpha(t,:));
        alpha(t,:)=alpha(t,:)/scale(t);
    end
    
    %compute rescaled backward messages
    beta(NT,:)=ones(1,K)/scale(NT);
    for t=NT-1:-1:1
        cur_yused = y_used + dX(cur_state_seq(t))*length(delta_y) + dY(cur_state_seq(t));
        cur_xused = x_used + dX(cur_state_seq(t+1))*length(delta_y) + dY(cur_state_seq(t+1));
        beta(t,:)=(beta(t+1,:).*B(t+1,:))*(upA_padded(cur_yused,cur_xused)')/scale(t);
% beta(t,:)=(beta(t+1,:).*B(t+1,:))*(squeeze(cur_Amat(t+1,:,:))')/scale(t);
    end
    
    %compute posteriors over hidden states
    gamma=(alpha.*beta);
    gamma = gamma./repmat(sum(gamma,2),1,K);
    
    %     %compute chi (posterior of two consecutive hidden states)
    %     chi=zeros(NT-1,K*K);
    %     for t=1:NT-1
    %         temp=A.*(alpha(t,:)' * (beta(t+1,:).*B(t+1,:)));
    %         chi(t,:)=temp(:)'/sum(temp(:));
    %     end
    
    lscale = log(scale);
    loglik(n) = sum(lscale); %rabiner eq(103), scale is defined as inverse here though
    
    %NOW FOR VITERBI SEQUENCE ESTIMATION
    fprintf('Computing viterbi sequence\n');
    %initialize variables
    delta=zeros(NT,K); %delta is the maximized log probability as a function of time
    psi=zeros(NT,K); %Psi stores the most likely preceeding state
    %initialization
    delta(1,:) = log(Pi)+log(B(1,:));    % Eq. 105(a) Rabiner
    
    for t=2:NT
        %     for k = 1:K
        %         temp = delta(t-1,:) + log(A(:,k))';
        %         [delta(t,k),psi(t,k)] = max(temp);
        %     end
     
%         temp = repmat(delta(t-1,:),K,1) + squeeze(log(cur_Amat(t,:,:)));
        
        cur_yused = y_used + dX(cur_state_seq(t-1))*length(delta_y) + dY(cur_state_seq(t-1));
        cur_xused = x_used + dX(cur_state_seq(t))*length(delta_y) + dY(cur_state_seq(t));
        temp = repmat(delta(t-1,:),K,1) + log(upA_padded(cur_yused,cur_xused));
%         temp(isinf(upA_padded(cur_yused,cur_xused)) = -Inf;
        %         temp = repmat(delta(t-1,:),K,1) + log(A)';
        [delta(t,:),psi(t,:)] = max(temp,[],2);
        delta(t,:) = delta(t,:) + log(B(t,:));
    end
    
    % Backtracking Viterbi
    state_seq = zeros(size(cur_state_seq));
    [llik_best,state_seq(NT)] = max(delta(NT,:));
    for t=NT-1:-1:1,
        state_seq(t) = psi(t+1,state_seq(t+1));
    end
    
    delta_X_seqs(n:end,:) = bsxfun(@plus,delta_X_seqs(n:end,:),dX(state_seq)');
    delta_Y_seqs(n:end,:) = bsxfun(@plus,delta_Y_seqs(n:end,:),dY(state_seq)');
    cur_state_seq = cur_state_seq + dX(state_seq)*length(delta_y) + dY(state_seq); %transform into current state sequence   
    disp('DONE')
    
    %%
    Xmat_resh = reshape(cur_Xmat,size(cur_Xmat,1),SDIM,SDIM);
    for i = 1:size(cur_Xmat,1)
        Xmat_resh(i,:,:) = dist_shift2d(squeeze(Xmat_resh(i,:,:)),-dX(state_seq(i)),2,0);
        Xmat_resh(i,:,:) = dist_shift2d(squeeze(Xmat_resh(i,:,:)),-dY(state_seq(i)),1,0);
    end
    Xmat_new = reshape(Xmat_resh,size(cur_Xmat,1),size(cur_Xmat,2));
    
    %%
    for cellid = 1:10
        fprintf('Fitting cell %d of %d\n',cellid,10);
        mod_fit(n,cellid) = fitGLM_lexp(mod_fit(n-1,cellid),Xmat_new,spikebins{cellid},'tots');
    end
    cur_Xmat = Xmat_new;
end

cd /Users/James/James_scripts/bruce
save modfits_avgeye_25_ds4 mod_fit dsfrac Fsd delta_X_seqs delta_Y_seqs gamma loglik all_sac_* all_fix_* blockids all_stim_*


%%
figure
it_num = 6;

xax = linspace(-SDIM/2,SDIM/2,SDIM)/Fsd;
yax = linspace(-SDIM/2,SDIM/2,SDIM)/Fsd;

% xl = [-2 2];
% yl = [-2 2];
xl = [-1.2 1.2];
yl = [-1.2 1.2];

for cellid = 1:5
    k_vec = get_k_mat(ml_sta_orig(cellid));
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
    k_vec = get_k_mat(ml_sta_orig(cellid));
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
% nmods = 3;
% cellid = 3;
% init_kerns = randn(kern_l,nmods);
% init_signs = [1 1 -1];
% init_betas = 2*ones(nmods,1);
% % glm0 = createGLM_lexp(init_kerns,init_signs,init_betas,defmod,basis,[],[]);
% glm0 = createGLM_quad(init_kerns,init_signs,defmod,basis,[],[]);
% glm0.image_type = '2d';
% 
% glm0 = normalizeRFs_full(glm0,cur_Xmat);
% ml_multifilt = fitGLM_lexp(glm0,cur_Xmat,spikebins{cellid},'tots');

%%
