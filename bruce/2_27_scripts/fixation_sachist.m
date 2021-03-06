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
all_sac_locs = [];
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
    %     sname = sprintf('image_patch_block%d_reye_50_dsf4',blockid);
    %     sname = sprintf('image_patch_block%d_leye_50_dsf4',blockid);
    sname = sprintf('image_patch_block%d_leye_25_dsf4',blockid);
    load(sname)
    
    % interpolate eye signal onto stimulus time axis
    %     eye_interp = interp1(eyets(1:end-1),reye_pos,recon_t);
    eye_interp = interp1(eyets(1:end-1),leye_pos,recon_t);
    %         eye_interp = interp1(eyets(1:end-1),leye_pos+reye_pos,recon_t)/2;
    
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
    cur_fix_locs = zeros(n_used_fixs,2);
    for i = 1:n_used_fixs
        %        cur_set = fix_start_inds(used_fix_start_inds(i)):fix_stop_inds(used_fix_start_inds(i)); %use entire fixation
        cur_set = fix_start_inds(used_fix_start_inds(i)) : (fix_start_inds(used_fix_start_inds(i))+round(min_fix_dur/stimres)); %use only window at begenning of fixation
        cur_avg_image = nanmean(ov_im_patch(cur_set,:,:));
        cur_fix_locs(i,:) = mean(eye_interp(cur_set,:));
        
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
    all_sac_locs = [all_sac_locs; cur_fix_locs];
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

backwin = 0.05;
forwardwin = 0.75;
T = backwin+forwardwin+dt;
Nbases = 25;
dt = 0.025;
kb_taxis = -backwin:dt:forwardwin;
% ymat = keat_basis(T,1,Nbases,dt,T);


time_since_fix = [];
fix_inds = [];
spikebins = [];
fix_trig_avg = zeros(10,length(kb_taxis));
fix_trig_pts = zeros(10,length(kb_taxis));
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    %compute binned spike count vectors at high temporal resolution
    recon_dt = recon_t(1):dt:recon_t(end);
    time_bin_edges = [(recon_dt(1)-dt/2) (recon_dt+dt/2)];
    spike_bin_vecs = zeros(10,length(recon_dt));
    for cellid = 1:10
        spike_times = Blocks{blockid}.spktimes{cellid};
        
        spikes_binned = histc(spike_times,time_bin_edges);
        spikes_binned(end) = [];
        spike_bin_vecs(cellid,:) = spikes_binned;
    end
    
    %construct x matrix for time since fixation
    cur_set = find(blockids == blockid);
    for i = 1:length(cur_set)
        cur_fix_inds = find(recon_dt > all_fix_start_times(cur_set(i))-backwin & recon_dt <= all_fix_stop_times(cur_set(i)));
        time_since_fix = [time_since_fix; (1:length(cur_fix_inds))'];
        fix_inds = [fix_inds; ones(length(cur_fix_inds),1)*cur_set(i)];
        
        for cellid = 1:10
            
            ep = min(length(kb_taxis),length(cur_fix_inds));
            fix_trig_avg(cellid,1:ep) = fix_trig_avg(cellid,1:ep) + spike_bin_vecs(cellid,cur_fix_inds(1:ep));
            fix_trig_pts(cellid,1:ep) = fix_trig_pts(cellid,1:ep) + ones(1,ep);
        end
        spikebins = [spikebins; spike_bin_vecs(:,cur_fix_inds)'];
        
    end
    
    %     time_since_fix(time_since_fix > length(kb_taxis)) = nan; %saturate at T
    
end
bad_pts = find(time_since_fix > length(kb_taxis));
time_since_fix(bad_pts) = [];
fix_inds(bad_pts) = [];
spikebins(bad_pts,:) = [];

fix_trig_rate = fix_trig_avg./fix_trig_pts/dt;

%%
clear bestk bestk_nospace
Xmat = tbrep(time_since_fix,1:length(kb_taxis));

for cellid = 1:10
    
    cur_mod = ml_sta_orig(cellid);
    k_mat = get_k_mat(cur_mod);
    %compute spatial model output.
    spatial_out = reshape(X,size(X,1),size(X,2)*size(X,3))*k_mat(:);
    %pass through internal NL
    min_x = min(cur_mod.mods(1).nlx);
    if strcmp(cur_mod.mods(1).nltype,'lexp')
        spatial_out = 1/cur_mod.mods(1).beta*log(1+exp(cur_mod.mods(1).beta*spatial_out)) - ...
            1/cur_mod.mods(1).beta*log(1+exp(cur_mod.mods(1).beta*min_x));
    end
    
    spatial_component = spatial_out(fix_inds);
    
    % Xmat = zeros(length(time_since_fix),Nbases);
    % for i = 1:Nbases
    %    Xmat(:,i) = ymat(i,time_since_fix);
    % end
    
    
    %%
    un_spk_cnts = unique(spikebins(:,cellid));
    cur_spikebins = [];
    for i = 1:length(un_spk_cnts)
        cur_set = find(spikebins(:,cellid) == un_spk_cnts(i));
        cur_spikebins = [cur_spikebins; repmat(cur_set(:),un_spk_cnts(i),1)];
    end
    spkbs = sort(cur_spikebins);
    
    %%
    
    K0 = zeros(length(kb_taxis),1);
    K0 = [K0; cur_mod.const]; %tack on the constant term
    lamrange = [0 1 length(kb_taxis)];
%     lamrange = [];
    lamrange2 = [50 1 length(kb_taxis)];
    Pcon = [];
    Pmono = [];
    hold_const = [];
    NLtype = 0;
    llist = [];
    silent = 0;
    % [fitp,grad] = GLMsolve_jmm( Xmat', spkbs, K0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype);
    options.Display = 'iter';
    SX = Xmat(spkbs,:);
    [bestk(cellid,:),LL(cellid),exitflag,output] = minFunc(@(K) LLelog_tempmod(K,Xmat,SX,spatial_component,spatial_component(spkbs),lamrange,lamrange2,llist,hold_const,kb_taxis),K0,options);
    
    spatial_component = zeros(size(spatial_component));
    mrate = length(spkbs)/size(Xmat,1);
    C0 = log(1+exp(mrate));
    K0 = [log(exp(fix_trig_rate(cellid,:)*dt-1))-C0 C0]';
    [bestk_nospace(cellid,:),LL_nospace(cellid),exitflag,output] = minFunc(@(K) LLelog_tempmod(K,Xmat,SX,spatial_component,spatial_component(spkbs),lamrange,lamrange2,llist,hold_const,kb_taxis),K0,options);
end

%%
cellid = 5;
close all

figure
plot(kb_taxis,fix_trig_rate(cellid,:),'r')

figure
 plot(kb_taxis,bestk(cellid,1:end-1)-mean(bestk(cellid,1:end-1)),'r')
hold on
plot(kb_taxis,bestk_nospace(cellid,1:end-1)-mean(bestk_nospace(cellid,1:end-1)),'b')

%%
% test_Xmat = tbrep(1:length(kb_taxis),1:length(kb_taxis));
% testK = bestk_nospace(cellid,:);
% k = testK(1:end-1)';
% b = testK(end);
% 
% kx = test_Xmat*k+b+zeros(length(kb_taxis),1);
% 
% ekx = exp(kx);
% r = log(1+ekx);
