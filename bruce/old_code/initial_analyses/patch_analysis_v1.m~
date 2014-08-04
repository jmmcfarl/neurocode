clear all
close all

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-5 5;-5 5];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
min_fix_dur = 0.2;
nlags = 4;

X = [];
spikebins = cell(10,1);
full_t = [];
full_isfiltered = [];
full_blockid = [];

norm_fac = 55; %approximate std dev of pixel intensities

% blockid = 1;
for blockid = 1:4
    
    fprintf('Processing block %d of %d...\n',blockid,4);
    
    block_times = Blocks{blockid}.blocktimes;
    stim_times = Blocks{blockid}.stimtime;
    n_stims = length(stim_times);
    
    cd ~/Data/bruce/2_27_12/saccades/
    load(sprintf('lemM232.5%d.em.sac.mat',blockid))
    
    
    % identify saccade start and stop times
    EyeStartT = Expt.Trials.Start/10000; % time of first eye sample
    EyeEndT = Expt.Trials.End/10000; % time of last eye sample
    Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
    eyets = EyeStartT:Eyedt:EyeEndT; %eye tracking time axis (sec)
    
    reye_pos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
    leye_pos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];
    
    %if we haven't already computed saccade times do so now
    sname = sprintf('sac_times_cor_block%d',blockid);
    if ~exist([sname '.mat'],'file')
        fprintf('Computing saccade times\n');
        avg_eyepos = (reye_pos + leye_pos)/2;
        clear sm_avg_eyepos eye_vel
        sm_avg_eyepos(:,1) = smooth(avg_eyepos(:,1),3);
        sm_avg_eyepos(:,2) = smooth(avg_eyepos(:,2),3);
        eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/Eyedt;
        eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/Eyedt;
        eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
        
        sac_inds = find(eye_speed(1:end-1) < sac_eyespeed & eye_speed(2:end) > sac_eyespeed);
        
        sac_start_inds = nan(size(sac_inds));
        sac_stop_inds = nan(size(sac_inds));
        for i = 1:length(sac_inds)
            temp = find(eye_speed(1:sac_inds(i)) < thresh_eyespeed,1,'last');
            if ~isempty(temp)
                sac_start_inds(i) = temp;
            end
            temp = find(eye_speed(sac_inds(i):end) < thresh_eyespeed,1,'first');
            if ~isempty(temp)
                sac_stop_inds(i) = sac_inds(i)+temp-1;
            end
        end
        
        sac_vec = zeros(size(reye_pos,1),1);
        for i = 1:length(sac_start_inds)
            if ~isnan(sac_start_inds(i)) & ~isnan(sac_stop_inds(i))
                sac_vec(sac_start_inds(i):sac_stop_inds(i)) = 1;
            end
        end
        sac_vec(1) = 0; sac_vec(end) = 0;
        sac_start_indsn = find(sac_vec(1:end-1) == 0 & sac_vec(2:end) == 1);
        sac_stop_indsn = find(sac_vec(1:end-1) == 1 & sac_vec(2:end) == 0);
        if length(sac_start_indsn) ~= length(sac_stop_indsn)
            error('saccade mis-alignment');
        end
        
        sac_start_times = eyets(sac_start_indsn);
        sac_stop_times = eyets(sac_stop_indsn);
        
        cd ~/Data/bruce/2_27_12/saccades/
        save(sname,'sac_start_times','sac_stop_times','sac_start_inds','sac_stop_inds','sac_vec')
    else
        load(sname)
    end
    
    % load patch video data
    cd ~/Data/bruce/2_27_12/stimrecon/
%     sname = sprintf('image_patch_block%d_reye_50_dsf4',blockid);
%     sname = sprintf('image_patch_block%d_leye_50_dsf4',blockid);
    sname = sprintf('image_patch_block%d_leye_50_dsf4',blockid);
    load(sname)
    
    % interpolate eye signal onto stimulus time axis
%     eye_interp = interp1(eyets(1:end-1),reye_pos,recon_t);
    eye_interp = interp1(eyets(1:end-1),leye_pos,recon_t);
    
    sac_vec_interp = zeros(size(recon_t));
    for i = 1:length(sac_start_times)
        cur_set = find(recon_t >= sac_start_times(i) & recon_t <= sac_stop_times(i));
        sac_vec_interp(cur_set) = 1;
    end
    
    fix_start_inds = 1+[0 find(sac_vec_interp(1:end-1) == 1 & sac_vec_interp(2:end) == 0)];
    fix_stop_inds = 1+[find(sac_vec_interp(1:end-1) == 0 & sac_vec_interp(2:end) == 1) length(recon_t)-1];
    if length(fix_start_inds) ~= length(fix_stop_inds)
        error('Fixation mis-alignment');
    end
    fix_durs = (fix_stop_inds - fix_start_inds)*stimres;
    too_short = find(fix_durs < min_fix_dur);
    
    fprintf('%d of %d fixations too short\n',length(too_short),length(fix_durs));
    fix_start_inds(too_short) = []; fix_stop_inds(too_short) = [];
    
    fix_vec = zeros(size(recon_t));
    for i = 1:length(fix_start_inds)
        fix_vec(fix_start_inds(i):fix_stop_inds(i)-1) = 1;
    end
    
    % create vector determining whether current stimulus was high-pass filtered
    filt_stims = mod(0:500,4) + 1 >2;
    is_stim_filtered = filt_stims(stim_num);
    
    % find times where eye-position is within central window
    in_window = (eye_interp(:,1) >= accept_window(1,1) & eye_interp(:,1) <= accept_window(1,2) & ...
        eye_interp(:,2) >= accept_window(2,1) & eye_interp(:,2) <= accept_window(2,2));
    
    if max(isnan(reshape(ov_im_patch(in_window,:,:),1,numel(ov_im_patch(in_window,:,:)))))
        error('RF patch goes off-screen');
    end
    
    % these are the time samples we are going to use in analysis
    used_times = find(in_window == 1 & fix_vec' == 1);
%     used_times = find(in_window == 1 & fix_vec' == 1 & is_stim_filtered'==1);
%     used_times = find(in_window == 1 & fix_vec' == 1 & is_stim_filtered'==0);
    fprintf('%d of %d time bins included\n',length(used_times),length(recon_t));
    
    %
    full_t = [full_t; recon_t(used_times)'];
    full_isfiltered = [full_isfiltered; is_stim_filtered(used_times)'];
    full_blockid = [full_blockid; ones(length(used_times),1)*blockid];
    %
    time_bin_edges = [(recon_t(1)-stimres/2) (recon_t+stimres/2)];
    offset = size(X,1);
    for cellid = 1:10
        spike_times = Blocks{blockid}.spktimes{cellid};
        
        spikes_binned = histc(spike_times,time_bin_edges);
        spikes_binned(end) = [];
        spikes_binned = spikes_binned(used_times);
        
        un_spk_cnts = unique(spikes_binned);
        cur_spikebins = [];
        for i = 1:length(un_spk_cnts)
            cur_set = find(spikes_binned == un_spk_cnts(i));
            cur_spikebins = [cur_spikebins; repmat(cur_set(:),un_spk_cnts(i),1)];
        end
        cur_spikebins = sort(cur_spikebins);
        spikebins{cellid} = [spikebins{cellid}; [(cur_spikebins+offset) ones(length(cur_spikebins),1)*blockid]];
        
        %
    end
    
    % create time-embedded stimulus matrix
    S = makeStimRows(ov_im_patch(used_times,:,:)/norm_fac,nlags,0);
    X = [X; S];

end

%%
cellid = 7;
sta = mean(X(spikebins{cellid}(:,1),:)) - mean(X);
plotfilterbank(sta',33,1:33^2)

%%
utvcv = cov(X);
stvcv = cov(X(spikebins{cellid}(:,1),:));
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
npos=3; nneg=3;
stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
nposdims = 3; nnegdims =3;
posdims = 1:nposdims; negdims = 1:nnegdims;
STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
plotfilterbank(STCbvs,33,1:33^2)

%%
SDIM = 33;
kern_l = SDIM^2*nlags;
nmods = 1;
init_kerns = sta(:);
init_signs = 1;
init_betas = 2;
init_kerns = randn(kern_l,nmods);
% init_signs = ones(nmods,1);


cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat
defmod.fsdim = SDIM^2;
defmod.pids = 1:defmod.fsdim;
defmod.h = 1;
defmod.SDIM = SDIM;
defmod.locLambda = 0;
defmod.lambda_dX = 100;
defmod.lambda_dT = 0;
defmod.lambda_L1x = 0;
basis = 'pix';

% glm0 = createGLM_quad(init_kerns,init_signs,defmod,basis,[],[]);
glm0 = createGLM_lexp(init_kerns,init_signs,init_betas,defmod,basis,[],[]);
glm0.image_type = '2d';
glm0 = normalizeRFs_full(glm0,X(ismember(full_blockid,[1 2 3]),:));
ml_sta = fitGLM_lexp(glm0,X(ismember(full_blockid,[1 2 3]),:),spikebins{cellid}(ismember(spikebins{cellid}(:,2),[1 2 3]),1),'tots')

% glm0 = normalizeRFs_full(glm0,X);
% ml_sta = fitGLM_lexp(glm0,X,spikebins{cellid}(:,1),'tots');
% 
% ml_sta2 = fitGLM_lexp(ml_sta2,X,spikebins{cellid}(:,1),'tots');

%%
cur_spikebins = spikebins{cellid}(ismember(spikebins{cellid}(:,2),[1 2 3]));
n_used_samps = sum(ismember(full_blockid,[1 2 3]));
avg_rate = length(cur_spikebins)/n_used_samps;

n_xv_samps = sum(ismember(full_blockid,4));
xv_spikebins = spikebins{cellid}(ismember(spikebins{cellid}(:,2),4)) - n_used_samps;
null_pred = ones(n_xv_samps,1)*avg_rate;
null_xvLL = -(sum(log(null_pred(xv_spikebins))) - sum(null_pred))/length(xv_spikebins);

mod_xvLL = getLLGLM_lexp(ml_sta,X(ismember(full_blockid,[4]),:),xv_spikebins,'tots')

