
clear all close all
addpath(genpath('~/Data/bruce/2_27_12'))
cd ~/Data/bruce/2_27_12/stimrecon/
addpath('~/James_scripts/bruce/');

%%
load Blocks.mat
% accept_window = [-5 5;-5 5];
accept_window = [-5 5;-4 4];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;

% original image resolution
% Pix2Deg = 1.1279 / 60;
Pix2Deg = 0.018837;

% down-sampling fraction for image
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
Nyd = Nyp/dsfrac;
Nxd = Nxp/dsfrac;

% desired temporal resolution
stimres = 0.025; %in s

Fsd = 1/Pix2Deg/dsfrac; %new sampling frequency (pix/deg)

edge_sig = 0.25*Fsd;
edgeGauss = fspecial('gaussian',round(10*edge_sig),edge_sig);

% RF_patch = [3 6.5; -4.5 -0];
% RF_patch = [-6.5 6.5; -5 5];
RF_patch = [0 6.5; -5 1];
used_x = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
used_y = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

gaussfilt = fspecial('gaussian',5,0.25);


kern_len = 50;
[X,Y] = meshgrid(-kern_len/2:kern_len/2,-kern_len/2:kern_len/2);
Pix2Deg = 0.018837;
% down-sampling fraction for image
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;

orientations = linspace(0,pi,9);
orientations(end) = [];
lambdas = 1./[0.5 1 1.5 2 2.5 3]*Fsd;
phases = [0 pi/2];
gabor_bank = zeros(length(orientations),length(lambdas),length(phases),size(X,1),size(X,2));
for i = 1:length(orientations)
    for j = 1:length(lambdas)
        for k = 1:length(phases)
        gabor_bank(i,j,k,:,:) = get_gabor_template(X,Y,0,0,orientations(i),lambdas(j),phases(k));
        end
    end
end


%%
% blockid = 1;
for blockid = 1:4
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
    
    %         if we haven't already computed saccade times do so now
    sname = sprintf('sac_times_corf_block%d',blockid);
    %     if ~exist([sname '.mat'],'file')
    %         fprintf('Computing saccade times\n');
    %         avg_eyepos = (reye_pos + leye_pos)/2;
    %         clear sm_avg_eyepos eye_vel
    %         sm_avg_eyepos(:,1) = smooth(avg_eyepos(:,1),3);
    %         sm_avg_eyepos(:,2) = smooth(avg_eyepos(:,2),3);
    %         eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/Eyedt;
    %         eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/Eyedt;
    %         eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
    %
    %         vergence = reye_pos - leye_pos;
    %         vervel = [0 0;diff(vergence)]/Eyedt;
    %         %         verspeed = abs(vervel);
    %         verspeed = sqrt(vervel(:,1).^2+vervel(:,2).^2);
    %         %         verspeed = max(vervel,[],2); %take the larger of the hor and ver vergence speeds (Read and Cummin 2003).
    %         blink_on = find(verspeed(2:end) > blink_thresh & verspeed(1:end-1) <= blink_thresh);
    %         blink_off = find(verspeed(2:end) <= blink_thresh & verspeed(1:end-1) > blink_thresh);
    %         blink_dur = eyets(blink_off) - eyets(blink_on);
    %         is_blink = find(blink_dur > min_blink_dur);
    %         blink_times = eyets(round((blink_on(is_blink)+blink_off(is_blink))/2));
    %
    %         %find saccade start and stop indices
    %         sac_inds = find(eye_speed(1:end-1) < sac_eyespeed & eye_speed(2:end) > sac_eyespeed);
    %
    %         sac_start_inds = nan(size(sac_inds));
    %         sac_stop_inds = nan(size(sac_inds));
    %         for i = 1:length(sac_inds)
    %             temp = find(eye_speed(1:sac_inds(i)) < thresh_eyespeed,1,'last');
    %             if ~isempty(temp)
    %                 sac_start_inds(i) = temp;
    %             end
    %             temp = find(eye_speed(sac_inds(i)+1:end) < thresh_eyespeed,1,'first');
    %             if ~isempty(temp)
    %                 sac_stop_inds(i) = sac_inds(i)+temp-1;
    %             end
    %         end
    %
    %         %identify start and stop times of unique saccades
    %         sac_vec = zeros(size(reye_pos,1),1);
    %         for i = 1:length(sac_start_inds)
    %             if ~isnan(sac_start_inds(i)) & ~isnan(sac_stop_inds(i))
    %                 sac_vec(sac_start_inds(i):sac_stop_inds(i)+sac_buffer_inds) = 1;
    %             end
    %         end
    %         sac_vec(length(eyets)+1:end) = [];
    %         sac_vec([1 end]) = 0;
    %         sac_start_indsn = find(sac_vec(1:end-1) == 0 & sac_vec(2:end) == 1);
    %         sac_stop_indsn = find(sac_vec(1:end-1) == 1 & sac_vec(2:end) == 0);
    %         if length(sac_start_indsn) ~= length(sac_stop_indsn)
    %             error('saccade mis-alignment');
    %         end
    %
    %         sac_start_times = eyets(sac_start_indsn);
    %         sac_stop_times = eyets(sac_stop_indsn);
    %
    %         %compute saccade amplitudes, velocities, and durations
    %         sac_dx = avg_eyepos(sac_stop_indsn,1) - avg_eyepos(sac_start_indsn,1);
    %         sac_dy = avg_eyepos(sac_stop_indsn,2) - avg_eyepos(sac_start_indsn,2);
    %         sac_amps = sqrt(sac_dx.^2+sac_dy.^2);
    %         sac_peakvel = zeros(size(sac_start_indsn));
    %         for i = 1:length(sac_start_indsn)
    %             sac_peakvel(i) = max(eye_speed(sac_start_indsn(i):sac_stop_indsn(i)));
    %         end
    %         sac_durs = sac_stop_times-sac_start_times;
    %
    %         cd ~/Data/bruce/2_27_12/saccades/
    %         save(sname,'sac_start_times','sac_stop_times','sac_start_inds','sac_stop_inds','sac_vec',...
    %             'sac_amps','sac_peakvel','sac_durs','blink_times')
    %     else
    load(sname)
    %     end
    
    recon_t = []; %time axis for reconstructed stim matrix (sec)
    stim_num = []; %corresponding vector of image indices
    ov_eye_pos = [];
    spike_bin_vecs = [];
    ov_used_inds = [];
    
    spk_trg_avg = zeros(10,length(orientations),length(lambdas),length(used_y),length(used_x));
    ov_avg = zeros(length(orientations),length(lambdas),length(used_y),length(used_x));
    

    filt_stims = mod(0:500,4) + 1 >2;
    Nimage = length(stim_times);
     avg_imset = zeros(Nimage,length(orientations),length(lambdas),length(used_y),length(used_x));;
    var_imset = zeros(Nimage,length(orientations),length(lambdas));
   for i = 1:Nimage
        onsetT = stim_times(i); %stimulus onset time (s)
        if i < Nimage %for the last image, take the end time as the block end time
            endT = stim_times(i+1)-stimres; %set end time as one time bin before the following stimulus presentation
            endT = endT - 0.2; %last 200ms doesn't have any image presentation!
        else
            endT = Blocks{blockid}.blocktimes(2,end);
        end
        
        % relevant eye signal
        eyeindx = round((onsetT-EyeStartT)/Eyedt):round((endT-EyeStartT)/Eyedt); %eye-tracking inds during current stimulus
        
        %make sure we aren't outside the range of the eye-tracking data
        eyeindx(eyeindx > length(eyets) | eyeindx > size(leye_pos,1)) = [];
        
        %time axis for reconstructing current image stimulus
        Tticks = onsetT: stimres: endT;
        Tticks(Tticks > eyets(eyeindx(end))) = [];
        
        %interpolate eye-tracking signal onto current time axis (linear interp)
        lEyelow = interp1(eyets(eyeindx), leye_pos(eyeindx,:), Tticks);
        %             lEyelow = interp1(eyets(eyeindx), reye_pos(eyeindx,:), Tticks);
        
        %allow for extrapolation of eye signal to first or last sample if
        %interpolation failed
        if any(isnan(lEyelow(end,:)))
            lEyelow(end,:) = lEyelow(end-1,:);
        end
        if any(isnan(lEyelow(1,:)))
            lEyelow(1,:) = lEyelow(2,:);
        end
        
        %compute binned spike count vectors at stimulus resolution
        time_bin_edges = [(Tticks(1)-stimres/2) (Tticks+stimres/2)];
        im_spike_bin_vecs = zeros(10,length(Tticks));
        for cellid = 1:10
            spikes_binned = histc(Blocks{blockid}.spktimes{cellid},time_bin_edges);
            spikes_binned(end) = [];
            im_spike_bin_vecs(cellid,:) = spikes_binned;
        end
        
        %create an interpolated 0-1 vector of saccade times
        sac_vec_interp = zeros(size(Tticks));
        cur_sacs = find(sac_start_times >= Tticks(1) & sac_stop_times <= Tticks(end));
        for ii = 1:length(cur_sacs)
            cur_set = find(Tticks >= sac_start_times(cur_sacs(ii)) & Tticks <= sac_stop_times(cur_sacs(ii)));
            sac_vec_interp(cur_set) = 1;
        end
        sac_vec_interp([1 end]) = 0;
        
        % create vector determining whether current stimulus was high-pass filtered
        is_stim_filtered(i) = filt_stims(stimID(i))==1;
        
        % find times where eye-position is within central window
        in_window = (lEyelow(:,1) >= accept_window(1,1) & lEyelow(:,1) <= accept_window(1,2) & ...
            lEyelow(:,2) >= accept_window(2,1) & lEyelow(:,2) <= accept_window(2,2));
        
        used_inds = find(in_window==1 & sac_vec_interp'==0);
        
        eye_pos_grid = round(lEyelow*Fsd);      
        
        fprintf('Image %d of %d\n',i,Nimage);
        if stimID(i)<10
            filename = ['ExptA0000' num2str(stimID(i)) '.png'];
        elseif stimID(i)>=10 && stimID(i)<100
            filename = ['ExptA000' num2str(stimID(i)) '.png'];
        end
        
        IMAGEorg = imread(filename);
        IMAGEorg = double(IMAGEorg); % convert to double format
        IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
        IMAGEorg = flipud(IMAGEorg);
        
        [Ny Nx] = size(IMAGEorg);
        
        IMAGE = imfilter(IMAGEorg,gaussfilt,'replicate');%slight smoothing filter (anti-aliasing)
        IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
        
        IMAGE = edgetaper(IMAGE,edgeGauss); %taper image edges
        
        im_filt_bank = zeros(length(orientations),length(lambdas),length(phases),Nyd,Nxd);
        for ii = 1:length(orientations)
            for jj = 1:length(lambdas)
                for kk = 1:length(phases)
                    im_filt_bank(ii,jj,kk,:,:) = conv2(IMAGE,squeeze(gabor_bank(ii,jj,kk,:,:)),'same');
                end
            end
        end
        im_filt_bank = sqrt(squeeze(im_filt_bank(:,:,1,:,:).^2 + im_filt_bank(:,:,2,:,:).^2));       

        cur_imset = zeros(length(used_inds),length(orientations),length(lambdas),length(used_y),length(used_x));
        for ii = 1:length(used_inds)
            cur_imset(ii,:,:,:,:) = im_filt_bank(:,:,used_y + eye_pos_grid(used_inds(ii),2),...
                used_x + eye_pos_grid(used_inds(ii),1));
            ov_avg = ov_avg + squeeze(cur_imset(ii,:,:,:,:));
%             spk_trg_avg = spk_trg_avg + bsxfun(@times,reshape(cur_set(ii,:,:,:,:),[1 length(orientations) length(lambdas) ...
%                 length(used_y) length(used_x)]),im_spike_bin_vecs(:,used_inds(ii)));
            spk_trg_avg = spk_trg_avg + bsxfun(@times,cur_imset(ii,:,:,:,:),im_spike_bin_vecs(:,used_inds(ii)));
        end
        avg_imset(i,:,:,:,:) = mean(cur_imset);
        cur_imset = shiftdim(cur_imset,3);
        oc = size(cur_imset);
        cur_imset = reshape(cur_imset,[prod(oc(1:3)) oc(4:5)]);
        var_imset(i,:,:) = var(cur_imset);
        
        ov_used_inds = [ov_used_inds; used_inds];
        ov_eye_pos = [ov_eye_pos; lEyelow];
        recon_t = [recon_t; Tticks'];
        stim_num = [stim_num; i*ones(length(Tticks),1)];
        spike_bin_vecs = [spike_bin_vecs; im_spike_bin_vecs'];
    end
    
    avg_within_im_var(blockid,:,:) = mean(var_imset);
    bet_im_var = squeeze(var(avg_imset));
    bet_im_var = shiftdim(bet_im_var,2);
    sb = size(bet_im_var);
    bet_im_var = reshape(bet_im_var,[sb(1)*sb(2) sb(3) sb(4)]);
    avg_between_im_var(blockid,:,:) = mean(bet_im_var);
    
    ov_avg = ov_avg/length(ov_used_inds);
    spk_trg_avg = bsxfun(@rdivide,spk_trg_avg,sum(spike_bin_vecs(ov_used_inds,:))');
    so = size(ov_avg);
     
    block_ov_avg(blockid,:,:,:,:) = ov_avg;
    block_spk_trg_avg(blockid,:,:,:,:,:) = spk_trg_avg;
    spk_trg_avg_sub(blockid,:,:,:,:,:) = bsxfun(@minus,spk_trg_avg,...
        reshape(ov_avg,[1 so]));
    
end

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
save gabor_models_v2 block_* spk_* avg_* orientations lambdas used_* phases RF_patch

%%
tot_var = squeeze(mean(avg_between_im_var) + mean(avg_within_im_var));
tot_std = sqrt(tot_var);
st = size(tot_var);

spk_trg_avg_norm_sub = bsxfun(@rdivide,spk_trg_avg_sub,reshape(tot_std,[1 1 st 1 1]));
avg_spktrg_norm_sub = squeeze(mean(spk_trg_avg_norm_sub));

%%
close all
cellid = 8;
cur_set = squeeze(avg_spktrg_norm_sub(cellid,:,:,:,:));
minval = min(cur_set(:));
maxval = max(cur_set(:));

for i = 1:length(orientations)
    for j = 1:length(lambdas)
        subplot(length(orientations),length(lambdas),(i-1)*length(lambdas)+j)
        imagesc(xax(used_x),yax(used_y),squeeze(cur_set(i,j,:,:)));
        set(gca,'ydir','normal');
        caxis([minval maxval]);
    end
end




