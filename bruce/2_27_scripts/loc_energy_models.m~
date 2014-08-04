
clear all close all
addpath(genpath('~/Data/bruce/2_27_12'))
cd ~/Data/bruce/2_27_12/stimrecon/
addpath('~/James_scripts/bruce/');

%%
load Blocks.mat
accept_window = [-5 5;-5 5];
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

% desired temporal resolution
stimres = 0.025; %in s

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
    
    cd ~/Data/bruce/2_27_12/stimrecon/
    load(sprintf('all_image_Energyarrays_v2_block%d',blockid));
    Nimage = size(all_Energy_array,5);
    % load images
    
    %normalize
    all_Energy_array = sqrt(all_Energy_array);
    all_Energy_array = bsxfun(@minus,all_Energy_array,mean(all_Energy_array,5));
    all_Energy_array = bsxfun(@rdivide,all_Energy_array,std(all_Energy_array,[],5));
    
    %zero pad images
    nzpad = 60;
    Esize = size(all_Energy_array);
    used_y = 1:Esize(1); used_x = 1:Esize(2);
    Energyp = cat(1,zeros([nzpad Esize(2:end)]),all_Energy_array);
    Energyp = cat(1,Energyp,zeros([nzpad Esize(2:end)]));
    Esize = size(Energyp);
    Energyp = cat(2,zeros([Esize(1) nzpad Esize(3:end)]),Energyp);
    Energyp = cat(2,Energyp,zeros([Esize(1) nzpad Esize(3:end)]));
    Esize = size(Energyp);
        
    recon_t = []; %time axis for reconstructed stim matrix (sec)
    stim_num = []; %corresponding vector of image indices
    ov_eye_pos = [];
    spike_bin_vecs = [];
    
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
        spike_bin_vecs = [spike_bin_vecs; im_spike_bin_vecs'];
        
        ov_eye_pos = [ov_eye_pos; lEyelow];
        recon_t = [recon_t; Tticks'];
        stim_num = [stim_num; i*ones(length(Tticks),1)];
    end
    ov_eye_pos_grid = round(ov_eye_pos/(grid_dx/Fsd));
    
    %create an interpolated 0-1 vector of saccade times
    sac_vec_interp = zeros(size(recon_t));
    for i = 1:length(sac_start_times)
        cur_set = find(recon_t >= sac_start_times(i) & recon_t <= sac_stop_times(i));
        sac_vec_interp(cur_set) = 1;
    end
    sac_vec_interp([1 end]) = 0;
    
    % create vector determining whether current stimulus was high-pass filtered
    filt_stims = mod(0:500,4) + 1 >2;
    is_stim_filtered = filt_stims(stimID(stim_num));
    
    % find times where eye-position is within central window
    in_window = (ov_eye_pos(:,1) >= accept_window(1,1) & ov_eye_pos(:,1) <= accept_window(1,2) & ...
        ov_eye_pos(:,2) >= accept_window(2,1) & ov_eye_pos(:,2) <= accept_window(2,2));
    
    used_inds = find(in_window==1 & sac_vec_interp==0);
    
    %%
    spk_trg_avg = zeros(10,length(used_y),length(used_x),Esize(3),Esize(4));
    ov_avg = zeros(length(used_y),length(used_x),Esize(3),Esize(4));
    for i = 1:length(used_inds)
        fprintf('%d of %d\n',i,length(used_inds));
        cur_Energy_set = squeeze(Energyp(nzpad+used_y+ov_eye_pos_grid(used_inds(i),2),...
            nzpad + used_x+ov_eye_pos_grid(used_inds(i),1),:,:,stim_num(used_inds(i))));
        ov_avg = ov_avg + cur_Energy_set;
        temp = size(cur_Energy_set);
        spk_trg_avg = spk_trg_avg + bsxfun(@times,reshape(cur_Energy_set,...
            [1 temp]),spike_bin_vecs(used_inds(i),:)');
    end
    ov_avg = ov_avg/length(used_inds);
    spk_trg_avg = bsxfun(@rdivide,spk_trg_avg,sum(spike_bin_vecs)');
    o = size(ov_avg);
    spk_trg_avg_sub(blockid,:,:,:,:,:) = bsxfun(@minus,spk_trg_avg,reshape(ov_avg,[1 o]));
    
end

%%
y_ax = (1:154)-77;
x_ax = (1:192)-96;
x_axd = x_ax*grid_dx/Fsd;
y_axd = y_ax*grid_dx/Fsd;
RF_patch = [3.5 6; -3.5 -1];
RF_patch_width = diff(RF_patch,[],2)
close all
freq_range = 3;
cellid =7;

for i = 1:4
    subplot(2,2,i)
    imagesc(y_axd,x_axd,squeeze(nanmean(spk_trg_avg_sub(:,cellid,:,:,freq_range,i))));
    set(gca,'ydir','normal');
    caxis([-0.2 0.2])
    r= rectangle('Position',[RF_patch(1,1) RF_patch(2,1) RF_patch_width(1) RF_patch_width(2)],...
        'EdgeColor','r','LineWidth',2);
end

%%
cellid = 4;
i = 1;
figure
    imagesc(y_axd,x_axd,squeeze(nanmean(spk_trg_avg_sub(:,cellid,:,:,freq_range,i))));
    set(gca,'ydir','normal');
    caxis([-0.1 0.2])
    r= rectangle('Position',[RF_patch(1,1) RF_patch(2,1) RF_patch_width(1) RF_patch_width(2)],...
        'EdgeColor','r','LineWidth',2);
xlim([0 8])
ylim([-8 5])

%%
load loc_energy_lefteye2 spk_trg_avg_sub
%%
cellid = 10;
% cellid = 4;
% cellid = 6;
usex = find(x_axd > 3 & x_axd < 6);
usey = find(y_axd > -3.5 & y_axd < 0);
freq_range = 3;
clear peakx yslice
figure
for i = 1:4
    subplot(2,2,i)
    temp = squeeze(spk_trg_avg_sub(i,cellid,:,:,freq_range,1));
    imagesc(y_axd,x_axd,temp);
    temp_slice = mean(temp(usey,usex));
    [~,peakx(i)] = max(temp_slice);
    yslice(i,:) = temp(usey,usex(peakx(i)));
    set(gca,'ydir','normal');
% caxis([0. 0.3])
xlim([2.5 6.5])
ylim([-5 1.5])
    r= rectangle('Position',[RF_patch(1,1) RF_patch(2,1) RF_patch_width(1) RF_patch_width(2)],...
        'EdgeColor','r','LineWidth',2);
end

%%


%%
freq_range = 3;
for i = 1:4
    subplot(2,2,i)
    imagesc(y_axd,x_axd,squeeze(nanmean(nanmean(spk_trg_avg_sub(:,:,:,:,freq_range,i)))));
    set(gca,'ydir','normal');
    caxis([-0.15 0.15])
    r= rectangle('Position',[RF_patch(1,1) RF_patch(2,1) RF_patch_width(1) RF_patch_width(2)],...
        'EdgeColor','r','LineWidth',2);
end

