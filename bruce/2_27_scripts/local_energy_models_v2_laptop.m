
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
Nyd = Nyp/dsfrac;
Nxd = Nxp/dsfrac;

% desired temporal resolution
stimres = 0.025; %in s


Fsd = 1/Pix2Deg/dsfrac; %new sampling frequency (pix/deg)
sigx = 0.5*Fsd; %sigma for Gaussian
sigy = 0.5*Fsd;

win_x = round(-3*sigx):round(3*sigx);
win_y = round(-3*sigy):round(3*sigy);
[X,Y] = meshgrid(win_x,win_y);
gaussfun = exp(-(X.^2/(2*sigx^2)+Y.^2/(2*sigy^2)));
gaussfun = gaussfun/sum(gaussfun(:));

edge_sig = 0.25*Fsd;
edgeGauss = fspecial('gaussian',round(10*edge_sig),edge_sig);

eps = 1e-10;
bandpass_freqs = [0.5 1;0.5 2; 1 2];
[freqx,freqy] = meshgrid(linspace(-Fsd/2,Fsd/2,length(win_x)),linspace(-Fsd/2,Fsd/2,length(win_x)));
freqnorm = sqrt(freqx.^2+freqy.^2);
freq = linspace(-Fsd/2,Fsd/2,length(win_x));
max_freq = 3.5;
used_freqs = find(abs(freqx) < max_freq & abs(freqy) < max_freq & freqy >= freqx);
freq_u = find(abs(freq) < max_freq);

grid_dx = sigx/6;
grid_dy = sigy/6;
gridX = round(grid_dx/2:grid_dx:(Nxd-grid_dx/2));
gridY = round(grid_dy/2:grid_dy:(Nyd-grid_dy/2));
gridX_p = (gridX-max(gridX)/2-grid_dx/2)/Fsd;
gridY_p = (gridY-max(gridY)/2-grid_dy/2)/Fsd;

RF_patch = [3 6.5; -4.5 -0];
used_gridx = find(gridX_p >= RF_patch(1,1) & gridX_p <= RF_patch(1,2));
used_gridy = find(gridY_p >= RF_patch(2,1) & gridY_p <= RF_patch(2,2));


gaussfilt = fspecial('gaussian',5,0.25);

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
    
    spk_trg_avg = zeros(10,length(used_gridy),length(used_gridx),length(used_freqs));
    ov_avg = zeros(length(used_gridy),length(used_gridx),length(used_freqs));
    
    filt_stims = mod(0:500,4) + 1 >2;
    Nimage = length(stim_times);
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
        
        eye_pos_grid = round(lEyelow/(grid_dx/Fsd));
        
        
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
        
        nzpad = 100;
        IMAGEp = cat(1,zeros(nzpad,size(IMAGE,2)),IMAGE);
        IMAGEp = cat(1,IMAGEp,zeros(nzpad,size(IMAGEp,2)));
        IMAGEp = cat(2,zeros(size(IMAGEp,1),nzpad),IMAGEp);
        IMAGEp = cat(2,IMAGEp,zeros(size(IMAGEp,1),nzpad));
        
        cd /Users/James/Data/bruce/2_27_12/stimrecon
        fname = [sprintf('Block%d_Image%d_FFT',blockid,i) '.mat'];
        %         if ~exist(fname,'file')
        disp('Computing Image FFT');
        cur_im_decomp = zeros(length(gridY),length(gridX),length(used_freqs));
        for xx = 1:length(gridX)
            for yy = 1:length(gridY)
                IMchunk = IMAGEp(win_y+gridY(yy)+nzpad,win_x+gridX(xx)+nzpad);
                IMchunk = IMchunk.*gaussfun;
                IMchunk = IMchunk - mean(IMchunk(:));
                F = fftshift(fft2(IMchunk));
                cur_im_decomp(yy,xx,:) = abs(F(used_freqs));
            end
        end
        save(fname,'cur_im_decomp');
        %         else
        %             load(fname);
        %         end
        for ii = 1:length(used_inds)
            cur_set = cur_im_decomp(used_gridy+eye_pos_grid(used_inds(ii),2),...
                used_gridx+eye_pos_grid(used_inds(ii),1),:);
            ov_avg = ov_avg + cur_set;
            spk_trg_avg = spk_trg_avg + bsxfun(@times,reshape(cur_set,[1 length(used_gridy) length(used_gridx) length(used_freqs)]),...
                im_spike_bin_vecs(:,used_inds(ii)));
        end
        
        ov_used_inds = [ov_used_inds; used_inds];
        ov_eye_pos = [ov_eye_pos; lEyelow];
        recon_t = [recon_t; Tticks'];
        stim_num = [stim_num; i*ones(length(Tticks),1)];
        spike_bin_vecs = [spike_bin_vecs; im_spike_bin_vecs'];
    end
    
    
    ov_avg = ov_avg/length(ov_used_inds);
    spk_trg_avg = bsxfun(@rdivide,spk_trg_avg,sum(spike_bin_vecs(ov_used_inds,:))');
    o = size(ov_avg);
    
    block_ov_avg(blockid,:,:,:) = ov_avg;
    block_spk_trg_avg(blockid,:,:,:,:) = spk_trg_avg;
    spk_trg_avg_sub(blockid,:,:,:,:) = bsxfun(@minus,spk_trg_avg,...
        reshape(ov_avg,[1 o]));
    
end

%%
zmat = zeros(length(freqx),length(freqy));
[indj,indi] = ind2sub(size(zmat),used_freqs);
xrange = min(indi):max(indi);
yrange = min(indj):max(indj);

cur_set = squeeze(spk_trg_avg_sub(1,4,:,:,:));
ub = mean(cur_set(:))+5*std(cur_set(:));
lb = mean(cur_set(:))-5*std(cur_set(:));
zerofreq = find(freqx(used_freqs)==0 & freqy(used_freqs)==0);
cur_set(:,:,zerofreq) = nan;

%%
block_ov_avg_norm = bsxfun(@rdivide,block_ov_avg,...
    reshape(mean(block_avg_pow),[1 1 1 length(used_freqs)]));
block_spktrg_avg_norm = bsxfun(@rdivide,block_spk_trg_avg,...
    reshape(mean(block_avg_pow),[1 1 1 1 length(used_freqs)]));
block_spktrg_avg_norm_sub = bsxfun(@minus,block_spktrg_avg_norm,...
    reshape(block_ov_avg_norm,[4 1 length(used_gridy) length(used_gridx) length(used_freqs)]));

%%
% blockid = 1;
zerofreq = find(abs(freqx(used_freqs))<0.1 & abs(freqy(used_freqs))<0.1);
temp = reshape(mean(block_ov_avg),size(block_ov_avg,2)*size(block_ov_avg,3),size(block_ov_avg,4));
ov_means = mean(temp);
ov_stds = std(temp);
% spk_trg_avg_sub_norm = bsxfun(@rdivide,spk_trg_avg_sub,reshape(ov_stds,[1,1,1,1,length(ov_stds)]));
% spk_trg_avg_sub_norm = spk_trg_avg_sub;
block_spktrg_avg_norm_sub(:,:,:,:,zerofreq) = nan;
% freqnorm = sqrt(freqx.^2+freqy.^2);
% low_freqs = find(freqnorm(used_freqs) < 0.8);
% block_spktrg_avg_norm_sub(:,:,:,:,low_freqs) = nan;

[tempx,tempy] = meshgrid(freq(xrange),freq(yrange));
% unused_freqs = find(tempx.^2+tempy.^2 < 1 | tempx.^2+tempy.^2 > 3.5);
unused_freqs = find(tempx.^2+tempy.^2 < 1);
% unused_freqs = [];

used_blocks = [1:4];
clear avg_spec_array
for cellid = 1:10
    % cellid = 7;
    cur_set = squeeze(nanmean(block_spktrg_avg_norm_sub(used_blocks,cellid,:,:,:),1));
    avg_spec_array(cellid,:,:,:,:) = zeros(size(cur_set,1),size(cur_set,2),length(yrange),length(xrange));
    cur_set_raw = squeeze(nanmean(spk_trg_avg_sub(used_blocks,cellid,:,:,:),1));
    avg_spec_array_raw(cellid,:,:,:,:) = zeros(size(cur_set,1),size(cur_set,2),length(yrange),length(xrange));
    temp = zmat;
    for ii = 1:size(cur_set,1)
        for jj = 1:size(cur_set,2)
            temp(used_freqs) = cur_set(ii,jj,:);
            temp = flipud(fliplr(temp));
            temp(used_freqs) = cur_set(ii,jj,:);
            avg_spec_array(cellid,ii,jj,:,:) = temp(yrange,xrange);
            
            temp(used_freqs) = cur_set_raw(ii,jj,:);
            temp = flipud(fliplr(temp));
            temp(used_freqs) = cur_set_raw(ii,jj,:);
            avg_spec_array_raw(cellid,ii,jj,:,:) = temp(yrange,xrange);
        end
    end
    
    temp = squeeze(avg_spec_array(cellid,:,:,:,:));
    temp(:,:,unused_freqs) = -Inf;
    [a,b] = nanmax(temp(:));
    [ypeakloc(cellid),xpeakloc(cellid),yfreqloc(cellid),xfreqloc(cellid)] = ind2sub(size(temp),b);
end

peak_xfreqs = freq(xrange(xfreqloc));
peak_yfreqs = freq(yrange(yfreqloc));

%%
clear ind_spec_array
for blockid = 1:4
    for cellid = 1:10
        % cellid = 7;
        cur_set = squeeze(block_spktrg_avg_norm_sub(blockid,cellid,:,:,:));
        %         cur_set(:,:,unused_freqs) = -Inf;
        ind_spec_array(blockid,cellid,:,:,:,:) = zeros(size(cur_set,1),size(cur_set,2),length(yrange),length(xrange));
        temp = zmat;
        for ii = 1:size(cur_set,1)
            for jj = 1:size(cur_set,2)
                temp(used_freqs) = cur_set(ii,jj,:);
                temp = flipud(fliplr(temp));
                temp(used_freqs) = cur_set(ii,jj,:);
                ind_spec_array(blockid,cellid,ii,jj,:,:) = temp(yrange,xrange);
            end
        end
        
        temp = squeeze(ind_spec_array(blockid,cellid,:,:,:,:));
        temp(:,:,unused_freqs) = -Inf;
        [a,b] = nanmax(temp(:));
        [ind_ypeakloc(blockid,cellid),ind_xpeakloc(blockid,cellid),ind_yfreqloc(blockid,cellid),ind_xfreqloc(blockid,cellid)] = ind2sub(size(temp),b);
    end
end

%%
close all
cellid = 3;
figure
imagesc(gridX_p(used_gridx),gridY_p(used_gridy),squeeze(avg_spec_array(cellid,:,:,yfreqloc(cellid),xfreqloc(cellid))));
set(gca,'ydir','normal')
xlabel('Horizontal position (degrees)','fontsize',14)
ylabel('Vertical position (degrees)','fontsize',14)
% caxis([0 3])

figure
subplot(2,1,1)
imagesc(freqx(used_freqs),freqy(used_freqs),squeeze(avg_spec_array(cellid,ypeakloc(cellid),xpeakloc(cellid),:,:)));
set(gca,'ydir','normal')
xlabel('Horizontal frequency (cyc/deg)','fontsize',14)
ylabel('Vertical frequency (cyc/deg)','fontsize',14)
subplot(2,1,2)
imagesc(freqx(used_freqs),freqy(used_freqs),squeeze(avg_spec_array_raw(cellid,ypeakloc(cellid),xpeakloc(cellid),:,:)));
set(gca,'ydir','normal')
xlabel('Horizontal frequency (cyc/deg)','fontsize',14)
ylabel('Vertical frequency (cyc/deg)','fontsize',14)

figure
for blockid = 1:4
    subplot(2,2,blockid)
    imagesc(gridX_p(used_gridx),gridY_p(used_gridy),squeeze(ind_spec_array(blockid,cellid,:,:,yfreqloc(cellid),xfreqloc(cellid))));
    set(gca,'ydir','normal');
    xlabel('Horizontal position (degrees)','fontsize',14)
    ylabel('Vertical position (degrees)','fontsize',14)
end

%%
close all
for cellid = 1:10;
    fname = sprintf('Avg_energy_map_cell%d_lefteye',cellid);
    figure
    imagesc(gridX_p(used_gridx),gridY_p(used_gridy),squeeze(avg_spec_array(cellid,:,:,yfreqloc(cellid),xfreqloc(cellid))));
    set(gca,'ydir','normal')
    xlabel('Horizontal position (degrees)','fontsize',14)
    ylabel('Vertical position (degrees)','fontsize',14)
    axis square
%     print('-dpdf',fname);
    print('-dpng',fname);
    close
    
    fname = sprintf('Peakpos_powerspectra_cell%d_lefteye',cellid);
    figure
    subplot(2,1,1)
    imagesc(freqx(used_freqs),freqy(used_freqs),squeeze(avg_spec_array(cellid,ypeakloc(cellid),xpeakloc(cellid),:,:)));
    set(gca,'ydir','normal')
    xlabel('Horizontal frequency (cyc/deg)','fontsize',14)
    ylabel('Vertical frequency (cyc/deg)','fontsize',14)
         axis square
   subplot(2,1,2)
    imagesc(freqx(used_freqs),freqy(used_freqs),squeeze(avg_spec_array_raw(cellid,ypeakloc(cellid),xpeakloc(cellid),:,:)));
    set(gca,'ydir','normal')
    xlabel('Horizontal frequency (cyc/deg)','fontsize',14)
    ylabel('Vertical frequency (cyc/deg)','fontsize',14)
        axis square
%     print('-dpdf',fname);
    print('-dpng',fname);
    close

    fname = sprintf('Block_Energydist_cell%d_lefteye',cellid);
    figure
    for blockid = 1:4
        subplot(2,2,blockid)
        imagesc(gridX_p(used_gridx),gridY_p(used_gridy),squeeze(ind_spec_array(blockid,cellid,:,:,yfreqloc(cellid),xfreqloc(cellid))));
        set(gca,'ydir','normal');
     axis square
     title(sprintf('Block %d',blockid));
       xlabel('Horizontal position (degrees)','fontsize',14)
        ylabel('Vertical position (degrees)','fontsize',14)
    end
%     print('-dpdf',fname);
    print('-dpng',fname);
    close

end

%%
hor_freqrange = find(abs(freqx(used_freqs)) >= 1 & abs(freqx(used_freqs)) <= 2 & freqy(used_freqs) == 0);
ver_freqrange = find(abs(freqy(used_freqs)) >= 1 & abs(freqy(used_freqs)) <= 2 & freqx(used_freqs) == 0);
clear hor_freqmap ver_freqmap
for cellid = 1:10
    %     cur_set = squeeze(mean(spk_trg_avg_sub_norm(used_blocks,cellid,:,:,:)));
    cur_set = squeeze(spk_trg_avg_sub_norm(:,cellid,:,:,:));
    hor_freqmap(:,cellid,:,:) = mean(cur_set(:,:,:,hor_freqrange),4);
    ver_freqmap(:,cellid,:,:) = mean(cur_set(:,:,:,ver_freqrange),4);
    
end



