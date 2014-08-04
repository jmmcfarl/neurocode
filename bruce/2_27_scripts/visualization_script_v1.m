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


Pix2Deg = 0.018837;

% down-sampling fraction for image
dsfrac = 2;
Fsd_sp = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd_sp; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd_sp;

RF_patch = [2.75 6.25; -3.5 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd_sp*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

fov_patch = [-1 1; -1 1]; %location of RFs in degrees [x1 x2;y1 y2]
fov_xpatch_inds = find(xax >= fov_patch(1,1) & xax <= fov_patch(1,2));
fov_ypatch_inds = find(yax >= fov_patch(2,1) & yax <= fov_patch(2,2));
fov_patch_width = diff(fov_patch,[],2);

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

cd /Users/James/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_d2
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

X_mean = cur_X_seq(end,:)/Fsd_sp;
X_std = post_stdX(end,:)/Fsd_sp;
Y_mean = cur_Y_seq(end,:)/Fsd_sp;
Y_std = post_stdY(end,:)/Fsd_sp;
est_pos_times = 0.5*all_fix_start_times(used_fixs)+0.5*all_fix_stop_times(used_fixs);
est_pos_blockids = blockids(used_fixs);
pos_error = [X_mean(:) Y_mean(:)];

SDIM = length(xpatch_inds);


su_rf_y0 = RF_patch_cent(2) + gabor_params_fin(1:9,2)/Fsd_sp;
su_rf_x0 = RF_patch_cent(1) + gabor_params_fin(1:9,1)/Fsd_sp;
su_sigma = gabor_params_fin(1:9,5)/Fsd_sp;
su_ecc = gabor_params_fin(1:9,6);
ra = su_sigma*4;
rb = ra./su_ecc;
ang = gabor_params_fin(1:9,3);

used_cells = [2 3 4 5 6 7 8];

lfp_backwin = 0.15;
lfp_forwin = 0.45;


used_lfp_ch = [4 16];
lfp_colors = [1 0 0; 0 0 1];
use_su_phase = 6;

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_all
    all_fix_durs = all_fix_stop_times - all_fix_start_times;

%% CYCLE THROUGH BLOCKS
for blockid = 2;
    
    fprintf('Processing block %d of %d...\n',blockid,4);
    
    block_times = Blocks{blockid}.blocktimes;
    stim_times = Blocks{blockid}.stimtime;
    stimID = Blocks{blockid}.stimids;
    
    n_stims = length(stim_times);
    
    %% LOAD AND PROCESS EYE-TRACKING DATA
    cd ~/Data/bruce/2_27_12/saccades/
    load(sprintf('lemM232.5%d.em.sac.mat',blockid))
    %     load(sprintf('lemM232.5%d.em.hor.mat',blockid))
    
    % identify saccade start and stop times
    EyeStartT = Expt.Trials.Start/10000; % time of first eye sample
    EyeEndT = Expt.Trials.End/10000; % time of last eye sample
    Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
    eyets = EyeStartT:Eyedt:(EyeEndT-Eyedt); %eye tracking time axis (sec)
    sac_buffer_inds = round(sac_buffer/Eyedt);
    
    reye_pos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
    leye_pos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];
    
    avg_eyepos = (reye_pos + leye_pos)/2;
    clear sm_avg_eyepos eye_vel
    sm_avg_eyepos(:,1) = smooth(avg_eyepos(:,1),3);
    sm_avg_eyepos(:,2) = smooth(avg_eyepos(:,2),3);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/Eyedt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/Eyedt;
    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
    vel_thresh = 5;
    eye_vel(eye_vel(:,1) > vel_thresh,1) = vel_thresh;
    eye_vel(eye_vel(:,1) < -vel_thresh,1) = -vel_thresh;
    eye_vel(eye_vel(:,2) > vel_thresh,2) = vel_thresh;
    eye_vel(eye_vel(:,2) < -vel_thresh,2) = -vel_thresh;
    
    %% LOAD AND PROCESS LFP DATA
    cd /Users/James/Data/bruce/2_27_12
    load(sprintf('lemM232A.5%d.lfp.mat',blockid));
    %get start times of each LFP trial
    n_lfp_trials = length(LFP.Trials);
    lfp_trial_start = nan(n_lfp_trials,1);
    for i = 1:n_lfp_trials
        lfp_trial_start(i) = LFP.Trials(i).Start/1e4; %originally in tenths of ms
        lfp_trial_stop(i) = LFP.Trials(i).End/1e4;
        lfp_dur(i) = size(LFP.Trials(i).LFP,1)/1000;
    end
    % lfp_trial_stop = lfp_trial_start+lfp_dur';
    block_trial_times = Blocks{blockid}.blocktimes;
    
    lfp_time = [];
    lfp_samps = [];
    for i = 1:n_lfp_trials
        lfp_time = [lfp_time linspace(lfp_trial_start(i),lfp_trial_start(i)+lfp_dur(i),size(LFP.Trials(i).LFP,1))];
        lfp_samps = [lfp_samps; LFP.Trials(i).LFP(:,used_lfp_ch)];
    end
    % lfp_samps = zscore(lfp_samps);
    
    Fs = 1000;
    niqf = Fs/2;
    [b,a] = butter(2,[1 200]/niqf);
    lfp_samps = filtfilt(b,a,lfp_samps);
    [b,a] = butter(2,[40 60]/niqf);
    lfp_samps_hf = filtfilt(b,a,lfp_samps);
    
    dsf = 2;
    lfp_sampsd = downsample(lfp_samps,dsf);
    lfp_sampsd = zscore(lfp_sampsd);
    lfp_sampsd_hf = downsample(lfp_samps_hf,dsf);
    lfp_sampsd_hf = zscore(lfp_sampsd_hf);
    lfp_sampsd_hfamp = abs(hilbert(lfp_sampsd_hf));
    lfp_timed = downsample(lfp_time,dsf);
    Fsd = Fs/dsf;
    
    %% indices of fixation start/stop relative to eye-tracking data
    cur_block_fixs = find(blockids == blockid);
    all_fix_start_inds = round(interp1(eyets,1:length(eyets),all_fix_start_times(cur_block_fixs)));
    all_fix_stop_inds = round(interp1(eyets,1:length(eyets),all_fix_stop_times(cur_block_fixs)));

    %% COMPUTE GAMMA POWER WITHIN EACH LATE FIXATION WINDOW
    gamma_pow = nan(length(cur_block_fixs),length(used_lfp_ch));
    long_enough = find(all_fix_durs(cur_block_fixs) >= 0.3);
    for i = 1:length(long_enough)
       cur_set = all_fix_start_inds(i):all_fix_stop_inds(i); 
        gamma_pow(long_enough(i),:) = nanmean(lfp_sampsd_hfamp(cur_set,:));
    end
    [sorted_gamma_pow,gamma_ord] = sort(gamma_pow(long_enough));
    gamma_ord = long_enough(gamma_ord);
    
    %% BIN SPIKING DATA AT LFP RESOLUTION
    spikes_binned = zeros(24,length(lfp_timed));
    for c = 1:10
        temp = round(interp1(lfp_timed,1:length(lfp_timed),Blocks{blockid}.spktimes{c}));
        temp(isnan(temp)) = [];
        all_spike_lfp_inds{c} = temp;
        
        cur_hist = hist(Blocks{blockid}.spktimes{c},lfp_timed);
        cur_hist(end) = 0;
        spikes_binned(c,:) = smooth(cur_hist,10);
    end
    for c = 1:14
        temp = round(interp1(lfp_timed,1:length(lfp_timed),Blocks{blockid}.mutimes{c}));
        temp(isnan(temp)) = [];
        all_mua_lfp_inds{c} = temp;
        
        cur_hist = hist(Blocks{blockid}.mutimes{c},lfp_timed);
        cur_hist(end) = 0;
        spikes_binned(c+10,:) = smooth(cur_hist,10);
    end
    
    spike_rates = spikes_binned*Fsd;
    spike_rates = zscore(spike_rates')';
    net_rate = mean(spike_rates);
    
    suprobes = Blocks{1}.suprobes;
    muprobes = Blocks{1}.muprobes;
    allprobes = [suprobes muprobes];
    set1 = find(allprobes < 10);
    set2 = find(allprobes >= 10 & allprobes <= 20);
    net_rate_1 = mean(spike_rates(set1,:));
    net_rate_2 = mean(spike_rates(set2,:));
    
    %% APPLY EYE-TRACKING CORRECTION
    cur_interp_points = find(est_pos_blockids==blockid);
    cur_interp_times = est_pos_times(cur_interp_points);
    cur_interp_errors = pos_error(cur_interp_points,:);
    cur_interp_times = [eyets(1); cur_interp_times; eyets(end)];
    cur_interp_errors = [cur_interp_errors(1,:); cur_interp_errors; cur_interp_errors(end,:)];
    interp_correction = interp1(cur_interp_times,cur_interp_errors,eyets);
    corrected_eyepos = avg_eyepos + interp_correction;
    corrected_eyepos_pix = round(corrected_eyepos/Pix2Deg/dsfrac);
        
    %% INITIATE LFP/SPIKE PLOT
    close all
    figure(1)
    for ll = 1:length(used_lfp_ch)
        plot(lfp_timed,lfp_sampsd(:,ll),'color',lfp_colors(ll,:))
        hold on
        plot(lfp_timed,lfp_sampsd_hf(:,ll)-5-4*(ll-1),'color',lfp_colors(ll,:))
%         plot(lfp_timed,lfp_sampsd_hfamp(:,ll)-5,'color',lfp_colors(ll,:))
        plot(lfp_timed,lfp_sampsd_hfamp(:,ll)-5-4*(ll-1),'k')
    end
    for c = 1:length(set1)
        if set1(c) > 10
            plot(Blocks{blockid}.mutimes{set1(c)-10},-ones(size(Blocks{blockid}.mutimes{set1(c)-10}))*0.3*(c-1)-12,'r.');
        else
            plot(Blocks{blockid}.spktimes{set1(c)},-ones(size(Blocks{blockid}.spktimes{set1(c)}))*0.3*(c-1)-12,'r.');
        end
    end
    for c = 1:length(set2)
        if set2(c) > 10
            plot(Blocks{blockid}.mutimes{set2(c)-10},-ones(size(Blocks{blockid}.mutimes{set2(c)-10}))*0.3*(c+length(set1)-1)-12,'b.');
        else
            plot(Blocks{blockid}.spktimes{set2(c)},-ones(size(Blocks{blockid}.spktimes{set2(c)}))*0.3*(c+length(set1)-1)-12,'b.');
        end
    end
    plot(lfp_timed(all_spike_lfp_inds{use_su_phase}),lfp_sampsd_hf(all_spike_lfp_inds{use_su_phase},2)-9,'k.')
    
    ylim([-20 4])
    lfp_yl = ylim();
    for i = 1:length(cur_block_fixs)
        line([all_fix_start_times(cur_block_fixs(i)) all_fix_start_times(cur_block_fixs(i))],lfp_yl,'color','k')
        line([all_fix_stop_times(cur_block_fixs(i)) all_fix_stop_times(cur_block_fixs(i))],lfp_yl,'color','r')
        if all_fix_durs(cur_block_fixs(i)) > 0.15
           line([all_fix_start_times(cur_block_fixs(i)) all_fix_start_times(cur_block_fixs(i))]+0.15,lfp_yl,'color','g') 
        end
        if all_fix_durs(cur_block_fixs(i)) > 0.25
           line([all_fix_start_times(cur_block_fixs(i)) all_fix_start_times(cur_block_fixs(i))]+0.25,lfp_yl,'color','g')             
        end
    end
    plot(lfp_timed,net_rate_1*2-12,'r')
    plot(lfp_timed,net_rate_2*2-16,'b')
    set(gcf,'Position',[1100 600 900 900])
    
    
    %% CYCLE THROUGH IMAGES IN CURRENT BLOCK
    cd ~/Data/bruce/2_27_12/ExptA
    for i = 1:n_stims
        
        %% load and preprocess image
        cur_im_fixs = find(all_stim_nums(cur_block_fixs) == i); %find fixations with current stimulus
        fprintf('Image %d of %d\n',i,n_stims);
        if stimID(i) <10
            filename = ['ExptA0000' num2str(stimID(i)) '.png'];
        elseif stimID(i) >=10 && stimID(i) <100
            filename = ['ExptA000' num2str(stimID(i)) '.png'];
        end
        
        IMAGEorg = imread(filename);
        IMAGEorg = double(IMAGEorg); % convert to double format
        IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
        
        [Ny Nx] = size(IMAGEorg);
        
        IMAGE = IMAGEorg;
        %             IMAGE = imfilter(IMAGEorg,gaussfilt,'replicate');%slight smoothing filter (anti-aliasing)
        IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
        [Nyd Nxd] = size(IMAGE); %down-sampled image size
        
        IMAGE = flipud(IMAGE);
        
        %zero pad image
        nzpad = 300;
        IMAGEp = cat(1,nan(nzpad,size(IMAGE,2)),IMAGE);
        IMAGEp = cat(1,IMAGEp,nan(nzpad,size(IMAGEp,2)));
        IMAGEp = cat(2,nan(size(IMAGEp,1),nzpad),IMAGEp);
        IMAGEp = cat(2,IMAGEp,nan(size(IMAGEp,1),nzpad));
 
        f2 = figure(2)
        imagesc(xax,yax,IMAGE); title('Static Image','fontsize',16); set(gca,'ydir','normal')
        colormap(gray)
        set(f2,'Position',[20 600 560*1.5 420*1.5])
        title('Static Image','fontsize',16); 
            hold on
            xlabel('Horizontal Position (degrees)','fontsize',14)
            ylabel('Vertical Position (degrees)','fontsize',14)
        
            f3 = figure(3)
         set(f3,'Position',[20 40 500 400])
        f4 = figure(4)
         set(f4,'Position',[500 40 500 400])

        %% cycle through the individual fixations for this image
        for j = 1:length(cur_im_fixs)
            
            cur_eye_samps = all_fix_start_inds(cur_im_fixs(j)):all_fix_stop_inds(cur_im_fixs(j));
            cur_fix_eyepos = median(corrected_eyepos(cur_eye_samps,:));
            cur_fix_eyepos_pix = median(corrected_eyepos_pix(cur_eye_samps,:));
            ypatch_inds_adj = round(ypatch_inds + cur_fix_eyepos_pix(2) + nzpad);
            xpatch_inds_adj = round(xpatch_inds + cur_fix_eyepos_pix(1) + nzpad);
            STstim_patch = IMAGEp(ypatch_inds_adj,xpatch_inds_adj);
            fov_ypatch_inds_adj = round(fov_ypatch_inds + cur_fix_eyepos_pix(2) + nzpad);
            fov_xpatch_inds_adj = round(fov_xpatch_inds + cur_fix_eyepos_pix(1) + nzpad);
            fov_stim_patch = IMAGEp(fov_ypatch_inds_adj,fov_xpatch_inds_adj);
            
            %% plot image
            figure(2)
            imagesc(xax,yax,IMAGE); 
            colormap(gray)
            r= rectangle('Position',[fov_patch(1,1)+cur_fix_eyepos(1) fov_patch(2,1)+cur_fix_eyepos(2) fov_patch_width(1) fov_patch_width(2)],...
                'EdgeColor','r','LineWidth',4);
            r= rectangle('Position',[RF_patch(1,1)+cur_fix_eyepos(1) RF_patch(2,1)+cur_fix_eyepos(2) RF_patch_width(1) RF_patch_width(2)],...
                'EdgeColor','k','LineWidth',4);
            %             for cc = 1:length(used_cells)
            %                 ellipse(ra(used_cells(cc)),rb(used_cells(cc)),ang(used_cells(cc)),su_rf_x0(used_cells(cc))+...
            %                     cur_fix_eyepos(1),su_rf_y0(used_cells(cc))+cur_fix_eyepos(2),'k','linewidth',2)
            %             end
            ellipse(mean(ra(used_cells)),mean(rb(used_cells)),circ_mean(ang(used_cells)),mean(su_rf_x0(used_cells))+cur_fix_eyepos(1),....
                mean(su_rf_y0(used_cells))+cur_fix_eyepos(2),'g','linewidth',2)
            
            figure(3)
            imagesc(fov_stim_patch);colormap(gray);
            set(gca,'ydir','normal');
            figure(4)
            imagesc(STstim_patch);colormap(gray);
            set(gca,'ydir','normal');
         
            figure(1)
            xlim([all_fix_start_times(cur_block_fixs(cur_im_fixs(j)))-lfp_backwin all_fix_start_times(cur_block_fixs(cur_im_fixs(j)))+lfp_forwin])
            input('');
        end
            
    end
    
end