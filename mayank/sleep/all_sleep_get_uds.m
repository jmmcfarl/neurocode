clear all
% close all
cd ~/Analysis/Mayank/sleep/
load sleep_dirs

fig_dir = '/Users/james/Analysis/Mayank/sleep/sleep_figs2/';
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))
addpath('~/James_scripts/mayank/hsmm_uds_code');
addpath('~/James_scripts/mayank/hsmm_state_detection/')
addpath('~/James_scripts/hsmm_uds_toolbox/')
addpath('~/James_scripts/mayank/sleep/');

%% DEFINE PARAMS
params.uds_range = [0.5 3]; %range of frequencies for defining UDS
params.lf_range = [0 0.2];  %range of low-freq power to characterize LF artifacts
params.use_range = [0.5 80]; %overall range of usable power (excluding LF artifacts)

params.mp_feature_frange = [0 5];  %freq range for defining MP features used to classify MP UDS [4]
params.mp_feature_frange_bb = [0 20]; %freq range for computing MP features [10]

params.n_dens_bins = 100; %number of pts to use for density estimation
% params.dens_smth_sig = 3; %sigma for visualizing density estimates (in units of bins)

params.mp_lcf = 0.05;  %lcf for hp-filter on MP signal (just gets rid of drift) [.01]
params.lfp_lcf = 0.2; %lcf for hp-filter on LFP signals (gets rid of LF artifact)
params.desired_Fs = 200;  %want a uniform Fs

params.min_seg_dur = 20; %minimum duration of UDS segments used for analysis (in sec)

params.dip_Nboot = 5000; %number of bootstrap samples for computing dip stat significance
params.min_dip_pvalue = 0.05; %alpha on dip test [0.01]

params.state_buffer = 0.15; %[0.15] buffer around state transitions (sec) for counting cortical blips as skipped

params.amp_prc_bins = [0:10:90]; %amplitude-binning (in percentile) for LFP blips

%defining criteria for UDS selection (for normalized LFP power)
params.lfp_uds_thresh = -0.6; %min log power concentrated in the UDS band [-0.5]
params.lfp_lf_max = 2; %max log power in the LF band (artifact detection)

params.compute_state_seqs = true; %compute MP state sequences?
params.hmm_dsf = 2; %down-sample-factor for HMM signal featurs
params.hmm_Fs = params.desired_Fs/params.hmm_dsf; %sample-frequency of HMM signal features (Hz)

%parameters for despiking MP data
params.spk_thresh = 20; %threshold (in robust SD units) for spikes
params.spk_back = 0.001; %window for spike rmoval
params.spk_for = 0.004;

params.shift_lfp_peaks = true; %account for delay between LFP and MP

params.min_n_events = 10; %minimum number of peaks in order to construct trig dists

%parameters for computing spectrograms
params.chron_params.tapers = [2 3];
params.chron_params.fpass = [0 80];
params.movingwin = [20 5];

%time windows for computing trig avgs (sec)
params.back_Twin = 0.75;
params.for_Twin = 0.75;

to_print = false; %for printing specgram fig

%% loop over all recs
% for dd = 14
for dd = 1:length(data)
    %     close all
    %%
    cd(data(dd).dir)
    pwd
    load('procData.mat','mp_*','ipsi_*','heka*','csc_time');
    
    all_data(dd).has_ipsi = true; %initialize
    %range of channels that might be used as our cortical LFP
    if length(data(dd).ipsiLFPs) == 8
        all_data(dd).poss_ctx_chs = 1:4;
    elseif length(data(dd).ipsiLFPs) == 16
        all_data(dd).poss_ctx_chs = 1:8;
    elseif length(data(dd).ipsiLFPs) == 7
        all_data(dd).poss_ctx_chs = 4:7;
    else
        all_data(dd).has_ipsi = false;
    end
    
    %% get sample freqs and determine additional dsf's
    csc_Fs = 1/nanmedian(diff(csc_time));
    mp_Fs = 1/nanmedian(diff(mp_t));
    
    if csc_Fs > 500
        mp_dsf = 10;
        csc_dsf = 10;
        heka_dsf = 10;
    else
        mp_dsf = 1;
        csc_dsf = 1;
        heka_dsf = 5;
    end
    
    %%
    if all_data(dd).has_ipsi && mp_Fs >= 150 %if there are any ipsilateral LFPs
        
        
        %% get MP signal (use heka/spike 2 if we have it)
        if ~isempty(heka_data)
            mp_data = heka_data(:);
            mp_time = heka_time(:);
            cur_dsf = heka_dsf;
            all_data(dd).used_heka = true;
        else
            mp_data = mp_d(:);
            mp_time = mp_t(:);
            cur_dsf = mp_dsf;
            all_data(dd).used_heka = false;
        end
        mp_Fs = 1/nanmedian(diff(mp_time));
        
        %high-pass filter and find spike times
        [spk_b,spk_a] = butter(2,100/(mp_Fs/2),'high');
        mp_high = filtfilt(spk_b,spk_a,mp_data);
        mp_high = mp_high/robust_std_dev(mp_high);
        spike_inds = find(diff(sign([0; diff(mp_high)])) < 0);
        spike_inds(mp_high(spike_inds) < params.spk_thresh) = [];
        spike_times = mp_time(spike_inds);
        clear mp_high
        
        %remove spikes from DC MP signal by interpolation
        cur_spk_inds = round(interp1(mp_time,1:length(mp_time),spike_times)); %indices of spikes
        spike_times(isnan(cur_spk_inds)) = [];
        cur_spk_inds(isnan(cur_spk_inds)) = [];
        spike_error = abs(mp_time(cur_spk_inds) - spike_times); %use only spikes that interpolate onto actual DC data
        bad = find(spike_error > 2/mp_Fs); %if interpolated time is too far from the real time get rid of these
        cur_spk_inds(bad) = [];
        blk_inds = -round(params.spk_back*mp_Fs):round(params.spk_for*mp_Fs); %set of times around each spike to interpolate out
        blocked_inds = bsxfun(@plus,cur_spk_inds,blk_inds);
        blocked_inds = blocked_inds(:);
        used_inds = setdiff(1:length(mp_data),blocked_inds); %non-spike samples
        mp_data = interp1(used_inds,mp_data(used_inds),1:length(mp_data)); %de-spiked data
                
        %% identify any missing data at the beginning and end of the recording
        mp_sp = find(abs(diff(mp_data)) > 0,1); %first pt where we actually have MP data
        mp_ep = 1+find(abs(diff(mp_data)) > 0,1,'last'); %last pt with MP data
        if data(dd).good_bounds(1) < mp_time(mp_sp)
            data(dd).good_bounds(1) = mp_time(mp_sp);
        end
        if data(dd).good_bounds(2) > mp_time(mp_ep)
            data(dd).good_bounds(2) = mp_time(mp_ep);
        end
        
        % find usable MP times
        use_mp_inds = find(mp_time >= data(dd).good_bounds(1) & mp_time <= data(dd).good_bounds(2));
        %         if ~isempty(findstr(data(dd).dir,'2015-08-16_13-06-25')) %this rec had a brief period where cell resealed
        %             bad_trange = [340 375]; %MP rec resealed here
        %             use_mp_inds(mp_time(use_mp_inds) >= bad_trange(1) & mp_time(use_mp_inds) <= bad_trange(2)) = [];
        %         end
        mp_time = mp_time(use_mp_inds); mp_data = mp_data(use_mp_inds);
     
        %%
        %now apply high-pass filtering, and down-sample MP data
        mp_Fsd = mp_Fs/cur_dsf;
        mp_data = decimate(mp_data,cur_dsf);
         if params.mp_lcf > 0
       [b,a] = butter(2,params.mp_lcf/(mp_Fsd/2),'high');
        mp_data = filtfilt(b,a,mp_data);
        end
        mp_time = downsample(mp_time,cur_dsf);

        %% process LFP data
        %find usable LFP times
        csc_time = downsample(csc_time,csc_dsf);
        csc_inds = find(csc_time >= data(dd).good_bounds(1) & csc_time <= data(dd).good_bounds(2));
        csc_time = csc_time(csc_inds);
        
        %filter and downsample LFPs
        csc_Fsd = csc_Fs/csc_dsf;
        [bb,aa] = butter(2,params.lfp_lcf/(csc_Fsd/2),'high');
        for ii = 1:length(ipsi_csc)
            ipsi_csc{ii} = decimate(ipsi_csc{ii},csc_dsf);
            ipsi_csc{ii} = ipsi_csc{ii}(csc_inds);
            ipsi_csc_hp{ii} = filtfilt(bb,aa,ipsi_csc{ii});
            ipsi_csc{ii} = ipsi_csc{ii}/robust_std_dev(ipsi_csc_hp{ii}); %normalize by robust SD est
            ipsi_csc_hp{ii} = ipsi_csc_hp{ii}/robust_std_dev(ipsi_csc_hp{ii}); %normalize by robust SD est
        end
        
        %% Get MP signal interpolated onto CSC time axis
        uu = find(~isnan(mp_time));
        mp_interp = interp1(mp_time(uu),mp_data(uu),csc_time);
        
        %handle any interpolation boundary issues
        bad_mp_inds = find(isnan(mp_interp));
        mp_interp(bad_mp_inds) = 0;
        
        if sum(isnan(mp_interp)) > 0.01*length(mp_interp)
            warning('More than 1% NAN components!');
        end
        mp_interp(isnan(mp_interp)) = 0;
        
        %% compute LFP and MP spectrograms
        
        params.chron_params.Fs = csc_Fsd;
        [b,a] = butter(2,[1/params.movingwin(1)]/(csc_Fs/2),'high'); %filter out stuff that's slower than the window size
        
        %compute spectrogram for MP signal
        [S_mp,t,f] = mtspecgramc(mp_interp(:),params.movingwin,params.chron_params); %MP spectrogram
        
        %compute specgrams for each LFP channel
        clear S_lfp
        use_cscs = 1:length(data(dd).ipsiLFPs);
        for cc = 1:length(use_cscs)
            cur_lfp = filtfilt(b,a,ipsi_csc{use_cscs(cc)}); %high-pass filter
            [S_lfp{cc},t,f]=mtspecgramc(cur_lfp(:),params.movingwin,params.chron_params);
        end
        
        %get freq ranges of interest
        uds_freqs = find(f >= params.uds_range(1) & f <= params.uds_range(2));
        lf_freqs = find(f >= params.lf_range(1) & f <= params.lf_range(2));
        use_freqs = find(f >= params.use_range(1) & f <= params.use_range(2));
        
        %get LFP power in each band
        [lfp_tot_pow,lfp_lf_pow,lfp_uds_pow] = deal(nan(length(use_cscs),length(t)));
        for cc = 1:length(use_cscs)
            lfp_uds_pow(cc,:) = trapz(f(uds_freqs),(S_lfp{cc}(:,uds_freqs)),2);
            lfp_lf_pow(cc,:) = trapz(f(lf_freqs),(S_lfp{cc}(:,lf_freqs)),2);
            lfp_tot_pow(cc,:) = trapz(f(use_freqs),S_lfp{cc}(:,use_freqs),2);
        end
        
        %%
        %average log-powers across time in each freq range (the log helps
        %normalize the var
        avg_uds_pow = nanmean(log10(lfp_uds_pow),2);
        avg_lf_pow = nanmean(log10(lfp_lf_pow),2);
        %         avg_hf_pow = nanmean(log10(lfp_hf_pow),2);
        
        %find channel with maximum UDS/HF power ratio (difference in
        %log-space)
        %         [~,ctx_ch] = max(avg_uds_pow(poss_ctx_chs)-0*avg_lf_pow(poss_ctx_chs));
        [~,ctx_ch] = max(avg_uds_pow(all_data(dd).poss_ctx_chs)); %just find channel maximizing UDS power
        all_data(dd).ctx_ch = all_data(dd).poss_ctx_chs(ctx_ch); %this is the channel used for cortical LFP
        
        uds_epochs = log10(lfp_uds_pow(all_data(dd).ctx_ch,:)) >= params.lfp_uds_thresh; %these epochs have sufficient LFP UDS
        lf_epochs = log10(lfp_lf_pow(all_data(dd).ctx_ch,:)) >= params.lfp_lf_max; %these epochs have too much LF power
        
        all_data(dd).lfp_pow_spec = mean(S_lfp{all_data(dd).ctx_ch});
        all_data(dd).mp_pow_spec = mean(S_mp);
        
        %% handle the fact that old recordings have slightly different FS
        new_taxis = csc_time(1):1/params.desired_Fs:csc_time(end);
        new_mp_data = interp1(csc_time,mp_interp,new_taxis);
        mp_Fsd = params.desired_Fs;
        mp_interp = new_mp_data;
        ctx_lfp = interp1(csc_time,ipsi_csc{ctx_ch},new_taxis);
        
        %% print spectrogram fig if desired
        if to_print
            f1 = figure();
            subplot(2,2,1); hold on
            imagesc(t,1:length(use_cscs),log10(lfp_uds_pow)); colorbar;
            plot(t,uds_epochs*4,'w','linewidth',2);
            line([0 t(end)],all_data(dd).ctx_ch + [0 0],'color','k');
            axis tight
            xlim([0 t(end)]);
            xlabel('Time (s)');
            ylabel('Channel num');
            title('LFP UDS power')
            
            subplot(2,2,3); hold on
            imagesc(t,1:length(use_cscs),log10(lfp_lf_pow)); colorbar;
            plot(t,lf_epochs*4,'w','linewidth',2);
            line([0 t(end)],all_data(dd).ctx_ch + [0 0],'color','k');
            axis tight
            xlim([0 t(end)]);
            xlabel('Time (s)');
            ylabel('Channel num');
            title('LFP low-freq power')
            
            good_uds = uds_epochs & ~lf_epochs;
            subplot(2,2,2); hold on
            pcolor(t,f,log10(S_lfp{all_data(dd).ctx_ch})'); shading flat; colorbar;
            plot(t,good_uds*2+0.1,'w','linewidth',2);
            ylim([0.1 10]);
            xlim([0 t(end)]);
            line([0 t(end)],[params.uds_range(1) params.uds_range(1)],'color','r','linestyle','--')
            line([0 t(end)],[params.uds_range(2) params.uds_range(2)],'color','r','linestyle','--')
            %         caxis([-2.5 1])
            set(gca,'yscale','log');
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            title('LFP best-channel specgram')
            
            subplot(2,2,4); hold on
            pcolor(t,f,log10(S_mp)'); shading flat; colorbar;
            plot(t,good_uds*2+0.1,'w','linewidth',2);
            ylim([0.1 10]);
            xlim([0 t(end)]);
            line([0 t(end)],[params.uds_range(1) params.uds_range(1)],'color','r','linestyle','--')
            line([0 t(end)],[params.uds_range(2) params.uds_range(2)],'color','r','linestyle','--')
            %         caxis([-2.5 1])
            set(gca,'yscale','log');
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            title('MP specgram')
            
            %             dname = find(data(dd).dir == '/',1,'last');
            %             dname = data(dd).dir(dname+1:end);
            %             fname = strcat(fig_dir,'NUDS_epochs_',data(dd).MPloc,dname);
            %             figufy(f1);
            %             set(f1,'PaperUnits','inches','PaperSize',[10 10],'PaperPosition',[0 0 20 10]);
            %             print(f1,'-dpng',fname);
            %             close(f1);
            %
        end
        %%
        removal_window = 0; %buffer window around desynch epochs to remove
        
        if ~isempty(findstr(data(dd).dir,'2015-08-16_13-06-25')) %this rec had a brief period where cell resealed
            bad_trange = [340 380]; %MP rec resealed here
            to_nan = find(t >= bad_trange(1) & t <= bad_trange(2)); 
            uds_epochs(to_nan) = 0; %force data during this time to be desynch
        end

        %mark boundary points of data epochs were going to use
        desynch_indicator = ~uds_epochs | lf_epochs; %this is data we don't want to use
        desynch_start_ids = 1+find(desynch_indicator(1:end-1) == 0 & desynch_indicator(2:end) == 1);
        desynch_stop_ids = 1+find(desynch_indicator(1:end-1) == 1 & desynch_indicator(2:end) == 0);
        
        if isempty(desynch_start_ids) && ~isempty(desynch_stop_ids)
            desynch_start_ids = [1];
        end
        
        %make sure you start and stop in desynchronized epochs in correct order
        if ~isempty(desynch_start_ids)
            if isempty(desynch_stop_ids)
                desynch_stop_ids = length(t);
            else
                if desynch_start_ids(1) > desynch_stop_ids(1)
                    desynch_start_ids = [1 desynch_start_ids];
                end
                if desynch_start_ids(end) > desynch_stop_ids(end)
                    desynch_stop_ids = [desynch_stop_ids length(t)];
                end
            end
        end
        
        if length(desynch_start_ids) ~= length(desynch_stop_ids)
            disp('error start and stop IDs not equal!')
        end
        
        %compute the duration of each putative desynchronized epoch (in seconds)
        desynch_durs = (desynch_stop_ids-desynch_start_ids)*params.movingwin(2);
        
        %now make a window around desynchronized times for data exclusion
        for w = 1:length(desynch_start_ids)
            if desynch_start_ids(w) <= removal_window
                desynch_start_ids(w) = 1;
            else
                desynch_start_ids(w) = desynch_start_ids(w)-removal_window;
            end
            if length(t)-desynch_stop_ids(w) <= removal_window
                desynch_stop_ids(w) = length(t);
            else
                desynch_stop_ids(w) = desynch_stop_ids(w)+removal_window;
            end
        end
        
        %now make sure there are no overlapping windows
        bad_desynch_start = [];
        for w = 2:length(desynch_start_ids)
            if desynch_start_ids(w) < desynch_stop_ids(w-1)
                bad_desynch_start = [bad_desynch_start w];
            end
        end
        desynch_start_ids(bad_desynch_start) = [];
        desynch_stop_ids(bad_desynch_start-1) = [];
        
        desynch_start_times = t(desynch_start_ids);
        desynch_stop_times = t(desynch_stop_ids);
        desynch_times = [desynch_start_times(:) desynch_stop_times(:)];
        desynch_ids = round(desynch_times*csc_Fsd);
        
        
        %% Compute signal features to use for UDS classification
        [mp_features,down_taxis] = hsmm_uds_get_lf_features(mp_interp,params.desired_Fs,params.hmm_Fs,params.mp_feature_frange);
        
        %get a broader band MP signal
        mp_features_bb = hsmm_uds_get_lf_features(mp_interp,params.desired_Fs,params.desired_Fs,params.mp_feature_frange_bb);
        
        % LOCATE THE SEGMENTS CONTAINING UDS
        UDS_segs = hsmm_uds_get_uds_segments(desynch_times,params.hmm_Fs,length(mp_features),params.min_seg_dur); %Nx2 matrix containing the index values of the beginning and end of each UDS segment
        is_uds = false(size(mp_features));
        for ii = 1:size(UDS_segs,1)
            is_uds(UDS_segs(ii,1):UDS_segs(ii,2)) = true;
        end
        
        %get UDS segs for higher-sample rate MP
        new_tax = (1:length(mp_features_bb))/params.desired_Fs;
        new_seg_inds = round(interp1(new_tax,1:length(new_tax),down_taxis(UDS_segs)));
        new_seg_inds(isnan(new_seg_inds(:,1)),1) = 1;
        new_seg_inds(isnan(new_seg_inds(:,2)),2) = length(new_tax);
        
        %%
        if sum(is_uds) > 0
            
            mp_dist_xx = linspace(prctile(mp_features_bb(is_uds),0.1),prctile(mp_features_bb(is_uds),99.9),params.n_dens_bins); %MP amplitude axis for density estimation
%             mp_dist_xx = linspace(prctile(mp_features(is_uds),0.1),prctile(mp_features(is_uds),99.9),params.n_dens_bins); %MP amplitude axis for density estimation
            
            %initializations
            [seg_dip,seg_p_value] = deal(nan(size(UDS_segs,1),1));
            mp_middle_crosspts = nan(size(UDS_segs,1),1);
            mp_bandwidths = nan(size(UDS_segs,1),1);
            mp_seg_pdf = nan(size(UDS_segs,1),length(mp_dist_xx));
            for ii = 1:size(UDS_segs,1) %for each UDS segment
%                 cur_inds = UDS_segs(ii,1):UDS_segs(ii,2);
%                 [mp_kpdf,~,mp_bandwidths(ii)] = ksdensity(mp_features(cur_inds),mp_dist_xx); %get MP density estimate for this segment
                cur_inds = new_seg_inds(ii,1):new_seg_inds(ii,2);
                [mp_kpdf,~,mp_bandwidths(ii)] = ksdensity(mp_features_bb(cur_inds),mp_dist_xx); %get MP density estimate for this segment
                mp_seg_pdf(ii,:) = mp_kpdf;
                
                %check the significance of bimodality
                [seg_dip(ii), seg_p_value(ii), xlow,xup]=hartigansdipsigniftest(downsample(mp_features_bb(cur_inds),2),params.dip_Nboot);
%                 [seg_dip(ii), seg_p_value(ii), xlow,xup]=hartigansdipsigniftest(mp_features(cur_inds),params.dip_Nboot);
                
                %if there is bimodality, compute the mean of up and down
                %state amps
                if seg_p_value(ii) < params.min_dip_pvalue
                    [~,dens_peaks] = findpeaks(mp_kpdf,'npeaks',2,'sortstr','descend','minpeakdistance',ceil(0.25/median(diff(mp_dist_xx))));
                    dens_peaks = mp_dist_xx(dens_peaks);
                    if length(dens_peaks) == 2
                        mp_middle_crosspts(ii) = nanmean(dens_peaks); %set the separation point to be the midpoint between the two modes
                    else
                        mp_middle_crosspts(ii) = nanmedian(mp_features_bb(cur_inds)); %if unimodal, just set this as the median
%                         mp_middle_crosspts(ii) = nanmedian(mp_features(cur_inds)); %if unimodal, just set this as the median
                    end
                end
            end
            avg_mp_bandwidth = nanmean(mp_bandwidths(seg_p_value < params.min_dip_pvalue)); %avg density estimation bandwidth across UDS segments with bimodality
            all_data(dd).mp_seg_pdf = mp_seg_pdf;
            all_data(dd).mp_dist_xx = mp_dist_xx;
            all_data(dd).avg_mp_bandwidth = avg_mp_bandwidth;
            all_data(dd).mp_seg_dipStat = seg_dip;
            all_data(dd).mp_seg_dip_p = seg_p_value;
            
            %% recompute usable data segments, using only segments with bimodal MP
            all_data(dd).lfp_uds_dur = sum(is_uds)/params.desired_Fs;
            all_data(dd).tot_dur = length(is_uds)/params.desired_Fs;
            
            bad_segs = find(seg_p_value > params.min_dip_pvalue);
            fprintf('%d/%d bad segs\n',length(bad_segs),length(seg_p_value));
            UDS_segs(bad_segs,:) = []; %get rid of the UDS segments without bimodal MP
            
            %% if we have some UDS
            if size(UDS_segs,1) > 0
                
                
                %%
                if params.compute_state_seqs
                    %fit an HMM to MP
                    clear hmm_params
                    hmm_params.meantype = 'variable'; %use 'variable' for time-varying state means. otherwise use 'fixed'
                    hmm_params.UDS_segs = UDS_segs;
                    hmm_params.movingwin = [20 5]; %moving window parameters [windowLength windowSlide](in seconds) for computing time-varying state means
                    [hmm] = hsmm_uds_initialize(mp_features(:),params.hmm_Fs,hmm_params);
                    
                    hmm.min_mllik = -8; %set minimum value of the maximal log likelihood for either state before switching to robust estimator.  If you don't want to use this feature, set this to -Inf
                    [hmm,gamma] = hsmm_uds_train_hmm(hmm,mp_features(:));
                    
                    [mp_state_seq,llik_best] = hsmm_uds_viterbi_hmm(hmm,mp_features(:));
                    
                    hmm.min_state_dur = 1;
                    hmm.max_state_dur = round(params.hmm_Fs*20);
                    hmm.dur_range = (1:hmm.max_state_dur)/hmm.Fs;
                    hsmm = hmm;
                    smoothed_seq = thresh_state_smooth_seg(mp_state_seq,params.hmm_Fs,100,100); %smooth out excessively short state durations to improve model fits
                    [state_durations] = compute_state_durations_seg(smoothed_seq,params.hmm_Fs); %compute the set of up and down state durations
                    for i = 1:hmm.K
                        %estimate the empirical state duration pmf
                        emp_pmf = hist(state_durations{i},hmm.dur_range);
                        emp_pmf = emp_pmf/sum(emp_pmf);
                        [mu,lambda] = inverse_gauss_mlfit(state_durations{i});%compute the ML parameters of the inverse gaussian dist
                        ig_pmf = inverse_gaussian_pmf(hmm.dur_range,mu,lambda);%estimate the IG PMF
                        hsmm.state(i).dur_type = 'inv_gauss'; %use inverse gaussian model in all cases
                        hsmm.state(i).dur_pars = [mu lambda];
                        hsmm.P(i,:) = ig_pmf;
                    end
                    % enforce any minimum state duration and renormalize
                    for i = 1:hmm.K
                        if hsmm.min_state_dur > 1
                            hsmm.P(i,1:hsmm.min_state_dur-1) = zeros(1,hsmm.min_state_dur-1);
                            hsmm.P(i,:) = hsmm.P(i,:)/sum(hsmm.P(i,:));
                        end
                    end
                    
                    % estimate HSMM params, and viterbi seq
                    [hsmm,gamma,hmm_window_post]=hsmm_uds_train_hsmm(hsmm,mp_features);
                    [hsmm_state_seq,max_lik] = hsmm_uds_viterbi_hsmm(hsmm,mp_features);
                    hsmm.posteriors = hmm_window_post;
                    hsmm.max_lik = max_lik;
                    pert_range_bb = [-0.15 0.15];
                    [mp_state_seq_bb] = hsmm_uds_pert_optimize_transitions_v2(mp_features_bb,hsmm,hsmm_state_seq,gamma,params.hmm_Fs,params.desired_Fs,pert_range_bb);
                    
                    %                     hmm.state(1).dur_type = 'geometric'; hmm.state(2).dur_type = 'geometric';
                    %                     pert_range_bb = [-0.15 0.15];
                    %                     mp_state_seq_bb = hsmm_uds_pert_optimize_transitions(mp_features_bb,hmm,mp_state_seq,gamma,params.hmm_Fs,params.desired_Fs,pert_range_bb);
                    %
                    
                    [new_seg_inds] = resample_uds_seg_inds_v2(UDS_segs,params.hmm_Fs,params.desired_Fs,mp_state_seq_bb);
                    [mp_state_durations] = compute_state_durations_seg(mp_state_seq_bb,params.desired_Fs);
                    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq_bb);
                    
                    mp_state_vec = nan(size(mp_features_bb));
                    for i = 1:size(new_seg_inds,1)
                        mp_state_vec(new_seg_inds(i,1):new_seg_inds(i,2)) = mp_state_seq_bb{i};
                    end
                    
                    save hsmm_state_seq mp_state_vec new_seg_inds hsmm params hmm
                else
                    new_tax = (1:length(mp_features_bb))/params.desired_Fs;
                    new_seg_inds = round(interp1(new_tax,1:length(new_tax),down_taxis(UDS_segs)));
                    new_seg_inds(isnan(new_seg_inds(:,1)),1) = 1;
                    new_seg_inds(isnan(new_seg_inds(:,2)),2) = length(new_tax);
                end
                
                is_uds = false(size(mp_features_bb));
                for i = 1:size(new_seg_inds,1)
                    is_uds(new_seg_inds(i,1):new_seg_inds(i,2)) = true;
                end
                
                all_data(dd).uds_dur = sum(is_uds)/params.desired_Fs;
                
                %% DETECT LOCAL EXTREMA OF FILTERED LFP SIGNAL
                % compute ctx LFP features
                [lfp_features,t_axis] = hsmm_uds_get_lf_features(ctx_lfp,params.desired_Fs,params.desired_Fs,params.uds_range); %'low-frequency amplitude'
                lfp_features = lfp_features/robust_std_dev(lfp_features(is_uds));
                
                %get local minima and maxima of lfp signal
                lfp_peaks = find(diff(sign(diff(lfp_features))) < 0)+1;
                lfp_valleys = find(diff(sign(diff(lfp_features))) > 0)+1;
                lfp_peaks(~is_uds(lfp_peaks)) = [];
                lfp_valleys(~is_uds(lfp_valleys)) = [];
                
                %get amplitudes of extrema
                lfp_peak_amps = lfp_features(lfp_peaks);
                lfp_valley_amps = -lfp_features(lfp_valleys);
                
                peak_checks = prctile(lfp_peak_amps,params.amp_prc_bins);
                valley_checks = prctile(lfp_valley_amps,params.amp_prc_bins);
                all_data(dd).peak_checks = peak_checks;
                all_data(dd).valley_checks = valley_checks;
                
                %% compute MP-LFP xcorr function
                [xcfun,xclags] = xcov(lfp_features(is_uds),mp_features_bb(is_uds),round(0.5*params.desired_Fs),'coeff');
                [all_data(dd).xc_peakval,all_data(dd).xc_peakloc] = max(xcfun);
                all_data(dd).xc_peakloc = xclags(all_data(dd).xc_peakloc);
                if params.shift_lfp_peaks %if accounting for temporal delay with LFP
                    lfp_peak_amps_shift = lfp_peak_amps;
                    lfp_valley_amps_shift = lfp_valley_amps;
                    lfp_peak_shift = lfp_peaks - all_data(dd).xc_peakloc;
                    lfp_valley_shift = lfp_valleys - all_data(dd).xc_peakloc;
                    bad_peaks = find(lfp_peak_shift < 1 | lfp_peak_shift > length(lfp_features));
                    bad_valleys = find(lfp_valley_shift < 1 | lfp_valley_shift > length(lfp_features));
                    lfp_peak_shift(bad_peaks) = []; lfp_peak_amps_shift(bad_peaks) = [];
                    lfp_valley_shift(bad_valleys) = []; lfp_valley_amps_shift(bad_valleys) = [];
                    cur_peak_amps = lfp_peak_amps_shift;
                    cur_valley_amps = lfp_valley_amps_shift;
                    cur_peaks = lfp_peak_shift;
                    cur_valleys = lfp_valley_shift;
                else
                    cur_peak_amps = lfp_peak_amps;
                    cur_valley_amps = lfp_valley_amps;
                    cur_peaks = lfp_peaks;
                    cur_valleys = lfp_valleys;
                end
                
                %%
                if params.compute_state_seqs
                    %find the set of down states that might have skipped ctx up
                    %blips
                    min_down_dur = 2*params.state_buffer;
                    poss_pers_downs = find(mp_state_durations{1} >= min_down_dur);
                    
                    %find the cortical up-blips skipped in each MP down
                    skipped_blips = cell(length(poss_pers_downs),1);
                    for ii = 1:length(poss_pers_downs)
                        cur_down_inds = (down_trans_inds(poss_pers_downs(ii)) + ...
                            round(params.desired_Fs*params.state_buffer)):(up_trans_inds(poss_pers_downs(ii)+1) - ...
                            round(params.desired_Fs*params.state_buffer));
                        skipped_blips{ii} = cur_peak_amps(ismember(cur_peaks,cur_down_inds));
                    end
                    
                    %calculate stats on skipped ctx up blips of different sizes
                    [peak_prob,skip_cnt,skip_prob] = deal(nan(length(peak_checks),1));
                    for ii = 1:length(peak_checks)
                        peak_prob(ii) = sum(cur_peak_amps >= peak_checks(ii))/sum(~isnan(cur_peak_amps));
                        skip_cnt(ii) = sum(cellfun(@(x) sum(x >= peak_checks(ii)),skipped_blips));
                        skip_prob(ii) = sum(cellfun(@(x) any(x >= peak_checks(ii)),skipped_blips))/length(poss_pers_downs);
                    end
                    
                    all_data(dd).lfp_uppeak_prob = peak_prob;
                    all_data(dd).pdown_cnt = skip_cnt;
                    all_data(dd).pdown_prob = skip_prob;
                    %% now repeat for pers ups
                    poss_pers_ups = find(mp_state_durations{2} >= min_down_dur); %set of MP up states that might skip down-blips
                    
                    %find ctx down-blips skipped by each mp up
                    skipped_blips = cell(length(poss_pers_ups),1);
                    for ii = 1:length(poss_pers_ups)
                        cur_up_inds = (up_trans_inds(poss_pers_ups(ii)) + ...
                            round(params.desired_Fs*params.state_buffer)):(down_trans_inds(poss_pers_ups(ii)) - ...
                            round(params.desired_Fs*params.state_buffer));
                        skipped_blips{ii} = cur_valley_amps(ismember(cur_valleys,cur_up_inds));
                    end
                    
                    %compute stats on skipped ctx down blips
                    [peak_prob,skip_cnt,skip_prob] = deal(nan(length(valley_checks),1));
                    for ii = 1:length(valley_checks)
                        peak_prob(ii) = sum(cur_valley_amps >= valley_checks(ii))/sum(~isnan(cur_valley_amps));
                        skip_cnt(ii) = sum(cellfun(@(x) sum(x >= valley_checks(ii)),skipped_blips));
                        skip_prob(ii) = sum(cellfun(@(x) any(x >= valley_checks(ii)),skipped_blips))/length(poss_pers_ups);
                    end
                    
                    all_data(dd).lfp_downpeak_prob = peak_prob;
                    all_data(dd).pup_cnt = skip_cnt;
                    all_data(dd).pup_prob = skip_prob;
                    
                end
                %%
                %temporal window around each LFP extremum to pull data snippets
                back_win = round(params.back_Twin*params.desired_Fs);
                for_win = round(params.for_Twin*params.desired_Fs);
                
                avg_cross_pt = nanmean(mp_middle_crosspts); %avg definition of state trheshold
                cross_pt_mp = find(mp_dist_xx >= avg_cross_pt,1,'first'); %index of threshold MP value
                
                %compute triggered average MP data
                all_up_tas = nan(length(peak_checks),back_win + for_win+1);
                all_up_dists = nan(length(peak_checks),back_win+for_win+1,length(mp_dist_xx));
                all_down_tas = nan(length(peak_checks),back_win + for_win+1);
                all_down_dists = nan(length(peak_checks),back_win+for_win+1,length(mp_dist_xx));
                
                [all_down_stas,all_up_stas] = deal(nan(length(peak_checks),back_win+for_win+1));
                for jj = 1:length(peak_checks)
                    cur_events = lfp_valleys(lfp_valley_amps >= valley_checks(jj)); %all LFP minima in the current set
                    if length(cur_events) >= params.min_n_events
                        [cur_ta,cur_lags,~,~,cur_mat] = get_event_trig_avg_v3(mp_features_bb,cur_events,back_win,for_win,2); %trig avg of MP amp
                        all_down_tas(jj,:) = cur_ta;
                        for ii = 1:length(cur_lags)
                            all_down_dists(jj,ii,:) = ksdensity(cur_mat(:,ii),mp_dist_xx,'bandwidth',all_data(dd).avg_mp_bandwidth);
%                             all_down_dists(jj,ii,:) = jmm_smooth_1d_cor(squeeze(all_down_dists(jj,ii,:)),params.dens_smth_sig); %gaussian smoothing
                        end
                        
                        if params.compute_state_seqs
                            [cur_sta] = get_event_trig_avg_v3(mp_state_vec,cur_events,back_win,for_win); %trig avg of MP state
                            all_down_stas(jj,:) = cur_sta;
                        end
                    end
                    
                    cur_events = lfp_peaks(lfp_peak_amps >= peak_checks(jj));
                    if length(cur_events) >= params.min_n_events
                        [cur_ta,cur_lags,~,~,cur_mat] = get_event_trig_avg_v3(mp_features_bb,cur_events,back_win,for_win,2);
                        all_up_tas(jj,:) = cur_ta;
                        for ii = 1:length(cur_lags)
                            all_up_dists(jj,ii,:) = ksdensity(cur_mat(:,ii),mp_dist_xx,'bandwidth',all_data(dd).avg_mp_bandwidth);
%                             all_up_dists(jj,ii,:) = jmm_smooth_1d_cor(squeeze(all_up_dists(jj,ii,:)),params.dens_smth_sig);
                        end
                        
                        if params.compute_state_seqs
                            [cur_sta] = get_event_trig_avg_v3(mp_state_vec,cur_events,back_win,for_win);
                            all_up_stas(jj,:) = cur_sta;
                        end
                    end
                end
                
                all_data(dd).utrig_Pmp_down = trapz(mp_dist_xx(1:(cross_pt_mp-1)),all_up_dists(:,:,1:(cross_pt_mp-1)),3);
                all_data(dd).utrig_Pmp_up = trapz(mp_dist_xx((cross_pt_mp+1):end),all_up_dists(:,:,(cross_pt_mp+1):end),3);
                all_data(dd).dtrig_Pmp_down = trapz(mp_dist_xx(1:(cross_pt_mp-1)),all_down_dists(:,:,1:(cross_pt_mp-1)),3);
                all_data(dd).dtrig_Pmp_up = trapz(mp_dist_xx((cross_pt_mp+1):end),all_down_dists(:,:,(cross_pt_mp+1):end),3);
                
                all_data(dd).utrig_mpstate = all_up_stas;
                all_data(dd).dtrig_mpstate = all_down_stas;
                all_data(dd).utrig_mpavg = all_up_tas;
                all_data(dd).dtrig_mpavg = all_down_tas;
                all_data(dd).uptrig_dists = all_up_dists;
                all_data(dd).downtrig_dists = all_down_dists;
                
                if params.compute_state_seqs
                    all_data(dd).avg_stateprob = nanmean(mp_state_vec(is_uds));
                end
                
                %%
            end
        end
        
    end
end

%%
cd ~/Analysis/Mayank/sleep/
save sleep_uds_analysis5 all_data params

%% select usable recs
min_uds_dur = 50; %minimum duration of data with LFP UDS and MP bimodality
poss_use = find(arrayfun(@(x) length(x.uds_dur),all_data)); %datasets that might be used
poss_use([all_data(poss_use).uds_dur] == 0) = []; %get rid of ones with no usable UDS
use_uds = poss_use([all_data(poss_use).uds_dur] >= min_uds_dur); %recs with at least minimal UDS dur
use_cells = use_uds;

mec_neurons = find(arrayfun(@(x) strcmp(x.MPloc,'MEC'),data));
ctx_neurons = find(arrayfun(@(x) strcmp(x.MPloc,'Ctx'),data));

use_mec_neurons = intersect(mec_neurons,use_cells);
use_ctx_neurons = intersect(ctx_neurons,use_cells);

%% plot MEC pers prob as function of blip amp
if params.compute_state_seqs
    mec_pdown_prob = cell2mat(arrayfun(@(x) x.pdown_prob,all_data(use_mec_neurons),'uniformoutput',0));
    ctx_pdown_prob = cell2mat(arrayfun(@(x) x.pdown_prob,all_data(use_ctx_neurons),'uniformoutput',0));
    
    mec_pup_prob = cell2mat(arrayfun(@(x) x.pup_prob,all_data(use_mec_neurons),'uniformoutput',0));
    ctx_pup_prob = cell2mat(arrayfun(@(x) x.pup_prob,all_data(use_ctx_neurons),'uniformoutput',0));
    
    mec_peak_checks = cell2mat(arrayfun(@(x) x.peak_checks',all_data(use_mec_neurons),'uniformoutput',0));
    ctx_peak_checks = cell2mat(arrayfun(@(x) x.peak_checks',all_data(use_ctx_neurons),'uniformoutput',0));
    mec_valley_checks = cell2mat(arrayfun(@(x) x.valley_checks',all_data(use_mec_neurons),'uniformoutput',0));
    ctx_valley_checks = cell2mat(arrayfun(@(x) x.valley_checks',all_data(use_ctx_neurons),'uniformoutput',0));
    
    %test for median difference in persistence prob at each amplitude
    %threshold
    [pup_p,pdown_p] = deal(nan(length(params.amp_prc_bins),1));
    for ii = 1:length(params.amp_prc_bins)
        pup_p(ii) = ranksum(mec_pup_prob(ii,:),ctx_pup_prob(ii,:));
        pdown_p(ii) = ranksum(mec_pdown_prob(ii,:),ctx_pdown_prob(ii,:));
    end
    
    f1 = figure();
    subplot(3,1,1); hold on
    plot(mec_peak_checks,mec_pdown_prob,'ko-','linewidth',1);
    plot(ctx_peak_checks,ctx_pdown_prob,'ro-','linewidth',1);
    ylim([0 1]);
    ylabel('Prob MP persistence');
    xlabel('LFP amplitude (z)');
    title('LFP UP blips');
    
    subplot(3,1,2); hold on
    plot(-mec_peak_checks,mec_pup_prob,'ko-','linewidth',1);
    plot(-ctx_peak_checks,ctx_pup_prob,'ro-','linewidth',1);
    ylim([0 1]);
    ylabel('Prob MP persistence');
    xlabel('LFP amplitude (z)');
    title('LFP DOWN blips');
    
    subplot(3,1,3);hold on
    plot(params.amp_prc_bins,pup_p,'o-')
    plot(params.amp_prc_bins,pdown_p,'go-')
    legend('Pers UPs','Pers DOWNs');
    xlabel('LFP blip amplitude (percentile)');
    ylabel('P-value');
    set(gca,'yscale','log');
    line([0 100],[0.05 0.05],'color','r','linestyle','--');
%         fig_width = 6; rel_height = 2;
%     figufy(f1);
%     fname = [fig_dir 'cortical_blip_pers2.pdf'];
%     exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
%     close(f1);
    
end

% %% plot blip-trig avg state probs
% xr = [-0.75 0.75];
% if compute_state_seqs
%     clear *_up_ts *_down_ts
%     for ii = 1:length(use_mec_neurons)
%         mec_up_ts(ii,:,:) = all_data(use_mec_neurons(ii)).utrig_mpstate;
%         mec_down_ts(ii,:,:) = all_data(use_mec_neurons(ii)).dtrig_mpstate;
%     end
%     for ii = 1:length(use_ctx_neurons)
%         ctx_up_ts(ii,:,:) = all_data(use_ctx_neurons(ii)).utrig_mpstate;
%         ctx_down_ts(ii,:,:) = all_data(use_ctx_neurons(ii)).dtrig_mpstate;
%     end
%     
%     %     mec_up_ts = bsxfun(@minus,mec_up_ts,[all_data(use_mec_neurons).avg_stateprob]');
%     %     ctx_up_ts = bsxfun(@minus,ctx_up_ts,[all_data(use_ctx_neurons).avg_stateprob]');
%     %     mec_down_ts = bsxfun(@minus,mec_down_ts,[all_data(use_mec_neurons).avg_stateprob]');
%     %     ctx_down_ts = bsxfun(@minus,ctx_down_ts,[all_data(use_ctx_neurons).avg_stateprob]');
%     
%     f1 = figure();
%     for ii = 1:length(peak_checks)
%         subplot(3,2,ii); hold on
%         plot(cur_lags/desired_Fs,squeeze(mec_up_ts(:,ii,:)),'k-');
%         plot(cur_lags/desired_Fs,squeeze(ctx_up_ts(:,ii,:)),'r-');
%         plot(cur_lags/desired_Fs,squeeze(nanmean(mec_up_ts(:,ii,:))),'k','linewidth',3);
%         plot(cur_lags/desired_Fs,squeeze(nanmean(ctx_up_ts(:,ii,:))),'r','linewidth',3);
%         xlim(xr);
%     end
%     
%     f1 = figure();
%     for ii = 1:length(peak_checks)
%         subplot(3,2,ii); hold on
%         plot(cur_lags/desired_Fs,squeeze(mec_down_ts(:,ii,:)),'k-');
%         plot(cur_lags/desired_Fs,squeeze(ctx_down_ts(:,ii,:)),'r-');
%         plot(cur_lags/desired_Fs,squeeze(nanmean(mec_down_ts(:,ii,:))),'k','linewidth',3);
%         plot(cur_lags/desired_Fs,squeeze(nanmean(ctx_down_ts(:,ii,:))),'r','linewidth',3);
%         xlim(xr);
%     end
% end

%% plot blip-trig avg MP amps
xr = [-params.back_Twin params.for_Twin];
clear *_up_ta *_down_ta
for ii = 1:length(use_mec_neurons)
    mec_up_ta(ii,:,:) = all_data(use_mec_neurons(ii)).utrig_mpavg;
    mec_down_ta(ii,:,:) = all_data(use_mec_neurons(ii)).dtrig_mpavg;
end
for ii = 1:length(use_ctx_neurons)
    ctx_up_ta(ii,:,:) = all_data(use_ctx_neurons(ii)).utrig_mpavg;
    ctx_down_ta(ii,:,:) = all_data(use_ctx_neurons(ii)).dtrig_mpavg;
end
use_amp_prc_bins = [20 40 60 80];
amp_bin_inds = find(ismember(params.amp_prc_bins,use_amp_prc_bins));

cur_lags = -round(params.desired_Fs*params.back_Twin):round(params.desired_Fs*params.for_Twin);
f1 = figure();
for ii = 1:length(amp_bin_inds)
    subplot(2,2,ii); hold on
    %     plot(cur_lags/params.desired_Fs,squeeze(nanmean(mec_up_ta(:,ii,:))),'k','linewidth',3);
    %     plot(cur_lags/params.desired_Fs,squeeze(nanmean(ctx_up_ta(:,ii,:))),'r','linewidth',3);
    shadedErrorBar(cur_lags/params.desired_Fs,squeeze(nanmean(mec_up_ta(:,amp_bin_inds(ii),:))),...
        squeeze(nanstd(mec_up_ta(:,amp_bin_inds(ii),:)))/sqrt(length(use_mec_neurons)));
    shadedErrorBar(cur_lags/params.desired_Fs,squeeze(nanmean(ctx_up_ta(:,amp_bin_inds(ii),:))),...
        squeeze(nanstd(ctx_up_ta(:,amp_bin_inds(ii),:)))/sqrt(length(use_ctx_neurons)),{'color','r'});
%     plot(cur_lags/params.desired_Fs,squeeze(mec_up_ta(:,amp_bin_inds(ii),:)),'k-');
%     plot(cur_lags/params.desired_Fs,squeeze(ctx_up_ta(:,amp_bin_inds(ii),:)),'r-');
    xlim(xr); ylim([-1 1.2]);
    grid on
    xlabel('Time (s)');
    ylabel('MP amplitude (z)');
    title(sprintf('%d percentile',use_amp_prc_bins(ii)));
end

f2 = figure();
for ii = 1:length(amp_bin_inds)
    subplot(2,2,ii); hold on
%     plot(cur_lags/params.desired_Fs,squeeze(nanmean(mec_down_ta(:,ii,:))),'k','linewidth',3);
%     plot(cur_lags/params.desired_Fs,squeeze(nanmean(ctx_down_ta(:,ii,:))),'r','linewidth',3);
    shadedErrorBar(cur_lags/params.desired_Fs,squeeze(nanmean(mec_down_ta(:,amp_bin_inds(ii),:))),...
        squeeze(nanstd(mec_down_ta(:,amp_bin_inds(ii),:)))/sqrt(length(use_mec_neurons)));
    shadedErrorBar(cur_lags/params.desired_Fs,squeeze(nanmean(ctx_down_ta(:,amp_bin_inds(ii),:))),...
        squeeze(nanstd(ctx_down_ta(:,amp_bin_inds(ii),:)))/sqrt(length(use_ctx_neurons)),{'color','r'});
%     plot(cur_lags/params.desired_Fs,squeeze(mec_down_ta(:,amp_bin_inds(ii),:)),'k-');
%     plot(cur_lags/params.desired_Fs,squeeze(ctx_down_ta(:,amp_bin_inds(ii),:)),'r-');
    xlim(xr); ylim([-1.5 1]); 
    grid on
    xlabel('Time (s)');
    ylabel('MP amplitude (z)');
    title(sprintf('%d percentile',use_amp_prc_bins(ii)));
end

    fig_width = 9; rel_height = 0.9;
figufy(f1);
fname = [fig_dir 'UP_bliptrig_avg_MP2.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

figufy(f2);
fname = [fig_dir 'DOWN_bliptrig_avg_MP2.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);