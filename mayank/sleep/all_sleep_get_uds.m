clear all
% close all
cd ~/Analysis/Mayank/sleep/
load sleep_dirs

fig_dir = '/Users/james/Analysis/Mayank/sleep/sleep_figs2/';
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))
addpath('~/James_scripts/hsmm_uds_toolbox/')

%% DEFINE FREQUENCY RANGES
uds_range = [0.5 3]; %range of frequencies for defining UDS
% hf_range = [30 80]; 
lf_range = [0 0.2]; %range of low-freq power to characterize LF artifacts
use_range = [0.5 80]; %overall range of usable power (excluding LF artifacts)

mp_feature_frange = [0 5]; %freq range for defining MP features used to classify MP UDS
mp_feature_frange_bb = [0 20]; %freq range for computing MP features

n_dens_bins = 100; %number of pts to use for density estimation
dens_smth_sig = 3; %sigma for visualizing density estimates (in units of bins)

mp_lcf = 0.01; %lcf for hp-filter on MP signal (just gets rid of drift)
lfp_lcf = 0.2; %lcf for hp-filter on LFP signals (gets rid of LF artifact)
desired_Fs = 200; %want a uniform Fs

min_seg_dur = 20; %minimum duration of UDS segments used for analysis (in sec)

dip_Nboot = 500; %number of bootstrap samples for computing dip stat significance
min_dip_pvalue = 0.05; %alpha on dip test

state_buffer = 0.25; %buffer around state transitions (sec) for counting cortical blips as skipped

peak_checks = [0 0.5 1 1.5 2 3]; %amplitude ranges of cortical LFP blips (z)

to_print = false;

%for normalized LFP power
lfp_uds_thresh = -0.75; %min log power concentrated in the UDS band
lfp_lf_max = 2; %max log power in the LF band (artifact detection)

%%
% for dd = 14
for dd = 1:length(data)
    %     close all
    %%
    cd(data(dd).dir)
    pwd
    load('procData.mat','mp_*','ipsi_*','heka*','csc_time');
    
    has_ipsi = true; %initialize
    %range of channels that might be used as our cortical LFP
    if length(data(dd).ipsiLFPs) == 8
        poss_ctx_chs = 1:4;
    elseif length(data(dd).ipsiLFPs) == 16
        poss_ctx_chs = 1:8;
    elseif length(data(dd).ipsiLFPs) == 7
        poss_ctx_chs = 4:7;
    else
        has_ipsi = false;
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
    if has_ipsi && mp_Fs >= 150 %if there are any ipsilateral LFPs
        
        %%
        mp_Fsd = mp_Fs/mp_dsf;
        [b,a] = butter(2,mp_lcf/(mp_Fsd/2),'high');
        mp_d = decimate(mp_d,mp_dsf);
        mp_data = filtfilt(b,a,mp_d); %filter out very slow drift in MP signal
        mp_time = downsample(mp_t,mp_dsf);
        
        if ~isempty(heka_data)
            heka_Fs = 1/nanmedian(diff(heka_time));
            heka_data = decimate(heka_data,heka_dsf);
            heka_Fsd = heka_Fs/heka_dsf;
            [b,a] = butter(2,mp_lcf/(heka_Fsd/2),'high'); %filter out very slow drift in MP signal
            mp_data = filtfilt(b,a,heka_data);
            mp_time = downsample(heka_time,heka_dsf);
        end
        
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
        mp_time = mp_time(use_mp_inds); mp_data = mp_data(use_mp_inds);
        
        %% process LFP data
        %find usable LFP times
        csc_time = downsample(csc_time,csc_dsf);
        csc_inds = find(csc_time >= data(dd).good_bounds(1) & csc_time <= data(dd).good_bounds(2));
        csc_time = csc_time(csc_inds);
        
        %filter LFPs
        csc_Fsd = csc_Fs/csc_dsf;
        [bb,aa] = butter(2,lfp_lcf/(csc_Fsd/2),'high');
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
        params.Fs = csc_Fsd;
        params.tapers = [2 3];
        params.fpass = [0 80];
        movingwin = [20 5];
        
        [b,a] = butter(2,[1/movingwin(1)]/(csc_Fs/2),'high'); %filter out stuff that's slower than the window size
                
        [S_mp,t,f] = mtspecgramc(mp_interp(:),movingwin,params); %MP spectrogram
        
        clear S_lfp
        use_cscs = 1:length(data(dd).ipsiLFPs);
        for cc = 1:length(use_cscs) 
            cur_lfp = filtfilt(b,a,ipsi_csc{use_cscs(cc)}); %high-pass filter
            [S_lfp{cc}]=mtspecgramc(cur_lfp(:),movingwin,params);
        end
        
        uds_freqs = find(f >= uds_range(1) & f <= uds_range(2));
%         hf_freqs = find(f >= hf_range(1) & f <= hf_range(2));
        lf_freqs = find(f >= lf_range(1) & f <= lf_range(2));
        use_freqs = find(f >= use_range(1) & f <= use_range(2));
        
        %get LFP power in each band
        lfp_uds_pow = nan(length(use_cscs),length(t));
%         lfp_hf_pow = nan(length(use_cscs),length(t));
        lfp_lf_pow = nan(length(use_cscs),length(t));
        lfp_tot_pow = nan(length(use_cscs),length(t));
        for cc = 1:length(use_cscs)
            lfp_uds_pow(cc,:) = trapz(f(uds_freqs),(S_lfp{cc}(:,uds_freqs)),2);
%             lfp_hf_pow(cc,:) = trapz(f(hf_freqs),(S_lfp{cc}(:,hf_freqs)),2);
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
        [~,ctx_ch] = max(avg_uds_pow(poss_ctx_chs)); %just find channel maximizing UDS power
        ctx_ch = poss_ctx_chs(ctx_ch); %this is the channel used for cortical LFP
        
        uds_epochs = log10(lfp_uds_pow(ctx_ch,:)) >= lfp_uds_thresh; %these epochs have sufficient LFP UDS
        lf_epochs = log10(lfp_lf_pow(ctx_ch,:)) >= lfp_lf_max; %these epochs have too much LF power
        
        all_data(dd).lfp_pow_spec = mean(S_lfp{ctx_ch});
        all_data(dd).mp_pow_spec = mean(S_mp);
        
        %% handle the fact that old recordings have slightly different FS
        %         if mp_Fsd > 201 %these are the old recs
        new_taxis = csc_time(1):1/desired_Fs:csc_time(end);
        new_mp_data = interp1(csc_time,mp_interp,new_taxis);
        mp_Fsd = desired_Fs;
        mp_interp = new_mp_data;
        ctx_lfp = interp1(csc_time,ipsi_csc{ctx_ch},new_taxis);
        %         end
        
        all_data(dd).Fsd = desired_Fs;
        all_data(dd).ctx_ch = ctx_ch;
        
        %%
        if to_print
            f1 = figure();
            subplot(2,2,1); hold on
            imagesc(t,1:length(use_cscs),log10(lfp_uds_pow)); colorbar;
            plot(t,uds_epochs*4,'w','linewidth',2);
            line([0 t(end)],ctx_ch + [0 0],'color','k');
            axis tight
            xlim([0 t(end)]);
            xlabel('Time (s)');
            ylabel('Channel num');
            title('LFP UDS power')
            
            subplot(2,2,3); hold on
            imagesc(t,1:length(use_cscs),log10(lfp_lf_pow)); colorbar;
            plot(t,lf_epochs*4,'w','linewidth',2);
            line([0 t(end)],ctx_ch + [0 0],'color','k');
            axis tight
            xlim([0 t(end)]);
            xlabel('Time (s)');
            ylabel('Channel num');
            title('LFP low-freq power')
            
            good_uds = uds_epochs & ~lf_epochs;
            subplot(2,2,2); hold on
            pcolor(t,f,log10(S_lfp{ctx_ch})'); shading flat; colorbar;
            plot(t,good_uds*2+0.1,'w','linewidth',2);
            ylim([0.1 10]);
            xlim([0 t(end)]);
            line([0 t(end)],[uds_range(1) uds_range(1)],'color','r','linestyle','--')
            line([0 t(end)],[uds_range(2) uds_range(2)],'color','r','linestyle','--')
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
            line([0 t(end)],[uds_range(1) uds_range(1)],'color','r','linestyle','--')
            line([0 t(end)],[uds_range(2) uds_range(2)],'color','r','linestyle','--')
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
        desynch_durs = (desynch_stop_ids-desynch_start_ids)*movingwin(2);
        
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
        hmm_dsf = 2; %down-sample-factor for HMM signal featurs
        hmm_Fs = desired_Fs/hmm_dsf; %sample-frequency of HMM signal features (Hz)
        mp_features = hsmm_uds_get_lf_features(mp_interp,desired_Fs,hmm_Fs,mp_feature_frange);

        %get a broader band MP signal
        mp_features_bb = hsmm_uds_get_lf_features(mp_interp,desired_Fs,desired_Fs,mp_feature_frange_bb);

        % LOCATE THE SEGMENTS CONTAINING UDS
        UDS_segs = hsmm_uds_get_uds_segments(desynch_times,hmm_Fs,length(mp_features),min_seg_dur); %Nx2 matrix containing the index values of the beginning and end of each UDS segment
        is_uds = false(size(mp_features));
        for ii = 1:size(UDS_segs,1)
            is_uds(UDS_segs(ii,1):UDS_segs(ii,2)) = true;
        end
        
        %%
        if sum(is_uds) > 0
            
            mp_dist_xx = linspace(prctile(mp_features_bb(is_uds),0.5),prctile(mp_features_bb(is_uds),99.5),n_dens_bins); %MP amplitude axis for density estimation
            
            %initializations
            [seg_dip,seg_p_value] = deal(nan(size(UDS_segs,1),1));
            mp_middle_crosspts = nan(size(UDS_segs,1),1);
            mp_bandwidths = nan(size(UDS_segs,1),1);
            mp_seg_pdf = nan(size(UDS_segs,1),length(mp_dist_xx));
            for ii = 1:size(UDS_segs,1) %for each UDS segment
                cur_inds = UDS_segs(ii,1):UDS_segs(ii,2);
%                 mp_pdf = hist(mp_features_bb(cur_inds),mp_dist_xx);
%                 mp_pdf([1 end]) = nan;
                [mp_kpdf,~,mp_bandwidths(ii)] = ksdensity(mp_features_bb(cur_inds),mp_dist_xx);
%                 mp_seg_pdf(ii,:) = mp_pdf/nanmean(mp_pdf)*nanmean(mp_kpdf);
                mp_seg_pdf(ii,:) = mp_kpdf;
                
                %check the significance of bimodality
                [seg_dip(ii), seg_p_value(ii), xlow,xup]=hartigansdipsigniftest(mp_features_bb(cur_inds),dip_Nboot);
                
                %if there is bimodality, compute the mean of up and down
                %state amps
                if seg_p_value(ii) < min_dip_pvalue
                    [~,dens_peaks] = findpeaks(mp_kpdf,'npeaks',2,'sortstr','descend','minpeakdistance',ceil(0.25/median(diff(mp_dist_xx))));
                    dens_peaks = mp_dist_xx(dens_peaks);
                    if length(dens_peaks) == 2
                        mp_middle_crosspts(ii) = nanmean(dens_peaks); %set the separation point to be the midpoint between the two modes
                    else
                        mp_middle_crosspts(ii) = nanmedian(mp_features_bb(cur_inds)); %if unimodal, just set this as the median
                    end
                end
            end
            avg_mp_bandwidth = nanmean(mp_bandwidths(seg_p_value < min_dip_pvalue)); %avg density estimation bandwidth across UDS segments with bimodality
            
            %% recompute usable data segments, using only segments with bimodal MP
            all_data(dd).lfp_uds_dur = sum(is_uds)/desired_Fs;
            all_data(dd).tot_dur = length(is_uds)/desired_Fs;
            
            bad_segs = find(seg_p_value > min_dip_pvalue);
            fprintf('%d/%d bad segs\n',length(bad_segs),length(seg_p_value));
            UDS_segs(bad_segs,:) = []; %get rid of the UDS segments without bimodal MP
                                                
            %% if we have some UDS
            if size(UDS_segs,1) > 0
                %fit an HMM to MP
                clear params
                params.meantype = 'variable'; %use 'variable' for time-varying state means. otherwise use 'fixed'
                params.UDS_segs = UDS_segs;
                params.movingwin = [20 5]; %moving window parameters [windowLength windowSlide](in seconds) for computing time-varying state means
                [hmm] = hsmm_uds_initialize(mp_features(:),hmm_Fs,params);
                
                hmm.min_mllik = -8; %set minimum value of the maximal log likelihood for either state before switching to robust estimator.  If you don't want to use this feature, set this to -Inf
                [hmm,gamma] = hsmm_uds_train_hmm(hmm,mp_features(:));
                
                [mp_state_seq,llik_best] = hsmm_uds_viterbi_hmm(hmm,mp_features(:));
                cur_Fs = desired_Fs;
                hmm.min_state_dur = 1;
                hmm.max_state_dur = round(hmm_Fs*20);
                hmm.dur_range = (1:hmm.max_state_dur)/hmm.Fs;
                hmm.state(1).dur_type = 'geometric'; hmm.state(2).dur_type = 'geometric';
                pert_range_bb = [-0.15 0.15];
                mp_state_seq_bb = hsmm_uds_pert_optimize_transitions(mp_features_bb,hmm,mp_state_seq,gamma,hmm_Fs,desired_Fs,pert_range_bb);
                
%                 new_seg_inds = UDS_segs;
                [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,desired_Fs,mp_state_seq_bb);
                [mp_state_durations] = compute_state_durations_seg(mp_state_seq_bb,cur_Fs);
                [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq_bb);
                
                mp_state_vec = nan(size(mp_features_bb));
                is_uds = false(size(mp_features_bb));
%                 mp_state_vec = nan(size(mp_features));
%                 is_uds = false(size(mp_features));
                for i = 1:size(new_seg_inds,1)
                    mp_state_vec(new_seg_inds(i,1):new_seg_inds(i,2)) = mp_state_seq_bb{i};
%                     mp_state_vec(new_seg_inds(i,1):new_seg_inds(i,2)) = mp_state_seq{i};
                    is_uds(new_seg_inds(i,1):new_seg_inds(i,2)) = true;
                end
                
                all_data(dd).uds_dur = sum(is_uds)/cur_Fs;

                %% DETECT LOCAL EXTREMA OF FILTERED LFP SIGNAL
                % compute ctx LFP features
                [lfp_features,t_axis] = hsmm_uds_get_lf_features(ctx_lfp,desired_Fs,cur_Fs,uds_range); %'low-frequency amplitude'
                lfp_features = lfp_features/robust_std_dev(lfp_features(is_uds));
                
                %get local minima and maxima of lfp signal
                lfp_peaks = find(diff(sign(diff(lfp_features))) < 0)+1;
                lfp_valleys = find(diff(sign(diff(lfp_features))) > 0)+1;
                lfp_peaks(~is_uds(lfp_peaks)) = [];
                lfp_valleys(~is_uds(lfp_valleys)) = [];
                
                %get amplitudes of extrema
                lfp_peak_amps = lfp_features(lfp_peaks);
                lfp_valley_amps = -lfp_features(lfp_valleys);
                
                %find the set of down states that might have skipped ctx up
                %blips
                min_down_dur = 2*state_buffer;
                poss_pers_downs = find(mp_state_durations{1} >= min_down_dur);
                
                %find the cortical up-blips skipped in each MP down
                skipped_blips = cell(length(poss_pers_downs),1);
                for ii = 1:length(poss_pers_downs)
                    cur_down_inds = (down_trans_inds(poss_pers_downs(ii)) + round(cur_Fs*state_buffer)):(up_trans_inds(poss_pers_downs(ii)+1) - round(cur_Fs*state_buffer));
                    skipped_blips{ii} = lfp_peak_amps(ismember(lfp_peaks,cur_down_inds));
                end
                
                %calculate stats on skipped ctx up blips of different sizes
                [peak_prob,skip_cnt,skip_prob] = deal(nan(length(peak_checks),1));
                for ii = 1:length(peak_checks)
                    peak_prob(ii) = sum(lfp_peak_amps >= peak_checks(ii))/sum(~isnan(lfp_peak_amps));
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
                    cur_up_inds = (up_trans_inds(poss_pers_ups(ii)) + round(cur_Fs*state_buffer)):(down_trans_inds(poss_pers_ups(ii)) - round(cur_Fs*state_buffer));
                    skipped_blips{ii} = lfp_valley_amps(ismember(lfp_valleys,cur_up_inds));
                end
                
                %compute stats on skipped ctx down blips
                [peak_prob,skip_cnt,skip_prob] = deal(nan(length(peak_checks),1));
                for ii = 1:length(peak_checks)
                    peak_prob(ii) = sum(lfp_valley_amps >= peak_checks(ii))/sum(~isnan(lfp_valley_amps));
                    skip_cnt(ii) = sum(cellfun(@(x) sum(x >= peak_checks(ii)),skipped_blips));
                    skip_prob(ii) = sum(cellfun(@(x) any(x >= peak_checks(ii)),skipped_blips))/length(poss_pers_ups);
                end
                
                all_data(dd).lfp_downpeak_prob = peak_prob;
                all_data(dd).pup_cnt = skip_cnt;
                all_data(dd).pup_prob = skip_prob;
                
                %%
                %temporal window around each LFP extremum to pull data snippets
                back_win = round(0.75*desired_Fs);
                for_win = round(0.75*desired_Fs);
                
                min_n_events = 10; %minimum number of peaks in order to construct trig dists
                
                avg_cross_pt = nanmean(mp_middle_crosspts); %avg definition of state trheshold
                cross_pt_mp = find(mp_dist_xx >= avg_cross_pt,1,'first'); %index of threshold MP value
                
                %compute triggered average MP data
                all_up_tas = nan(length(peak_checks),back_win + for_win+1);
                all_up_dists = nan(length(peak_checks),back_win+for_win+1,length(mp_dist_xx));
                all_down_tas = nan(length(peak_checks),back_win + for_win+1);
                all_down_dists = nan(length(peak_checks),back_win+for_win+1,length(mp_dist_xx));
                
                [all_down_stas,all_up_stas] = deal(nan(length(peak_checks),back_win+for_win+1));
                for jj = 1:length(peak_checks)
                    cur_events = lfp_valleys(lfp_valley_amps >= peak_checks(jj)); %all LFP minima in the current set
                    if length(cur_events) >= min_n_events
                        [cur_ta,cur_lags,~,~,cur_mat] = get_event_trig_avg_v3(mp_features_bb,cur_events,back_win,for_win,2); %trig avg of MP amp
                        all_down_tas(jj,:) = cur_ta;
                        for ii = 1:length(cur_lags)
                            all_down_dists(jj,ii,:) = ksdensity(cur_mat(:,ii),mp_dist_xx,'bandwidth',avg_mp_bandwidth);
                            all_down_dists(jj,ii,:) = jmm_smooth_1d_cor(squeeze(all_down_dists(jj,ii,:)),dens_smth_sig);
                        end
                        
                        [cur_sta] = get_event_trig_avg_v3(mp_state_vec,cur_events,back_win,for_win); %trig avg of MP state
                        all_down_stas(jj,:) = cur_sta;
                    end
                    
                    cur_events = lfp_peaks(lfp_peak_amps >= peak_checks(jj));
                    if length(cur_events) >= min_n_events
                        [cur_ta,cur_lags,~,~,cur_mat] = get_event_trig_avg_v3(mp_features_bb,cur_events,back_win,for_win,2);
                        all_up_tas(jj,:) = cur_ta;
                        for ii = 1:length(cur_lags)
                            all_up_dists(jj,ii,:) = ksdensity(cur_mat(:,ii),mp_dist_xx,'bandwidth',avg_mp_bandwidth);
                            all_up_dists(jj,ii,:) = jmm_smooth_1d_cor(squeeze(all_up_dists(jj,ii,:)),dens_smth_sig);
                        end
                        
                        [cur_sta] = get_event_trig_avg_v3(mp_state_vec,cur_events,back_win,for_win);
                        all_up_stas(jj,:) = cur_sta;
                    end
                end
                
                utrig_Pmp_down = trapz(mp_dist_xx(1:(cross_pt_mp-1)),all_up_dists(:,:,1:(cross_pt_mp-1)),3);
                utrig_Pmp_up = trapz(mp_dist_xx((cross_pt_mp+1):end),all_up_dists(:,:,(cross_pt_mp+1):end),3);
                dtrig_Pmp_down = trapz(mp_dist_xx(1:(cross_pt_mp-1)),all_down_dists(:,:,1:(cross_pt_mp-1)),3);
                dtrig_Pmp_up = trapz(mp_dist_xx((cross_pt_mp+1):end),all_down_dists(:,:,(cross_pt_mp+1):end),3);
                
                all_data(dd).utrig_Pmp_down = utrig_Pmp_down;
                all_data(dd).utrig_Pmp_up = utrig_Pmp_up;
                all_data(dd).dtrig_Pmp_down = dtrig_Pmp_down;
                all_data(dd).dtrig_Pmp_up = dtrig_Pmp_up;
                all_data(dd).utrig_mpstate = all_up_stas;
                all_data(dd).dtrig_mpstate = all_down_stas;
                all_data(dd).utrig_mpavg = all_up_tas;
                all_data(dd).dtrig_mpavg = all_down_tas;
                
                all_data(dd).avg_stateprob = nanmean(mp_state_vec(is_uds));
            end
        end
        
    end
end

%%
min_uds_dur = 50;
poss_use = find(arrayfun(@(x) length(x.uds_dur),all_data));
poss_use([all_data(poss_use).uds_dur] == 0) = [];
use_uds = poss_use([all_data(poss_use).uds_dur] >= min_uds_dur);
% use_mp = poss_use([all_data(poss_use).dipp] < 0.05);
% use_cells = intersect(use_mp,use_uds);
use_cells = use_uds;

mec_neurons = find(arrayfun(@(x) strcmp(x.MPloc,'MEC'),data));
ctx_neurons = find(arrayfun(@(x) strcmp(x.MPloc,'Ctx'),data));

use_mec_neurons = intersect(mec_neurons,use_cells);
use_ctx_neurons = intersect(ctx_neurons,use_cells);

%%
mec_pdown_prob = cell2mat(arrayfun(@(x) x.pdown_prob,all_data(use_mec_neurons),'uniformoutput',0));
ctx_pdown_prob = cell2mat(arrayfun(@(x) x.pdown_prob,all_data(use_ctx_neurons),'uniformoutput',0));

mec_pup_prob = cell2mat(arrayfun(@(x) x.pup_prob,all_data(use_mec_neurons),'uniformoutput',0));
ctx_pup_prob = cell2mat(arrayfun(@(x) x.pup_prob,all_data(use_ctx_neurons),'uniformoutput',0));

f1 = figure();
subplot(2,1,1); hold on
plot(peak_checks,mec_pdown_prob,'ko-');
plot(peak_checks,ctx_pdown_prob,'ro-');
ylim([0 1]);
subplot(2,1,2); hold on
plot(peak_checks,mec_pup_prob,'ko-');
plot(peak_checks,ctx_pup_prob,'ro-');
ylim([0 1]);

%%
clear *_up_ts *_down_ts
for ii = 1:length(use_mec_neurons)
    mec_up_ts(ii,:,:) = all_data(use_mec_neurons(ii)).utrig_mpstate;
    mec_down_ts(ii,:,:) = all_data(use_mec_neurons(ii)).dtrig_mpstate;
end
for ii = 1:length(use_ctx_neurons)
    ctx_up_ts(ii,:,:) = all_data(use_ctx_neurons(ii)).utrig_mpstate;
    ctx_down_ts(ii,:,:) = all_data(use_ctx_neurons(ii)).dtrig_mpstate;
end

mec_up_ts = bsxfun(@minus,mec_up_ts,[all_data(use_mec_neurons).avg_stateprob]');
ctx_up_ts = bsxfun(@minus,ctx_up_ts,[all_data(use_ctx_neurons).avg_stateprob]');
mec_down_ts = bsxfun(@minus,mec_down_ts,[all_data(use_mec_neurons).avg_stateprob]');
ctx_down_ts = bsxfun(@minus,ctx_down_ts,[all_data(use_ctx_neurons).avg_stateprob]');

f1 = figure();
for ii = 1:length(peak_checks)
    subplot(3,2,ii); hold on
    plot(cur_lags,squeeze(mec_up_ts(:,ii,:)),'k-');
    plot(cur_lags,squeeze(ctx_up_ts(:,ii,:)),'r-');    
    plot(cur_lags,squeeze(nanmean(mec_up_ts(:,ii,:))),'k','linewidth',3);
    plot(cur_lags,squeeze(nanmean(ctx_up_ts(:,ii,:))),'r','linewidth',3);
end

f1 = figure();
for ii = 1:length(peak_checks)
    subplot(3,2,ii); hold on
    plot(cur_lags,squeeze(mec_down_ts(:,ii,:)),'k-');
    plot(cur_lags,squeeze(ctx_down_ts(:,ii,:)),'r-');    
    plot(cur_lags,squeeze(nanmean(mec_down_ts(:,ii,:))),'k','linewidth',3);
    plot(cur_lags,squeeze(nanmean(ctx_down_ts(:,ii,:))),'r','linewidth',3);
end

%%
clear *_up_ta *_down_ta
for ii = 1:length(use_mec_neurons)
    mec_up_ta(ii,:,:) = all_data(use_mec_neurons(ii)).utrig_mpavg;
    mec_down_ta(ii,:,:) = all_data(use_mec_neurons(ii)).dtrig_mpavg;
end
for ii = 1:length(use_ctx_neurons)
    ctx_up_ta(ii,:,:) = all_data(use_ctx_neurons(ii)).utrig_mpavg;
    ctx_down_ta(ii,:,:) = all_data(use_ctx_neurons(ii)).dtrig_mpavg;
end

f1 = figure();
for ii = 1:length(peak_checks)
    subplot(3,2,ii); hold on
    plot(cur_lags,squeeze(mec_up_ta(:,ii,:)),'k-');
    plot(cur_lags,squeeze(ctx_up_ta(:,ii,:)),'r-');    
    plot(cur_lags,squeeze(nanmean(mec_up_ta(:,ii,:))),'k','linewidth',3);
    plot(cur_lags,squeeze(nanmean(ctx_up_ta(:,ii,:))),'r','linewidth',3);
end

f1 = figure();
for ii = 1:length(peak_checks)
    subplot(3,2,ii); hold on
    plot(cur_lags,squeeze(mec_down_ta(:,ii,:)),'k-');
    plot(cur_lags,squeeze(ctx_down_ta(:,ii,:)),'r-');    
    plot(cur_lags,squeeze(nanmean(mec_down_ta(:,ii,:))),'k','linewidth',3);
    plot(cur_lags,squeeze(nanmean(ctx_down_ta(:,ii,:))),'r','linewidth',3);
end