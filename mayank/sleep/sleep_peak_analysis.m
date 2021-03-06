clear all
% close all
cd ~/Analysis/Mayank/sleep/
load sleep_dirs

fig_dir = '/Users/james/Analysis/Mayank/sleep/sleep_figs2/';
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))

%% DEFINE FREQUENCY RANGES
uds_range = [0.5 3];
hf_range = [30 80];
lf_range = [0 0.2];
use_range = [0.5 80];
dens_smth_sig = 3;

mp_lcf = 0.01; %lcf for hp-filter on nlx mp signal
lfp_lcf = 0.2; %lcf for hp-filter on LFP signals

mua_smooth_sigma = 0.025; %smoothing sigma applied to MUA

to_print = true;

%for normalized LFP power
% lfp_uds_thresh = -2.5;
lfp_uds_thresh = -0.5;
% lfp_lf_max = 0.5;
lfp_lf_max = 2;

%%
% for dd = 29:length(data);
for dd = 14
    
    %%
    cd(data(dd).dir)
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
    
    %%
    if has_ipsi %if there are any ipsilateral LFPs
        
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
        mp_Fsd = mp_Fs/mp_dsf;
        [b,a] = butter(2,mp_lcf/(mp_Fsd/2),'high');
        mp_d = decimate(mp_d,mp_dsf);
        mp_d = (filtfilt(b,a,mp_d)); %filter out very slow drift in MP signal
        mp_t = downsample(mp_t,mp_dsf);
        
        if ~isempty(heka_data)
            heka_Fs = 1/nanmedian(diff(heka_time));
            heka_data = decimate(heka_data,heka_dsf);
            heka_Fsd = heka_Fs/heka_dsf;
            [b,a] = butter(2,mp_lcf/(heka_Fsd/2),'high');
            heka_data = filtfilt(b,a,heka_data);
            heka_time = downsample(heka_time,heka_dsf);
        end
        
        %% find usable MP times
        use_mp_inds = find(mp_t >= data(dd).good_bounds(1) & mp_t <= data(dd).good_bounds(2));
        use_heka_inds = find(heka_time >= data(dd).good_bounds(1) & heka_time <= data(dd).good_bounds(2));
        
        mp_t = mp_t(use_mp_inds); mp_d = mp_d(use_mp_inds);
        heka_time = heka_time(use_heka_inds); heka_data = heka_data(use_heka_inds);
        
        %%
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
            ipsi_csc{ii} = ipsi_csc{ii}/robust_std_dev(ipsi_csc_hp{ii});
            ipsi_csc_hp{ii} = ipsi_csc_hp{ii}/robust_std_dev(ipsi_csc_hp{ii});
        end
        
        %% process MUA rates
        mua_sm_sig = round(mua_smooth_sigma*csc_Fsd);
        ipsi_binned_spks = nan(length(csc_time),length(ipsi_mua_times));
        ipsi_sm_rate = nan(length(csc_time),length(ipsi_mua_times));
        for ii = 1:length(ipsi_mua_times)
            ipsi_binned_spks(:,ii) = hist(ipsi_mua_times{ii},csc_time);
            ipsi_sm_rate(:,ii) = zscore(jmm_smooth_1d_cor(ipsi_binned_spks(:,ii),mua_sm_sig));
        end
        
        %% Get MP signal interpolated onto CSC time axis
        uu = find(~isnan(mp_t));
        mp_interp = interp1(mp_t(uu),mp_d(uu),csc_time);
        
        if ~isempty(heka_time)
            fprintf('Using heka MP data\n');
            
            uu = find(~isnan(heka_time));
            mp_interp = interp1(heka_time(uu),heka_data(uu),csc_time);
        end
        
        %handle any interpolation boundary issues
        bad_mp_inds = find(isnan(mp_interp));
        mp_interp(bad_mp_inds) = 0;
        
        %% compute LFP and MP spectrograms
        params.Fs = csc_Fsd;
        params.tapers = [2 3];
        params.fpass = [0 80];
        movingwin = [10 1];
        clear C
        
        if sum(isnan(mp_interp)) > 0.01*length(mp_interp)
            warning('More than 1% NAN components!');
        end
        mp_interp(isnan(mp_interp)) = 0;
        
        [b,a] = butter(2,[0.05]/(csc_Fs/2),'high');
        % mp_interp = filtfilt(b,a,mp_interp);
        
        use_cscs = 1:1:length(data(dd).ipsiLFPs);
        
        [S_mp,t,f] = mtspecgramc(mp_interp(:),movingwin,params);
        
        clear S_lfp
        for cc = 1:length(use_cscs)
            cur_lfp = filtfilt(b,a,ipsi_csc{use_cscs(cc)});
            %     cur_lfp = ipsi_csc{use_cscs(cc)};
            %             [C{cc},phi,S12,S1,S2{cc},t,f]=cohgramc(mp_interp(:),cur_lfp(:),movingwin,params);
            [S_lfp{cc}]=mtspecgramc(cur_lfp(:),movingwin,params);
        end
        
        %         %high-resolution specgram
        %         mp_mwin = [2 1];
        %         [S_mp,mp_t,mp_f]=mtspecgramc(mp_interp(:),mp_mwin,params);
        
        uds_freqs = find(f >= uds_range(1) & f <= uds_range(2));
        hf_freqs = find(f >= hf_range(1) & f <= hf_range(2));
        lf_freqs = find(f >= lf_range(1) & f <= lf_range(2));
        use_freqs = find(f >= use_range(1) & f <= use_range(2));
        %         mp_use_freqs = find(mp_f >= use_range(1) & mp_f <= use_range(2));
        
        %get MP power in each band
        mp_uds_pow = trapz(f(uds_freqs),(S_mp(:,uds_freqs)),2);
        mp_hf_pow = trapz(f(hf_freqs),(S_mp(:,hf_freqs)),2);
        mp_lf_pow = trapz(f(lf_freqs),(S_mp(:,lf_freqs)),2);
                
        %get LFP power in each band
        lfp_uds_pow = nan(length(use_cscs),length(t));
        lfp_hf_pow = nan(length(use_cscs),length(t));
        lfp_lf_pow = nan(length(use_cscs),length(t));
        lfp_tot_pow = nan(length(use_cscs),length(t));
        for cc = 1:length(use_cscs)
            lfp_uds_pow(cc,:) = trapz(f(uds_freqs),(S_lfp{cc}(:,uds_freqs)),2);
            lfp_hf_pow(cc,:) = trapz(f(hf_freqs),(S_lfp{cc}(:,hf_freqs)),2);
            lfp_lf_pow(cc,:) = trapz(f(lf_freqs),(S_lfp{cc}(:,lf_freqs)),2);
            lfp_tot_pow(cc,:) = trapz(f(use_freqs),S_lfp{cc}(:,use_freqs),2);
        end
        
        %%
        %average log-powers across time in each freq range (the log helps
        %normalize the var
        avg_uds_pow = nanmean(log10(lfp_uds_pow),2);
        avg_lf_pow = nanmean(log10(lfp_lf_pow),2);
        avg_hf_pow = nanmean(log10(lfp_hf_pow),2);
        
        %find channel with maximum UDS/HF power ratio (difference in
        %log-space)
        %         [~,ctx_ch] = max(avg_uds_pow(poss_ctx_chs)-0*avg_lf_pow(poss_ctx_chs));
        [~,ctx_ch] = max(avg_uds_pow(poss_ctx_chs) - avg_lf_pow(poss_ctx_chs));
        ctx_ch = poss_ctx_chs(ctx_ch);
        
%         ctx_uds_hf_ratio = log10(lfp_uds_pow(ctx_ch,:)) - log10(lfp_hf_pow(ctx_ch,:));
        
        uds_epochs = log10(lfp_uds_pow(ctx_ch,:)) >= lfp_uds_thresh;
        lf_epochs = log10(lfp_lf_pow(ctx_ch,:)) >= lfp_lf_max;
        
        all_data(dd).lfp_pow_spec = mean(S_lfp{ctx_ch});
        all_data(dd).mp_pow_spec = mean(S_mp);
        
        % has_uds = true(length(csc_time),1);
        % for ii = 1:length(t)
        %    if ~ismember(ii,uds_epochs)
        %       curset = find(csc_time >= t(ii)-movingwin(1)/2 & csc_time <= t(ii) + movingwin(1)/2);
        %       has_uds(curset) = false;
        %    end
        % end
        
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
        addpath('~/James_scripts/hsmm_uds_toolbox/')
        removal_window = 0; %buffer window around desynch epochs
        
        %mark points where the SO power crosses below threshold
        %         desynch_indicator = ~uds_epochs | lf_epochs | mp_tot_pow' < -3;
        desynch_indicator = ~uds_epochs | lf_epochs;
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
        hmm_Fs = mp_Fsd/hmm_dsf; %sample-frequency of HMM signal features (Hz)
        lf_cut_off_freqs = [0.05 4]; %for "low-frequency amplitude" features, these are the cutoff-frequencies for the bandpass filter
        [bb,aa] = butter(2,lf_cut_off_freqs/(mp_Fsd/2));
        mp_features = filtfilt(bb,aa,mp_interp);
        mp_features = downsample(mp_features,hmm_dsf);
        
                        %% LOCATE THE SEGMENTS CONTAINING UDS
        T = length(mp_features); %number of samples
        min_seg_dur = 20; %minimum duration of UDS segments used for analysis (in sec)
        UDS_segs = hsmm_uds_get_uds_segments(desynch_times,hmm_Fs,T,min_seg_dur); %Nx2 matrix containing the index values of the beginning and end of each UDS segment
        
        is_uds = false(size(mp_features));
        for ii = 1:size(UDS_segs,1)
            is_uds(UDS_segs(ii,1):UDS_segs(ii,2)) = true;
        end
        
        all_data(dd).uds_dur = sum(is_uds)/csc_Fsd;
        all_data(dd).tot_dur = length(is_uds)/csc_Fsd;

        %%
        clear params
        params.meantype = 'variable'; %use 'variable' for time-varying state means. otherwise use 'fixed'
        params.UDS_segs = UDS_segs;
        params.movingwin = [15 5]; %moving window parameters [windowLength windowSlide](in seconds) for computing time-varying state means
        [hmm] = hsmm_uds_initialize(mp_features(:),hmm_Fs,params);
        
        % FIT AN HMM
        hmm.min_mllik = -8; %set minimum value of the maximal log likelihood for either state before switching to robust estimator.  If you don't want to use this feature, set this to -Inf
        hmm = hsmm_uds_train_hmm(hmm,mp_features(:));
        
        [mp_state_seq,llik_best] = hsmm_uds_viterbi_hmm(hmm,mp_features(:));
        mp_state_vec = nan(size(mp_features));
        for i = 1:size(UDS_segs,1)
            mp_state_vec(UDS_segs(i,1):UDS_segs(i,2)) = mp_state_seq{i};
        end
   
        [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,hmm_Fs,mp_state_seq);
        [mp_state_durations] = compute_state_durations_seg(mp_state_seq,hmm_Fs);
        [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
%%
        
        lf_cut_off_freqs = [0.5 3]; %for "low-frequency amplitude" features, these are the cutoff-frequencies for the bandpass filter
        [lfp_features,t_axis] = hsmm_uds_get_lf_features(ipsi_csc{ctx_ch},csc_Fsd,hmm_Fs,lf_cut_off_freqs); %'low-frequency amplitude'
        lfp_features = lfp_features/robust_std_dev(lfp_features);
        
        lf_cut_off_freqs = [0.5 20]; %for "low-frequency amplitude" features, these are the cutoff-frequencies for the bandpass filter
        [lfp_features_bb,t_axis] = hsmm_uds_get_lf_features(ipsi_csc{ctx_ch},csc_Fsd,hmm_Fs,lf_cut_off_freqs); %'low-frequency amplitude'
        lfp_features_bb = lfp_features_bb/robust_std_dev(lfp_features_bb);
        
        
        %%
        addpath('~/Analysis/Mayank/sleep/');
        %         xx = linspace(-3.5,3.5,100);
        mp_dist_xx = linspace(prctile(mp_features(is_uds),0.5),prctile(mp_features(is_uds),99.5),100);
        if sum(is_uds) > 0
            mp_pdf = hist(mp_features(is_uds),mp_dist_xx);
            mp_pdf([1 end]) = nan;
            [mp_kpdf,~,mp_bandwidth] = ksdensity(mp_features(is_uds),mp_dist_xx);
            mp_pdf = mp_pdf/nanmean(mp_pdf)*nanmean(mp_kpdf);
            
            [~,dens_peaks] = findpeaks(mp_kpdf,'npeaks',2,'sortstr','descend','minpeakdistance',ceil(0.25/median(diff(mp_dist_xx))));
            dens_peaks = mp_dist_xx(dens_peaks);
            if length(dens_peaks) == 2
                cross_pt = nanmean(dens_peaks); %set the separation point to be the midpoint between the two modes
            else
                cross_pt = nanmedian(mp_features(is_uds)); %if unimodal, just set this as the median
            end
            
            mp_ovPup = sum(mp_features(is_uds) > cross_pt)/sum(is_uds);
            mp_ovPdown = sum(mp_features(is_uds) < cross_pt)/sum(is_uds);
            all_data(dd).mp_ov_upProb = mp_ovPup;
            all_data(dd).mp_ov_downProb = mp_ovPdown;
            
            % gmfit = gmdistribution.fit(mp_features(is_uds),2,'SharedCov',true);
            % init_Sigma(:,:,1) = gmfit.Sigma;
            % init_Sigma(:,:,2) = gmfit.Sigma;
            % gminit = struct('mu',gmfit.mu,'Sigma',init_Sigma,'PComponents',gmfit.PComponents);
            % % gminit(2) = struct('mu',gmfit.mu(2),'Sigma',gmfit.Sigma,'PComponents',gmfit.PComponents(2));
            % gmfit = gmdistribution.fit(mp_features(is_uds),2,'Start',gminit);
            % g1 = gmfit.PComponents(1)*normpdf(xx,gmfit.mu(1),gmfit.Sigma(1));
            % g2 = gmfit.PComponents(2)*normpdf(xx,gmfit.mu(2),gmfit.Sigma(2));
            % % g2 = gmfit.PComponents(2)*normpdf(xx,gmfit.mu(2),gmfit.Sigma(1));
            
            % gm_posterior = posterior(gmfit,xx');
            % p1 = find(xx >= gmfit.mu(1),1); p2 = find(xx >= gmfit.mu(2),1);
            % peak_locs = [p1; p2]; peak_locs = sort(peak_locs);
            % between_peaks = peak_locs(1):peak_locs(2);
            % [~,cross_pt] = min(abs(gm_posterior(between_peaks,1) - 0.5));
            % cross_pt = xx(between_peaks(cross_pt));
            
            if sum(is_uds) > 0
                uu = find(is_uds & ~isnan(mp_features));
                [dip, p_value, xlow,xup]=hartigansdipsigniftest(mp_features(uu),100);
            else
                dip = nan; p_value = nan;
            end
            all_data(dd).dipstate = dip;
            all_data(dd).dipp = p_value;
            
            
            if to_print
                f1 = figure(); hold on
                plot(mp_dist_xx,mp_pdf,mp_dist_xx,mp_kpdf,'r');
                % plot(xx,g1,'k',xx,g2,'g');
                yl = ylim();
                line(cross_pt + [0 0],yl,'color','k');
                title(sprintf('MP Bimodality test: DipStat:%.3f  P:%.3g',dip,p_value));
                
%                 dname = find(data(dd).dir == '/',1,'last');
%                 dname = data(dd).dir(dname+1:end);
%                 fname = strcat(fig_dir,'MP_dist_',data(dd).MPloc,dname);
%                 figufy(f1);
%                 set(f1,'PaperUnits','inches','PaperSize',[6 5],'PaperPosition',[0 0 6 5]);
%                 print(f1,'-dpng',fname);
%                 close(f1);
            end
        end
        %% DETECT LOCAL EXTREMA OF FILTERED LFP SIGNAL
        lfp_peaks = find(diff(sign(diff(lfp_features))) < 0)+1;
        lfp_valleys = find(diff(sign(diff(lfp_features))) > 0)+1;
        
        lfp_peaks(~is_uds(lfp_peaks)) = [];
        lfp_valleys(~is_uds(lfp_valleys)) = [];
        
        lfp_peak_amps = lfp_features(lfp_peaks);
        lfp_valley_amps = -lfp_features(lfp_valleys);
        
        % poss_amps = [0 0.5 1:3];
        poss_amp_prc = [0 25 50 75 90];
        poss_amps_peaks = prctile(lfp_peak_amps,poss_amp_prc);
        poss_amps_valleys = prctile(lfp_valley_amps,poss_amp_prc);
        
        %temporal window around each LFP extremum to pull data snippets
        back_win = round(0.75*csc_Fsd);
        for_win = round(0.75*csc_Fsd);
        
        min_n_events = 10; %minimum number of peaks in order to construct trig dists
        %         lfp_bandwidth = 0.15; %smoothing bandwidth for LFP
        %         mp_bandwidth = 0.3; %smoothing bandwidth for MP
        
        %         mp_bndwidth = diff(prctile(mp_features(is_uds),[25 75]))/15;
        
        %         mp_dist_xx = linspace(-3,3,100);
        % gm_posterior = posterior(gmfit,xx');
        % p1 = find(xx >= gmfit.mu(1),1); p2 = find(xx >= gmfit.mu(2),1);
        % peak_locs = [p1; p2]; peak_locs = sort(peak_locs);
        % between_peaks = peak_locs(1):peak_locs(2);
        % [~,cross_pt] = min(abs(gm_posterior(between_peaks,1) - 0.5));
        % cross_pt = between_peaks(cross_pt);
        cross_pt_mp = find(mp_dist_xx >= cross_pt,1,'first'); %index of threshold MP value
        
        %compute triggered average MP data
        all_up_tas = nan(length(poss_amp_prc),back_win + for_win+1);
        all_up_dists = nan(length(poss_amp_prc),back_win+for_win+1,length(mp_dist_xx));
        all_down_tas = nan(length(poss_amp_prc),back_win + for_win+1);
        all_down_dists = nan(length(poss_amp_prc),back_win+for_win+1,length(mp_dist_xx));
        for jj = 1:length(poss_amp_prc)
            %     cur_events = lfp_valleys(lfp_features(lfp_valleys) < -poss_amps(jj));
            %     cur_events(~is_uds(cur_events)) = [];
            cur_events = lfp_valleys(lfp_valley_amps >= poss_amps_valleys(jj)); %all LFP minima in the current set
            if length(cur_events) >= min_n_events
                [cur_ta,cur_lags,~,~,cur_mat] = get_event_trig_avg_v3(mp_features,cur_events,back_win,for_win,2);
                all_down_tas(jj,:) = cur_ta;
                
                for ii = 1:length(cur_lags)
                    all_down_dists(jj,ii,:) = ksdensity(cur_mat(:,ii),mp_dist_xx,'bandwidth',mp_bandwidth);
                    all_down_dists(jj,ii,:) = jmm_smooth_1d_cor(squeeze(all_down_dists(jj,ii,:)),dens_smth_sig);
                end
            end
            
            %     cur_events = lfp_peaks(lfp_features(lfp_peaks) > poss_amps(jj));
            %     cur_events(~is_uds(cur_events)) = [];
            cur_events = lfp_peaks(lfp_peak_amps >= poss_amps_peaks(jj));
            if length(cur_events) >= min_n_events
                [cur_ta,cur_lags,~,~,cur_mat] = get_event_trig_avg_v3(mp_features,cur_events,back_win,for_win,2);
                all_up_tas(jj,:) = cur_ta;
                
                for ii = 1:length(cur_lags)
                    all_up_dists(jj,ii,:) = ksdensity(cur_mat(:,ii),mp_dist_xx,'bandwidth',mp_bandwidth);
                    all_up_dists(jj,ii,:) = jmm_smooth_1d_cor(squeeze(all_up_dists(jj,ii,:)),dens_smth_sig);
                end
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
        
        %%
        state_buffer = 0.25;
        min_down_dur = 2*state_buffer;
        poss_pers_downs = find(mp_state_durations{1} >= min_down_dur);
        
        skipped_blips = cell(length(poss_pers_downs),1);
        for ii = 1:length(poss_pers_downs)
           cur_down_inds = (down_trans_inds(poss_pers_downs(ii)) + round(hmm_Fs*state_buffer)):(up_trans_inds(poss_pers_downs(ii)+1) - round(hmm_Fs*state_buffer));
           skipped_blips{ii} = lfp_peak_amps(ismember(lfp_peaks,cur_down_inds));
        end
        
        peak_checks = [0.5 1 1.5 2 3 4];
        for ii = 1:length(peak_checks)
           peak_prob(ii) = sum(lfp_peak_amps >= peak_checks(ii))/sum(~isnan(lfp_peak_amps));
           skip_cnt(ii) = sum(cellfun(@(x) sum(x >= peak_checks(ii)),skipped_blips));
           skip_prob(ii) = sum(cellfun(@(x) any(x >= peak_checks(ii)),skipped_blips))/length(poss_pers_downs);
        end
        
        
        %%
        lfp_dist_xx = linspace(-5,3,100);
        if sum(is_uds) > 0
            [~,~,lfp_bandwidth] = ksdensity(lfp_features_bb(is_uds),lfp_dist_xx);
            all_up_tas_lfp = nan(length(poss_amp_prc),back_win + for_win+1);
            all_up_dists_lfp = nan(length(poss_amp_prc),back_win+for_win+1,length(mp_dist_xx));
            all_down_tas_lfp = nan(length(poss_amp_prc),back_win + for_win+1);
            all_down_dists_lfp = nan(length(poss_amp_prc),back_win+for_win+1,length(mp_dist_xx));
            for jj = 1:length(poss_amp_prc)
                %     cur_events = lfp_valleys(lfp_features(lfp_valleys) < -poss_amps(jj));
                %     cur_events(~is_uds(cur_events)) = [];
                cur_events = lfp_valleys(lfp_valley_amps >= poss_amps_valleys(jj));
                if length(cur_events) >= min_n_events
                    [cur_ta,cur_lags,~,~,cur_mat] = get_event_trig_avg_v3(lfp_features_bb,cur_events,back_win,for_win,2);
                    all_down_tas_lfp(jj,:) = cur_ta;
                    
                    for ii = 1:length(cur_lags)
                        all_down_dists_lfp(jj,ii,:) = ksdensity(cur_mat(:,ii),lfp_dist_xx,'bandwidth',lfp_bandwidth);
                        all_down_dists_lfp(jj,ii,:) = jmm_smooth_1d_cor(squeeze(all_down_dists_lfp(jj,ii,:)),dens_smth_sig);
                    end
                end
                
                %     cur_events = lfp_peaks(lfp_features(lfp_peaks) > poss_amps(jj));
                %     cur_events(~is_uds(cur_events)) = [];
                cur_events = lfp_peaks(lfp_peak_amps >= poss_amps_peaks(jj));
                if length(cur_events) >= min_n_events
                    [cur_ta,cur_lags,~,~,cur_mat] = get_event_trig_avg_v3(lfp_features_bb,cur_events,back_win,for_win,2);
                    all_up_tas_lfp(jj,:) = cur_ta;
                    
                    for ii = 1:length(cur_lags)
                        all_up_dists_lfp(jj,ii,:) = ksdensity(cur_mat(:,ii),lfp_dist_xx,'bandwidth',lfp_bandwidth);
                        all_up_dists_lfp(jj,ii,:) = jmm_smooth_1d_cor(squeeze(all_up_dists_lfp(jj,ii,:)),dens_smth_sig);
                    end
                end
                
            end
        end
        
        %%
        if sum(is_uds) > 0
            max_utrig_Pdown = max(utrig_Pmp_down(:));
            max_dtrig_Pup = max(dtrig_Pmp_up(:));
            if isnan(max_utrig_Pdown) || max_utrig_Pdown == 0
                max_utrig_Pdown = 0.1;
            end
            if isnan(max_dtrig_Pup) || max_dtrig_Pup == 0
                max_dtrig_Pup = 0.1;
            end
            lfp_dx = median(diff(lfp_dist_xx));
            mp_dx = median(diff(mp_dist_xx));
            
            if to_print
                f1 = figure();
                for ii = 1:length(poss_amp_prc)
                    subplot(5,3,3*(ii-1)+1)
                    imagesc(cur_lags/csc_Fsd,lfp_dist_xx,log(squeeze(all_up_dists_lfp(ii,:,:)*lfp_dx))');
                    caxis([-3 -0.5]-2.5);
                    set(gca,'ydir','normal');
                    xlabel('Lag (s)');
                    ylabel('Amplitude (z)');
                    title(sprintf('LFP %d pcnt',poss_amp_prc(ii)));
                    
                    subplot(5,3,3*(ii-1)+2)
                    imagesc(cur_lags/csc_Fsd,mp_dist_xx,log(squeeze(all_up_dists(ii,:,:)*mp_dx))'); hold on
                    xl = xlim(); line(xl,mp_dist_xx(cross_pt_mp) + [0 0],'color','w');
                    caxis([-5 -1]-2);
                    set(gca,'ydir','normal');
                    xlabel('Lag (s)');
                    ylabel('Amplitude (z)');
                    title(sprintf('MP %d pcnt',poss_amp_prc(ii)));
                    
                    subplot(5,3,3*(ii-1)+3); hold on
                    plot(cur_lags/csc_Fsd,squeeze(utrig_Pmp_down(ii,:)));
                    % plot(cur_lags/csc_Fsd,squeeze(utrig_Pmp_up(ii,:)),'r');
                    xlabel('Lag (s)');
                    ylabel('Probability');
                    xlim(cur_lags([1 end])/csc_Fsd);
                    ylim([0 max_utrig_Pdown])
                    xlim(xl); line(xl,mp_ovPdown+[0 0],'color','k');
                    title(sprintf('MP %d pcnt',poss_amp_prc(ii)));
                    
                end
                
                
                f2 = figure();
                for ii = 1:length(poss_amp_prc)
                    subplot(5,3,3*(ii-1)+1)
                    imagesc(cur_lags/csc_Fsd,lfp_dist_xx,log(squeeze(all_down_dists_lfp(ii,:,:)*lfp_dx))');
                    caxis([-3 -0.5]-2.5);
                    set(gca,'ydir','normal');
                    xlabel('Lag (s)');
                    ylabel('Amplitude (z)');
                    title(sprintf('LFP %d pcnt',poss_amp_prc(ii)));
                    
                    subplot(5,3,3*(ii-1)+2)
                    imagesc(cur_lags/csc_Fsd,mp_dist_xx,log(squeeze(all_down_dists(ii,:,:)*mp_dx))');
                    xl = xlim(); line(xl,mp_dist_xx(cross_pt_mp) + [0 0],'color','w');
                    caxis([-5 -1]-2);
                    set(gca,'ydir','normal');
                    xlabel('Lag (s)');
                    ylabel('Amplitude (z)');
                    title(sprintf('MP %d pcnt',poss_amp_prc(ii)));
                    
                    subplot(5,3,3*(ii-1)+3); hold on
                    plot(cur_lags/csc_Fsd,squeeze(dtrig_Pmp_up(ii,:)));
                    % plot(cur_lags/csc_Fsd,squeeze(dtrig_Pmp_down(ii,:)),'r');
                    xlabel('Lag (s)');
                    ylabel('Probability');
                    xlim(cur_lags([1 end])/csc_Fsd);
                    ylim([0 max_dtrig_Pup])
                    xlim(xl); line(xl,mp_ovPup+[0 0],'color','k');
                    title(sprintf('MP %d pcnt',poss_amp_prc(ii)));
                    
                end
                
                
                %%
                dname = find(data(dd).dir == '/',1,'last');
                dname = data(dd).dir(dname+1:end);
                fname = strcat(fig_dir,'NUPTRIG_',data(dd).MPloc,dname);
                figufy(f1);
                set(f1,'PaperUnits','inches','PaperSize',[10 14],'PaperPosition',[0 0 10 14]);
                print(f1,'-dpng',fname);
                close(f1);
                
                
                dname = find(data(dd).dir == '/',1,'last');
                dname = data(dd).dir(dname+1:end);
                fname = strcat(fig_dir,'NDOWNTRIG_',data(dd).MPloc,dname);
                figufy(f2);
                set(f2,'PaperUnits','inches','PaperSize',[10 14],'PaperPosition',[0 0 10 14]);
                print(f2,'-dpng',fname);
                close(f2);
            end
        end
    end
end