clear all
close all
cd ~/Analysis/Mayank/sleep/
% load sleep_dirs
load sleep_dirs_old
fig_dir = '/Users/james/Analysis/Mayank/sleep/sleep_figs2/';

addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))

%%
for dd = 1:length(data);
    
    %%
    cd(data(dd).dir)
    load('procData.mat','mp_*','ipsi_*','heka*','csc_time');
    
    has_ipsi = true;
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
    if has_ipsi
        % ipsi_csc = contra_csc;
        
        %%
        
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
        [b,a] = butter(2,[0.01]/(mp_Fsd/2),'high');
        mp_d = decimate(mp_d,mp_dsf);
        mp_d = (filtfilt(b,a,mp_d));
        mp_t = downsample(mp_t,mp_dsf);
        
        if ~isempty(heka_data)
            heka_Fs = 1/nanmedian(diff(heka_time));
            heka_data = decimate(heka_data,heka_dsf);
            heka_Fsd = heka_Fs/heka_dsf;
            [b,a] = butter(2,0.05/(heka_Fsd/2),'high');
            heka_data = filtfilt(b,a,heka_data);
            heka_time = downsample(heka_time,heka_dsf);
        end
        
        %%
        use_mp_inds = find(mp_t >= data(dd).good_bounds(1) & mp_t <= data(dd).good_bounds(2));
        use_heka_inds = find(heka_time >= data(dd).good_bounds(1) & heka_time <= data(dd).good_bounds(2));
        
        mp_t = mp_t(use_mp_inds); mp_d = mp_d(use_mp_inds);
        heka_time = heka_time(use_heka_inds); heka_data = heka_data(use_heka_inds);
        
        %%
        csc_time = downsample(csc_time,csc_dsf);
        csc_inds = find(csc_time >= data(dd).good_bounds(1) & csc_time <= data(dd).good_bounds(2));
        csc_time = csc_time(csc_inds);
        
        csc_Fsd = csc_Fs/csc_dsf;
        [bb,aa] = butter(4,0.2/(csc_Fsd/2),'high');
        for ii = 1:length(ipsi_csc)
            ipsi_csc{ii} = decimate(ipsi_csc{ii},csc_dsf);
            ipsi_csc{ii} = ipsi_csc{ii}(csc_inds);
            ipsi_csc_hp{ii} = filtfilt(bb,aa,ipsi_csc{ii});
            ipsi_csc{ii} = ipsi_csc{ii}/robust_std_dev(ipsi_csc_hp{ii});
            ipsi_csc_hp{ii} = ipsi_csc_hp{ii}/robust_std_dev(ipsi_csc_hp{ii});
        end
        
        mua_sm_sig = round(0.025*csc_Fsd);
        ipsi_binned_spks = nan(length(csc_time),length(ipsi_mua_times));
        ipsi_sm_rate = nan(length(csc_time),length(ipsi_mua_times));
        for ii = 1:length(ipsi_mua_times)
            ipsi_binned_spks(:,ii) = hist(ipsi_mua_times{ii},csc_time);
            ipsi_sm_rate(:,ii) = zscore(jmm_smooth_1d_cor(ipsi_binned_spks(:,ii),mua_sm_sig));
        end
        
        %%
        uu = find(~isnan(mp_t));
        mp_interp = interp1(mp_t(uu),mp_d(uu),csc_time);
        
        if ~isempty(heka_time)
            fprintf('Using heka MP data\n');
            
            uu = find(~isnan(heka_time));
            mp_interp = interp1(heka_time(uu),heka_data(uu),csc_time);
            mp_interp(isnan(mp_interp)) = 0;
        end
        
        params.Fs = csc_Fsd;
        params.tapers = [2 3];
        params.fpass = [0 40];
        movingwin = [20 10];
        clear C
        
        [b,a] = butter(2,[0.05]/(csc_Fs/2),'high');
        mp_interp(isnan(mp_interp)) = 0;
        % mp_interp = filtfilt(b,a,mp_interp);
        
        use_cscs = 1:1:length(data(dd).ipsiLFPs);
        
        clear C S2 S1
        for cc = 1:length(use_cscs)
            cur_lfp = filtfilt(b,a,ipsi_csc{use_cscs(cc)});
            %     cur_lfp = ipsi_csc{use_cscs(cc)};
            [C{cc},phi,S12,S1,S2{cc},t,f]=cohgramc(mp_interp(:),cur_lfp(:),movingwin,params);
        end
        
        mp_mwin = [2 1];
        [S_mp,mp_t,mp_f]=mtspecgramc(mp_interp(:),mp_mwin,params);
        
        uds_range = [0.4 3];
        hf_range = [20 80];
        lf_range = [0 0.25];
        use_range = [2 50];
        
        uds_freqs = find(f >= uds_range(1) & f <= uds_range(2));
        hf_freqs = find(f >= hf_range(1) & f <= hf_range(2));
        lf_freqs = find(f >= lf_range(1) & f <= lf_range(2));
        use_freqs = find(f >= use_range(1) & f <= use_range(2));
        mp_use_freqs = find(mp_f >= use_range(1) & mp_f <= use_range(2));
        
        mp_uds_pow = log10(trapz(f(uds_freqs),(S1(:,uds_freqs)),2));
        mp_hf_pow = log10(trapz(f(hf_freqs),(S1(:,hf_freqs)),2));
        mp_lf_pow = log10(trapz(f(lf_freqs),(S1(:,lf_freqs)),2));
        % mp_tot_pow = log10(trapz(f(use_freqs),(S1(:,use_freqs)),2));
        mp_tot_pow = log10(trapz(mp_f(mp_use_freqs),(S_mp(:,mp_use_freqs)),2));
        
        mp_tot_pow = (mp_tot_pow - median(mp_tot_pow))/robust_std_dev(mp_tot_pow);
        
        clear lfp*pow
        for cc = 1:length(use_cscs)
            %    lfp_uds_pow(cc,:) = trapz(f(uds_freqs),log10(S2{cc}(:,uds_freqs)),2);
            %    lfp_hf_pow(cc,:) = trapz(f(hf_freqs),log10(S2{cc}(:,hf_freqs)),2);
            %    lfp_lf_pow(cc,:) = trapz(f(lf_freqs),log10(S2{cc}(:,lf_freqs)),2);
            lfp_uds_pow(cc,:) = log10(trapz(f(uds_freqs),(S2{cc}(:,uds_freqs)),2));
            lfp_hf_pow(cc,:) = log10(trapz(f(hf_freqs),(S2{cc}(:,hf_freqs)),2));
            lfp_lf_pow(cc,:) = log10(trapz(f(lf_freqs),(S2{cc}(:,lf_freqs)),2));
        end
        
        lfp_uds_pow = interp1(t,lfp_uds_pow',mp_t)';
        lfp_hf_pow = interp1(t,lfp_hf_pow',mp_t)';
        lfp_lf_pow = interp1(t,lfp_lf_pow',mp_t)';
        %%
        avg_uds_pow = nanmean(lfp_uds_pow,2);
        avg_lf_pow = nanmean(lfp_lf_pow,2);
        [~,ctx_ch] = max(avg_uds_pow(poss_ctx_chs)-0*avg_lf_pow(poss_ctx_chs));
        ctx_ch = poss_ctx_chs(ctx_ch);
        
        % lfp_uds_thresh = -2.5;
        lfp_uds_thresh = -0.5;
        % lfp_lf_max = 0.5;
        lfp_lf_max = 1;
        uds_epochs = (lfp_uds_pow(ctx_ch,:) >= lfp_uds_thresh);
        lf_epochs = (lfp_lf_pow(ctx_ch,:) >= lfp_lf_max);
        
        % has_uds = true(length(csc_time),1);
        % for ii = 1:length(t)
        %    if ~ismember(ii,uds_epochs)
        %       curset = find(csc_time >= t(ii)-movingwin(1)/2 & csc_time <= t(ii) + movingwin(1)/2);
        %       has_uds(curset) = false;
        %    end
        % end
        
        %%
        f1 = figure();
        subplot(3,1,1); hold on
        imagesc(t,1:length(use_cscs),lfp_uds_pow); colorbar;
        plot(mp_t,uds_epochs*4,'w','linewidth',2);
        line([0 t(end)],ctx_ch + [0 0],'color','k');
        axis tight
        xlim([0 t(end)]);
        xlabel('Time (s)');
        ylabel('Channel num');
        
        subplot(3,1,2); hold on
        imagesc(t,1:length(use_cscs),lfp_lf_pow); colorbar;
        plot(mp_t,lf_epochs*4,'w','linewidth',2);
        line([0 t(end)],ctx_ch + [0 0],'color','k');
        axis tight
        xlim([0 t(end)]);
        xlabel('Time (s)');
        ylabel('Channel num');
        
        good_uds = uds_epochs & ~lf_epochs;
        subplot(3,1,3); hold on
        pcolor(t,f,log10(S2{ctx_ch})'); shading flat; colorbar;
        plot(mp_t,good_uds*3+0.1,'w','linewidth',2);
        ylim([0.1 4]);
        caxis([-2 1])
        set(gca,'yscale','log');
        xlim([0 t(end)]);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        
        dname = find(data(dd).dir == '/',1,'last');
        dname = data(dd).dir(dname+1:end);
        fname = strcat(fig_dir,'NUDS_epochs_',data(dd).MPloc,dname);
        figufy(f1);
        set(f1,'PaperUnits','inches','PaperSize',[8 10],'PaperPosition',[0 0 8 10]);
        print(f1,'-dpng',fname);
        close(f1);
        
        %%
        addpath('~/James_scripts/hsmm_uds_toolbox/')
        removal_window = 0; %buffer window around desynch epochs
        
        %mark points where the SO power crosses below threshold
        desynch_indicator = ~uds_epochs | lf_epochs | mp_tot_pow' < -3;
        desynch_start_ids = 1+find(desynch_indicator(1:end-1) == 0 & desynch_indicator(2:end) == 1);
        desynch_stop_ids = 1+find(desynch_indicator(1:end-1) == 1 & desynch_indicator(2:end) == 0);
        
        if isempty(desynch_start_ids) && ~isempty(desynch_stop_ids)
            desynch_start_ids = [1];
        end
        
        %make sure you start and stop in desynchronized epochs in correct order
        if ~isempty(desynch_start_ids)
            if isempty(desynch_stop_ids)
                desynch_stop_ids = length(mp_t);
            else
                if desynch_start_ids(1) > desynch_stop_ids(1)
                    desynch_start_ids = [1 desynch_start_ids];
                end
                if desynch_start_ids(end) > desynch_stop_ids(end)
                    desynch_stop_ids = [desynch_stop_ids length(mp_t)];
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
            if length(mp_t)-desynch_stop_ids(w) <= removal_window
                desynch_stop_ids(w) = length(mp_t);
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
        
        desynch_start_times = mp_t(desynch_start_ids);
        desynch_stop_times = mp_t(desynch_stop_ids);
        desynch_times = [desynch_start_times(:) desynch_stop_times(:)];
        desynch_ids = round(desynch_times*csc_Fsd);
        
        %% Compute signal features to use for UDS classification
        hmm_dsf = 1; %down-sample-factor for HMM signal featurs
        hmm_Fs = csc_Fsd/hmm_dsf; %sample-frequency of HMM signal features (Hz)
        lf_cut_off_freqs = [0.025 10]; %for "low-frequency amplitude" features, these are the cutoff-frequencies for the bandpass filter
        [bb,aa] = butter(2,lf_cut_off_freqs/(csc_Fsd/2));
        mp_features = filtfilt(bb,aa,mp_interp);
        mp_features = downsample(mp_features,hmm_dsf);
        
%         [mp_features,t_axis] = hsmm_uds_get_lf_features(mp_interp,csc_Fsd,hmm_Fs,lf_cut_off_freqs); %'low-frequency amplitude'
        % mp_features = mp_features/robust_std_dev(mp_features);
        
        lf_cut_off_freqs = [0.2 3]; %for "low-frequency amplitude" features, these are the cutoff-frequencies for the bandpass filter
        [lfp_features,t_axis] = hsmm_uds_get_lf_features(ipsi_csc{ctx_ch},csc_Fsd,hmm_Fs,lf_cut_off_freqs); %'low-frequency amplitude'
        % hf_cut_off_freqs = [20 80]; %for "high-frequency power" features, these are the cutoff-frequencies for the bandpass filter
        % smooth_sigma = 0.15; %smoothing sigma for computing HF power in seconds
        % [lfp_features,t_axis] = hsmm_uds_get_hf_features(ipsi_csc{ctx_ch},csc_Fs,hmm_Fs,hf_cut_off_freqs,smooth_sigma); %'high-frequency power'
        lfp_features = lfp_features/robust_std_dev(lfp_features);
        % lfp_features(lfp_features > 10) = 10;
        
        lf_cut_off_freqs = [0.2 20]; %for "low-frequency amplitude" features, these are the cutoff-frequencies for the bandpass filter
        [lfp_features_bb,t_axis] = hsmm_uds_get_lf_features(ipsi_csc{ctx_ch},csc_Fsd,hmm_Fs,lf_cut_off_freqs); %'low-frequency amplitude'
        lfp_features_bb = lfp_features_bb/robust_std_dev(lfp_features_bb);
        
        
        %% LOCATE THE SEGMENTS CONTAINING UDS
        T = length(mp_features); %number of samples
        min_seg_dur = 20; %minimum duration of UDS segments used for analysis (in sec)
        UDS_segs = hsmm_uds_get_uds_segments(desynch_times,hmm_Fs,T,min_seg_dur); %Nx2 matrix containing the index values of the beginning and end of each UDS segment
        
        is_uds = false(size(lfp_features));
        for ii = 1:size(UDS_segs,1)
            is_uds(UDS_segs(ii,1):UDS_segs(ii,2)) = true;
        end
        
        %%
        addpath('~/Analysis/Mayank/sleep/');
%         xx = linspace(-3.5,3.5,100);
        mp_dist_xx = linspace(prctile(mp_features(is_uds),0.1),prctile(mp_features(is_uds),99.9),100);
        if sum(is_uds) > 0
            mp_pdf = hist(mp_features(is_uds),mp_dist_xx);
            mp_kpdf = ksdensity(mp_features(is_uds),mp_dist_xx);
            mp_pdf = mp_pdf/mean(mp_pdf)*mean(mp_kpdf);
            
            [~,dens_peaks] = findpeaks(mp_kpdf,'npeaks',2,'sortstr','descend','minpeakdistance',round(0.25/median(diff(mp_dist_xx))));
            dens_peaks = mp_dist_xx(dens_peaks);
            if length(dens_peaks) == 2
                cross_pt = mean(dens_peaks);
            else
                cross_pt = median(mp_features(is_uds));
            end
            
            mp_ovPup = sum(mp_features(is_uds) > cross_pt)/sum(is_uds);
            mp_ovPdown = sum(mp_features(is_uds) < cross_pt)/sum(is_uds);
            
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
                [dip, p_value, xlow,xup]=hartigansdipsigniftest(mp_features(is_uds),100);
            else
                dip = nan; p_value = nan;
            end
            
            f1 = figure(); hold on
            plot(mp_dist_xx,mp_pdf,mp_dist_xx,mp_kpdf,'r');
            % plot(xx,g1,'k',xx,g2,'g');
            yl = ylim();
            line(cross_pt + [0 0],yl,'color','k');
            title(sprintf('Dip: %.3f  P: %.3f',dip,p_value));
            
            dname = find(data(dd).dir == '/',1,'last');
            dname = data(dd).dir(dname+1:end);
            fname = strcat(fig_dir,'MP_dist_',data(dd).MPloc,dname);
            figufy(f1);
            set(f1,'PaperUnits','inches','PaperSize',[6 5],'PaperPosition',[0 0 6 5]);
            print(f1,'-dpng',fname);
            close(f1);
        end
        %%
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
        
        back_win = round(0.75*csc_Fsd);
        for_win = round(0.75*csc_Fsd);
        
        min_n_events = 10;
        lfp_bndwidth = 0.15;
        mp_bndwidth = diff(prctile(mp_features(is_uds),[25 75]))/15;
        
%         mp_dist_xx = linspace(-3,3,100);
        % gm_posterior = posterior(gmfit,xx');
        % p1 = find(xx >= gmfit.mu(1),1); p2 = find(xx >= gmfit.mu(2),1);
        % peak_locs = [p1; p2]; peak_locs = sort(peak_locs);
        % between_peaks = peak_locs(1):peak_locs(2);
        % [~,cross_pt] = min(abs(gm_posterior(between_peaks,1) - 0.5));
        % cross_pt = between_peaks(cross_pt);
        cross_pt_mp = find(mp_dist_xx >= cross_pt,1,'first');
        
        all_up_tas = nan(length(poss_amp_prc),back_win + for_win+1);
        all_up_dists = nan(length(poss_amp_prc),back_win+for_win+1,length(mp_dist_xx));
        all_down_tas = nan(length(poss_amp_prc),back_win + for_win+1);
        all_down_dists = nan(length(poss_amp_prc),back_win+for_win+1,length(mp_dist_xx));
        for jj = 1:length(poss_amp_prc)
            %     cur_events = lfp_valleys(lfp_features(lfp_valleys) < -poss_amps(jj));
            %     cur_events(~is_uds(cur_events)) = [];
            cur_events = lfp_valleys(lfp_valley_amps >= poss_amps_valleys(jj));
            if length(cur_events) >= min_n_events
                [cur_ta,cur_lags,~,~,cur_mat] = get_event_trig_avg_v3(mp_features,cur_events,back_win,for_win,2);
                all_down_tas(jj,:) = cur_ta;
                
                for ii = 1:length(cur_lags)
                    all_down_dists(jj,ii,:) = ksdensity(cur_mat(:,ii),mp_dist_xx,'bandwidth',mp_bndwidth);
                end
            end
            
            %     cur_events = lfp_peaks(lfp_features(lfp_peaks) > poss_amps(jj));
            %     cur_events(~is_uds(cur_events)) = [];
            cur_events = lfp_peaks(lfp_peak_amps >= poss_amps_peaks(jj));
            if length(cur_events) >= min_n_events
                [cur_ta,cur_lags,~,~,cur_mat] = get_event_trig_avg_v3(mp_features,cur_events,back_win,for_win,2);
                all_up_tas(jj,:) = cur_ta;
                
                for ii = 1:length(cur_lags)
                    all_up_dists(jj,ii,:) = ksdensity(cur_mat(:,ii),mp_dist_xx,'bandwidth',mp_bndwidth);
                end
            end
            
        end
        
        utrig_Pmp_down = trapz(mp_dist_xx(1:(cross_pt_mp-1)),all_up_dists(:,:,1:(cross_pt_mp-1)),3);
        utrig_Pmp_up = trapz(mp_dist_xx((cross_pt_mp+1):end),all_up_dists(:,:,(cross_pt_mp+1):end),3);
        dtrig_Pmp_down = trapz(mp_dist_xx(1:(cross_pt_mp-1)),all_down_dists(:,:,1:(cross_pt_mp-1)),3);
        dtrig_Pmp_up = trapz(mp_dist_xx((cross_pt_mp+1):end),all_down_dists(:,:,(cross_pt_mp+1):end),3);
        
        %%
        lfp_dist_xx = linspace(-5,3,100);
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
                    all_down_dists_lfp(jj,ii,:) = ksdensity(cur_mat(:,ii),lfp_dist_xx,'bandwidth',lfp_bndwidth);
                end
            end
            
            %     cur_events = lfp_peaks(lfp_features(lfp_peaks) > poss_amps(jj));
            %     cur_events(~is_uds(cur_events)) = [];
            cur_events = lfp_peaks(lfp_peak_amps >= poss_amps_peaks(jj));
            if length(cur_events) >= min_n_events
                [cur_ta,cur_lags,~,~,cur_mat] = get_event_trig_avg_v3(lfp_features_bb,cur_events,back_win,for_win,2);
                all_up_tas_lfp(jj,:) = cur_ta;
                
                for ii = 1:length(cur_lags)
                    all_up_dists_lfp(jj,ii,:) = ksdensity(cur_mat(:,ii),lfp_dist_xx,'bandwidth',lfp_bndwidth);
                end
            end
            
        end
        
        %%
        if sum(is_uds) > 0
        max_utrig_Pdown = max(utrig_Pmp_down(:));
        max_dtrig_Pup = max(dtrig_Pmp_up(:));
        if isnan(max_utrig_Pdown) | max_utrig_Pdown == 0
            max_utrig_Pdown = 0.1;
        end
        if isnan(max_dtrig_Pup) | max_dtrig_Pup == 0
            max_dtrig_Pup = 0.1;
        end
        f1 = figure();
        for ii = 1:length(poss_amp_prc)
            subplot(5,3,3*(ii-1)+1)
            imagesc(cur_lags/csc_Fsd,lfp_dist_xx,log(squeeze(all_up_dists_lfp(ii,:,:)))');
            caxis([-3 -0.5]);
            set(gca,'ydir','normal');
            xlabel('Lag (s)');
            ylabel('Amplitude (z)');
            title(sprintf('LFP %d pcnt',poss_amp_prc(ii)));
            
            subplot(5,3,3*(ii-1)+2)
            imagesc(cur_lags/csc_Fsd,mp_dist_xx,log(squeeze(all_up_dists(ii,:,:)))'); hold on
            xl = xlim(); line(xl,mp_dist_xx(cross_pt_mp) + [0 0],'color','w');
            caxis([-5 -1]);
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
            imagesc(cur_lags/csc_Fsd,lfp_dist_xx,log(squeeze(all_down_dists_lfp(ii,:,:)))');
            caxis([-3 -0.5]);
            set(gca,'ydir','normal');
            xlabel('Lag (s)');
            ylabel('Amplitude (z)');
            title(sprintf('LFP %d pcnt',poss_amp_prc(ii)));
            
            subplot(5,3,3*(ii-1)+2)
            imagesc(cur_lags/csc_Fsd,mp_dist_xx,log(squeeze(all_down_dists(ii,:,:)))');
            xl = xlim(); line(xl,mp_dist_xx(cross_pt_mp) + [0 0],'color','w');
            caxis([-5 -1]);
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