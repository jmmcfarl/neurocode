%% LOAD PROCESSED DATA
clear all
close all

data_dir = '/media/NTlab_data3/Data/bruce/CPdata/';
fig_dir = '/home/james/Analysis/bruce/ChoiceProb/';
addpath('~/James_scripts/bruce/bruce_code');

addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))

data_sets = what(data_dir);data_sets = data_sets.mat;

for eee = 1:length(data_sets)
    % Expt_name = 'G050';
    Expt_name = data_sets{eee};
    
    load(strcat(data_dir, Expt_name));
    Expt_name = Expt_name(1:end-4);
    
    %% Assemble stim, Robs, LFP into trial structure
    Ntrials = length(AllExpt.Expt.Trials);
    dt = 0.005;  % assuming 10 ms frame time exactly
    trialDur = 2; %in sec.
    Nframes = 200; %number of video frames (@100Hz) per trial
    NT = ceil(trialDur/dt); %number of time bins per trial
    
    BlockStartid = AllExpt.Expt.Header.BlockStartid;
    n_blocks = length(BlockStartid);
    
    
    %% SU and MU data (combined)
    Nunits = length(AllExpt.Header);
    Ucellnumbers = [AllExpt.Header(:).cellnumber];
    Uchans = round([AllExpt.Header(:).probe]);
    
    %% relevant trial variables
    trialNums = [AllExpt.Expt.Trials(:).Trial];
    trialOR = [AllExpt.Expt.Trials(:).or];
    trialOB = [AllExpt.Expt.Trials(:).ob];
    trialrwDir = [AllExpt.Expt.Trials(:).rwdir];
    trialrespDir = [AllExpt.Expt.Trials(:).RespDir];
    trialSe = [AllExpt.Expt.Trials(:).se];
    trialId = [AllExpt.Expt.Trials(:).id];
    
    zsOB = max(abs(trialOB));
    msOB = min(abs(trialOB));
    %% identify the ids of the trials at the beginning and end of each block
    block_trial_boundaries = find(ismember(trialId,BlockStartid));
    if length(block_trial_boundaries) ~= n_blocks
        warning('Block alignment problem');
        block_trial_boundaries = nan(1,length(BlockStartid));
        for ii = 1:length(BlockStartid)
            block_trial_boundaries(ii) = find(trialId >= BlockStartid(ii),1);
        end
    end
    block_trial_boundaries = [block_trial_boundaries(:) [block_trial_boundaries(2:end)-1 Ntrials]'];
    
    %% get total trial-by-trial spike counts for each unit
    tot_spks_per_trial = zeros(Ntrials,Nunits);
    spk_delay_buff = 0.05;
    spk_end_exclude = 0;
    spk_beg_exclude = 0.5;
    for tr = 1:Ntrials
        for nn = 1:Nunits
            trindx = find( AllExpt.Spikes{nn}.Trial == trialNums(tr));
            if ~isempty(trindx)
                tot_spks_per_trial(tr,nn) = sum(double(AllExpt.Spikes{nn}.Spikes{trindx})*1e-4 <= (trialDur + spk_delay_buff - spk_end_exclude) & ...
                    (double(AllExpt.Spikes{nn}.Spikes{trindx})*1e-4 >= spk_delay_buff + spk_beg_exclude));
            end
        end
    end
    
    %find any blocks where the unit has no spikes and set these values to nans
    %(unit presumably not clustered in those blocks)
    spikes_per_block = nan(n_blocks,Nunits);
    for bb = 1:n_blocks
        spikes_per_block(bb,:) = sum(tot_spks_per_trial(block_trial_boundaries(bb,1):block_trial_boundaries(bb,2),:));
    end
    
    tot_spks_per_trial_norm = tot_spks_per_trial;
    for cc = 1:Nunits
        bad_blocks = find(spikes_per_block(:,cc) == 0);
        for ii = bad_blocks'
            tot_spks_per_trial_norm(block_trial_boundaries(ii,1):block_trial_boundaries(ii,2),cc) = nan;
        end
    end
    
    %% Get Binned Spikes
    fullRobs = zeros(NT,Ntrials,Nunits);
    for tr = 1:Ntrials
        for nn = 1:Nunits
            trindx = find( AllExpt.Spikes{nn}.Trial == trialNums(tr));
            if ~isempty(trindx)
                spks = double( AllExpt.Spikes{nn}.Spikes{trindx} ) * 1e-4;  % no adjustment from start of trial
                if ~isempty(spks)
                    cur_Robs = histc( spks, (0:(NT))*dt);
                    fullRobs(:,tr,nn) = cur_Robs(1:end-1);
                end
            end
        end
    end
    for cc = 1:Nunits
        fullRobs(:,isnan(tot_spks_per_trial_norm(:,cc)),cc) = nan;
    end
    
    trialSpkCnts = squeeze(nansum(fullRobs,1));
    totSpkCnts = squeeze(nansum(trialSpkCnts,1));
    avg_spk_rates = nanmean(reshape(fullRobs,[],Nunits));
    
    %% CALCULATE CHOICE PROBS
    tot_spks_per_trial_norm = tot_spks_per_trial;
    tot_spks_per_trial_norm(tot_spks_per_trial == 0) = nan;
    
    sig_trials = find(abs(trialOB) == msOB & trialrespDir ~= 0);
    stim_up = sig_trials(trialrwDir(sig_trials)==1);
    stim_down = sig_trials(trialrwDir(sig_trials)==-1);
    sig_prob = nan(Nunits,1);
    sp_pval = nan(Nunits,1);
    for cc = 1:Nunits
        %         cur_resp_ax = 0:(nanmax(tot_spks_per_trial_norm(sig_trials,cc)));
        %         if ~isnan(cur_resp_ax)
        %         up_hist = histc(tot_spks_per_trial_norm(stim_up,cc),cur_resp_ax);
        %         down_hist = histc(tot_spks_per_trial_norm(stim_down,cc),cur_resp_ax);
        %         fract_cor = cumsum(up_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_up,cc)));
        %         fract_incor = cumsum(down_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_down,cc)));
        %         sig_prob(cc) = trapz(fract_incor,fract_cor);
        %         [~,sp_pval(cc)] = ttest2(tot_spks_per_trial_norm(stim_up,cc),tot_spks_per_trial_norm(stim_down,cc));
        %         end
        if nansum(tot_spks_per_trial_norm(sig_trials,cc)) > 0
        [~,~,~,sig_prob(cc)]=perfcurve(trialrwDir(sig_trials), tot_spks_per_trial_norm(sig_trials,cc),1);
        end
    end
    
    test_trials = find(trialOB == zsOB & trialrespDir ~= 0);
    % test_trials = find(trialrespDir ~= 0);
    resp_up = test_trials(trialrespDir(test_trials) == 1);
    resp_down = test_trials(trialrespDir(test_trials) == -1);
    choice_prob = nan(Nunits,1);
    cp_pval = nan(Nunits,1);
    for cc = 1:Nunits
        %         cur_resp_ax = 0:(nanmax(tot_spks_per_trial_norm(test_trials,cc)));
        %         if ~isnan(cur_resp_ax)
        %         up_hist = histc(tot_spks_per_trial_norm(resp_up,cc),cur_resp_ax);
        %         down_hist = histc(tot_spks_per_trial_norm(resp_down,cc),cur_resp_ax);
        %         true_pos = cumsum(up_hist)/sum(~isnan(tot_spks_per_trial_norm(resp_up,cc)));
        %         false_pos = cumsum(down_hist)/sum(~isnan(tot_spks_per_trial_norm(resp_down,cc)));
        %         choice_prob(cc) = trapz(false_pos,true_pos);
        %         [~,~,~,CP_Neuron(cc)]=perfcurve(trialrespDir(test_trials), tot_spks_per_trial_norm(test_trials,cc),1);
        %         [~,cp_pval(cc)] = ttest2(tot_spks_per_trial_norm(resp_up,cc),tot_spks_per_trial_norm(resp_down,cc));
        %         end
        if nansum(tot_spks_per_trial_norm(test_trials,cc)) > 0
        [~,~,~,choice_prob(cc)]=perfcurve(trialrespDir(test_trials), tot_spks_per_trial_norm(test_trials,cc),1);
        end
    end
    
    %%
    tot_spks_z = nanzscore(tot_spks_per_trial_norm);
    avg_spks = nanmean(tot_spks_z,2);
    
    [~,~,~,avg_SP]=perfcurve(trialrwDir(sig_trials), avg_spks(sig_trials),1);
    
    [~,~,~,avg_CP]=perfcurve(trialrespDir(test_trials), avg_spks(test_trials),1);
    
    [avg_SP_res,avg_CP_res] = deal(nan(Nunits,1));
    for cc = 1:Nunits
        y = tot_spks_per_trial_norm(:,cc);
        X = [avg_spks ones(length(avg_spks),1)];
        r = regress(y,X);
        pred_y = X*r;
        res_y = y - pred_y;
        if any(~isnan(res_y(sig_trials)))
        [~,~,~,avg_SP_res(cc)]=perfcurve(trialrwDir(sig_trials), res_y(sig_trials),1);
        end
        if any(~isnan(res_y(test_trials)))
        [~,~,~,avg_CP_res(cc)]=perfcurve(trialrespDir(test_trials), res_y(test_trials),1);
        end
    end
    
    %%
    %     n_factors = 5;
    %     use_cells = find(~any(isnan(tot_spks_z)));
    % %     [lambda,psi,T,stats] = factoran(tot_spks_z(:,use_cells),n_factors);
    %
    %      [COEFF, SCORE, LATENT] = princomp(tot_spks_z(:,use_cells));
    %
    %     for ii = 1:n_factors
    % %         fact_out = tot_spks_z(:,use_cells)*lambda(:,ii);
    %         fact_out = tot_spks_z(:,use_cells)*COEFF(:,ii);
    %       [~,~,~,fact_SP(ii)]=perfcurve(trialrwDir(sig_trials), fact_out(sig_trials),1);
    %       [~,~,~,fact_CP(ii)]=perfcurve(trialrespDir(test_trials), fact_out(test_trials),1);
    %     end
    
    %%
    fullRobs_norm = bsxfun(@rdivide,fullRobs,reshape(avg_spk_rates,[1 1 Nunits]));
    avg_rup = squeeze(nanmean(fullRobs_norm(:,resp_up,:),2));
    avg_rdown = squeeze(nanmean(fullRobs_norm(:,resp_down,:),2));
    avg_sup = squeeze(nanmean(fullRobs_norm(:,stim_up,:),2));
    avg_sdown = squeeze(nanmean(fullRobs_norm(:,stim_down,:),2));
    
    sm_win = round(0.025/dt);
    for ii = 1:Nunits
        avg_rup(:,ii) = jmm_smooth_1d_cor(avg_rup(:,ii),sm_win);
        avg_rdown(:,ii) = jmm_smooth_1d_cor(avg_rdown(:,ii),sm_win);
        avg_sup(:,ii) = jmm_smooth_1d_cor(avg_sup(:,ii),sm_win);
        avg_sdown(:,ii) = jmm_smooth_1d_cor(avg_sdown(:,ii),sm_win);
    end
    
    
    tax = (1:NT)*dt;
    f1 = figure('visible','off');
    subplot(2,1,1);
    hold on;
    plot(tax,nanmean(avg_rup,2),'b',tax,nanmean(avg_rdown,2),'r',tax,nanmean(avg_sup,2),'k',tax,nanmean(avg_sdown,2),'g');
    legend('Chose up','Chose down','Stim up','Stim down','Location','Southeast');
    shadedErrorBar(tax,nanmean(avg_rup,2),nanstd(avg_rup,[],2)/sqrt(Nunits),{'color','b'});
    shadedErrorBar(tax,nanmean(avg_rdown,2),nanstd(avg_rdown,[],2)/sqrt(Nunits),{'color','r'});
    shadedErrorBar(tax,nanmean(avg_sup,2),nanstd(avg_sup,[],2)/sqrt(Nunits),{'color','k'});
    shadedErrorBar(tax,nanmean(avg_sdown,2),nanstd(avg_sdown,[],2)/sqrt(Nunits),{'color','g'});
    line(tax([1 end]),[1 1],'color','k');
    
    subplot(2,1,2);
    hold on
    plot(tax,nanmean(avg_rup-avg_rdown,2),'b',tax,nanmean(avg_sup-avg_sdown,2),'k');
    legend('Choice up-down','Stim up-down','Location','Southeast');
    shadedErrorBar(tax,nanmean(avg_rup-avg_rdown,2),nanstd(avg_rup-avg_rdown,[],2)/sqrt(Nunits),{'color','b'});
    shadedErrorBar(tax,nanmean(avg_sup-avg_sdown,2),nanstd(avg_sup-avg_sdown,[],2)/sqrt(Nunits),{'color','k'});
    line(tax([1 end]),[0 0],'color','k');
    
    fig_width = 6; rel_height = 1.6;
    fig_name = [fig_dir sprintf('%s_CPSP_PSTHcompare.pdf',Expt_name)];
    figufy(f1);
    exportfig(f1,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    close(f1);
    
    %%
    xr = [0.3 0.7]; yr = [0 1];
    h = scatterhist(choice_prob,sig_prob,'nbins',50);
    set(gcf,'visible','off');
    hold on
    plot(choice_prob(Ucellnumbers > 0),sig_prob(Ucellnumbers > 0),'r.');
    line(xr,[0.5 0.5],'color','k'); line([0.5 0.5],yr,'color','k');
    xlim(xr); ylim(yr);
    title(sprintf('ori: %d',unique(trialOR)));
    
    fig_width = 5; rel_height = 1;
    fig_name = [fig_dir sprintf('%s_CPSP_hist.pdf',Expt_name)];
    figufy(gcf);
    exportfig(gcf,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    close(gcf);
    
    %%
    xr = [0.3 0.7]; yr = [0 1];
    h = scatterhist(avg_CP_res,avg_SP_res,'nbins',50);
    set(gcf,'visible','off');
    hold on
    plot(avg_CP_res(Ucellnumbers > 0),avg_SP_res(Ucellnumbers > 0),'r.');
    line(xr,[0.5 0.5],'color','k'); line([0.5 0.5],yr,'color','k');
    xlim(xr); ylim(yr);
    title(sprintf('avg CP: %.4f  avg SP: %.4f',avg_CP,avg_SP));
    
    fig_width = 5; rel_height = 1;
    fig_name = [fig_dir sprintf('%s_CPSPRES_hist.pdf',Expt_name)];
    figufy(gcf);
    exportfig(gcf,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    close(gcf);
    
    %%
    un_OBs = [-60 -70 -80 zsOB 80 70 60];
    OB_x = [-3 -2 -1 0 1 2 3];
    trialOB_prime = trialOB;
    if min(trialOB) > 0
        trialOB_prime(trialrwDir==-1) = -trialOB_prime(trialrwDir==-1);
    end
    
    fract_one = nan(length(un_OBs),1);
    fract_one_uci = nan(length(un_OBs),1);
    fract_one_lci = nan(length(un_OBs),1);
    zval = 1.96;
    for ss = 1:length(un_OBs)
        cur_OB_neg = find(trialOB_prime == un_OBs(ss) & trialrespDir ~= 0);
        curN = length(cur_OB_neg);
        fract_one(ss) = sum(trialrespDir(cur_OB_neg) == 1)/curN;
        
        %uncertainty estimates based on Wilson score interval formula (this
        %gives better approx of 95% CI
        fract_one_uci(ss) = 1/(1+zval^2/curN)*(fract_one(ss) + zval^2/(2*curN) + zval*sqrt(fract_one(ss)/curN*(1-fract_one(ss)) + zval^2/(4*curN^2)));
        fract_one_lci(ss) = 1/(1+zval^2/curN)*(fract_one(ss) + zval^2/(2*curN) - zval*sqrt(fract_one(ss)/curN*(1-fract_one(ss)) + zval^2/(4*curN^2)));
    end
    
    sig_strengths = [zsOB 80 70 60];
    cmet_avgs = nan(Nunits,length(sig_strengths));
    nmet_avgs = nan(Nunits,length(sig_strengths));
    for ss = 1:length(sig_strengths)
        cur_OB_pos = find(abs(trialOB_prime) == sig_strengths(ss) & trialrwDir == 1 & trialrespDir ~= 0);
        cur_OB_neg = find(abs(trialOB_prime) == sig_strengths(ss) & trialrwDir == -1 & trialrespDir ~= 0);
        cmet_avgs(:,ss) = nanmean(tot_spks_per_trial_norm(cur_OB_pos,:)) - nanmean(tot_spks_per_trial_norm(cur_OB_neg,:));
        
        cur_OB_pos = find(abs(trialOB_prime) == sig_strengths(ss) & trialrespDir == 1);
        cur_OB_neg = find(abs(trialOB_prime) == sig_strengths(ss) & trialrespDir == -1);
        nmet_avgs(:,ss) = nanmean(tot_spks_per_trial_norm(cur_OB_pos,:)) - nanmean(tot_spks_per_trial_norm(cur_OB_neg,:));
    end
    % cmet_avgs = bsxfun(@rdivide,cmet_avgs,nanmean(tot_spks_per_trial_norm)');
    % nmet_avgs = bsxfun(@rdivide,nmet_avgs,nanmean(tot_spks_per_trial_norm)');
    cmet_avgs = bsxfun(@rdivide,cmet_avgs,nanstd(tot_spks_per_trial_norm)');
    nmet_avgs = bsxfun(@rdivide,nmet_avgs,nanstd(tot_spks_per_trial_norm)');
    
    f1 = figure('visible','off');
    subplot(2,2,1);
    errorbar(OB_x,fract_one,fract_one-fract_one_lci,fract_one_uci-fract_one);
    xlabel('Signal strength');
    ylabel('Fraction chose pos');
    ylim([0 1]);
    xl = xlim();
    line(xl,[0.5 0.5],'color','k');
    axis tight
    
    subplot(2,2,2);
    errorbar(0:(length(sig_strengths)-1),nanmean(nmet_avgs),nanstd(nmet_avgs)/sqrt(Nunits));
    hold on
    errorbar(0:(length(sig_strengths)-1),nanmean(cmet_avgs),nanstd(cmet_avgs)/sqrt(Nunits),'r');
    legend('Choice-based','Stim-based');
    axis tight
    ylm = max(abs(ylim()));ylim([-ylm ylm]);
    xl = xlim();
    line(xl,[0 0],'color','k');
    
    subplot(2,2,3);
    errorbar(0:(length(sig_strengths)-1),nanmean(abs(nmet_avgs)),nanstd(abs(nmet_avgs))/sqrt(Nunits));
    hold on
    errorbar(0:(length(sig_strengths)-1),nanmean(abs(cmet_avgs)),nanstd(abs(cmet_avgs))/sqrt(Nunits),'r');
    legend('Choice-based','Stim-based');
    axis tight
    
    resp_up_trials = find(trialrespDir == 1);
    resp_down_trials = find(trialrespDir == -1);
    ET_data_up = cat(3,AllExpt.Expt.Trials(resp_up_trials).EyeData);
    ET_data_down = cat(3,AllExpt.Expt.Trials(resp_down_trials).EyeData);
    
    ET_data_up = int2double(ET_data_up, AllExpt.Expt.Header.emscale);
    ET_data_down = int2double(ET_data_down, AllExpt.Expt.Header.emscale);
    
    trange = (size(ET_data_up,1)-100):(size(ET_data_up,1)-10);
    ET_data_up = permute(ET_data_up,[1 3 2]); ET_data_down = permute(ET_data_down,[1 3 2]);
    
    ET_data_up = bsxfun(@minus,ET_data_up,nanmedian(ET_data_up,2));
    ET_data_down = bsxfun(@minus,ET_data_down,nanmedian(ET_data_down,2));

    %     ET_data_up = bsxfun(@minus,ET_data_up,reshape(nanmedian(reshape(ET_data_up,[],4)),1,1,4));
%     ET_data_down = bsxfun(@minus,ET_data_down,reshape(nanmedian(reshape(ET_data_down,[],4)),1,1,4));
    
    em_r = [-8 8];
    subplot(2,2,4); hold on
    plot(squeeze(ET_data_up(trange,:,3)),squeeze(ET_data_up(trange,:,4)),'b.','markersize',1);
    plot(squeeze(ET_data_down(trange,:,3)),squeeze(ET_data_down(trange,:,4)),'r.','markersize',1);
    plot(squeeze(ET_data_up(trange,:,1)),squeeze(ET_data_up(trange,:,2)),'k.','markersize',1);
    plot(squeeze(ET_data_down(trange,:,1)),squeeze(ET_data_down(trange,:,2)),'g.','markersize',1);
    xlim(em_r); ylim(em_r); line(em_r,[0 0],'color','k'); line([0 0],em_r,'color','k');
    fig_width = 8; rel_height = 1;
    fig_name = [fig_dir sprintf('%s_ratecurves.pdf',Expt_name)];
    figufy(f1);
    exportfig(f1,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    close(f1);
    
    %% USE WAVELET ANALYSIS TO COMPUTE PHASE-LOCKING SPECTRA FOR EACH UNIT
    if isfield(AllExpt.Expt.Header,'LFPsamplerate')
        LFP_Fs = double(1/AllExpt.Expt.Header.LFPsamplerate);
        if Expt_name(1) == 'G'
        LFP_offset = AllExpt.Expt.Header.preperiod/1e4;
        LFP_trial_taxis = (1:length(AllExpt.Expt.Header.LFPtimes))/LFP_Fs - LFP_offset;
        else
            LFP_trial_taxis = AllExpt.Expt.Header.LFPtimes;
            if nanmax(LFP_trial_taxis) > 1e3
                LFP_trial_taxis = LFP_trial_taxis/1e4;
            end
        end
        LFP_dsf = 2;
        LFP_Fsd = LFP_Fs/LFP_dsf;
        %anti-aliasing filter and high-pass filter
        aa_hcf = LFP_Fsd/2*0.8;
        % [b_aa,a_aa] = butter(4,aa_hcf/(LFP_Fs/2),'low');
        aa_lcf = 0.5;
        [b_aa,a_aa] = butter(2,[aa_lcf aa_hcf]/(LFP_Fs/2));
        
        LFP_trial_taxis_ds = downsample(LFP_trial_taxis,LFP_dsf);
        TLEN = length(LFP_trial_taxis_ds);
        
        if Expt_name(1) == 'G'
            nprobes = 96;
            uprobes = 1:4:nprobes;
        else
            nprobes = 24;
            uprobes = 1:2:nprobes;
        end
        
        %     %wavelet parameters
        %     nwfreqs = 25;
        %     min_freq = 1.; max_freq = 80;
        %     min_scale = 1/max_freq*LFP_Fsd;
        %     max_scale = 1/min_freq*LFP_Fsd;
        %     wavetype = 'cmor1-1';
        %     scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
        %     wfreqs = scal2frq(scales,wavetype,1/LFP_Fsd);
        %
        %     % beg_buffer = round(0.15/dt);
        %     beg_buffer = round(0/dt);
        %     end_buffer = round(0/dt);
        %     trial_dur = round(2/dt);
        %     R_trial_taxis = (1+beg_buffer:(trial_dur - end_buffer))*dt;
        %     TLEN = length(R_trial_taxis);
        
        %     trial_LFP_real = nan(Ntrials,TLEN,length(wfreqs),length(uprobes));
        %     trial_LFP_imag = nan(Ntrials,TLEN,length(wfreqs),length(uprobes));
        trial_LFPs = nan(Ntrials,TLEN,length(uprobes));
        for tr = 1:Ntrials
%             fprintf('Trial %d of %d\n',tr,Ntrials);
            if size(AllExpt.Expt.Trials(tr).LFP,2) == nprobes
            cur_LFPs = double(int2double(AllExpt.Expt.Trials(tr).LFP(:,uprobes),AllExpt.Expt.Header.lfpscale));
            bad_LFPs = isnan(cur_LFPs(:,1));
            cur_LFPs(isnan(cur_LFPs)) = 0;
            cur_LFPs = filtfilt(b_aa,a_aa,cur_LFPs);
            cur_LFPs = downsample(cur_LFPs,LFP_dsf);
            bad_LFPs = downsample(bad_LFPs,LFP_dsf);
            
            %         cur_cwt = nan(size(cur_LFPs,1),length(wfreqs),length(uprobes));
            %         for cc = 1:length(uprobes)
            %             cur_cwt(:,:,cc) = cwt(cur_LFPs(:,cc),scales,wavetype)';
            %         end
            
            cur_LFPs(bad_LFPs,:) = nan;
            trial_LFPs(tr,:,:) = cur_LFPs;
            %         cur_cwt(bad_LFPs,:,:) = nan;
            %         trial_LFP_real(tr,:,:,:) = interp1(LFP_trial_taxis_ds,real(cur_cwt),R_trial_taxis);
            %         trial_LFP_imag(tr,:,:,:) = interp1(LFP_trial_taxis_ds,imag(cur_cwt),R_trial_taxis);
            %         trial_LFPs(tr,:,:) = interp1(LFP_trial_taxis_ds,cur_LFPs,R_trial_taxis);
            end
        end
        %     trial_LFP_real = permute(trial_LFP_real,[2 1 3 4]);
        %     trial_LFP_imag = permute(trial_LFP_imag,[2 1 3 4]);
        trial_LFPs = permute(trial_LFPs,[2 1 3]);
        
        %     trial_LFP_real = reshape(trial_LFP_real,[],length(wfreqs)*length(uprobes));
        %     trial_LFP_imag = reshape(trial_LFP_imag,[],length(wfreqs)*length(uprobes));
        
        %%
        trial_LFPs = nanzscore(reshape(trial_LFPs,[],length(uprobes)));
        trial_LFPs = reshape(trial_LFPs,[],Ntrials*length(uprobes));
        
        %%
        if sum(~isnan(trial_LFPs(:))) > 0
        beg_buff = 0.25;
        end_buff = 0.25;
        trial_dur = 2;
        uinds = find(LFP_trial_taxis_ds >= beg_buff & (trial_dur - LFP_trial_taxis_ds) >= end_buff);
        
        params.tapers = [10 19];
        params.Fs = LFP_Fsd;
        params.fpass = [aa_lcf aa_hcf];
        
        [S,f] = mtspectrumc(trial_LFPs(uinds,:),params);
        trial_Pspec = reshape(log10(S),length(f),Ntrials,length(uprobes));
        trial_Pspec = permute(trial_Pspec,[2 1 3]);
        trial_Pspec = nanzscore(trial_Pspec);
        
        avg_sigUp = squeeze(nanmean(trial_Pspec(stim_up,:,:)));
        avg_sigDown = squeeze(nanmean(trial_Pspec(stim_down,:,:)));
        avg_chUp = squeeze(nanmean(trial_Pspec(resp_up,:,:)));
        avg_chDown = squeeze(nanmean(trial_Pspec(resp_down,:,:)));
        
        sig_diffspec = avg_sigUp - avg_sigDown;
        ch_diffspec = avg_chUp - avg_chDown;
        %%
        f1 = figure('visible','off');hold on
        plot(f,mean(sig_diffspec,2),'b'); plot(f,mean(ch_diffspec,2),'r')
        legend('Signal based','Choice based');
        shadedErrorBar(f,mean(sig_diffspec,2),std(sig_diffspec,[],2)/sqrt(length(uprobes)),{'color','b'});
        shadedErrorBar(f,mean(avg_chUp-avg_chDown,2),std(avg_chUp-avg_chDown,[],2)/sqrt(length(uprobes)),{'color','r'});
        set(gca,'xscale','log')
        axis tight
        line(params.fpass,[0 0],'color','k');
        
        fig_width = 3; rel_height = 0.7;
        fig_name = [fig_dir sprintf('%s_LFP_diffspec.pdf',Expt_name)];
        figufy(f1);
        exportfig(f1,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
        close(f1);
        
        %%
        avg_trial_Pspec = squeeze(nanmean(trial_Pspec,3));
        
        Pspec_CP = nan(length(f),1);
        for ii = 1:length(f)
            [~,~,~,Pspec_CP(ii)]=perfcurve(trialrespDir(test_trials), avg_trial_Pspec(test_trials,ii),1);
        end
        
        %         trial_Pspec_resh = reshape(trial_Pspec,Ntrials,[]);
        C = corr(avg_trial_Pspec,avg_spks);
        f1 = figure('visible','off');
        [ax,h1,h2] = plotyy(f,Pspec_CP,f,C);
        set(get(ax(1),'ylabel'),'String','CP');
        set(get(ax(2),'ylabel'),'String','Correlation');
        set(ax(1),'xscale','log');
        set(ax(2),'xscale','log');
        set(ax(1),'ylim',[0.35 0.65]);
        set(ax(2),'ylim',[-0.7 0.7])
        set(get(ax(1),'xlabel'),'String','Frequency (Hz)');
        line(params.fpass,[0 0],'color','k');
        
        fig_width = 3; rel_height = 0.7;
        fig_name = [fig_dir sprintf('%s_LFP_CP.pdf',Expt_name)];
        figufy(f1);
        exportfig(f1,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
        close(f1);
        end
    end
end

