%% LOAD PROCESSED DATA
clear all
close all

data_dir = '/media/NTlab_data3/Data/bruce/CPdata/';
fig_dir = '/home/james/Analysis/bruce/ChoiceProb/';

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
    for tr = 1:Ntrials
        for nn = 1:Nunits
            trindx = find( AllExpt.Spikes{nn}.Trial == trialNums(tr));
            if ~isempty(trindx)
                tot_spks_per_trial(tr,nn) = sum(AllExpt.Spikes{nn}.Spikes{trindx}*1e-4 <= (trialDur + spk_delay_buff) & ...
                    (AllExpt.Spikes{nn}.Spikes{trindx}*1e-4 >= spk_delay_buff));
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
    
    sig_trials = find(abs(trialOB) == msOB);
    stim_up = sig_trials(trialrwDir(sig_trials)==1);
    stim_down = sig_trials(trialrwDir(sig_trials)==-1);
    sig_prob = nan(Nunits,1);
    sp_pval = nan(Nunits,1);
    for cc = 1:Nunits
        cur_resp_ax = 0:(nanmax(tot_spks_per_trial_norm(sig_trials,cc)));
        if ~isnan(cur_resp_ax)
        up_hist = histc(tot_spks_per_trial_norm(stim_up,cc),cur_resp_ax);
        down_hist = histc(tot_spks_per_trial_norm(stim_down,cc),cur_resp_ax);
        fract_cor = cumsum(up_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_up,cc)));
        fract_incor = cumsum(down_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_down,cc)));
        sig_prob(cc) = trapz(fract_incor,fract_cor);
        [~,sp_pval(cc)] = ttest2(tot_spks_per_trial_norm(stim_up,cc),tot_spks_per_trial_norm(stim_down,cc));
        end
    end
    
    test_trials = find(trialOB == zsOB & trialrespDir ~= 0);
    % test_trials = find(trialrespDir ~= 0);
    resp_up = test_trials(trialrespDir(test_trials) == 1);
    resp_down = test_trials(trialrespDir(test_trials) == -1);
    choice_prob = nan(Nunits,1);
    cp_pval = nan(Nunits,1);
    for cc = 1:Nunits
        cur_resp_ax = 0:(nanmax(tot_spks_per_trial_norm(test_trials,cc)));
        if ~isnan(cur_resp_ax)
        up_hist = histc(tot_spks_per_trial_norm(resp_up,cc),cur_resp_ax);
        down_hist = histc(tot_spks_per_trial_norm(resp_down,cc),cur_resp_ax);
        true_pos = cumsum(up_hist)/sum(~isnan(tot_spks_per_trial_norm(resp_up,cc)));
        false_pos = cumsum(down_hist)/sum(~isnan(tot_spks_per_trial_norm(resp_down,cc)));
        choice_prob(cc) = trapz(false_pos,true_pos);
        [~,~,~,CP_Neuron(cc)]=perfcurve(trialrespDir(test_trials), tot_spks_per_trial_norm(test_trials,cc),1);
        [~,cp_pval(cc)] = ttest2(tot_spks_per_trial_norm(resp_up,cc),tot_spks_per_trial_norm(resp_down,cc));
        end
    end
    
    %%
    xr = [0.3 0.7]; yr = [0 1];
    h = scatterhist(choice_prob,sig_prob,'nbins',50);
    hold on
    plot(choice_prob(Ucellnumbers > 0),sig_prob(Ucellnumbers > 0),'r.');
    line(xr,[0.5 0.5],'color','k'); line([0.5 0.5],yr,'color','k');
    xlim(xr); ylim(yr);
    title(sprintf('ori: %d',unique(trialOR)));
    
    fig_width = 6; rel_height = 1;
    fig_name = [fig_dir sprintf('%s_CPSP_hist.pdf',Expt_name)];
    figufy(gcf);
    exportfig(gcf,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    close(gcf);
    
    %%
    un_OBs = [-60 -70 -80 msOB 80 70 60];
    OB_x = [-3 -2 -1 0 1 2 3];
    trialOB_prime = trialOB;
    if min(trialOB) > 0
        trialOB_prime(trialrwDir==-1) = -trialOB_prime(trialrwDir==-1);
    end
    
    fract_one = nan(length(un_OBs),1);
    for ss = 1:length(un_OBs)
        cur_OB_neg = find(trialOB_prime == un_OBs(ss) & trialrespDir ~= 0);
        fract_one(ss) = sum(trialrespDir(cur_OB_neg) == 1)/length(cur_OB_neg);
    end
    
    sig_strengths = [msOB 80 70 60];
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
    
    f1 = figure();
    subplot(2,2,1);
    plot(OB_x,fract_one,'o-');
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
    ET_data_up = double(cat(3,AllExpt.Expt.Trials(resp_up_trials).EyeData));
    ET_data_down = double(cat(3,AllExpt.Expt.Trials(resp_down_trials).EyeData));
    trange = (size(ET_data_up,1)-110):(size(ET_data_up,1)-10);
    ET_data_up = permute(ET_data_up,[1 3 2]); ET_data_down = permute(ET_data_down,[1 3 2]);
    
    ET_data_up = bsxfun(@minus,ET_data_up,reshape(nanmedian(reshape(ET_data_up,[],4)),1,1,4));
    ET_data_down = bsxfun(@minus,ET_data_down,reshape(nanmedian(reshape(ET_data_down,[],4)),1,1,4));

    subplot(2,2,4); hold on
    plot(squeeze(ET_data_up(trange,:,3)),squeeze(ET_data_up(trange,:,4)),'b.','markersize',1);
    plot(squeeze(ET_data_down(trange,:,3)),squeeze(ET_data_down(trange,:,4)),'r.','markersize',1);
    xlim([-2.5 2.5]*1e4); ylim([-2.5 2.5]*1e4); line([-2.5 2.5]*1e4,[0 0],'color','k'); line([0 0],[-2.5 2.5]*1e4,'color','k');
    fig_width = 8; rel_height = 1;
    fig_name = [fig_dir sprintf('%s_ratecurves.pdf',Expt_name)];
    figufy(f1);
    exportfig(f1,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
    close(f1);
    
    %%
    
end

