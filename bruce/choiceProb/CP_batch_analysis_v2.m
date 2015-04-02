%% LOAD PROCESSED DATA
clear all
close all

data_dir = '/media/NTlab_data3/Data/bruce/CPdata/';
fig_dir = '/home/james/Analysis/bruce/ChoiceProb/';
addpath('~/James_scripts/bruce/bruce_code');

addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))

data_sets = what(data_dir);data_sets = data_sets.mat;

%%
for eee = 1:length(data_sets)
    % Expt_name = 'G050';
    Expt_name = data_sets{eee};
    
    load(strcat(data_dir, Expt_name));
    Expt_name = Expt_name(1:end-4);
    
    CPDATA(eee).expt_name = Expt_name;
    
    %% Assemble stim, Robs, LFP into trial structure
    Ntrials = length(AllExpt.Expt.Trials);
    dt = 0.05;  % assuming 10 ms frame time exactly
    trialDur = 2; %in sec.
    Nframes = 200; %number of video frames (@100Hz) per trial
    NT = ceil(trialDur/dt); %number of time bins per trial
    
    BlockStartid = AllExpt.Expt.Header.BlockStartid;
    n_blocks = length(BlockStartid);
    
    %% SU and MU data (combined)
    Nunits = length(AllExpt.Header);
    Ucellnumbers = [AllExpt.Header(:).cellnumber];
    Uchans = round([AllExpt.Header(:).probe]);
    
    CPDATA(eee).cellnums = Ucellnumbers;
    
    %% relevant trial variables
    trialNums = [AllExpt.Expt.Trials(:).Trial];
    trialOR = [AllExpt.Expt.Trials(:).or];
    trialOB = [AllExpt.Expt.Trials(:).ob];
    trialrwDir = [AllExpt.Expt.Trials(:).rwdir];
    trialrespDir = [AllExpt.Expt.Trials(:).RespDir];
    trialSe = [AllExpt.Expt.Trials(:).se];
    trialId = [AllExpt.Expt.Trials(:).id];
    
    zsOB = max(abs(trialOB)); %zero-signal condition
    msOB = min(abs(trialOB)); %max signal condition
    
    %%
        
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
    spk_delay_buff = 0.0;
    spk_end_exclude = 0;
    spk_beg_exclude = 0.;
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
    
    trialSpkCnts = squeeze(nanmean(fullRobs,1))*NT;
    totSpkCnts = squeeze(nansum(trialSpkCnts,1));
    avg_spk_rates = nanmean(reshape(fullRobs,[],Nunits));
    
    %% CALCULATE CHOICE PROBS
    sig_trials = find(abs(trialOB) == msOB & trialrespDir ~= 0); %max signal trials
    stim_up = sig_trials(trialrwDir(sig_trials)==1);
    stim_down = sig_trials(trialrwDir(sig_trials)==-1);
    test_trials = find(trialOB == zsOB & trialrespDir ~= 0); %zero signal trials
    resp_up = test_trials(trialrespDir(test_trials) == 1);
    resp_down = test_trials(trialrespDir(test_trials) == -1);
    
    
    %     sig_prob = nan(Nunits,1);
    %     choice_prob = nan(Nunits,1);
    %     for cc = 1:Nunits
    %         if nansum(trialSpkCnts(sig_trials,cc)) > 0
    %             [~,~,~,sig_prob(cc)]=perfcurve(trialrwDir(sig_trials), trialSpkCnts(sig_trials,cc),1);
    %         end
    %         if nansum(trialSpkCnts(test_trials,cc)) > 0
    %             [~,~,~,choice_prob(cc)]=perfcurve(trialrespDir(test_trials), trialSpkCnts(test_trials,cc),1);
    %         end
    %     end
    
    poss_win_exclude = [0 0.25 0.5 0.75 1 1.25 1.5];
    [sig_prob_EX,choice_prob_EX,sig_prob_BX,choice_prob_BX] = deal(nan(length(poss_win_exclude),Nunits));
    for ww = 1:length(poss_win_exclude)
        fprintf('Excluding window %d of %d\n',ww,length(poss_win_exclude));
        end_buff = round(poss_win_exclude(ww)/dt);
        cur_inds = 1:(NT-end_buff);
        cur_trial_cnts = squeeze(nanmean(fullRobs(cur_inds,:,:),1))*length(cur_inds);
        for cc = 1:Nunits
            if nansum(cur_trial_cnts(sig_trials,cc)) > 0
                [~,~,~,sig_prob_EX(ww,cc)]=perfcurve(trialrwDir(sig_trials), cur_trial_cnts(sig_trials,cc),1);
            end
            if nansum(cur_trial_cnts(test_trials,cc)) > 0
                [~,~,~,choice_prob_EX(ww,cc)]=perfcurve(trialrespDir(test_trials), cur_trial_cnts(test_trials,cc),1);
            end
        end
        beg_buff = round(poss_win_exclude(ww)/dt);
        cur_inds = (1+beg_buff):NT;
        cur_trial_cnts = squeeze(nanmean(fullRobs(cur_inds,:,:),1))*length(cur_inds);
        for cc = 1:Nunits
            if nansum(cur_trial_cnts(sig_trials,cc)) > 0
                [~,~,~,sig_prob_BX(ww,cc)]=perfcurve(trialrwDir(sig_trials), cur_trial_cnts(sig_trials,cc),1);
            end
            if nansum(cur_trial_cnts(test_trials,cc)) > 0
                [~,~,~,choice_prob_BX(ww,cc)]=perfcurve(trialrespDir(test_trials), cur_trial_cnts(test_trials,cc),1);
            end
        end
    end
    
    CPDATA(eee).sig_prob_EX = sig_prob_EX;
    CPDATA(eee).sig_prob_BX = sig_prob_BX;
    CPDATA(eee).choice_prob_EX = choice_prob_EX;
    CPDATA(eee).choice_prob_BX = choice_prob_BX;
    CPDATA(eee).poss_win_exclude = poss_win_exclude;
    
    %%
    norm_fullRobs = reshape(nanzscore(reshape(fullRobs,[],Nunits)),NT,Ntrials,Nunits);
    avg_Robs = squeeze(nanmean(norm_fullRobs,3));
    
    [timedep_sp,timedep_cp] = deal(nan(NT,1));
    for ii = 1:NT
        [~,~,~,timedep_sp(ii)]=perfcurve(trialrwDir(sig_trials), avg_Robs(ii,sig_trials)',1);
        [~,~,~,timedep_cp(ii)]=perfcurve(trialrespDir(test_trials), avg_Robs(ii,test_trials)',1);
    end
    
    avg_trial = nanmean(avg_Robs,1);
    [~,~,~,avg_SP]=perfcurve(trialrwDir(sig_trials), avg_trial(sig_trials),1);
    [~,~,~,avg_CP]=perfcurve(trialrespDir(test_trials), avg_trial(test_trials),1);
    
    CPDATA(eee).timedep_avg_CP = timedep_cp;
    CPDATA(eee).timedep_avg_SP = timedep_sp;
    CPDATA(eee).timedep_TAX = (1:NT)*dt;
    CPDATA(eee).avg_SP = avg_SP;
    CPDATA(eee).avg_CP = avg_CP;
    
    %%
    
%     f1 = figure('visible','off');
%     subplot(1,2,1);hold on
%     errorbar(poss_win_exclude,nanmean(sig_prob_EX,2),nanstd(sig_prob_EX,[],2)/sqrt(Nunits));
%     errorbar(poss_win_exclude,nanmean(sig_prob_BX,2),nanstd(sig_prob_BX,[],2)/sqrt(Nunits),'r');
%     errorbar(poss_win_exclude,nanmean(choice_prob_EX,2),nanstd(choice_prob_EX,[],2)/sqrt(Nunits),'k');
%     errorbar(poss_win_exclude,nanmean(choice_prob_BX,2),nanstd(choice_prob_BX,[],2)/sqrt(Nunits),'g');
%     legend('Sig EX','Sig BX','Choice EX','Choice BX');
%     
%     subplot(1,2,2);hold on
%     plot((1:NT)*dt,timedep_sp,'b',(1:NT)*dt,timedep_cp,'r');
%     legend('Sig prob','Choice prob');
%     line([0 2],[0.5 0.5],'color','k');
%     
%     fig_width = 10; rel_height = 0.4;
%     fig_name = [fig_dir sprintf('%s_CPtimedep.pdf',Expt_name)];
%     figufy(f1);
%     exportfig(f1,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
%     close(f1);
    
    %%
    un_OBs = [-60 -70 -80 zsOB 80 70 60];
    OB_x = [-3 -2 -1 0 1 2 3]; %labels of different signal conditions for plotting
    trialOB_prime = trialOB;
    if min(trialOB) > 0 %if OB values arent stored as negative, convert them
        trialOB_prime(trialrwDir==-1) = -trialOB_prime(trialrwDir==-1);
    end
    
    fract_one = nan(length(un_OBs),1);
    fract_one_ci = nan(length(un_OBs),2);
    zval = 1.96;
    for ss = 1:length(un_OBs)
        cur_OB_neg = find(trialOB_prime == un_OBs(ss) & trialrespDir ~= 0);
        curN = length(cur_OB_neg);
        fract_one(ss) = sum(trialrespDir(cur_OB_neg) == 1)/curN;
        
        %uncertainty estimates based on Wilson score interval formula (this
        %gives better approx of 95% CI
        fract_one_ci(ss,2) = 1/(1+zval^2/curN)*(fract_one(ss) + zval^2/(2*curN) + zval*sqrt(fract_one(ss)/curN*(1-fract_one(ss)) + zval^2/(4*curN^2)));
        fract_one_ci(ss,1) = 1/(1+zval^2/curN)*(fract_one(ss) + zval^2/(2*curN) - zval*sqrt(fract_one(ss)/curN*(1-fract_one(ss)) + zval^2/(4*curN^2)));
    end
    
    CPDATA(eee).psych_fract_one = fract_one;
    CPDATA(eee).psych_fract_one_ci = fract_one_ci;
    
end

%%
close all
ex_in = 5;

pool_SU_abs_CPs = [];
pool_SU_rel_CPs = [];
pool_SU_CPs = [];
pool_SU_SPs = [];
pool_SU_an = [];
for ii = 1:length(CPDATA)
    if data_sets{ii}(1) == 'G'
        dtype(ii) = 1;
    else
        dtype(ii) = 0;
    end
    
    sign_flip = -sign((CPDATA(ii).psych_fract_one(1)-0.5));
    all_sfs(ii) = sign_flip;
    
    cur_SUs = CPDATA(ii).cellnums > 0;
    cur_CPs = CPDATA(ii).choice_prob_BX(ex_in,:);
    cur_SPs = CPDATA(ii).sig_prob_BX(ex_in,:);
    
    cur_CPs = (cur_CPs-0.5)*sign_flip + 0.5;
    cur_SPs = (cur_SPs - 0.5)*sign_flip + 0.5;
    
    avg_CP = CPDATA(ii).avg_CP;
    avg_SP = CPDATA(ii).avg_SP;
    
    avg_CP = (avg_CP - 0.5)*sign_flip + 0.5;
    avg_SP = (avg_SP - 0.5)*sign_flip + 0.5;
    
    timedep_avg_CP = CPDATA(ii).timedep_avg_CP;
    timedep_avg_CP = (timedep_avg_CP - 0.5)*sign_flip + 0.5;
    all_timedep_CP(ii,:) = timedep_avg_CP;
    
    all_avg_CPs(ii) = avg_CP;
    all_avg_SPs(ii) = avg_SP;
    
    
    avg_abs_CP(ii) = mean(cur_CPs-0.5)*sign(avg_CP-0.5);
    avg_rel_CP(ii) = mean((cur_CPs-0.5).*sign(cur_SPs-0.5));
    avg_rel_CP(ii) = mean((cur_CPs-0.5).*sign(cur_SPs-0.5));
    
    cur_CPs_abs = (cur_CPs(cur_SUs)-0.5)*sign(avg_CP-0.5);
    cur_CPs_rel = (cur_CPs(cur_SUs)-0.5).*sign(cur_SPs(cur_SUs)-0.5);
        
    pool_SU_an = cat(1,pool_SU_an,ones(sum(cur_SUs),1)*dtype(ii));
    pool_SU_abs_CPs = cat(1,pool_SU_abs_CPs,cur_CPs_abs(cur_SUs)');
    pool_SU_rel_CPs = cat(1,pool_SU_rel_CPs,cur_CPs_rel(cur_SUs)');
    pool_SU_CPs = cat(1,pool_SU_CPs,cur_CPs(cur_SUs)');
    pool_SU_SPs = cat(1,pool_SU_SPs,cur_SPs(cur_SUs)');
end

f1 = figure(); hold on
plot(avg_abs_CP(dtype == 1)+0.5,avg_rel_CP(dtype==1)+0.5,'o');
plot(avg_abs_CP(dtype == 0)+0.5,avg_rel_CP(dtype==0)+0.5,'ro');
legend('JB','LEM');
line([-0.05 0.15]+0.5,[-0.05 0.15]+0.5,'color','k')
line([-0.05 0.15]+0.5,[0 0]+0.5,'color','k','linestyle','--');
line([0 0]+0.5,[-0.05 0.15]+0.5,'color','k','linestyle','--');
axis tight
xlabel('Absolute CP');
ylabel('Relative CP');



f2 = figure(); hold on
plot(all_avg_CPs(dtype==1),all_avg_SPs(dtype==1),'o');
plot(all_avg_CPs(dtype==0),all_avg_SPs(dtype==0),'ro');
legend('JB','LEM');
xlim([0 1]); xl = xlim();
ylim([0 1]); yl  = ylim();
line(xl,[0.5 0.5],'color','k');
line([0.5 0.5],yl,'color','k');
xlabel('Avg CP');
ylabel('Avg SP');


% f2 = figure();
% hold on
% plot(pool_SU_abs_CPs(pool_SU_an==1),pool_SU_rel_CPs(pool_SU_an==1),'o');
% plot(pool_SU_abs_CPs(pool_SU_an==0),pool_SU_rel_CPs(pool_SU_an==0),'ro');

% f3 = figure();
% hold on
% plot(pool_SU_CPs(pool_SU_an==1),pool_SU_SPs(pool_SU_an==1),'o');
% plot(pool_SU_CPs(pool_SU_an==0),pool_SU_SPs(pool_SU_an==0),'ro');


f3 = figure(); 
xe = linspace(0.2,0.8,100);
n1 = histc(pool_SU_rel_CPs+0.5,xe);
n2 = histc(pool_SU_abs_CPs+0.5,xe);
subplot(3,1,1);
hold on
plot(xe,cumsum(n1)/sum(n1))
plot(xe,cumsum(n2)/sum(n2),'r');
xlim([0.25 0.75]);
ylim([0 1]);
line([0.5 0.5],[0 1],'color','k');
line([0 1],[0.5 0.5],'color','k','linestyle','--');
title('Overall');
xlabel('CP');
ylabel('Cumulative prob');
legend('Relative','Absolute');

n1 = histc(pool_SU_rel_CPs(pool_SU_an == 1)+0.5,xe);
n2 = histc(pool_SU_abs_CPs(pool_SU_an==1)+0.5,xe);
subplot(3,1,2);
hold on
plot(xe,cumsum(n1)/sum(n1))
plot(xe,cumsum(n2)/sum(n2),'r');
xlim([0.25 0.75]);
ylim([0 1]);
line([0.5 0.5],[0 1],'color','k');
line([0 1],[0.5 0.5],'color','k','linestyle','--');
title('JB');
xlabel('CP');
ylabel('Cumulative prob');
legend('Relative','Absolute');

n1 = histc(pool_SU_rel_CPs(pool_SU_an == 0)+0.5,xe);
n2 = histc(pool_SU_abs_CPs(pool_SU_an== 0)+0.5,xe);
subplot(3,1,3);
hold on
plot(xe,cumsum(n1)/sum(n1))
plot(xe,cumsum(n2)/sum(n2),'r');
xlim([0.25 0.75]);
ylim([0 1]);
line([0.5 0.5],[0 1],'color','k');
line([0 1],[0.5 0.5],'color','k','linestyle','--');
title('LEM');
xlabel('CP');
ylabel('Cumulative prob');
legend('Relative','Absolute');


f4 = figure(); hold on
h1 = shadedErrorBar(CPDATA(1).timedep_TAX,mean(all_timedep_CP(dtype==1,:)),std(all_timedep_CP(dtype==1,:))/sqrt(sum(dtype==1)),{'color','b'});
h2 = shadedErrorBar(CPDATA(1).timedep_TAX,mean(all_timedep_CP(dtype==0,:)),std(all_timedep_CP(dtype==0,:))/sqrt(sum(dtype==0)),{'color','r'});
line([0 2],[0.5 0.5],'color','k');
xlabel('Time (s)');
ylabel('Avg CP');
legend([h1.mainLine h2.mainLine],{'JB','LEM'});

% fig_width = 5; rel_height = 0.8;
% fig_name = [fig_dir 'Avg_unit_CPSP.pdf'];
% figufy(f1);
% exportfig(f1,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% fig_width = 5; rel_height = 0.8;
% fig_name = [fig_dir 'Avg_rec_CPSP.pdf'];
% figufy(f2);
% exportfig(f2,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
% fig_width = 5; rel_height = 2.4;
% fig_name = [fig_dir 'Avg_unit_CP_dists.pdf'];
% figufy(f3);
% exportfig(f3,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);
% 
% fig_width = 5; rel_height = 0.8;
% fig_name = [fig_dir 'Avg_CP_timecourse.pdf'];
% figufy(f4);
% exportfig(f4,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);
