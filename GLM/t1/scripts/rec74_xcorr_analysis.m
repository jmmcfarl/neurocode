clear all; close all;
datdir = '~/Data/blanche/rec_74';
cd(datdir)
cd matlabdata/
dt   = 19.9920031987/1000;
sdim = 64;

load spksRec74.mat
load dstimpsRec74.mat
load spksegsRec74.mat
reps_per_seq = 24;
rep_len = 250;
rep_inds = floor(((1:reps_per_seq*rep_len)-1)/rep_len);
rep_times = rep_inds*rep_len*dt;
nrand = 100;

rstimfiles = {'ara-64p-50h-2m-mn125-ct015.s0', 'wra-64p-50h-2m-mn125-ct015.s1', 'pra-64p-50h-2m-mn125-ct015.s2', 'nra-64p-50h-2m-mn125-ct015.s3', ...
    'nra-64p-50h-2m-mn125-ct025.s4', 'ara-64p-50h-2m-mn125-ct025.s5', 'wra-64p-50h-2m-mn125-ct025.s6', 'pra-64p-50h-2m-mn125-ct025.s7',...
    'pra-64p-50h-2m-mn125-ct035.s8', 'ara-64p-50h-2m-mn125-ct035.s9', 'nra-64p-50h-2m-mn125-ct035.s10', 'wra-64p-50h-2m-mn125-ct035.s11',...
    'nra-64p-50h-2m-mn125-ct045.s12', 'ara-64p-50h-2m-mn125-ct045.s13', 'pra-64p-50h-2m-mn125-ct045.s14', 'wra-64p-50h-2m-mn125-ct045.s15'};

stimfiles  = cellfun(@(x)x(1:26),rstimfiles,'UniformOutput',0);
nstims = length(stimfiles);
used_conds = [5 8 9 11 13 15];

ncells = length(spks74);
cell_ids = 1:ncells;
cell_ids(9:end) = cell_ids(9:end) + 1; %make up for missing cell 9 in this data set relative to rec75/76

%% DEFINE XCORR PARAMETERS
tmax = 120; %2 min segs
rep_axis = 1:reps_per_seq;

%for low precision xcorrs
ltp_timeres = 0.0025; %in sec
ltp_tclen   = 0.125;
ltp_ncbins = ltp_tclen/ltp_timeres;
ltp_bin_ax = 0:ltp_timeres:tmax;
ltp_Nbins = length(ltp_bin_ax);
ltp_t_ax = (-ltp_ncbins:ltp_ncbins)*ltp_timeres;
ltp_cond_avg_ccs = cell(ncells,ncells);

%for medium precision xcorrs
mtp_timeres = 0.0005; %in sec
mtp_tclen   = 0.01;
mtp_ncbins = mtp_tclen/mtp_timeres;
mtp_bin_ax = 0:mtp_timeres:tmax;
mtp_Nbins = length(mtp_bin_ax);
mtp_t_ax = (-mtp_ncbins:mtp_ncbins)*mtp_timeres;
mtp_cond_avg_ccs = cell(ncells,ncells);

%for precise xcorrs
htp_timeres = 0.0001; %in sec
htp_tclen   = 0.005;
htp_ncbins = htp_tclen/htp_timeres;
htp_bin_ax = 0:htp_timeres:tmax;
htp_Nbins = length(htp_bin_ax);
htp_t_ax = (-htp_ncbins:htp_ncbins)*htp_timeres;
temp_t_ax = ((-htp_ncbins-1):(htp_ncbins+1))*htp_timeres;
htp_cond_avg_ccs = cell(ncells,ncells);

%% XCORR ANALYSIS
ncells = length(spks74);
cd ~/James_scripts/GLM/t1/rec74_xcorrs/
ov_rep_inds_ltp = ceil(((1:length(ltp_bin_ax)))/rep_len/dt*ltp_timeres);
ov_rep_inds_ltp(ov_rep_inds_ltp > 24) = 24;
trial_t_axis_ltp = 0:ltp_timeres:5;
ov_rep_inds_mtp = ceil(((1:length(mtp_bin_ax)))/rep_len/dt*mtp_timeres);
ov_rep_inds_mtp(ov_rep_inds_mtp > 24) = 24;
trial_t_axis_mtp = 0:mtp_timeres:5;
for icell1 = 1:(ncells-1);
    % icell1 = 11;
    for icell2 = (icell1+1):ncells;
        fprintf('C%d_C%d\n',icell1,icell2);
        % icell2 = 13;
        psth_xcorr_ltp = zeros(length(used_conds),2*ltp_ncbins+1);
        emp_xcorr_ltp = zeros(length(used_conds),2*ltp_ncbins+1);
        emp_xcorr_repcorr_ltp = zeros(length(used_conds),2*ltp_ncbins+1);
        psth_xcorr_mtp = zeros(length(used_conds),2*mtp_ncbins+1);
        emp_xcorr_mtp = zeros(length(used_conds),2*mtp_ncbins+1);
        emp_xcorr_repcorr_mtp = zeros(length(used_conds),2*mtp_ncbins+1);
        emp_xcorr_shuff_mtp = zeros(length(used_conds),nrand,2*mtp_ncbins+1);
        for n = 1:length(used_conds)
            fprintf('Cond %d of %d\n',n,length(used_conds));
            spktimes1 = spksegs74{icell1,used_conds(n)};
            spktimes2 = spksegs74{icell2,used_conds(n)};
            
            %compute xcorr of avg rates across repeats
            spk_rep_inds1 = ceil(spktimes1/rep_len/dt);
            spk_rep_inds2 = ceil(spktimes2/rep_len/dt);            
            rep_spk_rates1 = hist(spk_rep_inds1,rep_axis)/5;
            rep_spk_rates2 = hist(spk_rep_inds2,rep_axis)/5;            
            rep_xcorr{icell1,icell2}(n,:) = xcov(rep_spk_rates1,rep_spk_rates2,[],'coeff');
  
            
            %for Low temporal precision
            strain1 = histc(spktimes1,ltp_bin_ax)'/ltp_timeres;
            strain2 = histc(spktimes2,ltp_bin_ax)'/ltp_timeres;            
            %correct for changes in avg rate across repeats
            strain1_corr = strain1 - rep_spk_rates1(ov_rep_inds_ltp);
            strain2_corr = strain2 - rep_spk_rates2(ov_rep_inds_ltp);           
            emp_xcorr_ltp(n,:) = xcov(strain1,strain2,ltp_ncbins,'coeff'); %plot(allccs{ipair}); drawnow
            emp_xcorr_repcorr_ltp(n,:) = xcov(strain1_corr,strain2_corr,ltp_ncbins,'coeff'); %plot(allccs{ipair}); drawnow            

            spktimes1_trialsub = spktimes1 - (spk_rep_inds1-1)*rep_len*dt;
            spktimes2_trialsub = spktimes2 - (spk_rep_inds2-1)*rep_len*dt;
            psth1 = hist(spktimes1_trialsub,trial_t_axis_ltp)/ltp_timeres/reps_per_seq;
            psth2 = hist(spktimes2_trialsub,trial_t_axis_ltp)/ltp_timeres/reps_per_seq;
            psth_xcorr_ltp(n,:) = xcov(psth1,psth2,ltp_ncbins,'coeff');
   
            
            %for Med temporal precision
            strain1 = histc(spktimes1,mtp_bin_ax)'/mtp_timeres;
            strain2 = histc(spktimes2,mtp_bin_ax)'/mtp_timeres;            
            %correct for changes in avg rate across repeats
            strain1_corr = strain1 - rep_spk_rates1(ov_rep_inds_mtp);
            strain2_corr = strain2 - rep_spk_rates2(ov_rep_inds_mtp);           
            emp_xcorr_mtp(n,:) = xcov(strain1,strain2,mtp_ncbins,'coeff'); %plot(allccs{ipair}); drawnow
            emp_xcorr_repcorr_mtp(n,:) = xcov(strain1_corr,strain2_corr,mtp_ncbins,'coeff'); %plot(allccs{ipair}); drawnow            

            spktimes1_trialsub = spktimes1 - (spk_rep_inds1-1)*rep_len*dt;
            spktimes2_trialsub = spktimes2 - (spk_rep_inds2-1)*rep_len*dt;
            psth1 = hist(spktimes1_trialsub,trial_t_axis_mtp)/mtp_timeres/reps_per_seq;
            psth2 = hist(spktimes2_trialsub,trial_t_axis_mtp)/mtp_timeres/reps_per_seq;
            psth_xcorr_mtp(n,:) = xcov(psth1,psth2,mtp_ncbins,'coeff');

            for r = 1:nrand
                rand_trials = randperm(reps_per_seq);
                rand_rep_inds = rand_trials(ov_rep_inds_mtp);
                spktimes2_trialshuf = spktimes2_trialsub + (rand_trials(spk_rep_inds2)'-1)*rep_len*dt;
                strain2 = histc(spktimes2_trialshuf,mtp_bin_ax)'/mtp_timeres;
                strain2_corr = strain2 - rep_spk_rates2(rand_rep_inds);
                emp_xcorr_shuff_mtp(n,r,:) = xcov(strain1_corr,strain2_corr,mtp_ncbins,'coeff'); %plot(allccs{ipair}); drawnow
            end
           
        end
        
        htp_cond_avg_ccs{icell1,icell2} = zeros(1,length(htp_t_ax));
        ucond_cnt = 0;
        for n = 1:length(used_conds)
            spktimes1 = spksegs74{icell1,used_conds(n)};
            spktimes2 = spksegs74{icell2,used_conds(n)};
            
            %FASTER FOR BINNING < ~0.5ms
            if length(spktimes1) > 3 & length(spktimes2) > 3
            cur_xcov = zeros(1,length(temp_t_ax));
            for i = 1:length(spktimes1)
                cur_xcov = cur_xcov + hist(spktimes2-spktimes1(i),temp_t_ax);
            end
            cur_xcov([1 end]) = [];
            cur_xcov = (cur_xcov - length(spktimes1)*length(spktimes2)/...
                htp_Nbins)/sqrt((length(spktimes1) - length(spktimes1)^2/htp_Nbins)...
                *(length(spktimes2) - length(spktimes2)^2/htp_Nbins));
            
            htp_cond_avg_ccs{icell1,icell2} = htp_cond_avg_ccs{icell1,icell2} + cur_xcov;
            ucond_cnt = ucond_cnt + 1;
            end
        end
        htp_cond_avg_ccs{icell1,icell2} = htp_cond_avg_ccs{icell1,icell2}/ucond_cnt;      
        
        avg_rep_xcorr{icell1,icell2} = mean(rep_xcorr{icell1,icell2});
        emp_xcorr_shuff_avgs_mtp = squeeze(mean(emp_xcorr_shuff_mtp));       
        avg_psth_xcorr_ltp{icell1,icell2} = mean(psth_xcorr_ltp);
        avg_psth_xcorr_mtp{icell1,icell2} = mean(psth_xcorr_mtp);
        avg_xcorr_ltp{icell1,icell2} = mean(emp_xcorr_ltp);
        avg_xcorr_mtp{icell1,icell2} = mean(emp_xcorr_mtp);
        avg_xcorr_repcorr_ltp{icell1,icell2} = mean(emp_xcorr_repcorr_ltp);
        avg_xcorr_repcorr_mtp{icell1,icell2} = mean(emp_xcorr_repcorr_mtp);
         avg_xcorr_shuff_mtp{icell1,icell2} = mean(emp_xcorr_shuff_avgs_mtp);
        uci_xcorr_shuff_mtp{icell1,icell2} = prctile(emp_xcorr_shuff_avgs_mtp,95);
        lci_xcorr_shuff_mtp{icell1,icell2} = prctile(emp_xcorr_shuff_avgs_mtp,5);
        
        f2 = figure('visible','off');
        set(f2,'PaperUnits','centimeters');
        set(f2, 'PaperSize', [30 60]);
        set(f2,'PaperPosition',[0,0,(get(f2,'PaperSize'))])
        subplot(3,1,1)
        plot(1000*ltp_t_ax,psth_xcorr_ltp); hold on
        plot(1000*ltp_t_ax,avg_psth_xcorr_ltp{icell1,icell2},'k','linewidth',2)
        xlim(1000*[-0.125 0.125])
        xlabel('Time lag (ms)','fontsize',14)
        ylabel('Correlation','fontsize',14)
        title('PSTH xcorr (2.5ms bins)','fontsize',14)
        subplot(3,1,2)
        plot(1000*mtp_t_ax,emp_xcorr_repcorr_mtp); hold on
        plot(1000*mtp_t_ax,avg_xcorr_repcorr_mtp{icell1,icell2},'k','linewidth',2)
        plot(1000*mtp_t_ax,uci_xcorr_shuff_mtp{icell1,icell2},'k--','linewidth',2)
        plot(1000*mtp_t_ax,lci_xcorr_shuff_mtp{icell1,icell2},'k--','linewidth',2)
        xlim(1000*[-0.01 0.01])
        xlabel('Time lag (ms)','fontsize',14)
        ylabel('Correlation','fontsize',14)
        title('Spike time xcorr (500us bins)','fontsize',14)
        subplot(3,1,3)
        plot(1000*htp_t_ax,fliplr(htp_cond_avg_ccs{icell1,icell2}),'k')
        xlim(1000*[-0.003 0.003])
        xlabel('Time lag (ms)','fontsize',14)
        ylabel('Correlation','fontsize',14)
        title('Spike time xcorr (100us bins)','fontsize',14)
        f_name = sprintf('C%d_C%d',cell_ids(icell1),cell_ids(icell2));
        print(f_name,'-dpng'); close all
       
%         smooth_cov = smooth(avg_xcorr{icell1,icell2},7,'lowess');
%         smooth_cov_repcorr = smooth(avg_xcorr_repcorr{icell1,icell2},7,'lowess');
%         f2 = figure;
%         set(f2,'PaperUnits','centimeters');
%         set(f2, 'PaperSize', [70 50]);
%         set(f2,'PaperPosition',[0,0,(get(f2,'PaperSize'))])
%         subplot(2,2,1)
%         plot(lags,smooth_cov,'b','linewidth',2); hold on
%         plot(lags,smooth_cov_repcorr,'r','linewidth',2)
%         plot(lags,avg_xcorr_shuff{icell1,icell2},'k')
%         plot(lags,uci_xcorr_shuff{icell1,icell2},'k--')
%         plot(lags,lci_xcorr_shuff{icell1,icell2},'k--')
%         xlabel('Time lag (s)','fontsize',14)
%         ylabel('Correlation','fontsize',14)
%         axis tight
%         title('Spike time xcorr','fontsize',14)
%         subplot(2,2,2)
%         psth_xcorr_sm = zeros(size(psth_xcorr));
%         for nn = 1:length(used_conds)
%             psth_xcorr_sm(nn,:) = smooth(psth_xcorr(nn,:),14,'lowess');
%         end
%         plot(lags,psth_xcorr_sm); hold on
%         plot(lags,smooth(avg_psth_xcorr{icell1,icell2},7,'lowess'),'k','linewidth',2)
%          xlabel('Time lag (s)','fontsize',14)
%         ylabel('Correlation','fontsize',14)
%         axis tight
%         title('PSTH xcorr','fontsize',14)
%        subplot(2,2,3)
%         plot(lags,avg_xcorr{icell1,icell2}); hold on
%         plot(lags,avg_xcorr_repcorr{icell1,icell2},'r')
%         plot(lags,avg_xcorr_shuff{icell1,icell2},'k','linewidth',2)
%         plot(lags,uci_xcorr_shuff{icell1,icell2},'k--')
%         plot(lags,lci_xcorr_shuff{icell1,icell2},'k--')
%         xlabel('Time lag (s)','fontsize',14)
%         ylabel('Correlation','fontsize',14)
%         xlim([-0.015 0.015])
%          title('Spike time xcorr','fontsize',14)
%        subplot(2,2,4)
%         plot(htp_t_ax,fliplr(htp_cond_avg_ccs{icell1,icell2}),'k')
%         xlim([-0.003 0.003]);
%         xlabel('Time lag (s)','fontsize',14)
%         ylabel('Correlation','fontsize',14)
%          title('Spike time xcorr (100us bins)','fontsize',14)
% %         plot(-23:23,rep_xcorr{icell1,icell2})
% %         hold on
% %         plot(-23:23,avg_rep_xcorr{icell1,icell2},'k','linewidth',2)
% %         axis tight
% %         xlabel('Repeat lag','fontsize',14)
% %         ylabel('Correlation','fontsize',14)
%         f_name = sprintf('C%d_C%d',cell_ids(icell1),cell_ids(icell2));
%         print(f_name,'-dpng'); close all
    end
end

save rec74_xcorr_data ltp_* mtp_* htp_* *_xcorr*

%% FIRING RATE ANALYSIS
cd ~/James_scripts/GLM/t1/rec74_xcorrs/
for icell = 1:ncells
    for n = 1:16
        spktimes = spksegs74{icell,n};
        spk_rep_inds = ceil(spktimes/rep_len/dt);
        rep_spk_rates(icell,n,:) = hist(spk_rep_inds,rep_axis)/5;
    end
end
avg_cond_rates = mean(rep_spk_rates,3);
avg_rates = mean(avg_cond_rates,2);
avg_rep_rates = squeeze(mean(rep_spk_rates,1));
norm_spk_rates = bsxfun(@rdivide,rep_spk_rates,avg_rates);

% for n = 1:length(used_conds)
%     figure
%     set(gcf,'Position',[500 1000 700 1000])
%     subplot(2,1,1)
%     imagesc(squeeze(rep_spk_rates(:,n,:)));colorbar
%     ylabel('Cell number','fontsize',14)
%     xlabel('Repeat number','fontsize',14)
%     title('Absolute rates','fontsize',16)
%     subplot(2,1,2)
%     imagesc(squeeze(norm_spk_rates(:,n,:)));colorbar
%     ylabel('Cell number','fontsize',14)
%     xlabel('Repeat number','fontsize',14)
%     title('Relative rates','fontsize',16)
%     f_name = sprintf('Repeat_ratemat_%s',stimfiles{n});
%     print(f_name,'-dpng'); close all
% end
