clear all
close all
%%
load G:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('G:\Code\WC_anal\general\')
addpath('G:\WC_Germany\Overall_EC\')
addpath('G:\Code\Chronux\spectral_analysis\continuous\')
addpath('G:\WC_Germany\hsmm_state_detection\\')
addpath('G:\Code\GCCA_toolbox_jan2011')
drive_letter = 'G';

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

dsf = 40;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,[0.05/niqf 15/niqf]);

%%
% for d = 1:length(sess_data)
     d=1;
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    load used_data lf8 lf3 wcv_minus_spike
    
    wcv_d = filtfilt(b,a,wcv_minus_spike);
    wcv_d = downsample(wcv_d,dsf)/sess_data(d).gains(1);
    lf8_d = filtfilt(b,a,lf8);
    lf8_d = downsample(lf8_d,dsf)/sess_data(d).gains(8);
    lf3_d = filtfilt(b,a,lf3);
    lf3_d = downsample(lf3_d,dsf)/sess_data(d).gains(3);
    t_axis = (1:length(wcv_d))/Fsd;
    
    wcv_d = zscore(wcv_d);
    lf8_d = zscore(lf8_d);
    lf3_d = zscore(lf3_d);
    
    
    %     load ./pa_hsmm_state_seq8
    %     [new_seg_inds] = resample_uds_seg_inds(hmm8.UDS_segs,hmm8.Fs,Fsd,length(lf8_d));
    
%     tot_X = [wcv_d(:) lf8_d(:) lf3_d(:)]';
    tot_X = [lf8_d(:) lf3_d(:)]';
    tot_N = size(tot_X,2);
    nvar = size(tot_X,1);
    MINP = 10;
    MAXP = 45;
    PVAL    =   0.01;       % probability threshold for Granger causality significance
    NLAGS   =   -1;         % if -1, best model order is assessed automatically
    Fs      =   Fsd;        % sampling frequency  (for spectral analysis only)
    freqs   =   [0.02:0.02:10];    % frequency range to analyze (spectral analysis only)
    
%     winlen = round(Fsd*100);
%     winslide = round(Fsd*50);
winlen = length(wcv_d);
winslide = winlen;
    nwins = 1+floor((tot_N-winlen)/winslide);
    t = nan(nwins,1);
    for i = 1:nwins
        xs = (i-1)*winslide+1;
        xe = xs+winlen-1;
        t(i) = (xs+xe)/2/Fsd;
        X = tot_X(:,xs:xe);       
        N       =   size(X,2);       % number of observations
        
%         % check covariance stationarity
%         disp('checking for covariance stationarity ...');
%         uroot = cca_check_cov_stat(X,10);
%         inx = find(uroot);
%         if sum(uroot) == 0,
%             disp('OK, data is covariance stationary by ADF');
%         else
%             disp('WARNING, data is NOT covariance stationary by ADF');
%             disp(['unit roots found in variables: ',num2str(inx')]);
%         end
%         
%         % check covariance stationarity again using KPSS test
%         [kh,kpss] = cca_kpss(X);
%         inx = find(kh==0);
%         if isempty(inx),
%             disp('OK, data is covariance stationary by KPSS');
%         else
%             disp('WARNING, data is NOT covariance stationary by KPSS');
%             disp(['unit roots found in variables: ',num2str(inx')]);
%         end
        
        % find best model order
%         if NLAGS == -1,
%             disp('finding best model order ...');
            [bic,aic] = cca_find_model_order(X,MINP,MAXP);
%             disp(['best model order by Bayesian Information Criterion = ',num2str(bic)]);
%             disp(['best model order by Aikaike Information Criterion = ',num2str(aic)]);
            NLAGS = aic;
%         end
        used_lags(d,i) = NLAGS;
        %-------------------------------------------------------------------------
        % analyze time-domain granger
        
        % find time-domain conditional Granger causalities [THIS IS THE KEY FUNCTION]
%         disp('finding conditional Granger causalities ...');
ret = cca_granger_regress(X,NLAGS,1);   % STATFLAG = 1 i.e. compute stats
% part_ret = cca_partialgc(X,NLAGS,1);   % STATFLAG = 1 i.e. compute stats
        % reta = cca_autonomy_regress(X,NLAGS);
        
%         % check that residuals are white
%         dwthresh = 0.05/nvar;    % critical threshold, Bonferroni corrected
%         waut = zeros(1,nvar);
%         for ii=1:nvar,
%             if ret.waut<dwthresh,
%                 waut(ii)=1;
%             end
%         end
%         inx = find(waut==1);
%         if isempty(inx),
%             disp('All residuals are white by corrected Durbin-Watson test');
%         else
%             disp(['WARNING, autocorrelated residuals in variables: ',num2str(inx)]);
%         end
        
%         % check model consistency, ie. proportion of correlation structure of the
%         % data accounted for by the MVAR model
%         if ret.cons>=80,
%             disp(['Model consistency is OK (>80%), value=',num2str(ret.cons)]);
%         else
%             disp(['Model consistency is <80%, value=',num2str(ret.cons)]);
%         end
        
%         % analyze adjusted r-square to check that model accounts for the data (2nd
%         % check)
%         rss = ret.rss_adj;
%         inx = find(rss<0.3);
%         if isempty(inx)
%             disp(['Adjusted r-square is OK: >0.3 of variance is accounted for by model, val=',num2str(mean(rss))]);
%         else
%             disp(['WARNING, low (<0.3) adjusted r-square values for variables: ',num2str(inx)]);
%             disp(['corresponding values are ',num2str(rss(inx))]);
%             disp('try a different model order');
%         end
        
%         %     find significant Granger causality interactions (Bonferonni correction)
%         [PR,q] = cca_findsignificance(ret,PVAL,1);
%         % [PRa,qa] = cca_findsignificance_autonomy(reta,PVAL,1);
%         disp(['testing significance at P < ',num2str(PVAL), ', corrected P-val = ',num2str(q)]);
        
        % extract the significant causal interactions only
        %     GC = ret.gc;
        %     GC2(d,:,:) = GC.*PR;
        GC{d}(i,:,:) = ret.gc;
%         part_GC{d}(i,:,:) = part_ret.gc;
        %     GA(d,:,:) = reta.gaut;
        
        % % calculate causal connectivity statistics
        % disp('calculating causal connectivity statistics');
        % causd = cca_causaldensity(GC,PR);
        % causf = cca_causalflow(GC,PR);
        %
        % disp(['time-domain causal density = ',num2str(causd.cd)]);
        % disp(['time-domain causal density (weighted) = ',num2str(causd.cdw)]);
        %
        % % create Pajek readable file
        % cca_pajek(PR,GC,sfile);
        
        % %-------------------------------------------------------------------------
        % % plot time-domain granger results
        % figure(1); clf reset;
        % FSIZE = 8;
        % colormap(flipud(bone));
        %
        % % plot raw time series
        % for i=2:nvar,
        %     X(i,:) = X(i,:)+(10*(i-1));
        % end
        % subplot(231);
        % set(gca,'FontSize',FSIZE);
        % plot(X');
        % axis('square');
        % set(gca,'Box','off');
        % xlabel('time');
        % set(gca,'YTick',[]);
        % xlim([0 N]);
        % title('Causal Connectivity Toolbox v2.0');
        %
        % % plot granger causalities as matrix
        % subplot(232);
        % set(gca,'FontSize',FSIZE);
        % imagesc(GC2);
        % axis('square');
        % set(gca,'Box','off');
        % title(['Granger causality, p<',num2str(PVAL)]);
        % xlabel('from');
        % ylabel('to');
        % set(gca,'XTick',[1:N]);
        % set(gca,'XTickLabel',1:N);
        % set(gca,'YTick',[1:N]);
        % set(gca,'YTickLabel',1:N);
        %
        % % plot granger causalities as a network
        % subplot(233);
        % cca_plotcausality(GC2,[],5);
        %
        % % plot causal flow  (bar = unweighted, line = weighted)
        % subplot(234);
        % set(gca,'FontSize',FSIZE);
        % set(gca,'Box','off');
        % mval1 = max(abs(causf.flow));
        % mval2 = max(abs(causf.wflow));
        % mval = max([mval1 mval2]);
        % bar(1:nvar,causf.flow,'m');
        % ylim([-(mval+1) mval+1]);
        % xlim([0.5 nvar+0.5]);
        % set(gca,'XTick',[1:nvar]);
        % set(gca,'XTickLabel',1:nvar);
        % title('causal flow');
        % ylabel('out-in');
        % hold on;
        % plot(1:nvar,causf.wflow);
        % axis('square');
        %
        % % plot unit causal densities  (bar = unweighted, line = weighted)
        % subplot(235);
        % set(gca,'FontSize',FSIZE);
        % set(gca,'Box','off');
        % mval1 = max(abs(causd.ucd));
        % mval2 = max(abs(causd.ucdw));
        % mval = max([mval1 mval2]);
        % bar(1:nvar,causd.ucd,'m');
        % ylim([-0.25 mval+1]);
        % xlim([0.5 nvar+0.5]);
        % set(gca,'XTick',[1:nvar]);
        % set(gca,'XTickLabel',1:nvar);
        % title('unit causal density');
        % hold on;
        % plot(1:nvar,causd.ucdw);
        % axis('square');
        %
        % %-------------------------------------------------------------------------
        % % analyze frequency-domain granger
        %
        SPECTHRESH = 0.2.*ones(1,length(freqs));    % bootstrap not used in this demo
        %
        % find pairwise frequency-domain Granger causalities [KEY FUNCTION]
%         disp('finding pairwise frequency-domain Granger causalities ...');
        [GW,COH,pp]=cca_pwcausal(X,1,N,NLAGS,Fs,freqs,0);
        
        % % calculate freq domain causal connectivity statistics
        % disp('calculating causal connectivity statistics');
        % causd = cca_causaldensity_spectral(GW,SPECTHRESH);
        % causf = cca_causalflow_spectral(GW,SPECTHRESH);
        %
        % totalcd = sum(causd.scdw);
        % disp(['freq-domain causal density (weighted) = ',num2str(totalcd)]);
        %
        % %-------------------------------------------------------------------------
        % % plot frequency-domain granger results
        % figure(2); clf reset;
        % FSIZE = 8;
        % colormap(flipud(bone));
        %
        % % plot fft for each variable
        % ct = 1;
        % for i=1:nvar,
        %     subplot(3,nvar,ct);
        %     cca_spec(X(i,:),Fs,1);
        %     title(['v',num2str(i)]);
        %     if i==1,
        %         ylabel('fft: amplitude');
        %     end
        %     ct = ct+1;
        %     set(gca,'Box','off');
        % end
        %
        % % plot causal density spectrum for each variable
        % for i=1:nvar,
        %     subplot(3,nvar,ct);
        %     plot(causd.sucdw(i,:));
        %     if i==1,
        %         ylabel('unit cd');
        %     end
        %     ct = ct+1;
        %     set(gca,'Box','off');
        % end
        %
        % % plot causal flow spectrum for each variable
        % for i=1:nvar,
        %     subplot(3,nvar,ct);
        %     plot(causf.swflow(i,:));
        %     if i==1,
        %         ylabel('unit flow');
        %     end
        %     ct = ct+1;
        %     set(gca,'Box','off');
        % end
        %
        % % plot network causal density
        % figure(3); clf reset;
        % plot(causd.scdw);
        % set(gca,'Box','off');
        % title(['spectral cd, total=',num2str(totalcd),', thresh=',num2str(SPECTHRESH)]);
        % xlabel('Hz');
        % ylabel('weighted cd');
        %
%         % plot causal interactions
%         figure(4); clf reset;
%         cca_plotcausality_spectral(GW,freqs);
%         
        % % plot coherence
        % figure(5); clf reset;
        % cca_plotcoherence(COH,freqs);
        %     NLAGS = bic;
        %    ret = cca_granger_regress(X,NLAGS,1);
        %
        %     % end
        
        %%%%
    end
    
    mean_GC(d,:,:) = squeeze(nanmean(GC{d},1));
    median_GC(d,:,:) = squeeze(nanmedian(GC{d},1));
%     mean_pGC(d,:,:) = squeeze(nanmean(part_GC{d}));
%     median_pGC(d,:,:) = squeeze(nanmedian(part_GC{d}));

%     figure
%     plot(t,squeeze(GC{d}(:,3,1)))
%     hold on
%     plot(t,squeeze(GC{d}(:,1,3)),'r')
%     plot(t,squeeze(GC{d}(:,3,2)),'k')
% %     plot(t,squeeze(part_GC{d}(:,3,1)),'--')
% %     plot(t,squeeze(part_GC{d}(:,1,3)),'r--')
% %     plot(t,squeeze(part_GC{d}(:,3,2)),'k--')
%         t_names = ['G:\WC_Germany\persistent_9_27_2010\caus_check\' s_name];
% print('-dpng',t_names), close
% end
% 
%%
