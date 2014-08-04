clear all

load C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\UDS_dur_data_over_smooth
load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\correl_data
load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\lag_data

dsf = 4;
Fsd = 2016/dsf;
winsize = 30;
niqf = 2016/2;
hcf = 2/niqf;
[b,a] = butter(2,hcf,'low');
cnt = 0;
maxlag = 10*Fsd;

for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike

    lf8 = filtfilt(b,a,lf8);
    wcv = filtfilt(b,a,wcv_minus_spike);

    lf8_d = downsample(lf8,dsf);
    wcv_d = downsample(wcv,dsf);

    lf8_d = zscore(lf8_d);
    wcv_d = zscore(wcv_d);

    pers_check = (max_up_dur{d} > 5);
    pers_check8 = (max_up_dur8{d} > 5);

    if max(pers_check) == 1
        %% Get C, Sw, and S8 during pers and non-pers times
        pers_times = find(max_up_dur{d} > 5 & max_up_dur8{d} < 5);
        npers_times = find(max_up_dur{d} < 5 | max_up_dur8{d} > 5);
        time_d = size(l_w_corr{d},1);
        pers_times(pers_times>time_d) = [];
        npers_times(npers_times>time_d) = [];
%         pers_C{d} = mean(C{d}(pers_times,:));
%         npers_C{d} = mean(C{d}(npers_times,:));
%         pers_Sw{d} = mean(Sw{d}(pers_times,:));
%         npers_Sw{d} = mean(Sw{d}(npers_times,:));
%         pers_S8{d} = mean(S8{d}(pers_times,:));
%         npers_S8{d} = mean(S8{d}(npers_times,:));
        pers_xcorr{d} = mean(l_w_corr{d}(pers_times,:));
        npers_xcorr{d} = mean(l_w_corr{d}(npers_times,:));
        pers_wcorr{d} = mean(w_acorr{d}(pers_times,:));
        npers_wcorr{d} = mean(w_acorr{d}(npers_times,:));
        pers_lcorr{d} = mean(l_acorr{d}(pers_times,:));
        npers_lcorr{d} = mean(l_acorr{d}(npers_times,:));

%         %% PLot pers and non-pers spectral analyses
%         Fig = figure(1)
%         clf
%         set(Fig,'PaperUnits','centimeters');
%         set(gcf, 'PaperSize', [40 60]);% paper size is in [width height] format
%         set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%         subplot(3,1,1)
%         plot(f{d},pers_C{d},'linewidth',2)
%         hold on
%         plot(f{d},npers_C{d},'r','linewidth',2)
%         xlim([0 2]);ylim([0 1])
%         subplot(3,1,2)
%         plot(f{d},10*log10(pers_Sw{d}),'linewidth',2)
%         hold on
%         plot(f{d},10*log10(npers_Sw{d}),'r','linewidth',2)
%         xlim([0 2])
% 
%         subplot(3,1,3)
%         plot(f{d},10*log10(pers_S8{d}),'linewidth',2)
%         hold on
%         plot(f{d},10*log10(npers_S8{d}),'r','linewidth',2)
%         xlim([0 2])
% 
%         tname = ['C:\WC_Germany\JMM_Analysis_pyr\pers_state_anal\state_dep_spectra_' f_names{d}];
%         print('-dpng',tname);
%         close all

        Fig = figure(1)
        clf
        set(Fig,'PaperUnits','centimeters');
        set(gcf, 'PaperSize', [40 60]);% paper size is in [width height] format
        set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
        subplot(3,1,1)
        plot(lags,pers_xcorr{d},'linewidth',2)
        hold on
        plot(lags,npers_xcorr{d},'r','linewidth',2)
        subplot(3,1,2)
        plot(lags,pers_wcorr{d},'linewidth',2)
        hold on
        plot(lags,npers_wcorr{d},'r','linewidth',2)
        subplot(3,1,3)
        plot(lags,pers_lcorr{d},'linewidth',2)
        hold on
        plot(lags,npers_lcorr{d},'r','linewidth',2)
        tname = ['C:\WC_Germany\JMM_Analysis_pyr\pers_state_anal\state_dep_tdanal_' f_names{d}];
        print('-dpng',tname);
        close all


    end


    %get rid of instances where both are persistent
    pers_check = pers_check &  ~pers_check8;

    %find start and stop times of wcv pers (in seconds)
    dpc = diff(pers_check);
    pstart = find(dpc == 1);
    pstop = find(dpc == -1)+floor(winsize/2)-2;

    if pers_check(1) == 1
        pstart = [1 pstart];
    end
    if pers_check(end) == 1
        pstop = [pstop length(pers_check)];
    end

    %make sure all pers and non-pers states are at least 30 long
    pers_dur = pstop-pstart;
    npers_dur = pstart(2:end) - pstop(1:end-1);

    bad_npers = find(npers_dur < 30);
    pstop(bad_npers) = [];
    pstart(bad_npers+1) = [];

%     pstart
%     pstop

    %convert to sample index
    pstart = pstart*Fsd;
    pstop = pstop*Fsd;

    pers_wcv = [];
    npers_wcv = [];
    pers_lf8 = [];
    npers_lf8 = [];
%     %piece together data during pers states
%     if max(pers_check) == 1
%         for t = 1:length(pstart)
% 
%             pers_wcv = [pers_wcv;wcv_d(pstart(t):pstop(t))];
%             pers_lf8 = [pers_lf8;lf8_d(pstart(t):pstop(t))];
% 
%             [p_xcor,lags] = xcov(wcv_d(pstart(t):pstop(t)),lf8_d(pstart(t):pstop(t)),maxlag,'coeff');
% 
%             pers_xcor(t,:) = p_xcor;
% 
%         end
%     end

%     %piece together data during npers states
%     if min(pers_check) == 0
%         if max(pers_check) == 1
% 
%             for t = 1:length(pstart)-1
% 
%                 npers_wcv = [npers_wcv;wcv_d(pstop(t):pstart(t+1))];
%                 npers_lf8 = [npers_lf8;lf8_d(pstop(t):pstart(t+1))];
% 
%                 wcv_temp = wcv_d(pstop(t):pstart(t+1));
%                 lf8_temp = lf8_d(pstop(t):pstart(t+1));
%                 [n_xcor,lags] = xcov(wcv_temp,lf8_temp,maxlag,'coeff');
% 
%                 npers_xcor(t,:) = n_xcor(1:10081);
% 
%             end
% 
%             npers_wcv = [wcv_d(1:pstart(1));npers_wcv];
%             npers_lf8 = [lf8_d(1:pstart(1));npers_lf8];
%             npers_wcv = [npers_wcv;wcv_d(pstop(end):end)];
%             npers_lf8 = [npers_lf8;lf8_d(pstop(end):end)];
% 
%         else
%             npers_wcv = wcv_d;
%             npers_lf8 = lf8_d;
%         end
% 
%         if ~isempty(pers_wcv) & ~isempty(npers_wcv)
%             cnt = cnt+1;
%             perc_pers(cnt) = round(100*length(pers_wcv)/(length(pers_wcv)+length(npers_wcv)));
% %             Fig = figure(1)
% %             clf
% %             set(Fig,'PaperUnits','centimeters');
% %             set(gcf, 'PaperSize', [20 30]);% paper size is in [width height] format
% %             set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
% %             subplot(2,1,1)
%             [w_kde_p(cnt,:),w_kgrid_p] = gpkde(pers_wcv,-3,[-5 5 500]);
%             [w_kde_n(cnt,:),w_kgrid_n] = gpkde(npers_wcv,-3,[-5 5 500]);
% %             plot(w_kgrid_p,w_kde_p(cnt,:),'linewidth',2)
% %             hold on
% %             plot(w_kgrid_n,w_kde_n(cnt,:),'r','linewidth',2)
% %             legend('Pers','No Pers')
% %             title(num2str(perc_pers(cnt)))
% %             subplot(2,1,2)
%             [l_kde_p(cnt,:),l_kgrid_p] = gpkde(pers_lf8,-3,[-5 5 500]);
%             [l_kde_n(cnt,:),l_kgrid_n] = gpkde(npers_lf8,-3, [-5 5 500]);
% %             plot(l_kgrid_p,l_kde_p(cnt,:),'linewidth',2)
% %             hold on
% %             plot(l_kgrid_n,l_kde_n(cnt,:),'r','linewidth',2)
% %             legend('Pers','No Pers')
% %             tname = ['C:\WC_Germany\JMM_Analysis_pyr\pers_state_anal\state_dep_hist_' f_names{d}];
% %             print('-dpng',tname);
% %             close all
% 
%             %plot xcorr in both states
%             m_pers_xcorr(cnt,:) = mean(pers_xcor);
%             m_npers_xcorr(cnt,:) = mean(npers_xcor);
%             %            Fig4 = figure(4)
%             %            clf
%             %            set(Fig4,'PaperUnits','centimeters');
%             %            set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%             %            set(Fig4,'PaperPosition',[0,0,(get(Fig4,'PaperSize'))])
%             %             plot(lags/Fsd,m_pers_xcorr)
%             %             hold on
%             %             plot(lags/Fsd,m_npers_xcorr,'r')
%             %             legend('Pers','Npers')
%             %             title(num2str(perc_pers(cnt)))
%             % tname = ['C:\WC_Germany\JMM_Analysis_pyr\pers_state_anal\state_dep_xcorr_' f_names{d}];
%             %     print('-dpng',tname);
%             %             close all
% 
% 
%         else
%             if isempty(pers_wcv)
%                 disp('no pers')
%             else
%                 disp('only pers')
%             end
%         end

%     end
end


% %            Fig = figure(1)
% %            clf
% %            set(Fig,'PaperUnits','centimeters');
% %            set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
% %            set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
% % %            perc_pers = [perc_pers zeros(1,17-cnt)];
% %            hist(perc_pers,50)
% %
% %
% %
%            m_w_p = mean(w_kde_p);
%            u_w_p = m_w_p+2*std(w_kde_p)/sqrt(14);
%            l_w_p = m_w_p-2*std(w_kde_p)/sqrt(14);
%
%            m_w_n = mean(w_kde_n);
%            u_w_n = m_w_n+2*std(w_kde_n)/sqrt(14);
%            l_w_n = m_w_n-2*std(w_kde_n)/sqrt(14);
%
%            m_l_p = mean(l_kde_p);
%            u_l_p = m_l_p+2*std(l_kde_p)/sqrt(14);
%            l_l_p = m_l_p-2*std(l_kde_p)/sqrt(14);
%
%            m_l_n = mean(l_kde_n);
%            u_l_n = m_l_n+2*std(l_kde_n)/sqrt(14);
%            l_l_n = m_l_n-2*std(l_kde_n)/sqrt(14);
%
%                 Fig2 = figure(2)
%            clf
%            set(Fig2,'PaperUnits','centimeters');
%            set(Fig2, 'PaperSize', [30 30]);% paper size is in [width height] format
%            set(Fig2,'PaperPosition',[0,0,(get(Fig2,'PaperSize'))])
%             subplot(2,1,1)
%             plot(w_kgrid_p,m_w_p,'linewidth',2)
%             hold on
%             plot(w_kgrid_p,u_w_p,'--')
%             plot(w_kgrid_p,l_w_p,'--')
%             plot(w_kgrid_n,m_w_n,'r','linewidth',2)
%             plot(w_kgrid_n,u_w_n,'r--')
%             plot(w_kgrid_n,l_w_n,'r--')
%             subplot(2,1,2)
%                         plot(l_kgrid_p,m_l_p,'linewidth',2)
%             hold on
%             plot(l_kgrid_p,u_l_p,'--')
%             plot(l_kgrid_p,l_l_p,'--')
%             plot(l_kgrid_n,m_l_n,'r','linewidth',2)
%             plot(l_kgrid_n,u_l_n,'r--')
%             plot(l_kgrid_n,l_l_n,'r--')
%
%
% %%Repeat using only cells with less than 50% pers activity
%
% %get rid of 7_4_B
% w_kde_p(10,:) = [];
% w_kde_n(10,:) = [];
% l_kde_p(10,:) = [];
% l_kde_n(10,:) = [];
% perc_pers(10) = [];
%
% bad_cells = find(perc_pers > 50);
%
% w_kde_p(bad_cells,:) = [];
% w_kde_n(bad_cells,:) = [];
% l_kde_p(bad_cells,:) = [];
% l_kde_n(bad_cells,:) = [];
%
%            m_w_p = mean(w_kde_p);
%            u_w_p = m_w_p+2*std(w_kde_p)/sqrt(10);
%            l_w_p = m_w_p-2*std(w_kde_p)/sqrt(10);
%
%            m_w_n = mean(w_kde_n);
%            u_w_n = m_w_n+2*std(w_kde_n)/sqrt(10);
%            l_w_n = m_w_n-2*std(w_kde_n)/sqrt(10);
%
%            m_l_p = mean(l_kde_p);
%            u_l_p = m_l_p+2*std(l_kde_p)/sqrt(10);
%            l_l_p = m_l_p-2*std(l_kde_p)/sqrt(10);
%
%            m_l_n = mean(l_kde_n);
%            u_l_n = m_l_n+2*std(l_kde_n)/sqrt(10);
%            l_l_n = m_l_n-2*std(l_kde_n)/sqrt(10);
%
%                 Fig3 = figure(3)
%            clf
%            set(Fig3,'PaperUnits','centimeters');
%            set(Fig3, 'PaperSize', [30 30]);% paper size is in [width height] format
%            set(Fig3,'PaperPosition',[0,0,(get(Fig3,'PaperSize'))])
%             subplot(2,1,1)
%             plot(w_kgrid_p,m_w_p,'linewidth',2)
%             hold on
%             plot(w_kgrid_p,u_w_p,'--')
%             plot(w_kgrid_p,l_w_p,'--')
%             plot(w_kgrid_n,m_w_n,'r','linewidth',2)
%             plot(w_kgrid_n,u_w_n,'r--')
%             plot(w_kgrid_n,l_w_n,'r--')
%             subplot(2,1,2)
%                         plot(l_kgrid_p,m_l_p,'linewidth',2)
%             hold on
%             plot(l_kgrid_p,u_l_p,'--')
%             plot(l_kgrid_p,l_l_p,'--')
%             plot(l_kgrid_n,m_l_n,'r','linewidth',2)
%             plot(l_kgrid_n,u_l_n,'r--')
%             plot(l_kgrid_n,l_l_n,'r--')
%
%




% m_p_xco = mean(m_pers_xcorr);
% m_n_xco = mean(m_npers_xcorr);
% u_p_xc = m_p_xco+2*std(m_pers_xcorr)/sqrt(14);
% l_p_xc = m_p_xco-2*std(m_pers_xcorr)/sqrt(14);
% u_n_xc = m_n_xco+2*std(m_npers_xcorr)/sqrt(14);
% l_n_xc = m_n_xco-2*std(m_npers_xcorr)/sqrt(14);
% 
% plot(lags/Fsd,m_p_xco,'linewidth',2)
% hold on
% plot(lags/Fsd,u_p_xc,'--')
% plot(lags/Fsd,l_p_xc,'--')
% plot(lags/Fsd,m_n_xco,'r','linewidth',2)
% plot(lags/Fsd,u_n_xc,'r--')
% plot(lags/Fsd,l_n_xc,'r--')

