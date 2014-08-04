clear all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\overall_info_file_data

dsf = 8;
Fsd = 2016/dsf;
params.Fs = Fsd;
params.tapers = [2 3];
params.fpass = [0 10];
params.err = 0;
movingwin = [5 1];

niqf = 2016/2;
hcf1 = 100/niqf;
lcf1 = 0.5/niqf;
[b,a] = butter(2,[lcf1 hcf1]);

for d = 1:62
    disp(sprintf('session %d',d))
    cd(over_dir{d});
    load used_data wcv_minus_spike lf8 lf3 lf2


    %bandlimit signals
    down_w = filtfilt(b,a,wcv_minus_spike);
    down_8 = filtfilt(b,a,lf8);
    down_3 = filtfilt(b,a,lf3);
    down_2 = filtfilt(b,a,lf2);

    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
    down_3 = downsample(down_3,dsf);
    down_2 = downsample(down_2,dsf);

    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_3 = zscore(down_3);
    down_2 = zscore(down_2);

    %     % minimize cortical correlation with hipp lfp
    %     temp_cor = corrcoef(down_3,down_8);
    %     hc_cor = temp_cor(2,1)
    %     down3_sub = down_3-hc_cor*down_8;

    [Sw{d},t{d},f{d}]=mtspecgramc(down_w,movingwin,params);
%     [S8{d},t{d},f{d}]=mtspecgramc(down_8,movingwin,params);
    %         if ~isnan(lf3_gain(d))
%     [S3{d},t{d},f{d}]=mtspecgramc(down_3,movingwin,params);
 
% [Sw{d},f{d},t{d}] = SPECTROGRAM(down_w,round(1.5*Fsd),round(0.5*Fsd),[],Fsd);
% 
% 
% [S8{d},f{d},t{d}] = SPECTROGRAM(down_8,1.5,0.5,[],Fsd);
% [S3{d},f{d},t{d}] = SPECTROGRAM(down_3,1.5,0.5,[],Fsd);

    
%     f_beg_so = find(f{d} > 0.4,1,'first');
%     f_end_so = find(f{d} > 1.5,1,'first');
%     f_theta_beg = find(f{d} > 2.5,1,'first');
%     f_theta_end = find(f{d} > 5,1,'first');
%     df = median(diff(f{d}));
%     so_pow_w{d} = trapz(10*log10(Sw{d}(:,f_beg_so:f_end_so)),2)*df;
%     so_pow_8{d} = trapz(10*log10(S8{d}(:,f_beg_so:f_end_so)),2)*df;
%     so_pow_3{d} = trapz(10*log10(S3{d}(:,f_beg_so:f_end_so)),2)*df;
%     theta_pow_w{d} = trapz(10*log10(Sw{d}(:,f_theta_beg:f_theta_end)),2)*df;
%     theta_pow_8{d} = trapz(10*log10(S8{d}(:,f_theta_beg:f_theta_end)),2)*df;
%     theta_pow_3{d} = trapz(10*log10(S3{d}(:,f_theta_beg:f_theta_end)),2)*df;
  
    
    %                 [C3{d},phi3{d},S12,S1,S2,t{d},f{d},confC3,phistd,cerr]=cohgramc(down_w,down_3,movingwin,params);
    % phi3{d}(C3{d}<confC3) = nan;

    %                 C3{d}(C3{d}<confC3) = nan;
    %         end
    %         if ~isnan(lf2_gain(d))
%     [S2{d},t{d},f{d}]=mtspecgramc(down_2,movingwin,params);
    %         [C2{d},phi2{d},S12,S1,S2,t{d},f{d},confC2,phistd,cerr]=cohgramc(down_w,down_2,movingwin,params);
    % phi2{d}(C2{d}<confC2) = nan;

    %         C2{d}(C2{d}<confC2) = nan;

    %         end

    %     [C8{d},phi8{d},S12,S1,S2,t{d},f{d},confC8,phistd,cerr]=cohgramc(down_w,down_8,movingwin,params);
    % phi8{d}(C8{d}<confC8) = nan;

    %     C8{d}(C8{d}<confC8) = nan;


%     comod_w8{d} = corr(Sw{d},S8{d});

    %     if ~isnan(lf3_gain(d))
%     comod_w3{d} = corr(Sw{d},S3{d});
%     comod_83{d} = corr(S8{d},S3{d});
    %     end

%     pcolor(f{d},f{d},comod_w8{d});shading flat
%     caxis([-0.2 1]); colorbar
%     tname = ['C:\WC_Germany\overall_calcs\specgram\comod_w8_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close all
% 
%     % if ~isnan(lf3_gain(d))
%     pcolor(f{d},f{d},comod_w3{d});shading flat
%     caxis([-0.2 1]); colorbar
%     tname = ['C:\WC_Germany\overall_calcs\specgram\comod_w3_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close all
% 
%     pcolor(f{d},f{d},comod_83{d});shading flat
%     caxis([-0.2 1]); colorbar
%     tname = ['C:\WC_Germany\overall_calcs\specgram\comod_83_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close all

    % end

    
        Fig = figure(1);
        clf
        set(Fig,'PaperUnits','centimeters');
        set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
        set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
        pcolor(t{d},f{d},10*log10(Sw{d}'));shading flat;colorbar
        caxis([-30 0]);colorbar
        ylim([0 10])
        tname = ['C:\WC_Germany\overall_calcs\specgram\mp_wb_new_' num2str(cell_type(d)) '_' over_names{d}];
        print('-dpng',tname);
        close

%         Fig = figure(1);
%         clf
%         set(Fig,'PaperUnits','centimeters');
%         set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%         set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%         subplot(2,1,1)
%         pcolor(t{d},f{d},10*log10(S8{d}'));shading flat;colorbar
%         caxis([-45 0]);colorbar
%         ylim([0 10])
%         hold on
%                 plot(t{d},theta_pow_8{d}./so_pow_8{d},'k')
%         subplot(2,1,2)
%         plot(t{d},so_pow_8{d},'linewidth',2)
%         hold on
%         plot(t{d},theta_pow_8{d},'r','linewidth',2)
%         tname = ['C:\WC_Germany\overall_calcs\specgram_class\lf8_' num2str(cell_type(d)) '_' over_names{d}];
%         print('-dpng',tname);
%         close
%         
%         
%     %
%          if ~isnan(lf3_gain(d))
%         Fig = figure(1);
%         clf
%         set(Fig,'PaperUnits','centimeters');
%         set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%         set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%         subplot(2,1,1)
%         pcolor(t{d},f{d},10*log10(S3{d}'));shading flat;colorbar
%         caxis([-45 0]);colorbar
%         ylim([0 10])
%         hold on
%                 plot(t{d},theta_pow_3{d}./so_pow_3{d},'k')
%         subplot(2,1,2)
%         plot(t{d},so_pow_3{d},'linewidth',2)
%         hold on
%         plot(t{d},theta_pow_3{d},'r','linewidth',2)
%         tname = ['C:\WC_Germany\overall_calcs\specgram_class\lf3_' num2str(cell_type(d)) '_' over_names{d}];
%         print('-dpng',tname);
%         close       
%          end
%          
         
    %
    % %      if ~isnan(lf2_gain(d))
    %      Fig = figure(1);
    %     clf
    %     set(Fig,'PaperUnits','centimeters');
    %     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    %     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    %     pcolor(t{d},f{d},10*log10(S2{d}'));shading flat;colorbar
    %     caxis([-45 0]);colorbar
    %     ylim([0 10])
    %     tname = ['C:\WC_Germany\overall_calcs\specgram\lf2_' num2str(cell_type(d)) '_' over_names{d}];
    %     print('-dpng',tname);
    %     ylim([1 40])
    %        tname = ['C:\WC_Germany\overall_calcs\specgram\lf2_wb' num2str(cell_type(d)) '_' over_names{d}];
    %      print('-dpng',tname);
    % %      end
    %
    %
%         Fig = figure(1);
%         clf
%         set(Fig,'PaperUnits','centimeters');
%         set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%         set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%         subplot(2,1,1)
%         pcolor(t{d},f{d},10*log10(Sw{d}'));shading flat;colorbar
%         caxis([-45 0]);colorbar
%         ylim([0 10])
%         hold on
%         plot(t{d},theta_pow_w{d}./so_pow_w{d},'k')
%         subplot(2,1,2)
%         plot(t{d},so_pow_w{d},'linewidth',2)
%         hold on
%         plot(t{d},theta_pow_w{d},'r','linewidth',2)
%         tname = ['C:\WC_Germany\overall_calcs\specgram_class\mp_' num2str(cell_type(d)) '_' over_names{d}];
%         print('-dpng',tname);
% close

% plot(t{d},so_pow_w{d})
% hold on
% plot(t{d},so_pow_8{d},'r')
% plot(t{d},so_pow_3{d},'k')
% plot(t{d},theta_pow_w{d},'--')
% plot(t{d},theta_pow_8{d},'r--')
% plot(t{d},theta_pow_3{d},'k--')
% legend('MP SO','LF8 SO','LF3 SO','MP Theta','LF8 Theta','LF3 Theta')
%         tname = ['C:\WC_Germany\overall_calcs\specgram_class\compare_' num2str(cell_type(d)) '_' over_names{d}];
%         print('-dpng',tname);
% close

        
        
    %         Fig = figure(1)
    %     clf
    %     set(Fig,'PaperUnits','centimeters');
    %     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    %     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    %     subplot(2,1,1)
    %     pcolor(t{d},f{d},(C8{d}'));shading flat;colorbar
    %     ylim([0 10])
    %     subplot(2,1,2)
    %     pcolor(t{d},f{d},(phi8{d}'));shading flat;colorbar
    %         ylim([0 10])
    %
    %     tname = ['C:\WC_Germany\overall_calcs\cohgram\new_8_' num2str(cell_type(d)) '_' over_names{d}];
    %     print('-dpng',tname);
    %     close
    %
    %     if ~isnan(lf3_gain(d))
    %           Fig = figure(1)
    %     clf
    %     set(Fig,'PaperUnits','centimeters');
    %     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    %     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    %     subplot(2,1,1)
    %     pcolor(t{d},f{d},(C3{d}'));shading flat;colorbar
    %         ylim([0 10])
    %         subplot(2,1,2)
    %         pcolor(t{d},f{d},(phi3{d}'));shading flat;colorbar
    %         ylim([0 10])
    %
    %         tname = ['C:\WC_Germany\overall_calcs\cohgram\new_3_' num2str(cell_type(d)) '_' over_names{d}];
    %     print('-dpng',tname);
    %     close
    %     end
    %
    %     if ~isnan(lf2_gain(d))
    %     Fig = figure(1)
    %     clf
    %     set(Fig,'PaperUnits','centimeters');
    %     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    %     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    %     subplot(2,1,1)
    %     pcolor(t{d},f{d},(C2{d}'));shading flat;colorbar
    %     ylim([0 10])
    %     subplot(2,1,2)
    %     pcolor(t{d},f{d},(phi2{d}'));shading flat;colorbar
    %         ylim([0 10])
    %
    %     tname = ['C:\WC_Germany\overall_calcs\cohgram\new_2_' num2str(cell_type(d)) '_' over_names{d}];
    %     print('-dpng',tname);
    %     close
    %     end

end

% save C:\WC_Germany\overall_calcs\specgram\specgram_data_coh t f C8 C3 C2 phi*
