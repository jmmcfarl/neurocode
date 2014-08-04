clear all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\overall_info_file_data

dsf = 8;
Fsd = 2016/dsf;
params.Fs = Fsd;
params.tapers = [4 7];
params.fpass = [0 10];
params.err = [2 .05];
movingwin = [20 5];

niqf = 2016/2;
hcf1 = 100/niqf;
[b,a] = butter(2,hcf1,'low');

maxLag = 3*Fsd;
lags = -maxLag:maxLag;

lcf = 0.2/Fsd*2;
hcf = 10/Fsd*2;
[b2,a2] = butter(2,[lcf hcf]);

for d = 1:62

    disp(sprintf('session %d',d))
    cd(over_dir{d});
    load used_data wcv_minus_spike lf8 lf3 lf2


    %bandlimit signals
    down_w = filtfilt(b,a,wcv_minus_spike/mp_gain(d));
    down_8 = filtfilt(b,a,lf8/lf8_gain(d));
    down_3 = filtfilt(b,a,lf3/lf3_gain(d));
    down_2 = filtfilt(b,a,lf2/lf2_gain(d));

    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
    down_3 = downsample(down_3,dsf);
    down_2 = downsample(down_2,dsf);

    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_3 = zscore(down_3);
    down_2 = zscore(down_2);

    good3 = ~isnan(lf3_gain(d));
    
    %     % minimize cortical correlation with hipp lfp
    %     temp_cor = corrcoef(down_3,down_8);
    %     hc_cor = temp_cor(2,1)
    %     down3_sub = down_3-hc_cor*down_8;

    [Sw{d},t{d},f{d}]=mtspecgramc(down_w,movingwin,params);
    [S8{d},t{d},f{d}]=mtspecgramc(down_8,movingwin,params);
    %         if ~isnan(lf3_gain(d))
    [S3{d},t{d},f{d}]=mtspecgramc(down_3,movingwin,params);
    
    
    f_beg_so = find(f{d} > 0.4,1,'first');
    f_end_so = find(f{d} > 1.5,1,'first');
    f_theta_beg = find(f{d} > 2.5,1,'first');
    f_theta_end = find(f{d} > 5,1,'first');
    df = median(diff(f{d}));
    so_pow_w{d} = trapz(10*log10(Sw{d}(:,f_beg_so:f_end_so)),2)*df;
    so_pow_8{d} = trapz(10*log10(S8{d}(:,f_beg_so:f_end_so)),2)*df;
    so_pow_3{d} = trapz(10*log10(S3{d}(:,f_beg_so:f_end_so)),2)*df;
    theta_pow_w{d} = trapz(10*log10(Sw{d}(:,f_theta_beg:f_theta_end)),2)*df;
    theta_pow_8{d} = trapz(10*log10(S8{d}(:,f_theta_beg:f_theta_end)),2)*df;
    theta_pow_3{d} = trapz(10*log10(S3{d}(:,f_theta_beg:f_theta_end)),2)*df;
  
    filt_w = filtfilt(b2,a2,down_w);
    filt_8 = filtfilt(b2,a2,down_8);
    filt_3 = filtfilt(b2,a2,down_3);
    
    totalDur = length(down_w)/Fsd;
    numWins = round((totalDur-movingwin(1))/movingwin(2));
    Aw{d} = zeros(numWins,length(lags));
    A8{d} = zeros(numWins,length(lags));
    if good3
    A3{d} = zeros(numWins,length(lags));
    end
    for i = 1:numWins
        
       segBeg = 1+(i-1)*movingwin(2)*Fsd;
       segEnd = segBeg+movingwin(1)*Fsd;
       Aw{d}(i,:) = xcov(filt_w(segBeg:segEnd),maxLag,'coeff');
       A8{d}(i,:) = xcov(filt_8(segBeg:segEnd),maxLag,'coeff');
if good3
          A3{d}(i,:) = xcov(filt_3(segBeg:segEnd),maxLag,'coeff');

end
    end
    
    
    
    At{d} = ((1:numWins)-1)*movingwin(2) + movingwin(1)/2;
    
    
        Fig = figure(1);
        clf
        set(Fig,'PaperUnits','centimeters');
        set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
        set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
        subplot(2,1,1)
        pcolor(t{d},f{d},10*log10(S8{d}'));shading flat;colorbar
        caxis([-45 0]);colorbar
        ylim([0 10])
        hold on
        plot(t{d},theta_pow_8{d}./so_pow_8{d},'k')
        subplot(2,1,2)
        pcolor(At{d},lags/Fsd,A8{d}');shading flat; colorbar
        ylim([0 maxLag/Fsd])
         caxis([-0.4 0.4])
        tname = ['C:\WC_Germany\overall_calcs\specgram_auto\lf8_' num2str(cell_type(d)) '_' over_names{d}];
        print('-dpng',tname);
        close
        
        
    %
         if ~isnan(lf3_gain(d))
        Fig = figure(1);
        clf
        set(Fig,'PaperUnits','centimeters');
        set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
        set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
        subplot(2,1,1)
        pcolor(t{d},f{d},10*log10(S3{d}'));shading flat;colorbar
        caxis([-45 0]);colorbar
        ylim([0 10])
        hold on
                plot(t{d},theta_pow_3{d}./so_pow_3{d},'k')
        subplot(2,1,2)
        pcolor(At{d},lags/Fsd,A3{d}');shading flat; colorbar
        ylim([0 maxLag/Fsd])
         caxis([-0.4 0.4])
        tname = ['C:\WC_Germany\overall_calcs\specgram_auto\lf3_' num2str(cell_type(d)) '_' over_names{d}];
        print('-dpng',tname);
        close       
         end
         

        Fig = figure(1);
        clf
        set(Fig,'PaperUnits','centimeters');
        set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
        set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
        subplot(2,1,1)
        pcolor(t{d},f{d},10*log10(Sw{d}'));shading flat;colorbar
        caxis([-45 0]);colorbar
        ylim([0 10])
        hold on
        plot(t{d},theta_pow_w{d}./so_pow_w{d},'k')
        subplot(2,1,2)
        pcolor(At{d},lags/Fsd,Aw{d}');shading flat; colorbar
        ylim([0 maxLag/Fsd])
         caxis([-0.4 0.4])
        tname = ['C:\WC_Germany\overall_calcs\specgram_auto\mp_' num2str(cell_type(d)) '_' over_names{d}];
        print('-dpng',tname);
close


end

% save C:\WC_Germany\overall_calcs\specgram\specgram_data_coh t f C8 C3 C2 phi*
