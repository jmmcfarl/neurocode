clear all
close all

load('C:\WC_Germany\tetrode_recs\tetrode_dir.mat')

dsf = 8;
Fsd = 2016/dsf;

niqf = 2016/2;
hcf1 = 100/niqf;
[b,a] = butter(2,hcf1,'low');
[b2,a2] = butter(2,[0.1/niqf 10/niqf]);

params.Fs = Fsd;
params.err = 0;
params.tapers = [2 3];
params.fpass = [0 10];

maxLag = 5*Fsd;
lags = -maxLag:maxLag;

winlength = 20;
winslide = 5;
noverlap = winlength-winslide;
movingwin = [winlength winslide];

for d = 1:length(dir_array)

    disp(sprintf('session %d',d))
    cd(dir_array{d});
    load lfp_data lf1_data lf12_data
    
    wcv_minus_spike = lf1_data;
    lf8 = lf12_data;


    %bandlimit signals
    down_w = filtfilt(b,a,wcv_minus_spike);
    down_8 = filtfilt(b,a,lf8);
    down_wf = filtfilt(b2,a2,wcv_minus_spike);
    down_8f = filtfilt(b2,a2,lf8);

    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
    down_wf = downsample(down_wf,dsf);
    down_8f = downsample(down_8f,dsf);
    
    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_wf = zscore(down_wf);
    down_8f = zscore(down_8f);
    

%     [Sw,f{d},t{d},Pw{d}] = spectrogram(down_w,round(winlength*Fsd),round(noverlap*Fsd),[],Fsd);       
%     [S8,f{d},t{d},P8{d}] = spectrogram(down_8,round(winlength*Fsd),round(noverlap*Fsd),[],Fsd);  

     [Pw{d},t{d},f{d}]=mtspecgramc(down_w,movingwin,params);
    [P8{d},t{d},f{d}]=mtspecgramc(down_8,movingwin,params);
   
    f_beg_so = find(f{d} > 0.2,1,'first');
    f_end_so = find(f{d} > 1.5,1,'first');
    f_theta_beg = find(f{d} > 3.5,1,'first');
    f_theta_end = find(f{d} > 6,1,'first');

%     max_theta_8 = max(10*log10(Pw{d}(f_theta_beg:f_theta_end,:)));
%     max_so_8 = max(10*log10(Pw{d}(f_beg_so:f_end_so,:)));
    max_theta_8 = max(10*log10(P8{d}(:,f_theta_beg:f_theta_end)),[],2);
    max_so_8 = max(10*log10(P8{d}(:,f_beg_so:f_end_so)),[],2);
    
    
    totalDur = length(down_w)/Fsd;
    numWins = round((totalDur-winlength)/(winlength-noverlap));
    Aw{d} = zeros(numWins,length(lags));
    A8{d} = zeros(numWins,length(lags));
    for i = 1:numWins       
       segBeg = 1+(i-1)*(winlength-noverlap)*Fsd;
       segEnd = segBeg+winlength*Fsd;
       Aw{d}(i,:) = xcov(down_wf(segBeg:segEnd),maxLag,'coeff');
       A8{d}(i,:) = xcov(down_8f(segBeg:segEnd),maxLag,'coeff');
    end
    
    At{d} = ((1:numWins)-1)*(winlength-noverlap) + winlength/2;
    
    
    Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    pcolor(t{d},f{d},10*log10(P8{d}'));shading flat;
    caxis([-25 2]);
    ylim([0 7])
    subplot(2,1,2)
    pcolor(At{d},lags/Fsd,A8{d}');shading flat; colorbar
    ylim([0 maxLag/Fsd])
    caxis([-0.4 0.4])
% plot(t{d},max_theta_8,'r')
% hold on
% plot(t{d},max_so_8,'k')
% plot(t{d},max_theta_8-max_so_8)
% xlim([t{d}(1) t{d}(end)])
    tname = ['C:\WC_Germany\tetrode_recs\specgram_auto\lf8_' f_names{d}];
    print('-dpng',tname);
    close

     Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    pcolor(t{d},f{d},10*log10(Pw{d}'));shading flat;colorbar
    caxis([-25 5]);colorbar
    ylim([0 7])
    subplot(2,1,2)
    pcolor(At{d},lags/Fsd,Aw{d}');shading flat; colorbar
    ylim([0 maxLag/Fsd])
    caxis([-0.4 0.4])
    tname = ['C:\WC_Germany\tetrode_recs\specgram_auto\mp_' f_names{d}];
    print('-dpng',tname);
    close
       
    %
         

end

% save C:\WC_Germany\persistent_revised\specgram_auto\specgram_data_coh_2 t f Pw P8 At Aw A8
