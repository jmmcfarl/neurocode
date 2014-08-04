clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Cortical_analysis\UDS_synch_state_dur\UDS_synch_state_dur_data


dsf = 8;
Fsd = 2016/dsf;

niqf = 2016/2;
hcf1 = 100/niqf;
[b,a] = butter(2,hcf1,'low');

params.Fs = Fsd;
params.err = 0;
params.tapers = [2 3];
params.fpass = [0 15];

maxLag = 5*Fsd;
lags = -maxLag:maxLag;

winlength = 20;
winslide = 5;
noverlap = winlength-winslide;
movingwin = [winlength winslide];


for d = 1:length(sess_data)
    
    
    cd(sess_data(d).directory)
    disp(num2str(d))
    
    load used_data wcv_minus_spike lf8 lf5
    
        down_w = filtfilt(b,a,wcv_minus_spike);
    down_8 = filtfilt(b,a,lf8);
    down_5 = filtfilt(b,a,lf5);

        down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
    down_5 = downsample(down_5,dsf);

        down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_5 = zscore(down_5);
    
     [Pw{d},t{d},f{d}]=mtspecgramc(down_w,movingwin,params);
    [P8{d},t{d},f{d}]=mtspecgramc(down_8,movingwin,params);
    [P5{d},t{d},f{d}]=mtspecgramc(down_5,movingwin,params);

    logpow_w(d,:) = nanmean(10*log10(Pw{d}));
    logpow_8(d,:) = nanmean(10*log10(P8{d}));
    logpow_5(d,:) = nanmean(10*log10(P5{d}));
    
    
        Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    pcolor(t{d},f{d},10*log10(P8{d}'));shading flat;
    caxis([-35 2]);
    ylim([0 15])
             cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
tname = ['C:\WC_Germany\Cortical_analysis\specgram\lf8_' cell_name];
    print('-dpng',tname);
    close

          Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    pcolor(t{d},f{d},10*log10(P5{d}'));shading flat;
    caxis([-35 2]);
    ylim([0 15])
             cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
tname = ['C:\WC_Germany\Cortical_analysis\specgram\lf5_' cell_name];
    print('-dpng',tname);
    close
  
           Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    pcolor(t{d},f{d},10*log10(Pw{d}'));shading flat;
    caxis([-35 2]);
    ylim([0 15])
             cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
tname = ['C:\WC_Germany\Cortical_analysis\specgram\mp_' cell_name];
    print('-dpng',tname);
    close
 
end
