clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update

dsf = 8;
Fsd = 2016/dsf;
params.Fs = Fsd;
params.tapers = [4 7];
params.err = [0];
params.fpass = [10 100];
movingwin = [0.5 0.25];
niqf = 2016/2;
hcf1 = 150/niqf;
[b1,a1] = butter(2,hcf1,'low');

for d = 1:length(dir_array)
    d=14
    disp(sprintf('session %d',d))
    cd(dir_array{d});
    pwd
    
    load used_data wcv_minus_spike lf8
    
    %bandlimit signals
    down_w = filtfilt(b1,a1,wcv_minus_spike);
    down_8 = filtfilt(b1,a1,lf8);

    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
    
    
    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    
    [Sw,t,f]=mtspecgramc(down_w,movingwin,params);
    fbeg = find(f > 20,1,'first');
    fend = find(f > 80,1,'first');
    gamm_pow = trapz(10*log10(Sw(:,fbeg:fend)),2);
    
    t_axis = (1:length(down_w))/Fsd;
    tseg_axis = (1:length(gamm_pow))*movingwin(2);
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     plot_vector(Sw(d,:),f,'l',squeeze(Swerr(d,:,:)),'b')
%     hold on
%     plot_vector(S8(d,:),f,'l',squeeze(S8err(d,:,:)),'r')
%     xlim([0 45])
%     tname = ['C:\WC_Germany\Persistent_activity\power_spectra\wband_' f_names{d}];
%     print('-dpng',tname);
%     close
%     
    
    
    clear down_w down_8 wcv* lf8 
    
end


