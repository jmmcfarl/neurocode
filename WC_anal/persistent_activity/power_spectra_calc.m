clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\desynch_extract\desynch_points

dsf = 8;
Fsd = 2016/dsf;
params.Fs = Fsd;
params.tapers = [4 7];
params.err = [2 .05];
params.fpass = [0 45];
window = [100 100];
niqf = 2016/2;
hcf1 = 100/niqf;
[b1,a1] = butter(2,hcf1,'low');

for d = 1:length(dir_array)
    
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
    
    if ~isempty(desynch_start{d})
        sMarkers = [[1;desynch_stop{d}'] [desynch_start{d}';length(down_w)]]
    else
        sMarkers = [1 length(down_w)];
    end
    
    [ Sw(d,:), f, Swerr(d,:,:) ]= mtspectrumc_unequal_length_trials(down_w,window, params, sMarkers);
    [ S8(d,:), f, S8err(d,:,:) ]= mtspectrumc_unequal_length_trials(down_8,window, params, sMarkers);
    
%     [Sw(d,:),f,varS,C,Swerr(d,:,:)] = mtspectrumsegc(down_w,window,params);
%     [S8(d,:),f,varS,C,S8err(d,:,:)] = mtspectrumsegc(down_8,window,params);
    
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    plot_vector(Sw(d,:),f,'l',squeeze(Swerr(d,:,:)),'b')
    hold on
    plot_vector(S8(d,:),f,'l',squeeze(S8err(d,:,:)),'r')
    xlim([0 45])
    tname = ['C:\WC_Germany\Persistent_activity\power_spectra\wband_' f_names{d}];
    print('-dpng',tname);
    close
    
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    plot_vector(Sw(d,:),f,'l',squeeze(Swerr(d,:,:)),'b')
    hold on
    plot_vector(S8(d,:),f,'l',squeeze(S8err(d,:,:)),'r')
    xlim([0 1])
    tname = ['C:\WC_Germany\Persistent_activity\power_spectra\' f_names{d}];
    print('-dpng',tname);
    close

    
    
    clear down_w down_8 wcv* lf8 
    
end


save C:\WC_Germany\Persistent_activity\power_spec_data Sw S8 f Swerr S8err