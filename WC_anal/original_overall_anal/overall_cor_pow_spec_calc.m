clear all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\desynch_extract\desynch_points
load C:\WC_Germany\overall_calcs\overall_info_file_data

dsf = 8;
Fsd = 2016/dsf;
params.Fs = Fsd;
params.tapers = [4 7];
params.err = [2 .05];
params.fpass = [0 45];
movingwin = [100 100];
niqf = 2016/2;
hcf1 = 100/niqf;
[b1,a1] = butter(2,hcf1,'low');

for d = 1:71
    
    disp(sprintf('session %d',d))
    cd(over_dir{d});
    pwd
    
    load used_data wcv_minus_spike lf8 lf3 lf2 
    
    %bandlimit signals
    down_w = filtfilt(b1,a1,wcv_minus_spike);
    down_8 = filtfilt(b1,a1,lf8);
    down_3 = filtfilt(b1,a1,lf3);
    down_2 = filtfilt(b1,a1,lf2);
    
    %correct for variable gains
    down_w = down_w/mp_gain(d);
    down_8 = down_8/lf8_gain(d);
    down_3 = down_3/lf3_gain(d);
    down_2 = down_2/lf2_gain(d);
    
    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
    down_3 = downsample(down_3,dsf);
    down_2 = downsample(down_2,dsf);
    
     
    if ~isempty(desynch_start{d})
        sMarkers = [[1;desynch_stop{d}'] [desynch_start{d}';length(down_w)]]
    else
        sMarkers = [1 length(down_w)];
    end
    
    [Sw(d,:),f,Swerr(d,:,:)]= mtspectrumc_unequal_length_trials(down_w, movingwin, params, sMarkers); 
    [S8(d,:),f,S8err(d,:,:)]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, sMarkers);
    if ~bad_lf3(d)
        [S3(d,:),f,S3err(d,:,:)]= mtspectrumc_unequal_length_trials(down_3, movingwin, params, sMarkers); 
    end
    if ~bad_lf2(d)
        [S2(d,:),f,S2err(d,:,:)]= mtspectrumc_unequal_length_trials(down_2, movingwin, params, sMarkers); 
    end

    
    plot(f,10*log10(Sw(d,:)),'linewidth',2)
    hold on
    plot(f,10*log10(squeeze(Swerr(d,1,:))),'--')
    plot(f,10*log10(squeeze(Swerr(d,2,:))),'--')
    xlim([0 10])
    tname = ['C:\WC_Germany\overall_calcs\cor_pow_spec\w_spec_wb' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    xlim([0 1])
    tname = ['C:\WC_Germany\overall_calcs\cor_pow_spec\w_spec' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
   
    
        plot(f,10*log10(S8(d,:)),'linewidth',2)
    hold on
    plot(f,10*log10(squeeze(S8err(d,1,:))),'--')
    plot(f,10*log10(squeeze(S8err(d,2,:))),'--')
    xlim([0 10])
    tname = ['C:\WC_Germany\overall_calcs\cor_pow_spec\8_spec_wb' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    xlim([0 1])
    tname = ['C:\WC_Germany\overall_calcs\cor_pow_spec\8_spec' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

    if ~bad_lf3(d)
    plot(f,10*log10(S3(d,:)),'linewidth',2)
    hold on
    plot(f,10*log10(squeeze(S3err(d,1,:))),'--')
    plot(f,10*log10(squeeze(S3err(d,2,:))),'--')
    xlim([0 10])
    tname = ['C:\WC_Germany\overall_calcs\cor_pow_spec\3_spec_wb' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    xlim([0 1])
    tname = ['C:\WC_Germany\overall_calcs\cor_pow_spec\3_spec' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
    end
    
    if ~bad_lf2(d)
            plot(f,10*log10(S2(d,:)),'linewidth',2)
    hold on
    plot(f,10*log10(squeeze(S2err(d,1,:))),'--')
    plot(f,10*log10(squeeze(S2err(d,2,:))),'--')
    xlim([0 10])
    tname = ['C:\WC_Germany\overall_calcs\cor_pow_spec\2_spec_wb' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    xlim([0 1])
    tname = ['C:\WC_Germany\overall_calcs\cor_pow_spec\2_spec' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
    end
    
    if ~bad_lf2(d) & ~bad_lf3(d)
     plot(f,10*log10(Sw(d,:)),'linewidth',2)
    hold on
    plot(f,10*log10(S8(d,:)),'r','linewidth',2)
    plot(f,10*log10(S3(d,:)),'g','linewidth',2)
    plot(f,10*log10(S2(d,:)),'k','linewidth',2)
legend('MP','LF8','LF3','LF2')
    xlim([0 10])
    tname = ['C:\WC_Germany\overall_calcs\cor_pow_spec\compare_wb' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    xlim([0 1])
    tname = ['C:\WC_Germany\overall_calcs\cor_pow_spec\compare_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
    end
    
     clear down_w down_8 wcv* lf8 lf3 down_3*
    
end


save C:\WC_Germany\overall_calcs\cor_pow_spec\cor_pow_spec_data f S*