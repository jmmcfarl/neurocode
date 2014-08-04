clear all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\desynch_extract\desynch_points

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

for d = 1:length(over_dir)
    
    disp(sprintf('session %d',d))
    cd(over_dir{d});
    pwd
    
    load used_data wcv_minus_spike lf8 lf3 lf2 
    
    %bandlimit signals
    down_w = filtfilt(b1,a1,wcv_minus_spike);
    down_8 = filtfilt(b1,a1,lf8);
    down_3 = filtfilt(b1,a1,lf3);
    down_2 = filtfilt(b1,a1,lf2);
    
    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
    down_3 = downsample(down_3,dsf);
    down_2 = downsample(down_2,dsf);
    
    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_3 = zscore(down_3);
    down_2 = zscore(down_2);
    
       temp_cor = corrcoef(down_3,down_8);
    hc_cor = temp_cor(2,1)
    down3_sub = down_3-hc_cor*down_8;
 
    if ~isempty(desynch_start{d})
        sMarkers = [[1;desynch_stop{d}'] [desynch_start{d}';length(down_w)]]
    else
        sMarkers = [1 length(down_w)];
    end
    
%     [Cmn(d,:),Phimn(d,:),Smn,Smm,f,ConfC(d),PhiStd,Cerr(d,:,:)] = coherencyc_unequal_length_trials([down_w down_8],window, params, sMarkers );
    [Cmn3(d,:),Phimn3(d,:),Smn,Smm,f,ConfC3(d),PhiStd,Cerr3(d,:,:)] = coherencyc_unequal_length_trials([down_w down_3],window, params, sMarkers );
    [Cmn3s(d,:),Phimn3s(d,:),Smn,Smm,f,ConfC3s(d),PhiStd,Cerr3s(d,:,:)] = coherencyc_unequal_length_trials([down_w down3_sub],window, params, sMarkers );
       [Cmn2(d,:),Phimn2(d,:),Smn,Smm,f,ConfC2(d),PhiStd,Cerr2(d,:,:)] = coherencyc_unequal_length_trials([down_w down_2],window, params, sMarkers );
 

%     plot(f,Cmn(d,:),'linewidth',2)
%     hold on
%     line([0 45],[ConfC(d) ConfC(d)],'Color','k')
%     xlim([0 45])
%     tname = ['C:\WC_Germany\overall_calcs\coherence\w_band' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close
%     
%     plot(f,Cmn(d,:),'linewidth',2)
%     hold on
%     plot(f,squeeze(Cerr(d,1,:)),'--')
%     plot(f,squeeze(Cerr(d,2,:)),'--')
%     line([0 1],[ConfC(d) ConfC(d)],'Color','k')
%     xlim([0 1])
%     tname = ['C:\WC_Germany\overall_calcs\coherence\' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close

%      plot(f,Cmn3(d,:),'linewidth',2)
%     hold on
%     line([0 45],[ConfC3(d) ConfC3(d)],'Color','k')
%     xlim([0 45])
%     tname = ['C:\WC_Germany\overall_calcs\coherence\hipp_w_band' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close
%     
%     plot(f,Cmn3(d,:),'linewidth',2)
%     hold on
%     plot(f,squeeze(Cerr3(d,1,:)),'--')
%     plot(f,squeeze(Cerr3(d,2,:)),'--')
%     line([0 1],[ConfC3(d) ConfC3(d)],'Color','k')
%     xlim([0 1])
%     tname = ['C:\WC_Germany\overall_calcs\coherence\hipp_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close
%     
%       plot(f,Cmn3s(d,:),'linewidth',2)
%     hold on
%     line([0 45],[ConfC3s(d) ConfC3s(d)],'Color','k')
%     xlim([0 45])
%     tname = ['C:\WC_Germany\overall_calcs\coherence\orth_hipp_w_band' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close
%     
%     plot(f,Cmn3s(d,:),'linewidth',2)
%     hold on
%     plot(f,squeeze(Cerr3s(d,1,:)),'--')
%     plot(f,squeeze(Cerr3s(d,2,:)),'--')
%     line([0 1],[ConfC3s(d) ConfC3s(d)],'Color','k')
%     xlim([0 1])
%     tname = ['C:\WC_Germany\overall_calcs\coherence\orth_hipp_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close

      plot(f,Cmn2(d,:),'linewidth',2)
    hold on
    plot(f,squeeze(Cerr2(d,1,:)),'--')
    plot(f,squeeze(Cerr2(d,2,:)),'--')
    line([0 1],[ConfC2(d) ConfC2(d)],'Color','k')
    xlim([0 1])
    tname = ['C:\WC_Germany\overall_calcs\coherence\hipp_lf2_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

    clear down_w down_8 wcv* lf8 lf3 down_3*
    
end


save C:\WC_Germany\overall_calcs\coherence\coherence_data_hipp Cmn3* f Phimn3* Cerr3* ConfC3* Cmn2* Phimn2* Cerr2* ConfC2*