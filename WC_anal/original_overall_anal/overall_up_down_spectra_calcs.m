clear all
close all


load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\desynch_extract\desynch_points
load C:\WC_Germany\overall_calcs\overall_info_file_data
load C:\WC_Germany\overall_calcs\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data

dsf = 1;
Fsd = 2016/dsf;
params.Fs = Fsd;
params.tapers = [2 3];
params.err = [2 .05];
params.fpass = [0 20];
movingwin = [1.5 0.25];
% niqf = 2016/2;
% hcf1 = 100/niqf;
% [b1,a1] = butter(2,hcf1,'low');

for d = 1:62
d = 19;
    disp(sprintf('session %d',d))
    cd(over_dir{d});
    pwd

    load used_data wcv_minus_spike lf8 lf3 lf2

    %bandlimit signals
    %     down_w = filtfilt(b1,a1,wcv_minus_spike);
    %     down_8 = filtfilt(b1,a1,lf8);
    %     down_3 = filtfilt(b1,a1,lf3);
    %     down_2 = filtfilt(b1,a1,lf2);

    %correct for variable gains
    down_w = wcv_minus_spike/mp_gain(d);
    down_8 = lf8/lf8_gain(d);
%     down_3 = lf3/lf3_gain(d);
%     down_2 = lf2/lf2_gain(d);

    %     down_w = downsample(down_w,dsf);
    %     down_8 = downsample(down_8,dsf);
    %     down_3 = downsample(down_3,dsf);
    %     down_2 = downsample(down_2,dsf);


    up_markers8 = [up_trans8{d}(synch_ups8{d})'*8 down_trans8{d}(synch_ups8{d})'*8];
    up_markersw = [up_trans{d}(synch_ups{d})'*8 down_trans{d}(synch_ups{d})'*8];
    down_markers8 = [down_trans8{d}(synch_downs8{d}(1:end-1))'*8 up_trans8{d}(synch_downs8{d}(2:end))'*8];
    down_markersw = [down_trans{d}(synch_downs{d}(1:end-1))'*8 down_trans{d}(synch_downs{d}(2:end))'*8];


%% MP signal
    [Sw_up{d},f,Swerr_up{d}]= mtspectrumc_unequal_length_trials(down_w, movingwin, params, up_markersw);
    [Sw_down{d},f,Swerr_down{d}]= mtspectrumc_unequal_length_trials(down_w, movingwin, params, down_markersw);
%     [Sw_lup{d},f,Swerr_lup{d}]= mtspectrumc_unequal_length_trials(down_w, movingwin, params, up_markers8);
%     [Sw_ldown{d},f,Swerr_ldown{d}]= mtspectrumc_unequal_length_trials(down_w, movingwin, params, down_markers8);

    plot(f,10*log10(Sw_up{d}))
    hold on
    plot(f,10*log10(Sw_down{d}),'r')
    legend('MP up state','MP down state')
    plot(f,10*log10(Swerr_up{d}(1,:)),'--')
        plot(f,10*log10(Swerr_up{d}(2,:)),'--')
    plot(f,10*log10(Swerr_down{d}(1,:)),'r--')
        plot(f,10*log10(Swerr_down{d}(2,:)),'r--')
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\mp_mpstate_wb_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    xlim([0 20])
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\mp_mpstate_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
     
        plot(f,10*log10(Sw_lup{d}))
    hold on
    plot(f,10*log10(Sw_ldown{d}),'r')
    legend('LFP up state','LFP down state')
    plot(f,10*log10(Swerr_lup{d}(1,:)),'--')
        plot(f,10*log10(Swerr_lup{d}(2,:)),'--')
    plot(f,10*log10(Swerr_ldown{d}(1,:)),'r--')
        plot(f,10*log10(Swerr_ldown{d}(2,:)),'r--')
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\mp_lfpstate_wb_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    xlim([0 20])
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\mp_lfpstate_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
    
        
%% LFP signal
    [S8_up{d},f,S8err_up{d}]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, up_markersw);
    [S8_down{d},f,S8err_down{d}]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, down_markersw);
    [S8_lup{d},f,S8err_lup{d}]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, up_markers8);
    [S8_ldown{d},f,S8err_ldown{d}]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, down_markers8);

        plot(f,10*log10(S8_up{d}))
    hold on
    plot(f,10*log10(S8_down{d}),'r')
    legend('MP up state','MP down state')
    plot(f,10*log10(S8err_up{d}(1,:)),'--')
        plot(f,10*log10(S8err_up{d}(2,:)),'--')
    plot(f,10*log10(S8err_down{d}(1,:)),'r--')
        plot(f,10*log10(S8err_down{d}(2,:)),'r--')
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\lf8_mpstate_wb_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    xlim([0 20])
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\lf8_mpstate_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
     
        plot(f,10*log10(S8_lup{d}))
    hold on
    plot(f,10*log10(S8_ldown{d}),'r')
    legend('LFP up state','LFP down state')
    plot(f,10*log10(S8err_lup{d}(1,:)),'--')
        plot(f,10*log10(S8err_lup{d}(2,:)),'--')
    plot(f,10*log10(S8err_ldown{d}(1,:)),'r--')
        plot(f,10*log10(S8err_ldown{d}(2,:)),'r--')
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\lf8_lfpstate_wb_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    xlim([0 20])
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\lf8_lfpstate_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

%% MP LF8 coh
%     [Cmn8_up{d},Phimn8_up{d},Smn,Smm,f,ConfC8_up{d},PhiStd_up,Cerr8_up{d}] = coherencyc_unequal_length_trials([down_8 down_w], movingwin, params,up_markersw );
%     [Cmn8_down{d},Phimn8_down{d},Smn,Smm,f,ConfC8_down{d},PhiStd_down{d},Cerr8_down{d}] = coherencyc_unequal_length_trials([down_8 down_w], movingwin, params,down_markersw );
%     [Cmn8_lup{d},Phimn8_lup{d},Smn,Smm,f,ConfC8_lup{d},PhiStd_lup{d},Cerr8_lup{d}] = coherencyc_unequal_length_trials([down_8 down_w], movingwin, params,up_markers8 );
%     [Cmn8_ldown{d},Phimn8_ldown{d},Smn,Smm,f,ConfC8_ldown{d},PhiStd_ldown{d},Cerr8_ldown{d}] = coherencyc_unequal_length_trials([down_8 down_w], movingwin, params,down_markers8 );
%    
%     plot(f,Cmn8_up{d})
%     hold on
%     plot(f,Cmn8_down{d},'r')
%     legend('MP up state','MP down state')
%     plot(f,Cerr8_up{d}(1,:),'--')
%     plot(f,Cerr8_up{d}(2,:),'--')
%     plot(f,Cerr8_down{d}(1,:),'r--')
%     plot(f,Cerr8_down{d}(2,:),'r--')
%     xlim([0 40])
%     tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\coh8_mpstate_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close
%      
%     plot(f,Cmn8_lup{d})
%     hold on
%     plot(f,Cmn8_ldown{d},'r')
%     legend('LFP up state','LFP down state')
%     plot(f,Cerr8_lup{d}(1,:),'--')
%     plot(f,Cerr8_lup{d}(2,:),'--')
%     plot(f,Cerr8_ldown{d}(1,:),'r--')
%     plot(f,Cerr8_ldown{d}(2,:),'r--')
%     xlim([0 40])
%     tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\coh8_lfpstate_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close

%% HIpp LFP
    if ~bad_lf3(d)

        [S3_up{d},f,S3err_up{d}]= mtspectrumc_unequal_length_trials(down_3, movingwin, params, up_markersw);
        [S3_down{d},f,S3err_down{d}]= mtspectrumc_unequal_length_trials(down_3, movingwin, params, down_markersw);
        [S3_lup{d},f,S3err_lup{d}]= mtspectrumc_unequal_length_trials(down_3, movingwin, params, up_markers8);
        [S3_ldown{d},f,S3err_ldown{d}]= mtspectrumc_unequal_length_trials(down_3, movingwin, params, down_markers8);

    plot(f,10*log10(S3_up{d}))
    hold on
    plot(f,10*log10(S3_down{d}),'r')
    legend('MP up state','MP down state')
    plot(f,10*log10(S3err_up{d}(1,:)),'--')
        plot(f,10*log10(S3err_up{d}(2,:)),'--')
    plot(f,10*log10(S3err_down{d}(1,:)),'r--')
        plot(f,10*log10(S3err_down{d}(2,:)),'r--')
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\lf3_mpstate_wb_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    xlim([0 20])
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\lf3_mpstate_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
     
        plot(f,10*log10(S3_lup{d}))
    hold on
    plot(f,10*log10(S3_ldown{d}),'r')
    legend('LFP up state','LFP down state')
    plot(f,10*log10(S3err_lup{d}(1,:)),'--')
        plot(f,10*log10(S3err_lup{d}(2,:)),'--')
    plot(f,10*log10(S3err_ldown{d}(1,:)),'r--')
        plot(f,10*log10(S3err_ldown{d}(2,:)),'r--')
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\lf3_lfpstate_wb_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    xlim([0 20])
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\lf3_lfpstate_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

%         [Cmn3_up{d},Phimn3_up{d},Smn,Smm,f,ConfC3_up{d},PhiStd_up,Cerr3_up{d}] = coherencyc_unequal_length_trials([down_3 down_w], movingwin, params,up_markersw );
%         [Cmn3_down{d},Phimn3_down{d},Smn,Smm,f,ConfC3_down{d},PhiStd_down{d},Cerr3_down{d}] = coherencyc_unequal_length_trials([down_3 down_w], movingwin, params,down_markersw );
%         [Cmn3_lup{d},Phimn3_lup{d},Smn,Smm,f,ConfC3_lup{d},PhiStd_lup{d},Cerr3_lup{d}] = coherencyc_unequal_length_trials([down_3 down_w], movingwin, params,up_markers8 );
%         [Cmn3_ldown{d},Phimn3_ldown{d},Smn,Smm,f,ConfC3_ldown{d},PhiStd_ldown{d},Cerr3_ldown{d}] = coherencyc_unequal_length_trials([down_3 down_w], movingwin, params,down_markers8 );
% 
%             plot(f,Cmn3_up{d})
%     hold on
%     plot(f,Cmn3_down{d},'r')
%     legend('MP up state','MP down state')
%     plot(f,Cerr3_up{d}(1,:),'--')
%     plot(f,Cerr3_up{d}(2,:),'--')
%     plot(f,Cerr3_down{d}(1,:),'r--')
%     plot(f,Cerr3_down{d}(2,:),'r--')
%     xlim([0 40])
%     tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\coh3_mpstate_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close
%      
%     plot(f,Cmn3_lup{d})
%     hold on
%     plot(f,Cmn3_ldown{d},'r')
%     legend('LFP up state','LFP down state')
%     plot(f,Cerr3_lup{d}(1,:),'--')
%     plot(f,Cerr3_lup{d}(2,:),'--')
%     plot(f,Cerr3_ldown{d}(1,:),'r--')
%     plot(f,Cerr3_ldown{d}(2,:),'r--')
%     xlim([0 40])
%     tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\coh3_lfpstate_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close

        
    end

    
%% Deep hipp LFP

    if ~bad_lf2(d)

        [S2_up{d},f,S2err_up{d}]= mtspectrumc_unequal_length_trials(down_2, movingwin, params, up_markersw);
        [S2_down{d},f,S2err_down{d}]= mtspectrumc_unequal_length_trials(down_2, movingwin, params, down_markersw);
        [S2_lup{d},f,S2err_lup{d}]= mtspectrumc_unequal_length_trials(down_2, movingwin, params, up_markers8);
        [S2_ldown{d},f,S2err_ldown{d}]= mtspectrumc_unequal_length_trials(down_2, movingwin, params, down_markers8);

        
    plot(f,10*log10(S2_up{d}))
    hold on
    plot(f,10*log10(S2_down{d}),'r')
    legend('MP up state','MP down state')
    plot(f,10*log10(S2err_up{d}(1,:)),'--')
        plot(f,10*log10(S2err_up{d}(2,:)),'--')
    plot(f,10*log10(S2err_down{d}(1,:)),'r--')
        plot(f,10*log10(S2err_down{d}(2,:)),'r--')
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\lf2_mpstate_wb_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    xlim([0 20])
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\lf2_mpstate_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
     
        plot(f,10*log10(S2_lup{d}))
    hold on
    plot(f,10*log10(S2_ldown{d}),'r')
    legend('LFP up state','LFP down state')
    plot(f,10*log10(S2err_lup{d}(1,:)),'--')
        plot(f,10*log10(S2err_lup{d}(2,:)),'--')
    plot(f,10*log10(S2err_ldown{d}(1,:)),'r--')
        plot(f,10*log10(S2err_ldown{d}(2,:)),'r--')
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\lf2_lfpstate_wb_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    xlim([0 20])
    tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\lf2_lfpstate_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

%         [Cmn2_up{d},Phimn2_up{d},Smn,Smm,f,ConfC2_up{d},PhiStd_up,Cerr2_up{d}] = coherencyc_unequal_length_trials([down_2 down_w], movingwin, params,up_markersw );
%         [Cmn2_down{d},Phimn2_down{d},Smn,Smm,f,ConfC2_down{d},PhiStd_down{d},Cerr2_down{d}] = coherencyc_unequal_length_trials([down_2 down_w], movingwin, params,down_markersw );
%         [Cmn2_lup{d},Phimn2_lup{d},Smn,Smm,f,ConfC2_lup{d},PhiStd_lup{d},Cerr2_lup{d}] = coherencyc_unequal_length_trials([down_2 down_w], movingwin, params,up_markers8 );
%         [Cmn2_ldown{d},Phimn2_ldown{d},Smn,Smm,f,ConfC2_ldown{d},PhiStd_ldown{d},Cerr2_ldown{d}] = coherencyc_unequal_length_trials([down_2 down_w], movingwin, params,down_markers8 );
% 
%         
%                     plot(f,Cmn3_up{d})
%     hold on
%     plot(f,Cmn2_down{d},'r')
%     legend('MP up state','MP down state')
%     plot(f,Cerr2_up{d}(1,:),'--')
%     plot(f,Cerr2_up{d}(2,:),'--')
%     plot(f,Cerr2_down{d}(1,:),'r--')
%     plot(f,Cerr2_down{d}(2,:),'r--')
%     xlim([0 40])
%     tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\coh2_mpstate_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close
%      
%     plot(f,Cmn2_lup{d})
%     hold on
%     plot(f,Cmn2_ldown{d},'r')
%     legend('LFP up state','LFP down state')
%     plot(f,Cerr2_lup{d}(1,:),'--')
%     plot(f,Cerr2_lup{d}(2,:),'--')
%     plot(f,Cerr2_ldown{d}(1,:),'r--')
%     plot(f,Cerr2_ldown{d}(2,:),'r--')
%     xlim([0 40])
%     tname = ['C:\WC_Germany\overall_calcs\up_down_spectra\coh2_lfpstate_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',tname);
%     close


    end
    
    
    
end