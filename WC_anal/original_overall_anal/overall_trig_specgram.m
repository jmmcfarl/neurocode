clear all
close all


load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\desynch_extract\desynch_points
load C:\WC_Germany\overall_calcs\overall_info_file_data
load C:\WC_Germany\overall_calcs\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data

dsf = 8;
Fsd = 2016/dsf;
params.Fs = 2016;
params.tapers = [1 1];
params.err = 0;
params.fpass = [0 200];
params.fpad = 1;
params.trialave = 1;
win = [0.5 0.5];
movingwin = [0.1 0.02];
% niqf = 2016/2;
% hcf1 = 100/niqf;
% [b1,a1] = butter(2,hcf1,'low');

for d = 1:62

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
    down_3 = lf3/lf3_gain(d);
    down_2 = lf2/lf2_gain(d);

    %     down_w = downsample(down_w,dsf);
    %     down_8 = downsample(down_8,dsf);
    %     down_3 = downsample(down_3,dsf);
    %     down_2 = downsample(down_2,dsf);


    up_markers8 = up_trans8{d}(synch_ups8{d})'/Fsd;
    up_markersw = up_trans{d}(synch_ups{d})'/Fsd;
    down_markers8 = down_trans8{d}(synch_downs8{d}(1:end-1))'/Fsd;
    down_markersw = down_trans{d}(synch_downs{d}(1:end-1))'/Fsd;


%% MP signal
    [Sw_up{d},t,f]=mtspecgramtrigc(down_w,up_markersw,win,movingwin,params); 
    [Sw_down{d},t,f]=mtspecgramtrigc(down_w,down_markersw,win,movingwin,params);
    [Sw_lup{d},t,f]=mtspecgramtrigc(down_w,up_markers8,win,movingwin,params);
    [Sw_ldown{d},t,f]=mtspecgramtrigc(down_w,down_markers8,win,movingwin,params);

    pcolor(t-0.5,f,10*log10(Sw_up{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\mp_mpup_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
 
        pcolor(t-0.5,f,10*log10(Sw_down{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\mp_mpdown_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

        pcolor(t-0.5,f,10*log10(Sw_lup{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\mp_lfpup_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

        pcolor(t-0.5,f,10*log10(Sw_ldown{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\mp_lfpdown_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

%% LFP signal
    
        [S8_up{d},t,f]=mtspecgramtrigc(down_8,up_markersw,win,movingwin,params); 
    [S8_down{d},t,f]=mtspecgramtrigc(down_8,down_markersw,win,movingwin,params);
    [S8_lup{d},t,f]=mtspecgramtrigc(down_8,up_markers8,win,movingwin,params);
    [S8_ldown{d},t,f]=mtspecgramtrigc(down_8,down_markers8,win,movingwin,params);

      pcolor(t-0.5,f,10*log10(S8_up{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\lf8_mpup_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
 
        pcolor(t-0.5,f,10*log10(S8_down{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\lf8_mpdown_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

        pcolor(t-0.5,f,10*log10(S8_lup{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\lf8_lfpup_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

        pcolor(t-0.5,f,10*log10(S8_ldown{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\lf8_lfpdown_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
  
%% LF3 signal
    
        [S3_up{d},t,f]=mtspecgramtrigc(down_3,up_markersw,win,movingwin,params); 
    [S3_down{d},t,f]=mtspecgramtrigc(down_3,down_markersw,win,movingwin,params);
    [S3_lup{d},t,f]=mtspecgramtrigc(down_3,up_markers8,win,movingwin,params);
    [S3_ldown{d},t,f]=mtspecgramtrigc(down_3,down_markers8,win,movingwin,params);

      pcolor(t-0.5,f,10*log10(S8_up{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\lf3_mpup_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
 
        pcolor(t-0.5,f,10*log10(S8_down{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\lf3_mpdown_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

        pcolor(t-0.5,f,10*log10(S8_lup{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\lf3_lfpup_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

        pcolor(t-0.5,f,10*log10(S8_ldown{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\lf3_lfpdown_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
    
%% LF2 signal
    
        [S2_up{d},t,f]=mtspecgramtrigc(down_2,up_markersw,win,movingwin,params); 
    [S2_down{d},t,f]=mtspecgramtrigc(down_2,down_markersw,win,movingwin,params);
    [S2_lup{d},t,f]=mtspecgramtrigc(down_2,up_markers8,win,movingwin,params);
    [S2_ldown{d},t,f]=mtspecgramtrigc(down_2,down_markers8,win,movingwin,params);

      pcolor(t-0.5,f,10*log10(S8_up{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\lf2_mpup_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close
 
        pcolor(t-0.5,f,10*log10(S8_down{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\lf2_mpdown_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

        pcolor(t-0.5,f,10*log10(S8_lup{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\lf2_lfpup_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

        pcolor(t-0.5,f,10*log10(S8_ldown{d}'));shading flat
    tname = ['C:\WC_Germany\overall_calcs\trig_specgram\lf2_lfpdown_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

end

    save C:\WC_Germany\overall_calcs\trig_specgram\trig_specgram_data t f S*
