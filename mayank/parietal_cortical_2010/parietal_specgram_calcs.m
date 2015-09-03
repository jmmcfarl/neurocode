clear all
close all

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8

dsf = 16;
Fsd = 2016/dsf;

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

params.Fs = Fsd;
params.err = 0;
params.tapers = [2 3];
params.fpass = [0 15];

winlength = 20;
winslide = 5;
noverlap = winlength-winslide;
movingwin = [winlength winslide];

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_start_times(interneurons) = [];
desynch_stop_times(interneurons) = [];


for d = 1:length(sess_data)


    cd(sess_data(d).directory)
    disp(num2str(d))

    load used_data wcv_minus_spike lf8

    down_w = filtfilt(b,a,wcv_minus_spike);
    down_8 = filtfilt(b,a,lf8);

    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);

    down_w = zscore(down_w);
    down_8 = zscore(down_8);

    [Pw{d},t{d},f{d}]=mtspecgramc(down_w,movingwin,params);
    [P8{d},t{d},f{d}]=mtspecgramc(down_8,movingwin,params);

    logpow_w = 10*log10(Pw{d});
    logpow_8 = 10*log10(P8{d});

    theta_freqs = find(f{d} > 3 & f{d} < 8);
    uds_freqs = find(f{d} < 2);
    dw = f{d}(2)-f{d}(1);
    tdr_w{d} = trapz(logpow_w(:,theta_freqs),2)./trapz(logpow_w(:,uds_freqs),2);
    tdr_8{d} = trapz(logpow_8(:,theta_freqs),2)./trapz(logpow_8(:,uds_freqs),2);
    uds_pow_w{d} = trapz(logpow_w(:,uds_freqs),2)*dw;
    uds_pow_8{d} = trapz(logpow_8(:,uds_freqs),2)*dw;
    uds_comfreq_w{d} = sum(logpow_w(:,uds_freqs).*repmat(f{d}(uds_freqs),length(t{d}),1),2)./sum(logpow_w(:,uds_freqs),2);
    uds_comfreq_8{d} = sum(logpow_8(:,uds_freqs).*repmat(f{d}(uds_freqs),length(t{d}),1),2)./sum(logpow_8(:,uds_freqs),2);
    
    avg_comfreq_w(d) = mean(uds_comfreq_w{d});
    avg_comfreq_8(d) = mean(uds_comfreq_8{d});
    avg_uds_pow_w(d) = mean(uds_pow_w{d});
    avg_uds_pow_8(d) = mean(uds_pow_8{d});
    avg_tdr_w(d) = mean(tdr_w{d});
    avg_tdr_8(d) = mean(tdr_8{d});
    
    Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(3,1,1)
    pcolor(t{d},f{d},10*log10(Pw{d}'));shading flat;
    caxis([-35 2]);
    ylim([0 15]), set(gca,'yscale','log')
    subplot(3,1,2)
     pcolor(t{d},f{d},10*log10(P8{d}'));shading flat;
    caxis([-35 2]);
    ylim([0 15]), set(gca,'yscale','log')
    subplot(3,1,3)
   plot(t{d},tdr_w{d},t{d},tdr_8{d},'r')
   xlim([t{d}(1) t{d}(end)])
    cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
    tname = ['G:\WC_Germany\parietal_cortical_2010\specgram\' cell_name];
    print('-dpng',tname);
    close


end
cd G:\WC_Germany\parietal_cortical_2010\
save specgram_data uds_* tdr_* f t 