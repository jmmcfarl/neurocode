% clear all
% close all
%%
load F:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('F:\Code\WC_anal\general\')
addpath('F:\WC_Germany\Overall_EC\')

%%
dsf = 16;
Fsd = 2016/dsf;
niqf = 2016/2;
[b1,a1] = butter(2,[0.05/niqf 45/niqf]);
[b2,a2] = butter(2,[0.5/niqf 45/niqf]);

params_uds.Fs = Fsd;
params_uds.err = 0;
params_uds.tapers = [2 3];
params_uds.fpass = [0 10];
winlength_uds = 20;
winslide_uds = 2;
movingwin_uds = [winlength_uds winslide_uds];

params_the.Fs = Fsd;
params_the.err = 0;
params_the.tapers = [2 3];
params_the.fpass = [0 40];
winlength_the = 8;
winslide_the = 1;
movingwin_the = [winlength_uds winslide_uds];

for d = 1:length(sess_data)
    
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load used_data lf8 lf3 wcv_minus_spike

    lf8 = filtfilt(b1,a1,lf8);
    lf82 = filtfilt(b2,a2,lf8);
    lf8 = downsample(lf8,dsf);
    lf82 = downsample(lf82,dsf);
%         lf8 = zscore(lf8);
    lf8 = lf8/sess_data(d).gains(8);
    lf82 = lf82/sess_data(d).gains(8);

    lf3 = filtfilt(b1,a1,lf3);
    lf32 = filtfilt(b2,a2,lf3);
    lf3 = downsample(lf3,dsf);
    lf32 = downsample(lf32,dsf);
%         lf3 = zscore(lf3);
    lf3 = lf3/sess_data(d).gains(3);
    lf32 = lf32/sess_data(d).gains(3);

    wcv = filtfilt(b1,a1,wcv_minus_spike);
    wcv2 = filtfilt(b2,a2,wcv_minus_spike);
    wcv = downsample(wcv,dsf);
    wcv2 = downsample(wcv2,dsf);
    %     wcv = zscore(wcv);
    wcv = wcv/sess_data(d).gains(1);
    wcv2 = wcv2/sess_data(d).gains(1);
     
    [Pw_uds,t_uds,f_uds] = mtspecgramc(wcv,movingwin_uds,params_uds);
    [P8_uds,t_uds,f_uds] = mtspecgramc(lf8,movingwin_uds,params_uds);
    [P3_uds,t_uds,f_uds] = mtspecgramc(lf3,movingwin_uds,params_uds);
    [Pw_the,t_the,f_the] = mtspecgramc(wcv2,movingwin_the,params_the);
    [P8_the,t_the,f_the] = mtspecgramc(lf82,movingwin_the,params_the);
    [P3_the,t_the,f_the] = mtspecgramc(lf32,movingwin_the,params_the);
    
    [mp_upow_stats(d),mp_peak_upow,mp_tot_upow,mp_peak_ufreq] = ...
        get_power_stats(Pw_uds,f_uds,[0.05 1.5]);
    [lf8_upow_stats(d),lf8_peak_upow,lf8_tot_upow,lf8_peak_ufreq] = ...
        get_power_stats(P8_uds,f_uds,[0.05 1.5]);
    [lf3_upow_stats(d),lf3_peak_upow,lf3_tot_upow,lf3_peak_ufreq] = ...
        get_power_stats(P3_uds,f_uds,[0.05 1.5]);

    [mp_tpow_stats(d),mp_peak_tpow,mp_tot_tpow,mp_peak_tfreq] = ...
        get_power_stats(Pw_the,f_the,[1.5 6]);
    [lf8_tpow_stats(d),lf8_peak_tpow,lf8_tot_tpow,lf8_peak_tfreq] = ...
        get_power_stats(P8_the,f_the,[1.5 6]);
    [lf3_tpow_stats(d),lf3_peak_tpow,lf3_tot_tpow,lf3_peak_tfreq] = ...
        get_power_stats(P3_the,f_the,[1.5 6]);

    mp_peak_dtr(d) = mean(mp_peak_upow - mp_peak_tpow);
    mp_tot_dtr(d) = mean(mp_tot_upow./mp_tot_tpow);
    temp = corrcoef(mp_peak_ufreq,lf8_peak_ufreq);
    mp_lf8_uds_freq_cor(d) = temp(2,1);
    temp = corrcoef(mp_peak_upow,lf8_peak_upow);
    mp_lf8_uds_pow_cor(d) = temp(2,1);
    temp = corrcoef(mp_peak_tfreq,lf8_peak_tfreq);
    mp_lf8_theta_freq_cor(d) = temp(2,1);
    temp = corrcoef(mp_peak_tpow,lf8_peak_tpow);
    mp_lf8_theta_pow_cor(d) = temp(2,1);
    
    fprintf([sess_data(d).name ' specgram stats:\n']);
    fprintf('avg log peak UDS pow: %.2f: \n',mp_upow_stats(d).avg_lpow);   
    fprintf('avg total UDS pow: %.2f: \n',mp_upow_stats(d).avg_tpow);   
    fprintf('avg peak UDS freq: %.2f: \n',mp_upow_stats(d).avg_pfreq);   
    fprintf('avg log peak theta pow: %.2f: \n',mp_tpow_stats(d).avg_lpow);   
    fprintf('avg total theta pow: %.2f: \n',mp_tpow_stats(d).avg_tpow);   
    fprintf('avg peak theta freq: %.2f: \n',mp_tpow_stats(d).avg_pfreq);   
    fprintf('avg peak TDR: %.2f: \n',mp_peak_dtr(d));   
    fprintf('avg total TDR: %.2f: \n',mp_tot_dtr(d));   
    
    if isnan(sess_data(d).gains(3))
        P3_uds(:) = 0;
        P3_the(:) = 0;
        lf3_peak_upow(:) =0;
        lf3_peak_tpow(:) = 0;
        disp('LF3 not recorded')   
    end
    
    
%% generate figures
figure('visible','off');
subplot(3,4,[1 2])
pcolor(t_uds,f_uds,Pw_uds');shading flat;
ylim([0 4])
set(gca,'yscale','log')
caxis([0 prctile(mp_peak_upow,75)])
subplot(3,4,[5 6])
pcolor(t_uds,f_uds,P8_uds');shading flat;
ylim([0 4])
set(gca,'yscale','log')
caxis([0 prctile(lf8_peak_upow,75)])
subplot(3,4,[9 10])
pcolor(t_uds,f_uds,P3_uds');shading flat;
ylim([0 4])
set(gca,'yscale','log')
caxis([0 prctile(lf3_peak_upow,75)])
subplot(3,4,[3 4])
plot(t_uds,zscore(mp_peak_upow)), hold on
plot(t_uds,zscore(lf8_peak_upow),'r')
plot(t_uds,zscore(lf3_peak_upow),'k')
xlim([t_uds(1) t_uds(end)])
subplot(3,4,[7 8])
plot(t_uds,zscore(mp_tot_upow)), hold on
plot(t_uds,zscore(lf8_tot_upow),'r')
plot(t_uds,zscore(lf3_tot_upow),'k')
xlim([t_uds(1) t_uds(end)])
subplot(3,4,[11 12])
plot(t_uds,lf8_peak_ufreq,'r'), hold on
plot(t_uds,lf3_peak_ufreq,'k')
plot(t_uds,mp_peak_ufreq)
xlim([t_uds(1) t_uds(end)])
tname = ['F:\WC_Germany\overall_EC\specgrams\uds_spec_' s_name];
print('-dpng',tname);
close

figure('visible','off');
subplot(3,4,[1 2])
pcolor(t_the,f_the,log(Pw_the'));shading flat;
ylim([0 7])
caxis([-5 1])
subplot(3,4,[5 6])
pcolor(t_the,f_the,log(P8_the'));shading flat;
ylim([0 7])
caxis([-13 -6])
subplot(3,4,[9 10])
pcolor(t_the,f_the,log(P3_the'));shading flat;
ylim([0 7])
caxis([-13 -6])
subplot(3,4,[3 4])
plot(t_the,zscore(mp_peak_tpow)), hold on
plot(t_the,zscore(lf8_peak_tpow),'r')
plot(t_the,zscore(lf3_peak_tpow),'k')
xlim([t_the(1) t_the(end)])
subplot(3,4,[7 8])
plot(t_the,zscore(mp_tot_tpow)), hold on
plot(t_the,zscore(lf8_tot_tpow),'r')
plot(t_the,zscore(lf3_tot_tpow),'k')
xlim([t_the(1) t_the(end)])
subplot(3,4,[11 12])
plot(t_the,lf8_peak_tfreq,'r'), hold on
plot(t_the,lf3_peak_tfreq,'k')
plot(t_the,mp_peak_tfreq)
xlim([t_the(1) t_the(end)])
tname = ['F:\WC_Germany\overall_EC\specgrams\theta_spec_' s_name];
print('-dpng',tname);
close

end

cd F:\WC_Germany\overall_EC\
save overall_specgram_data mp_* lf8_* lf3_*
