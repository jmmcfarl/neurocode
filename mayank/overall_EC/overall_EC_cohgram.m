% clear all
% close all
%%
load G:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('G:\Code\WC_anal\general\')
addpath('G:\WC_Germany\Overall_EC\')

%%
dsf = 16;
Fsd = 2016/dsf;
niqf = 2016/2;
[b1,a1] = butter(2,[0.05/niqf 40/niqf]);
[b2,a2] = butter(2,[0.5/niqf 40/niqf]);

params_uds.Fs = Fsd;
params_uds.err = 0;
params_uds.tapers = [4 7];
params_uds.fpass = [0 10];
winlength_uds = 30;
winslide_uds = 3;
movingwin_uds = [winlength_uds winslide_uds];

params_the.Fs = Fsd;
params_the.err = 0;
params_the.tapers = [4 7];
params_the.fpass = [0 40];
winlength_the = 15;
winslide_the = 2;
movingwin_the = [winlength_uds winslide_uds];

for d = 1:length(sess_data)
    
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load used_data lf8 lf3 wcv_minus_spike

    lf8 = filtfilt(b1,a1,lf8);
    lf82 = filtfilt(b2,a2,lf8);
    lf8 = downsample(lf8,dsf);
    lf82 = downsample(lf82,dsf);

    lf3 = filtfilt(b1,a1,lf3);
    lf32 = filtfilt(b2,a2,lf3);
    lf3 = downsample(lf3,dsf);
    lf32 = downsample(lf32,dsf);

    wcv = filtfilt(b1,a1,wcv_minus_spike);
    wcv2 = filtfilt(b2,a2,wcv_minus_spike);
    wcv = downsample(wcv,dsf);
    wcv2 = downsample(wcv2,dsf);
     
    [Cw8_uds,phi,S12,S1,S2,t_uds,f_uds] = cohgramc(wcv,lf8,movingwin_uds,params_uds);
    [Cw3_uds,phi,S12,S1,S2,t_uds,f_uds] = cohgramc(wcv,lf3,movingwin_uds,params_uds);
    [C83_uds,phi,S12,S1,S2,t_uds,f_uds] = cohgramc(lf8,lf3,movingwin_uds,params_uds);
    [Cw8_the,phi,S12,S1,S2,t_the,f_the] = cohgramc(wcv,lf8,movingwin_the,params_the);
    [Cw3_the,phi,S12,S1,S2,t_the,f_the] = cohgramc(wcv,lf3,movingwin_the,params_the);
    [C83_the,phi,S12,S1,S2,t_the,f_the] = cohgramc(lf8,lf3,movingwin_the,params_the);
    
    theta_freqs = find(f_the > 1.5 & f_the < 6);
    uds_freqs = find(f_uds > 0.05 & f_uds < 1.5);
    max_Cw8_uds = max(Cw8_uds(:,uds_freqs),[],2);
    max_Cw8_the = max(Cw8_the(:,theta_freqs),[],2);
    max_Cw3_uds = max(Cw3_uds(:,uds_freqs),[],2);
    max_Cw3_the = max(Cw3_the(:,theta_freqs),[],2);
    max_C83_uds = max(C83_uds(:,uds_freqs),[],2);
    max_C83_the = max(C83_the(:,theta_freqs),[],2);
    
    avg_Cw8_uds(d) = tanh(mean(atanh(max_Cw8_uds)));
    avg_Cw3_uds(d) = tanh(mean(atanh(max_Cw3_uds)));
    avg_C83_uds(d) = tanh(mean(atanh(max_C83_uds)));
    avg_Cw8_theta(d) = tanh(mean(atanh(max_Cw8_the)));
    avg_Cw3_theta(d) = tanh(mean(atanh(max_Cw3_the)));
    avg_C83_theta(d) = tanh(mean(atanh(max_C83_the)));
        
    if isnan(sess_data(d).gains(3))
        Cw3_uds(:) = 0;
        C83_uds(:) = 0;
        max_Cw3_uds(:) =0;
        max_C83_uds(:) =0;
        max_Cw3_the(:) =0;
        max_C83_the(:) =0;
        disp('LF3 not recorded')   
    end
    
    
%% generate figures
subplot(2,2,1)
pcolor(t_uds,f_uds,Cw8_uds');shading flat;
ylim([0 4])
set(gca,'yscale','log')
caxis([0.2 0.9])
subplot(2,2,2)
pcolor(t_uds,f_uds,C83_uds');shading flat;
ylim([0 4])
set(gca,'yscale','log')
caxis([0.2 0.9])
subplot(2,2,3)
pcolor(t_uds,f_uds,Cw3_uds');shading flat;
ylim([0 4])
set(gca,'yscale','log')
caxis([0.2 0.9])
subplot(2,2,4)
plot(t_uds,max_Cw8_uds), hold on
plot(t_uds,max_Cw3_uds,'k')
plot(t_uds,max_C83_uds,'g')
xlim([t_uds(1) t_uds(end)])
tname = ['G:\WC_Germany\overall_EC\cohgrams\uds_coh_' s_name];
print('-dpng',tname);
close

subplot(2,2,1)
pcolor(t_the,f_the,Cw8_the');shading flat;
ylim([0 7])
caxis([0.2 0.8])
subplot(2,2,2)
pcolor(t_the,f_the,C83_the');shading flat;
ylim([0 7])
caxis([0.2 0.8])
subplot(2,2,3)
pcolor(t_the,f_the,Cw3_the');shading flat;
ylim([0 7])
caxis([0.2 0.8])
subplot(2,2,4)
plot(t_the,max_Cw8_the), hold on
plot(t_the,max_Cw3_the,'k')
plot(t_the,max_C83_the,'g')
xlim([t_the(1) t_the(end)])
tname = ['G:\WC_Germany\overall_EC\cohgrams\theta_coh_' s_name];
print('-dpng',tname);
close

end

cd G:\WC_Germany\overall_EC\
save overall_cohgram_data avg_*
