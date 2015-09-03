clear all
close all
%%
load G:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('G:\Code\WC_anal\general\')
addpath('G:\WC_Germany\Overall_EC\')
addpath('G:\Code\Chronux\spectral_analysis\continuous\')
addpath('G:\WC_Germany\hsmm_state_detection\\')
addpath('G:\Code\GCCA_toolbox_jan2011')
drive_letter = 'G';

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

% dsf = 15;
% Fsd = 2016/dsf;
% niqf = 2016/2;
% [b,a] = butter(2,[0.05/niqf 40/niqf]);
dsf = 100;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,[0.005/niqf 2.5/niqf]);

f_int = 0.01:.02:10;
%%
for d = 1:length(sess_data)
%      d=14;
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    load used_data lf8 lf3 wcv_minus_spike lf2
    lf3 = lf2;
    wcv_d = filtfilt(b,a,wcv_minus_spike);
    wcv_d = downsample(wcv_d,dsf)/sess_data(d).gains(1);
    lf8_d = filtfilt(b,a,lf8);
    lf8_d = downsample(lf8_d,dsf)/sess_data(d).gains(8);
    lf3_d = filtfilt(b,a,lf3);
    lf3_d = downsample(lf3_d,dsf)/sess_data(d).gains(3);
    t_axis = (1:length(wcv_d))/Fsd;
    
    wcv_d = zscore(wcv_d);
    lf8_d = zscore(lf8_d);
    lf3_d = zscore(lf3_d);
    
    
    %     load ./pa_hsmm_state_seq8
    %     [new_seg_inds] = resample_uds_seg_inds(hmm8.UDS_segs,hmm8.Fs,Fsd,length(lf8_d));
    
    tot_X = [wcv_d(:) lf8_d(:) lf3_d(:)]';
%     tot_X = [lf8_d(:) lf3_d(:)]';
    tot_N = size(tot_X,2);
    nvar = size(tot_X,1);
    MINP = 35;
    MAXP = 50;
    PVAL    =   0.01;       % probability threshold for Granger causality significance
    NLAGS   =   -1;         % if -1, best model order is assessed automatically
    Fs      =   Fsd;        % sampling frequency  (for spectral analysis only)
    freqs   =   [0.02:0.02:10];    % frequency range to analyze (spectral analysis only)
%     
%                 [bic,aic] = cca_find_model_order(tot_X,MINP,MAXP);
%             disp(['best model order by Bayesian Information Criterion = ',num2str(bic)]);
%             disp(['best model order by Aikaike Information Criterion = ',num2str(aic)]);

    
    clear parameters
    parameters.Fs = Fsd;
    parameters.h0Hz = 0.05;
    savefilename = '';
    parameters.analyzing = 1;
    parameters.latexing = 0;
    parameters.resultspath = 'G:\Code\FDMa';
    parameters.silent = 1;
    parameters.Ldo = 1;
    parameters.f0 = 0;
    parameters.siglev0 = 0.01;
    parameters.P0 = 40;
    parameters.FDMaversion = 1;
    parameters.plotting = 0;
    parameters.Pmax = 80;
    parameters.deltaT = 50;
    parameters.deltaF = 0.1;
    parameters.plotmaxfreq = 2;
    parameters.dynrange = 50;
    parameters.plot_diag = 1;
    parameters.plot_zoom = 0;
    parameters.plot_grid = 0;
    parameters.plot_fig = 0;
    parameters.plot_diagnosis=1;
    parameters.plot_minimal_labels = 1;
    parameters.plot_f0 = 0;
    parameters.graphing = 0;
    parameters.position_SIG= [100 720 960 360];
    parameters.position_SMP= [960 720 480 360];
    parameters.position_VAR= [1440 720 480 360];
    parameters.position_COH= [100 360 640 360];
    parameters.position_PDC= [640 360 640 360];
    parameters.node_labels = {'MP' 'LF8' 'LF3'};
    
    results= FDMa(tot_X', parameters);
    f = (1:ceil(tot_N/2))/tot_N*Fsd;
    
    for j = 1:3
        for k = 1:3
            PDC_mat(d,j,k,:) = interp1(f,squeeze(results.PDC(j,k,:)),f_int);
            PDCsig_mat(d,j,k,:) = interp1(f,squeeze(results.PDCsig(j,k,:)),f_int);
            COH_mat(d,j,k,:) = interp1(f,squeeze(results.COH(j,k,:)),f_int);
            SMP_mat(d,j,k,:) = interp1(f,squeeze(results.SMP(j,k,:)),f_int);
            VAR_mat(d,j,k,:) = interp1(f,squeeze(results.VAR(j,k,:)),f_int);
        end
    end
    Phat(d) = results.Phat;
%     [results, parameters]= FDMa_analyze(X0, parameters)

%%
figure
subplot(3,3,1)
plot(f_int,squeeze(10*log10(SMP_mat(d,1,1,:))));
xlim([0 2])
% ylim([0 26])
hold on
plot(f_int,squeeze(10*log10(VAR_mat(d,1,1,:))),'r');
subplot(3,3,5)
plot(f_int,squeeze(10*log10(SMP_mat(d,2,2,:))));
xlim([0 2])
% ylim([0 26])
hold on
plot(f_int,squeeze(10*log10(VAR_mat(d,2,2,:))),'r');
subplot(3,3,9)
plot(f_int,squeeze(10*log10(SMP_mat(d,3,3,:))));
xlim([0 2])
% ylim([0 26])
hold on
plot(f_int,squeeze(10*log10(VAR_mat(d,3,3,:))),'r');
subplot(3,3,2)
plot(f_int,squeeze(PDC_mat(d,1,2,:)));
hold on
plot(f_int,squeeze(COH_mat(d,2,1,:)),'r');
plot(f_int,squeeze(COH_mat(d,1,2,:)),'k');
xlim([0 2])
ylim([0 1])
subplot(3,3,4)
plot(f_int,squeeze(PDC_mat(d,2,1,:)));
hold on
plot(f_int,squeeze(COH_mat(d,2,1,:)),'r');
plot(f_int,squeeze(COH_mat(d,1,2,:)),'k');
xlim([0 2])
ylim([0 1])
subplot(3,3,3)
plot(f_int,squeeze(PDC_mat(d,1,3,:)));
hold on
plot(f_int,squeeze(COH_mat(d,3,1,:)),'r');
plot(f_int,squeeze(COH_mat(d,1,3,:)),'k');
xlim([0 2])
ylim([0 1])
subplot(3,3,7)
plot(f_int,squeeze(PDC_mat(d,3,1,:)));
hold on
plot(f_int,squeeze(COH_mat(d,3,1,:)),'r');
plot(f_int,squeeze(COH_mat(d,1,3,:)),'k');
xlim([0 2])
ylim([0 1])
subplot(3,3,6)
plot(f_int,squeeze(PDC_mat(d,2,3,:)));
hold on
plot(f_int,squeeze(COH_mat(d,3,2,:)),'r');
plot(f_int,squeeze(COH_mat(d,2,3,:)),'k');
xlim([0 2])
ylim([0 1])
subplot(3,3,8)
plot(f_int,squeeze(PDC_mat(d,3,2,:)));
hold on
plot(f_int,squeeze(COH_mat(d,3,2,:)),'r');
plot(f_int,squeeze(COH_mat(d,2,3,:)),'k');
xlim([0 2])
ylim([0 1])
t_name = ['G:\WC_Germany\persistent_9_27_2010\caus_check\pdc_lf2_' s_name];
print('-dpng',t_name),close

%%


end