clear all
close all
cd C:\WC_Germany\JMM_analysis_ste
load dir_tree_ste

dsf = 2;
Fsd = 2016/dsf;

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 300/niqf;
[b,a] = butter(2,[lcf hcf]);


params.Fs = Fsd;
params.err = 0;
params.tapers = [3 5];
params.fpass = [0 300];

winlength = 20;
winslide = 2;
movingwin = [winlength winslide];

for d = 1:length(dir_array)
    disp(sprintf('session %d',d))
    cd(dir_array{d});
    load used_data wcv_minus_spike lf8 lf3 lf2
% load used_data lf2 lf3 lf8

    %bandlimit signals
    down_w = filtfilt(b,a,wcv_minus_spike);
    down_8 = filtfilt(b,a,lf8);
    down_3 = filtfilt(b,a,lf3);
    down_2 = filtfilt(b,a,lf2);

down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
    down_3 = downsample(down_3,dsf);
    down_2 = downsample(down_2,dsf);

    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_3 = zscore(down_3);
    down_2 = zscore(down_2);

    [Pw,t,f]=mtspecgramc(down_w,movingwin,params);
    [P8,t,f]=mtspecgramc(down_8,movingwin,params);
    [P3,t,f]=mtspecgramc(down_3,movingwin,params);
    [P2,t,f]=mtspecgramc(down_2,movingwin,params);
   
    lPw = 10*log10(Pw)-repmat(mean(10*log10(Pw)),length(t),1);
    lP8 = 10*log10(P8)-repmat(mean(10*log10(P8)),length(t),1);
    lP3 = 10*log10(P3)-repmat(mean(10*log10(P3)),length(t),1);
    lP2 = 10*log10(P2)-repmat(mean(10*log10(P2)),length(t),1);

%     [Cw8,phi,S12,S1,S2,t,f]=cohgramc(down_w,down_8,movingwin,params);
%      [Cw3,phi,S12,S1,S2,t,f]=cohgramc(down_w,down_3,movingwin,params);
   
    
    Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(t,f,10*log10(Pw'));shading flat;
    pcolor(t,f,lPw');shading flat;
%     caxis([-25 2]);
caxis([-10 10])
    tname = ['C:\WC_Germany\JMM_analysis_ste\new_spec_anal\mp_ms_' f_names{d}];
    print('-dpng',tname);
    set(gca,'yscale','log')
        tname = ['C:\WC_Germany\JMM_analysis_ste\new_spec_anal\mp_ms_logy_' f_names{d}];
    print('-dpng',tname);
    close
         
    Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(t,f,10*log10(P8'));shading flat;
    pcolor(t,f,lP8');shading flat;
%     caxis([-25 2]);
caxis([-10 10])
    tname = ['C:\WC_Germany\JMM_analysis_ste\new_spec_anal\lf8_ms_' f_names{d}];
    print('-dpng',tname);
    set(gca,'yscale','log')
        tname = ['C:\WC_Germany\JMM_analysis_ste\new_spec_anal\lf8_ms_logy_' f_names{d}];
    print('-dpng',tname);
    close

    Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(t,f,10*log10(P3'));shading flat;
    pcolor(t,f,lP3');shading flat;
%     caxis([-25 2]);
caxis([-10 10])

    tname = ['C:\WC_Germany\JMM_analysis_ste\new_spec_anal\lf3_ms_' f_names{d}];
    print('-dpng',tname);
    set(gca,'yscale','log')
        tname = ['C:\WC_Germany\JMM_analysis_ste\new_spec_anal\lf3_ms_logy_' f_names{d}];
    print('-dpng',tname);
    close


    Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(t,f,10*log10(P2'));shading flat;
    pcolor(t,f,lP2');shading flat;
%     caxis([-25 2]);
caxis([-10 10])

    tname = ['C:\WC_Germany\JMM_analysis_ste\new_spec_anal\lf2_ms_' f_names{d}];
    print('-dpng',tname);
    set(gca,'yscale','log')
        tname = ['C:\WC_Germany\JMM_analysis_ste\new_spec_anal\lf2_ms_logy_' f_names{d}];
    print('-dpng',tname);
    close

%         Fig = figure(1);
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(t,f,Cw8');shading flat;
%     caxis([0 1]);
%     tname = ['C:\WC_Germany\persistent_revised\new_spec_anal\mp_lf8_coh_' f_names{d}];
%     print('-dpng',tname);
%     close
% 
%             Fig = figure(1);
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     pcolor(t,f,Cw3');shading flat;
%     caxis([0 1]);
%     tname = ['C:\WC_Germany\persistent_revised\new_spec_anal\mp_lf3_coh_' f_names{d}];
%     print('-dpng',tname);
%     close

clear P* lP*

end


