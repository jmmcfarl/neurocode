clear all
close all
%%
load F:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('F:\Code\WC_anal\general\')
addpath('F:\WC_Germany\Overall_EC\')

%%
dsf = 16;
Fsd = 2016/dsf;
niqf = 2016/2;
[b1,a1] = butter(2,[0.05/niqf 45/niqf]);

params.Fs = Fsd;
params.err = 0;
params.fpass = [0 40];
winlength = 8;
winslide = 2;
movingwin = [winlength winslide];

for d = 66:length(sess_data)

    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);

%     load used_data lf3 lf2 lf5 lf8 wcv_minus_spike
    load used_data wcv_minus_spike

%     lf2 = filtfilt(b1,b1,lf2);
%     lf2 = downsample(lf2,dsf)/sess_data(d).gains(2);
%     
%     lf3 = filtfilt(b1,a1,lf3);
%     lf3 = downsample(lf3,dsf)/sess_data(d).gains(3);
%     
%      lf5 = filtfilt(b1,a1,lf5);
%     lf5 = downsample(lf5,dsf)/sess_data(d).gains(5);
% 
%     lf8 = filtfilt(b1,a1,lf8);
%     lf8 = zscore(downsample(lf8,dsf));
% 
%     lf2_r = lf2-lf5;
    
    wcv = filtfilt(b1,a1,wcv_minus_spike);
    wcv = downsample(wcv,dsf)/sess_data(d).gains(1);

%     params.tapers = [4 7];
%     [Cw3,phi,S12,S1,S2,t,f] = cohgramc(wcv,lf3,movingwin,params);
%     [Cw2,phi,S12,S1,S2,t,f] = cohgramc(wcv,lf2,movingwin,params);
%     [Cw2r,phi,S12,S1,S2,t,f] = cohgramc(wcv,lf2_r,movingwin,params);
%     [Cw8,phi,S12,S1,S2,t,f] = cohgramc(wcv,lf8,movingwin,params);
%     [C83,phi,S12,S1,S2,t,f] = cohgramc(lf8,lf3,movingwin,params);
%     [C82,phi,S12,S1,S2,t,f] = cohgramc(lf8,lf2,movingwin,params);
%     [C82r,phi,S12,S1,S2,t,f] = cohgramc(lf8,lf2_r,movingwin,params);
    
    params.tapers = [2 3];
    [Sw,t2,f2]=mtspecgramc(wcv,movingwin,params);
%     [S8,t2,f2]=mtspecgramc(lf8,movingwin,params);
%     [S3,t2,f2]=mtspecgramc(lf3,movingwin,params);
%     [S2,t2,f2]=mtspecgramc(lf2,movingwin,params);
%     [S2r,t2,f2]=mtspecgramc(lf2_r,movingwin,params);
   
    
    
% generate figures
cax = [-5 0];
yl = [0 40];
l1 = [20 20];
l2 = [40 40];
% l3 = [60 60];

Fig = figure('visible','off');
set(Fig,'PaperUnits','centimeters');
set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(3,2,1)
    pcolor(t,f,log10(S1'));shading flat;
    ylim(yl)
    xl = xlim();
    line(xl,l1,'color','w')
    line(xl,l2,'color','w')
%     line(xl,l3,'color','w')
    set(gca,'yscale','log')
    caxis(cax)
    subplot(3,2,3)
    pcolor(t,f,log10(S2'));shading flat;
    ylim(yl)
    line(xl,l1,'color','w')
    line(xl,l2,'color','w')
%     line(xl,l3,'color','w')
    caxis(cax)
    set(gca,'yscale','log')
    subplot(3,2,5)
    pcolor(t,f,log10(S3'));shading flat;
    ylim(yl)
    caxis(cax)
    line(xl,l1,'color','w')
    line(xl,l2,'color','w')
%     line(xl,l3,'color','w')
    set(gca,'yscale','log')    
    subplot(3,2,2)
    pcolor(t,f,Cw3');shading flat;
    ylim(yl)
    caxis([0 1])
    line(xl,l1,'color','w')
    line(xl,l2,'color','w')
%     line(xl,l3,'color','w')
    set(gca,'yscale','log')
subplot(3,2,4)    
    pcolor(t,f,Cw8');shading flat;
    ylim(yl)
    caxis([0 1])
    set(gca,'yscale','log')
    line(xl,l1,'color','w')
    line(xl,l2,'color','w')
%     line(xl,l3,'color','w')
    subplot(3,2,6)
    pcolor(t,f,C83');shading flat;
    ylim(yl)
    caxis([0 1])
    line(xl,l1,'color','w')
    line(xl,l2,'color','w')
%     line(xl,l3,'color','w')
   set(gca,'yscale','log')
    tname = ['F:\WC_Germany\overall_EC\mt_theta_grams_4\' s_name];
    print('-dpng',tname);
    close
  
    %
    cax = [-3.5 0];
yl = [0 5];
Fig = figure('visible','off');
set(Fig,'PaperUnits','centimeters');
set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
Fig = figure('visible','off');
set(Fig,'PaperUnits','centimeters');
set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(3,2,1)
    pcolor(t,f,log10(S1'));shading flat;
    ylim(yl)
%     set(gca,'yscale','log')
    caxis(cax)
    subplot(3,2,3)
    pcolor(t,f,log10(S2'));shading flat;
    ylim(yl)
    caxis(cax)
%     set(gca,'yscale','log')
    subplot(3,2,5)
    pcolor(t,f,log10(S3'));shading flat;
    ylim(yl)
    caxis(cax)
%     set(gca,'yscale','log')    
    subplot(3,2,2)
    pcolor(t,f,Cw3');shading flat;
    ylim(yl)
    caxis([0 1])
%     set(gca,'yscale','log')
subplot(3,2,4)    
    pcolor(t,f,Cw8');shading flat;
    ylim(yl)
    caxis([0 1])
%     set(gca,'yscale','log')
    subplot(3,2,6)
    pcolor(t,f,C83');shading flat;
    ylim(yl)
    caxis([0 1])
%     set(gca,'yscale','log')
    tname = ['F:\WC_Germany\overall_EC\mt_theta_grams_3\' s_name];
    print('-dpng',tname);
    close

end

cd F:\WC_Germany\overall_EC\

