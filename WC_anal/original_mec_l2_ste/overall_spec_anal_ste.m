clear all

load C:\WC_Germany\JMM_Analysis_ste\dir_tree_ste
dsf = 8;
params.Fs = 2016/dsf;
% params.err = [2 0.05];
params.fpass = [0 20];
params.tapers = [4 7];
sm_win = 2;


win = [200];

niqf = 2016/2;
hcf = 100/niqf;
[b,a] = butter(2,hcf,'low');

lf = 5202;
Sw = zeros(length(dir_array),lf);
S8 = Sw;
Swerr = zeros(length(dir_array),2,lf);
S8err = Swerr;

for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike

    wcv_minus_spike = filtfilt(b,a,wcv_minus_spike);
    lf8 = filtfilt(b,a,lf8);   

    down_w = downsample(wcv_minus_spike,dsf);
    down_8 = downsample(lf8,dsf);
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
   
%     data_dur = log2(length(down_8));
%     if data_dur<17
%         params.pad = 1;
%     else
%         params.pad = 0;
%     end
    
    params.err = [2 .05];
    [C(d,:),phi(d,:),S12,S1,S2,f,confC(d)]=coherencysegc(down_w,down_8,win,params);
%     params.err = 0;
    [S8(d,:),f,varS,dummy,S8err(d,:,:)]=mtspectrumsegc(down_8,win,params);
    [Sw(d,:),f,varS,dummy,Swerr(d,:,:)] = mtspectrumsegc(down_w,win,params);
    
    S8(d,:) = jmm_smooth_1d(S8(d,:),sm_win);
    Sw(d,:) = jmm_smooth_1d(Sw(d,:),sm_win);
    S8err(d,1,:) = jmm_smooth_1d(squeeze(S8err(d,1,:)),sm_win);
    S8err(d,2,:) = jmm_smooth_1d(squeeze(S8err(d,2,:)),sm_win);
    Swerr(d,1,:) = jmm_smooth_1d(squeeze(Swerr(d,1,:)),sm_win);
    Swerr(d,2,:) = jmm_smooth_1d(squeeze(Swerr(d,2,:)),sm_win);
    
        
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [40 40]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    plot_vector(S8(d,:),f,'l',squeeze(S8err(d,:,:)),'r')
    hold on
    plot_vector(Sw(d,:),f,'l',squeeze(Swerr(d,:,:)))
    xlim([0 1])
    subplot(2,1,2)
    plot(f,C(d,:),'linewidth',2)
    line([0 2],[confC(d) confC(d)],'Color','k')
    xlim([0 1])
    tnames = ['C:\WC_Germany\JMM_Analysis_ste\spectra\' f_names{d}];
    print(tnames,'-dpng')
    close all
    
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [40 40]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    plot_vector(S8(d,:),f,'l',squeeze(S8err(d,:,:)),'r')
    hold on
    plot_vector(Sw(d,:),f,'l',squeeze(Swerr(d,:,:)))
    xlim([0 10])
    subplot(2,1,2)
    plot(f,C(d,:),'linewidth',2)
    line([0 10],[confC(d) confC(d)],'Color','k')
    xlim([0 10])
    tnames = ['C:\WC_Germany\JMM_Analysis_ste\spectra\wb_' f_names{d}];
    print(tnames,'-dpng')
    close all
    
end

save C:\WC_Germany\JMM_Analysis_ste\spectra C phi confC S8 Sw f

% %% overall figures
% 
Swm = mean(10*log10(Sw));
Swe = std(10*log10(Sw))/sqrt(17);
S8m = mean(10*log10(S8));
S8e = std(10*log10(S8))/sqrt(18);

Swl = Swm-2*Swe;
Swu = Swm+2*Swe;
S8l = S8m-2*S8e;
S8u = S8m+2*S8e;
% 
plot(f,Swm,'linewidth',2)
hold on
plot(f,Swu,'--')
plot(f,Swl,'--')
plot(f,S8m,'r','linewidth',2)
plot(f,S8u,'r--')
plot(f,S8l,'r--')
xlim([0 1])
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (dB)','FontSize',14)
% 
% 
Cm = mean(C);
Ce = std(C)/sqrt(17);
Cu = Cm+2*Ce;
Cl = Cm-2*Ce;
% 
plot(f,Cm,'linewidth',2)
hold on
plot(f,Cu,'--')
plot(f,Cl,'--')
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Coherency','FontSize',14)
xlim([0 2])
av_conf = mean(confC);
min_conf = min(confC);
max_conf = max(confC);
line([0 1],[av_conf av_conf],'Color','k')
line([0 1],[min_conf min_conf],'Color','k','linestyle','--')
line([0 1],[max_conf max_conf],'Color','k','linestyle','--')
ylim([0 1])