% clear all

load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
dsf = 16;
params.Fs = 2016/dsf;
params.err = [2 .05];
params.fpass = [0 10];
params.tapers = [4 7];
% W = 0.04;



win = [100 10];

Fsd = 2016/16;
niqf = 2016/2;
hcf = 60/niqf;
[b,a] = butter(2,hcf,'low');

for d = 1:length(dir_array)
    
    disp(sprintf('session %d',d))
    cd(dir_array{d});
    pwd
    
    
    load used_data wcv_minus_spike lf8
      
%     T = floor(length(lf8)/2016);
%     TW = floor(T*W);
%     
% %     params.tapers = [TW 2*TW-1-3];
%     
%     
%     disp(sprintf('T = %0.5g K = %0.5g',T,(2*TW-1)))
  
    
    
    
    wcv_minus_spike = filtfilt(b,a,wcv_minus_spike);
    lf8 = filtfilt(b,a,lf8);   

    down_w = downsample(wcv_minus_spike,dsf);
    down_8 = downsample(lf8,dsf);
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
   
    
    [C{d},phi{d},S12,S1,S2,f{d},confC{d},phistd{d},Cerr{d}]=coherencysegc(down_w,down_8,win,params)
    
    goodF{d} = find(Cerr{d}(1,:)>confC{d});
%     phi{d}(badF{d}) = nan;
    
%     [S8{d},f{d},S8err{d}]=mtspectrumc(down_8,params);
%     [Sw{d},f{d},Swerr{d}] = mtspectrumc(down_w,params);
    
   
% sm_win = round(W/mean(diff(f{d}))/4)
    
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [40 40]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     plot(f{d},jmm_smooth_1d(10*log10(Sw{d}),sm_win),'linewidth',2)
%     hold on
%     plot(f{d},jmm_smooth_1d(10*log10(Swerr{d}(1,:)),sm_win),'--','linewidth',2)
%     plot(f{d},jmm_smooth_1d(10*log10(Swerr{d}(2,:)),sm_win),'--','linewidth',2)
%     plot(f{d},jmm_smooth_1d(10*log10(S8{d}),sm_win),'r','linewidth',2)
%     plot(f{d},jmm_smooth_1d(10*log10(S8err{d}(1,:)),sm_win),'r--','linewidth',2)
%     plot(f{d},jmm_smooth_1d(10*log10(S8err{d}(2,:)),sm_win),'r--','linewidth',2)
%     xlim([0 1])
%     tname = ['C:\WC_Germany\JMM_Analysis_pyr\pow_spectra\' f_names{d}];
%     print('-dpng',tname);
%     close

    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [40 40]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    plot_vector(C{d},f{d},'n',Cerr{d},'b',2)
    xlim([0 1])
    line([0 1],[confC{d} confC{d}],'Color','r','linewidth',2)
    subplot(2,1,2)
    plot(f{d},phi{d},'k','linewidth',2)
    hold on
    plot(f{d}(goodF{d}),phi{d}(goodF{d}),'r.','linewidth',2,'MarkerSize',16)
    xlim([0 1])
    tname = ['C:\WC_Germany\JMM_Analysis_pyr\coherency\' f_names{d}];
    print('-dpng',tname);
    close



    save C:\WC_Germany\JMM_Analysis_pyr\coherency\coh_data C phi f confC Cerr phistd
    
  
    clear down_w down_8 wcv* lf8 
    
end