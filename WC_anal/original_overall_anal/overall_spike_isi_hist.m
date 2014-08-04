clear all

load('C:\WC_Germany\overall_calcs\HC_Cort\hc_cor_dir.mat')

Fs = 2016;
maxtime = 0.5;

params.fpass = [0 50];
params.err = [2 0.05];
params.Fs = Fs;
win = 100;

for d = 1:length(over_dir)
    
    cd(over_dir{d})
    pwd
    
    if exist('spike_time_br.mat')
        load spike_time_br
    else
        load spike_time
    end
    
    isis = diff(spkid)/Fs;
    
    if length(spkid) > 1000
        numBins = 500;
    else
        numBins = 100;
    end
    times = linspace(0,maxtime,numBins);
    isihist = hist(isis,times);
    
%     [S,f,R,varS,zerosp,C,Serr]=mtspectrumsegpt(spkid/Fs,win,params);
    
%     plot(f,S)
%     hold on
%     plot(f,Serr(1,:),'--')
%     plot(f,Serr(2,:),'--')
%     line([0 12],[R R],'Color','r')
%     xlim([0 12])
%         t_names = ['C:\WC_Germany\overall_calcs\HC_Cort\spike_isis\pow_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close

    
    bar(times,isihist)
    xlim([0 0.45])
    t_names = ['C:\WC_Germany\overall_calcs\HC_Cort\spike_isis\' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    
end