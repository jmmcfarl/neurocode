clear all
close all
drive_letter = 'C';

addpath(strcat(drive_letter,':\Code\smoothing\software'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))
addpath(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))
addpath(strcat(drive_letter,':\WC_Germany\persistent_9_27_2010\'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_uds_code'))
addpath(strcat(drive_letter,':\Code\maphmmbox\'))
addpath('C:\WC_Germany\hsmm_uds_code\')

cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir.mat

raw_Fs = 2016;

%%
for d = 1:length(combined_dir)
% for d = [56]
    cd(combined_dir{d})
    pwd
    
    last_slash = find(combined_dir{d} == '\',1,'last');
    f_name = combined_dir{d}(last_slash+1:end);
    
    load ./used_data wcv_minus_spike lf8 lf7
%     if ctx_lfp(d) == 7
%         lf8 = lf7;
%     end
    [desynch_times_lf7,desynch_inds,P_lf7,f,t] = locate_desynch_times_individual(lf7);
    
    desynch_ind = zeros(size(t));
    for i = 1:size(desynch_times_lf7,1)
        desynch_ind(find(t > desynch_times_lf7(i,1) & t < desynch_times_lf7(i,2))) = 1;
    end
    
    figure
    subplot(3,1,1)
    pcolor(t,f,P_lf7'); shading flat
    ylim([0 6]); set(gca,'yscale','log')
    hold on
    plot(t,desynch_ind+1,'w','linewidth',2)
    
     [desynch_times_lf7,desynch_inds,P_lf7,f,t] = locate_desynch_times_individual_v2(lf7);
%       [max_so,net_hf,t] = locate_desynch_times_individual_v2(lf7);
   
f_so = find(f > 0.2 & f < 1.5); %define the slow-oscillation band as 0.2-1.5Hz
f_hf = find(f > 4); %define the 'high-frequency' band as > 4Hz

%find max log power in theta and SO freq bands
df = f(2)-f(1);
max_so = max(10*log10(P_lf7(:,f_so)),[],2); %maximum log relative power in the SO range
net_hf = (trapz(log10(P_lf7(:,f_hf)),2)*df); %z-score of the integral HF power

    desynch_ind = zeros(size(t));
    for i = 1:size(desynch_times_lf7,1)
        desynch_ind(find(t > desynch_times_lf7(i,1) & t < desynch_times_lf7(i,2))) = 1;
    end

    plot(t,desynch_ind+0.75,'r','linewidth',2)
    subplot(3,1,2)
    plot(t,max_so,'r'), axis tight
    subplot(3,1,3)
    plot(t,net_hf), axis tight
    
    pname = ['C:\WC_Germany\sven_thomas_combined\desynch_epochs\v2_' f_name];
    print(pname,'-dpng')
    close

end