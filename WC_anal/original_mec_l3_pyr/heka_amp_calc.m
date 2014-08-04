%% for pyramidal cells
clear all
load pyr_heka_dir

for d = 1:length(f_loc)
    
    load(f_loc{d})
    dat_name = [f_loc{d} '_MP'];
    dat_name(1:24) = [];
    eval(['detrend(' dat_name ');'])
   eval(['[y,x] = gpkde(' dat_name ',-3,[-0.8 0.1 600]);'])
    
   [a,b] = findpeaks(y);
   
   plot(x,y,'linewidth',2)
   xlabel('Amplitude (Vx10-4)','FontSize',14)
   ylabel('Density','FontSize',14)
   
   if length(a) > 2
       
       %find highest two peaks
       [dummy,peakorder] = sort(a,'descend');
       twopeaks = x(b(peakorder(1:2)));
       title(sprintf('first peak %0.5g second peak %0.5g',twopeaks(1),twopeaks(2)))
       
       [dummy,uploc] = max(twopeaks);
       [dummy,downloc] = min(twopeaks);
       upstate_p(d) = x(b(uploc));
       downstate_p(d) = x(b(downloc));
       
   else
       disp('unimodal')
       title(num2str(x(b)))
       downstate_p(d) = x(b);
       upstate_p(d) = nan;
   end
   
   t_name = ['C:\WC_Germany\JMM_analysis_pyr\amplitude\' f_names{d}];
   print('-dpng',t_name)
   close all
   
   dist_est_p(d,:) = y;
   dist_grid_p(d,:) = x;
   
end

save C:\WC_Germany\JMM_analysis_pyr\heka_data_pyr dist*




