clear all
load C:\WC_Germany\JMM_analysis_pyr\heka_anal\pers_start_times
load C:\WC_Germany\JMM_analysis_pyr\heka_anal\sweep_times
load pyr_heka_dir
pers_start_times(11) = [];
pers_stop_times(11) = [];
npers_start_times(11) = [];
npers_stop_times(11) = [];

ot1 = 2;
ot2 = 2+3.2432;
sd1 = 33783/5e3;
sd2 = 49999/5e3;

sps = [49999 49999 49999 49999 49999 33783 49999 33783 33783 33783 33783 33783 33783 33783 33783 33783];
off_time = [ot1 ot1 ot1 ot1 ot1 ot2 ot1 ot2 ot2 ot2 ot2 ot2 ot2 ot2 ot2 ot2];
no_pa = [];
for d = 1:length(f_loc)

    load(f_loc{d})
    dat_name = [f_loc{d} '_MP'];
    dat_name(1:24) = [];
    eval(['detrend(' dat_name ');'])
    eval(['sig = ' dat_name ';']);
   eval(['h_length = length(' dat_name ');'])
   sig = zscore(sig);
   h_secs = h_length/5e3;
   num_sweeps = floor(h_secs/10);

   htaxis = zeros(1,h_length);
   
    for s = 1:num_sweeps
        
        htaxis((s-1)*sps(d)+1:s*sps(d)) = sweep_times{d}(s,1):1/5e3:sweep_times{d}(s,2)-1/5e3;
        
    end
    
    htaxis(s*sps(d)+1:end-1) = sweep_times{d}(s+1,1):1/5e3:sweep_times{d}(s+1,2);
    
   pers_sig = [];
   npers_sig = [];
   num_pers_states = length(pers_start_times{d});
   num_npers_states = length(npers_start_times{d});
   if num_pers_states > 0
      
       for p = 1:num_pers_states
          
           cur_points = find(htaxis > pers_start_times{d}(p) & htaxis < pers_stop_times{d}(p));
           pers_sig = [pers_sig;sig(cur_points)];
           
       end
       
       for p = 1:num_npers_states
          
           cur_points = find(htaxis > npers_start_times{d}(p) & htaxis < npers_stop_times{d}(p));
           npers_sig = [npers_sig;sig(cur_points)];
           
       end
      
%       [pers_dist(d,:),pers_grid] = gpkde(pers_sig,-3,[-0.8 0 600]);
%       [npers_dist(d,:),npers_grid] = gpkde(npers_sig,-3,[-0.8 0 600]);
%        
      [pers_dist(d,:),pers_grid] = gpkde(pers_sig,-3,[-4 4 600]);
      [npers_dist(d,:),npers_grid] = gpkde(npers_sig,-3,[-4 4 600]);

%       plot(pers_grid,pers_dist)
%       hold on
%       plot(npers_grid,npers_dist,'r')
%       xlabel('Voltage','FontSize',14)
%       ylabel('Prob Density','FontSize',14)
%       legend('Persistent States','Non-persistent States')
%       t_names = ['C:\WC_Germany\JMM_analysis_pyr\heka_anal\state_dep_dist_' f_names{d}];
%       print('-dpng',t_names)
%       close all
      
   else
       disp('no pers activity')
       no_pa = [no_pa d];
   end

   eval(['clear ' dat_name ])
   clear sig pers_sig npers_sig
   
end

pers_dist(no_pa,:) = [];
npers_dist(no_pa,:) = [];

m_p_d = mean(pers_dist);
u_p_d = m_p_d + 2*std(pers_dist)/sqrt(14);
l_p_d = m_p_d - 2*std(pers_dist)/sqrt(14);

m_n_d = mean(npers_dist);
u_n_d = m_n_d + 2*std(npers_dist)/sqrt(14);
l_n_d = m_n_d - 2*std(npers_dist)/sqrt(14);
