
%% for pyramidal cells
clear all
load C:\WC_Germany\JMM_analysis_ste\stellate_heka_dir

windowSize = 30; %in s
windowSlide = 5; %in s


for d = 1:length(f_loc)
    d
load(f_loc{d})
    dat_name = [f_loc{d} '_MP'];
        dat_name(1:24) = [];
        eval([dat_name '=' dat_name '*100;']);
        Fs = 5e3;
        dsf = 10;
    Fsd = Fs/dsf;
    eval(['new_data = ' dat_name ';'])
    new_data_down = downsample(new_data,dsf);
    
    total_dur = length(new_data_down)/Fsd;
    numWins = floor((total_dur-windowSize)/windowSlide);
    t_axis = (1:numWins)*windowSlide;

for w = 1:numWins
    
           begT = (w-1)*windowSlide;
        endT = begT + windowSize;
            begInd = begT*Fsd+1;
        endInd = min((begInd+windowSize*Fsd),length(new_data_down));
        data_seg = new_data_down(begInd:endInd);

  [data_dist(w,:),dat_range] = gpkde(data_seg,0.5,[-0.9*100 0.1*100 600]);
 
end

figure
    pcolor(t_axis,dat_range,data_dist');shading flat
    ylim([-90 -20])
    tname = ['C:\WC_Germany\persistent_revised\heka_amp\sliding_dist_' f_names{d}];
    print('-dpng',tname);
    close
clear dat* t_axis
end

% save C:\WC_Germany\persistent_revised\heka_amp\heka_amp_data x y

    