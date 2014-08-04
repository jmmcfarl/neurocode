%% for pyramidal cells
clear all
load C:\WC_Germany\persistent_revised\pers_revised_heka_dir

for d = 1:length(f_loc)
    d
load(f_loc{d})
if d <= 17
    dat_name = [f_loc{d} '_MP'];
        dat_name(1:24) = [];
        eval([dat_name '=' dat_name '*100;']);
elseif d >= 21 & d < 23
    dat_name = 'dataVec';
        eval([dat_name '= transpose(' dat_name ')*100;']);
elseif d == 23
    dat_name = 'dataVec';
        eval([dat_name '=' dat_name '*100;']);
else
    dat_name = 'data';
end
    eval(['new_data = ' dat_name ';'])
  [y{d},x{d}] = gpkde(new_data,1,[-0.9*100 0.1*100 600]);

plot(x{d},y{d})

    tname = ['C:\WC_Germany\persistent_revised\heka_amp\' f_names{d}];
    print('-dpng',tname);
    close


end

save C:\WC_Germany\persistent_revised\heka_amp\heka_amp_data x y

for d = 1:length(f_loc)
    
   [pks,pk_locs] = findpeaks(y{d});
   plot(x{d},y{d})
   x{d}(pk_locs)
   d
   pause
   close all
    
end
    

%%
avg_state_mp = ...
    [-66.5 -46.3;
    -65.0 -41.1;
    -69.3 -44.1;
    -74.1 -49.8;
    -66.3 -43.3;
    -69.3 -45.9;
    -66.0 -41.8;
    -63.5 -44.8;
    -68.5 -44.1;
    -72.8 -48.8;
    -61.8 -47.6;
    -69.1 -48.6;
    -63.1 -44.8;
    -56.4 -37.7;
    -61.6 -44.9;
    -63.3 -43.3;
    -65.6 -41.8;
    -74.1 -47.3;
    -68.0 -46.1;
    -70.8 -49.8;
    -76.0 nan;
    -68.8 -61.5;
    -66.3 -59.1;
%     -70.8 -49.8;
    -73.8 -51.9;
    -74.0 -51.8;
    -76.1 -56.8];
    