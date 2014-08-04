%% for pyramidal cells
clear all
load C:\WC_Germany\Persistent_activity\pyr_heka_dir

for d = 1:17

load(f_loc{d})
    dat_name = [f_loc{d} '_MP'];
    dat_name(1:24) = [];
    eval(['detrend(' dat_name ');'])
   eval(['[y{d},x{d}] = gpkde(100*' dat_name ',-3,[-0.8*100 0.1*100 600]);'])

plot(x{d},y{d})

    tname = ['C:\WC_Germany\Persistent_activity\heka_amp\' f_names{d}];
    print('-dpng',tname);
    close


end

save C:\WC_Germany\Persistent_activity\heka_amp\heka_amp_data x y

