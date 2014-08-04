Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
addpath('~/James_scripts/bruce/eye_tracking_finalcode/');

for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    bar_ori = 0;
    do_monoc_eyetracking_RLS_highres_v2(Expt_name,bar_ori);
    
        bar_ori = 90;
    do_monoc_eyetracking_RLS_highres_v2(Expt_name,bar_ori);
end