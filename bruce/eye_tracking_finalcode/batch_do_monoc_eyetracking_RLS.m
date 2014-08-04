
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
use_LOOXV = 1;
use_coils = [0 0];
n_fix_inf_it = 0;
n_drift_inf_it = 0;
recompute_init_mods = 1;

for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    
    bar_ori = 0;
    use_measured_pos = 1;
    do_monoc_eyetracking_RLS(Expt_name,bar_ori,use_LOOXV,use_coils,use_measured_pos,recompute_init_mods,n_fix_inf_it,n_drift_inf_it)
    use_measured_pos = 2;
    do_monoc_eyetracking_RLS(Expt_name,bar_ori,use_LOOXV,use_coils,use_measured_pos,recompute_init_mods,n_fix_inf_it,n_drift_inf_it)

    bar_ori = 90;
    use_measured_pos = 1;
    do_monoc_eyetracking_RLS(Expt_name,bar_ori,use_LOOXV,use_coils,use_measured_pos,recompute_init_mods,n_fix_inf_it,n_drift_inf_it)
    use_measured_pos = 2;
    do_monoc_eyetracking_RLS(Expt_name,bar_ori,use_LOOXV,use_coils,use_measured_pos,recompute_init_mods,n_fix_inf_it,n_drift_inf_it)
end