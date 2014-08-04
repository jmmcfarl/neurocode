clear all

cd /Users/James/Data/bruce/2_27_12/saccades
load higher_order_calibration


%% calibrate eye signals
for i = 1:5
    fprintf('Block %d of %d\n',i,5);
    cd /Users/James/Data/bruce/2_27_12/
    
    load(sprintf('lemM232.5%d.em.mat',i));
    n_trials = length(Expt.Trials);
    for t = 1:n_trials
        cur_lh = Expt.Trials(t).Eyevals.lh;
        cur_lv = Expt.Trials(t).Eyevals.lv;
        cur_rh = Expt.Trials(t).Eyevals.rh;
        cur_rv = Expt.Trials(t).Eyevals.rv;
        disph = cur_rh - cur_lh;
        dispv = cur_rv - cur_lv;
        EyeLdatascal = CalibrateEyeSignal([cur_lh cur_lv], EyeParaALL.gain(1:2,:), EyeParaALL.offset(1:2), EyeParaALL.highergain(1:2,:));
        EyeRdatascal = CalibrateEyeSignal([cur_rh cur_rv], EyeParaALL.gain(3:4,:), EyeParaALL.offset(3:4), EyeParaALL.highergain(3:4,:));
        disph_cal = EyeRdatascal(:,1)-EyeLdatascal(:,1);
        dispv_cal = EyeRdatascal(:,2)-EyeLdatascal(:,2);
        
        Expt.Trials(t).Eyevals.lh = EyeLdatascal(:,1);
        Expt.Trials(t).Eyevals.lv = EyeLdatascal(:,2);
        Expt.Trials(t).Eyevals.rh = EyeRdatascal(:,1);
        Expt.Trials(t).Eyevals.rv = EyeRdatascal(:,2);
    end
    cd /Users/James/Data/bruce/2_27_12/saccades/
    fname = sprintf('lemM232.5%d.em.hor.mat',i);
    save(fname,'Expt');
    clear Expt
end