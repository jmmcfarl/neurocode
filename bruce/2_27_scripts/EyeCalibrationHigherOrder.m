% eye calibration with higher order terms
clear all

% addpath('./Data/ExptA/EyeCalibration');
cd ~/Data/bruce/2_27_12/EyeCalibration/ 
load EyeCal17Feb2012
load lemE003.em.mat
Trials = Expt.Trials;
eyecalavg = EyeCal.Data;
load EyeCalParas.mat

trialN = length(Trials);
SamplingFreq = 1./Expt.Header.CRrates(1);

EyeCal = [];

%% Separate out fixations
MinDuration = 0.2;
StarT = Expt.Trials(1).Start/10000;
for i=1:trialN
    Target = [Expt.Trials(i).fx Expt.Trials(i).fy];
    L = min([length(Expt.Trials(i).Eyevals.rh) length(Expt.Trials(i).Eyevals.rv) ...
        length(Expt.Trials(i).Eyevals.lh) length(Expt.Trials(i).Eyevals.lv)]);
    lEyeXY = [Expt.Trials(i).Eyevals.lh(1:L) Expt.Trials(i).Eyevals.lv(1:L)];
    rEyeXY = [Expt.Trials(i).Eyevals.rh(1:L) Expt.Trials(i).Eyevals.rv(1:L)];
    

    [fixations saccades status] = BinocularStatus(lEyeXY, rEyeXY, SamplingFreq, 20, 3);
    
    if isempty(fixations) continue; end
    
    Nfix = size(fixations,1);
    Ldist = zeros(Nfix,1); Rdist = zeros(Nfix,1);
    Lfixcenter = zeros(Nfix,2); Rfixcenter = zeros(Nfix,2);
    Duration = (fixations(:,2) - fixations(:,1))/SamplingFreq;
    LMaxShift = zeros(Nfix,2);
    for fix=1:Nfix
        indxs = fixations(fix,1):fixations(fix,2);
        Lfixcenter(fix,:) = median(lEyeXY(indxs,:));
        Rfixcenter(fix,:) = median(rEyeXY(indxs,:));
        LMaxShift(fix,:) = max(lEyeXY(indxs,:)) - min(lEyeXY(indxs,:));
        Ldist(fix) = norm(Lfixcenter(fix,:) - Target);
        Rdist(fix) = norm(Rfixcenter(fix,:) - Target);
        Tfix = Expt.Trials(i).Start/10000+fixations(fix,1)/SamplingFreq-StarT;
        
        if Duration(fix)>MinDuration &&  max(LMaxShift(fix,1:2)) <.5
            currtrial = [Lfixcenter(fix,:) Rfixcenter(fix,:) Target ...
                i fixations(fix,1) fixations(fix,2) Tfix];
            EyeCal = [EyeCal;  currtrial];
        end
    end
end

fprintf(' %d fixations from %d trials \n', size(EyeCal,1), trialN);
%% fit calibration parameters using all data pts
eyecal = EyeCal';

[EyeParaALL DataPts] = EyeCalibrationFit(eyecal, 'highergain');
% EyeParaALL
% to display error plots
% [EyeParaALL DataPts] = EyeCalibrationFit(eyecal, 'highergain', 1);
cd /Users/James/Data/bruce/2_27_12/saccades
save higher_order_calibration EyeParaALL
% %% calibrate eye signals
% EyeLdatascal = CalibrateEyeSignal(eyecal(1:2,:)', EyeParaALL.gain(1:2,:), EyeParaALL.offset(1:2), EyeParaALL.highergain(1:2,:));
% EyeRdatascal = CalibrateEyeSignal(eyecal(3:4,:)', EyeParaALL.gain(3:4,:), EyeParaALL.offset(3:4), EyeParaALL.highergain(3:4,:));