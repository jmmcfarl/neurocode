% The script performs three tasks:
% (1) calibrate eye signal used saved gain/offset
% (2) identify saccade
% (3) label data as 'off-screen' 0 'fixating' 1 and 'saccading '2'
% Yuwei Cui, Last Modified Mar 7, 2012
clear all

datafolder = './Data/ExptB/';
savetofolder = './Data/ExptBproc/';

addpath(datafolder); 
addpath(savetofolder);
load EyeCalParas.mat

filenames = SearchFileNames(datafolder, '.em.mat')
% eye data
for id = 1:length(filenames)
    filename = filenames{id};
    load(filename)
    SamplingFreq = 1./Expt.Header.CRrates(1);
    Expt.Header
%     continue
    
    TrialN = length(Expt.Trials);
    for triali=1:TrialN
        Start = Expt.Trials.Start;
        End = Expt.Trials.End;
        
        Tlength = (End-Start)/10000;
        lEyeXY0 = [Expt.Trials(triali).Eyevals.lh Expt.Trials(triali).Eyevals.lv];
        rEyeXY0 = [Expt.Trials(triali).Eyevals.rh Expt.Trials(triali).Eyevals.rv];
        
%         EyePara.gain
        
        % use calibration parameters
        lEyeXY = CalibrateEyeSignal(lEyeXY0, EyePara.gain(1:2,:),EyePara.offset(1:2));
        rEyeXY = CalibrateEyeSignal(rEyeXY0, EyePara.gain(3:4,:),EyePara.offset(3:4));       
        
        EM = EyeM(lEyeXY,SamplingFreq);
        EM = EM.AssignStatus([-7 7 -7 7], 20);
        
        % process drifts
        lfixs = EM.fixations;
        FixT = [];sumdist=[];shiftx=[];shifty=[];
        Tstart = []; Tend=[];
        for fixi=1:length(lfixs)
            eyetrace = lfixs(fixi).fixtrace;
            Tstart(fixi) = lfixs(fixi).Tstart;
            Tend(fixi) = lfixs(fixi).Tend;
            Duration = size(eyetrace,1)/SamplingFreq;
  
            FixT(fixi) = Duration;
%             eyetrace = eyetrace(10:end-10,:);
            diff = eyetrace(2:end,:)- eyetrace(1:end-1,:);
            absdiff = abs(diff);
            sumdist(fixi) = sum(sqrt(absdiff(:,1).^2+absdiff(:,2).^2));
            shiftx(fixi)  = eyetrace(end,1) - eyetrace(1,1);
            shifty(fixi)  = eyetrace(end,2) - eyetrace(1,2);            
        end
        EyeMdrift = struct('Duration',FixT,'StartT',Tstart,'EndT',Tend,'SumDist',sumdist,'Shift',[shiftx shifty]);
        % process saccade
        sacN = length(EM.saccades);
        
        saccades = EM.saccades;
        sacT = [];
        direction = [];
        duration = [];
        amplitude = [];
        peakV = [];
        for i=1:sacN
            if strcmp( EM.saccades(i).label, 'good')
                sacT(end+1) = saccades(i).Tstart;
                direction(end+1) = saccades(i).direction;
                duration(end+1) = saccades(i).Tend-saccades(i).Tstart;
                amplitude(end+1) = saccades(i).Amp;
                peakV(end+1) = saccades(i).PeakVel;
            end
        end
        % peakV = peakV(peakV<1000);
        status = EM.status;
        sacT = Start/10000+sacT/SamplingFreq;
        EyeMsac = struct('status',status,'Start',Start, 'sacT',sacT,'direction',direction,'duration',duration,...
            'amplitude',amplitude,'peakV', peakV);
        
        % save('saccades','saccades');
        % save('saccadeTstart','sacT');
        % save('status','status');
        avgsacfreq = length(sacT)/(sum(EM.status~=0)*EM.tres);
        disp(['Exp Length (sec) : ', num2str(Tlength), ' Avg Saccade Frequency: ',num2str(avgsacfreq)]);
        Expt.Trials(triali).EyeM = EM;
        Expt.Trials(triali).EyeMsac=EyeMsac;
        Expt.Trials(triali).EyeMdrift=EyeMdrift;
        Expt.Trials(triali).Eyevals.rh = rEyeXY(:,1);
        Expt.Trials(triali).Eyevals.rv = rEyeXY(:,2);
        Expt.Trials(triali).Eyevals.lh = lEyeXY(:,1);
        Expt.Trials(triali).Eyevals.lv = lEyeXY(:,2);
    end
    save([savetofolder filename(1:end-3) 'sac.mat'],'Expt');
end

% EM.DisplayEyeM(15,20)
