classdef EyeM
    properties
        tres;       % temporal resolution (sec,1/samplingrate)
        eyeXY;      % eye location
        eyeV;       % eye velocity
        status;     % 0 off screen; 1: fixating; 2: during saccade
        saccades;
        fixations;
        screen;     % screen size [-X X -Y Y]
    end
    
    methods
        function EM = EyeM(eyeXY, samplingrate)
            EM.eyeXY = eyeXY;
            EM.tres = 1/samplingrate;
            eyet = (1:length(EM.eyeXY))*EM.tres;
            EM.eyeV = EyeM.GetEyeVelocity(eyeXY, eyet);
            
        end
        
        function EM = AssignStatus(EM, screen, vthreshsac, vthreshfix)
            if nargin<4
                vthreshfix = vthreshsac;
            end
            EM.screen = screen;
            NT = length(EM.eyeXY);
            EM.status=ones(NT,1);
            OffScreen = EM.eyeXY(:,1)<screen(1) | EM.eyeXY(:,1)>screen(2)...
                | EM.eyeXY(:,2)<screen(3) | EM.eyeXY(:,2)>screen(4);
            
            fprintf('Offscreen Ratio: %2.6f\n',sum(OffScreen/length(OffScreen)));
            fprintf('On Screen time (sec): %4.2f\n',sum(OffScreen==0)*EM.tres);
            EM.status(OffScreen==1) = 0; % 0: off screen
            EM.saccades = [];            
            
            speed = sqrt(EM.eyeV(:,1).^2+EM.eyeV(:,2).^2);
            
            % isolate saccade and fixation            
            startsac = 0; % 0: no sac 1: during sac
            startfix = 0; % 1: during fixation
            endsact = -1; startsact = -1;
            endfixt = -1; startfixt = -1;
            for t=1:NT
%                 if OffScreen(t)
%                     continue;
%                 end
                switch startfix
                    case 0 % not during fixation
                        if speed(t)<vthreshfix && OffScreen(t)==0
                            startfix = 1; startfixt = t;
                        end
                    case 1 % during fixation
                        if OffScreen(t) % end of bad fixation (offscreen)
                            startfix = 0; endfixt = t; 
                            fixlabel = 'offscreen';
                        elseif speed(t)>=vthreshfix % end of good fixation
                            startfix = 0; endfixt = t;
                            fixlabel = 'good';
                        end                
                end
                switch startsac
                    case 0 % not during saccade
                        if speed(t)>=vthreshsac && OffScreen(t)==0
                            startsac = 1; startsact = t;                                                    
                        end
                    case 1 % during sac
                        if OffScreen(t) % saccade goes off the screen
                            startsac = 0; endsact = t;                                 
                            saclabel = 'offscreen';
                        elseif speed(t)<vthreshsac % end of good sac
                            startsac = 0; endsact = t;
                            saclabel = 'good';
                        end
                end
                if endfixt==t % analyze fixation
                    if (endfixt-startfixt)*EM.tres>0.020                        
                        fixtrace = EM.eyeXY(startfixt:endfixt,:);
                        fixcenter = median(fixtrace);
                        fixtrace(:,1) = fixtrace(:,1)-fixcenter(1);
                        fixtrace(:,2) = fixtrace(:,2)-fixcenter(2);
                        lastfix = struct('label',fixlabel,'fixcenter',fixcenter,'fixtrace',fixtrace,...
                            'Tstart',startfixt,'Tend',endfixt);
                        if ~isempty(EM.fixations)
                            EM.fixations(end+1) = lastfix;
                        else
                            EM.fixations = lastfix;
                        end
                    end
                end
                        
                if endsact==t % analyze saccade
                    sactrace = EM.eyeXY(startsact:endsact,:);
                    startloc = sactrace(1,:);
                    endloc = sactrace(end,:);
                    amp = norm(endloc-startloc);
                    direction = atan2(endloc(2)-startloc(2),endloc(1)-startloc(1));
                    peakV = max(speed(startsact:endsact));
                    lastsac = struct('label',saclabel,'Tstart',startsact,'Tend',endsact,...
                        'PeakVel',peakV,...
                        'StartLoc',startloc,'EndLoc',endloc,'Amp',amp,'direction',direction);
                    if peakV<vthreshsac
                        keyboard;
                    elseif peakV>1000
                        continue;
                    end
                    EM.status(startsact:endsact) = 2;
                    if ~isempty(EM.saccades)
                        EM.saccades(end+1) = lastsac;
%                         fprintf('t = %5.4f PeakVel %d\n',endsact*EM.tres, lastsac.PeakVel);
                    else
                        EM.saccades = lastsac;
                    end
                end
            end            
        end
        
        function h=DisplayEyeM(EyeM, tstart, tend)
            h=figure(1);clf; rn=3;cn=1;
            if nargin>1
                ti = round(tstart/EyeM.tres):min([round(tend/EyeM.tres) length(EyeM.eyeXY)]);
            else
                ti = 1:length(EyeM.eyeXY);
            end
            taxis = ti*EyeM.tres;
            ShowSac = 0;
            % label microsaccade with triangle
            
            if ShowSac saci = find(EyeM.sact>ti(1) & EyeM.sact<ti(end));end
            
            subplot(rn,cn,1);
            plot(taxis,EyeM.eyeXY(ti,1),'b');hold on;
            plot(taxis,EyeM.eyeXY(ti,2),'r');hold on;
            Radi = sqrt(EyeM.eyeXY(:,1).^2+EyeM.eyeXY(:,2).^2);
%             plot(taxis,Radi(ti),'k--');hold on;
            if ShowSac
                plot(EyeM.sact(saci)*EyeM.tres, Radi(EyeM.sact(saci))*1.2,'vm');hold on;
            end
            set(gca,'XLim',[taxis(1) taxis(end)]); axis tight
            %             set(gca,'YLim',[-3 3]);
            
            xlabel('t (sec)'); ylabel('Deg');
            
            subplot(rn,cn,2);
            plot(taxis,EyeM.eyeV(ti,1),'b');hold on;
            plot(taxis,EyeM.eyeV(ti,2),'r');hold on;
            EyeSpds = sqrt(EyeM.eyeV(ti,1).^2+EyeM.eyeV(ti,2).^2);
            plot(taxis,EyeSpds,'k--');hold on;
            if ShowSac
                plot(EyeM.sact(saci)*EyeM.tres, EyeM.sacspds(saci)*1.2,'vm');hold on;
            end
            set(gca,'XLim',[taxis(1) taxis(end)]);
            xlabel('t (sec)'); ylabel('Deg/s'); axis tight
%                         set(gca,'YLim',[-40 40]);
            
            subplot(rn,cn,3);
            plot(taxis,EyeM.status(ti));axis tight;
            set(gca,'YLim',[-.1 2.1]);                   
        end
        
        function EM = SmoothEyeData(EM, cutoffFreq)            
            EM.eyeXY(:,1) = EyeM.CutOffHighfreq(EM.eyeXY(:,1), EM.tres, cutoffFreq);
            EM.eyeXY(:,2) = EyeM.CutOffHighfreq(EM.eyeXY(:,1), EM.tres, cutoffFreq);
            eyet = (1:length(EM.eyeXY))*EM.tres;
            EM.eyeV = EyeM.GetEyeVelocity(EM.eyeXY, eyet);            
%             EM.eyeV(:,1) = EyeM.CutOffHighfreq(EM.eyeV(:,1), EM.tres, cutoffFreq);
%             EM.eyeV(:,2) = EyeM.CutOffHighfreq(EM.eyeV(:,2), EM.tres, cutoffFreq);
        end
    end
    
    
    methods (Static)
        function eyev = GetEyeVelocity(eyeXY, eyet)
            % function: Calculate Eye velocity
            % Input:
            %       eyeXY(NT,2) Eye trajectory
            %       eyet(NT,1)  sampling times
            
            NT = length(eyeXY);
            eyex = eyeXY(:,1);
            eyey = eyeXY(:,2);
            
            % estimate sampling resolution (ms)
            eyedt = (eyet(end)-eyet(1))/(length(eyet)-1);
            neweyet = (0.5:NT-0.5)*eyedt;
            
            % resample eye trajectory at the average sampling freq.
            eyex = interp1(eyet,eyex,neweyet);
            eyex(isnan(eyex))=0;
            eyey = interp1(eyet,eyey,neweyet);
            eyey(isnan(eyey))=0;
            eyeXY(:,1) = eyex;
            eyeXY(:,2) = eyey;
            
            %             keyboard
            
            % calculate eye velocity
            eyev = EyeM.vecvel(eyeXY,1/eyedt,2);
            
            % eyex0=eyex;
            % eyey0=eyey;
            % smoothing eye-trajectory, cutoff high-freq noise
            %             eyex = CutOffHighfreq(eyex, eyedt, 60);
            %             eyey = CutOffHighfreq(eyey, eyedt, 60);
            
            %             calcualte eye velocity (1st-order derivative)
            %             EyeM.eyev = zeros(NT,2);
            %             EyeM.eyev(1:end-1,1) = diff(eyex)/eyedt;
            %             EyeM.eyev(1:end-1,2) = diff(eyey)/eyedt;
            
            
            % smoothing eye-trajectory, cutoff high-freq noise
            %             EyeM.eyev(:,1) = CutOffHighfreq(EyeM.eyev(:,1), eyedt, 60);
            %             EyeM.eyev(:,2) = CutOffHighfreq(EyeM.eyev(:,2), eyedt, 60);
            
        end
        
        
        function v = vecvel(xx,SAMPLING,TYPE)
            %------------------------------------------------------------
            %
            %  FUNCTION vecvel.m
            %  Calculation of eye velocity from position data
            %  Please cite: Engbert, R., & Kliegl, R. (2003) Microsaccades
            %  uncover the orientation of covert attention. Vision Research
            %  43: 1035-1045.
            %
            %  (Version 1.2, 01 JUL 05)
            %-------------------------------------------------------------
            %
            %  INPUT:
            %
            %  xy(1:N,1:2)     raw data, x- and y-components of the time series
            %  SAMPLING        sampling rate (number of samples per second)
            %  TYPE            velocity type: TYPE=2 recommended
            %
            %  OUTPUT:
            %
            %  v(1:N,1:2)      velocity, x- and y-components
            %
            %-------------------------------------------------------------
            N = length(xx);            % length of the time series
            v = zeros(N,2);
            
            switch TYPE
                case 1
                    v(2:N-1,:) = SAMPLING/2*[xx(3:end,:) - xx(1:end-2,:)];
                case 2
                    v(3:N-2,:) = SAMPLING/6*[xx(5:end,:) + xx(4:end-1,:) - xx(2:end-3,:) - xx(1:end-4,:)];
                    v(2,:) = SAMPLING/2*[xx(3,:) - xx(1,:)];
                    v(N-1,:) = SAMPLING/2*[xx(end,:) - xx(end-2,:)];
            end
        end
        
        
        
        function sig = CutOffHighfreq(sig0, dt, Fcutoff)
            % Function: Cutoff high frequency noise
            % Input:
            %       sig0: original signal
            %       dt: temporal resolution of sig0 (sec)
            %       Fcutoff: cutoff frequency  (Hz)
            % Output:
            %       sig: filtered signal
            % Yuwei Cui, Created by May 5, 2011
            
            N = length(sig0);
            
            T = N*dt;
            % Frequency resolution
            f0 = 1/T;
            
            Ncut = round(Fcutoff/f0);
            
            % Fourier transform
            four = fft(sig0);
            
            % Filter
            four(Ncut+2:N-Ncut) = 0;
            
            % Inverse transform
            sig = ifft(four);
            
            % keyboard
        end
    end
end