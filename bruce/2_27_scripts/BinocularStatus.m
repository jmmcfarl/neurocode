function [fixations saccades status] = BinocularStatus(LeyeXY, ReyeXY, SamplingFreq, vthreshsac, vthreshfix)
% binocular status analysis

assert(length(LeyeXY)==length(ReyeXY));
L = length(LeyeXY);

tres = 1/SamplingFreq; % sec
eyet = (0:L-1)*tres;
LeyeV = EyeM.GetEyeVelocity(LeyeXY, eyet);
ReyeV = EyeM.GetEyeVelocity(ReyeXY, eyet);

fixating = (abs(LeyeV(:,1))<vthreshfix & abs(LeyeV(:,2))<vthreshfix) & ...
    (abs(ReyeV(:,1))<vthreshfix & abs(ReyeV(:,2))<vthreshfix);

saccading = (abs(LeyeV(:,1))>vthreshsac | abs(LeyeV(:,2))>vthreshsac) | ...
    (abs(ReyeV(:,1))>vthreshsac | abs(ReyeV(:,2))>vthreshsac);

status = zeros(L, 1);
status(fixating==1) = 1;
status(saccading==1) = 2;
saccades = [];
fixations = [];

minfixdur = round(0.05/tres);
t=1;
currentstate = 0;
while t<L
    % update current state
    if status(t)==1 % fixating        
        nextstatet = find(fixating(t+1:end)~=1,1,'first') + t;
        if isempty(nextstatet) nextstatet=L; end
        fixdur = nextstatet-t ;  % fixation duration
        if fixdur > minfixdur
            fixations = [fixations; [t nextstatet-1]];
        end
        
        t = nextstatet;
    elseif status(t)==2
        nextstatet = find(saccading(t+1:end)~=1,1,'first') + t;
        if isempty(nextstatet) nextstatet=L; end
        sacdur = nextstatet-t ;  % fixation duration
    
        saccades = [saccades; [t nextstatet-1]];
    
        
        t = nextstatet;        
    end
    t = t+1;
end