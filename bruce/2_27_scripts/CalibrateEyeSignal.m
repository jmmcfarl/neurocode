function EyeCal = CalibrateEyeSignal(EyeSigRaw, Gain, Offset, HigherGain)
% Raw signal should be organized as
%    [Leye-H, Leye-V, Reye-H, Reye-V] for two eyes mode
% or [eye-H eye-V] for one eye mode
% Gain should be 4x2 or 2x2
% Offset should be 4x1 or 2x1
% Yuwei Cui Feb 28 2012, Last Modified Apr 3
[Ndata Neye] = size(EyeSigRaw);
Neye = round(Neye/2);
% fprintf('Signal for %d Eye(s)\n', Neye);

EyeLdatas = EyeSigRaw(:, [1 2]);
EyeLdatascal = zeros(size(EyeLdatas));
if Neye>1
    EyeRdatas = EyeSigRaw(:, [3 4]);
    EyeRdatascal = zeros(size(EyeRdatas));
end

for eyei = 1:Neye
    for hv=1:2
        gi = 2*(eyei-1)+hv;
        if eyei==1
            EyeLdatascal(:,hv)=EyeLdatas(:, 1:2)*Gain(gi,:)'+Offset(gi); % linear calibration
            if nargin>3
            EyeLdatascal(:,hv) = EyeLdatascal(:,hv)+EyeLdatas(:,1).^2*HigherGain(gi,1)+EyeLdatas(:,2).^2*HigherGain(gi,2)+EyeLdatas(:,1).*EyeLdatas(:,2)*HigherGain(gi,3);
            end
        elseif eyei==2
            EyeRdatascal(:,hv)=EyeRdatas(:, 1:2)*Gain(gi,:)'+Offset(gi); % linear calibration
            if nargin>3
            EyeRdatascal(:,hv) = EyeRdatascal(:,hv)+EyeRdatas(:,1).^2*HigherGain(gi,1)+EyeRdatas(:,2).^2*HigherGain(gi,2)+EyeRdatas(:,1).*EyeRdatas(:,2)*HigherGain(gi,3);
            end
        end
    end
end
if Neye>1
    EyeCal = [EyeLdatascal EyeRdatascal];
else
    EyeCal = EyeLdatascal;
end