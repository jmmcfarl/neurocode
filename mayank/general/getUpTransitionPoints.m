function [trans_points] = getUpTransitionPoints(Fs,lf8)


nyqf = Fs/2;
hif = 1/nyqf;
lif = 0.1/nyqf;
[b,a] = butter(2,[lif hif]);
loflf8 = filtfilt(b,a,lf8);

% hif = 50/nyqf;
% lif = 0.1/nyqf;
% [b,a] = butter(2,[lif hif]);
% medlf8 = filtfilt(b,a,lf8);

dsf = 20;

slope_loflf8 = [0;diff(loflf8)];

loflf8_z = loflf8 - mean(loflf8);
loflf8_z = loflf8_z/std(loflf8_z);

% medlf8 = medlf8 - mean(medlf8);
% medlf8 = medlf8/std(medlf8);

slope_loflf8_z = slope_loflf8- mean(slope_loflf8);
slope_loflf8_z = slope_loflf8_z/std(slope_loflf8_z);

down_slope = downsample(slope_loflf8_z,dsf);
[z,ds_up_points] = findpeaks(down_slope,'minpeakheight',1);

up_points = dsf*ds_up_points;



%get rid of transitions where loflf8 doesn't dip below mean inbetween
badPoints = [];
for i = 2:length(up_points)

    dipPoint = find(loflf8_z(up_points(i-1):up_points(i))<0,1,'first');

    if isempty(dipPoint)

        badPoints = [badPoints i];

    end

end

up_points(badPoints) = [];

%take as transition point the nearest 0 crossing around the up_point
for i = 1:length(up_points)

    %if up_point is below mean
    if lf8(up_points(i)) < 0
        temp_point = find(lf8(up_points(i):end)>0,1,'first')+up_points(i);
        if ~isempty(temp_point)
            trans_points(i) = temp_point;
        end
    else
        
        %if up point is above mean
        temp_point = find(lf8(1:up_points(i))<0,1,'last');
        if ~isempty(temp_point)
        trans_points(i) = temp_point;
        else
            trans_points(i) = 0;
        end
    end

end

trans_points(trans_points == 0) = [];


%enforce minimum transition separation
minSep = round(Fs*0.5);
diff_trans = [Fs*100 diff(trans_points)];
trans_points(diff_trans<minSep) = [];



