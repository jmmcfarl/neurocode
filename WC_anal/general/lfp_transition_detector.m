Fs = mean(CSC8_SampleFrequencies);
nyqf = Fs/2;
hif = 1/nyqf;
lif = 0.1/nyqf;
[b,a] = butter(2,[lif hif]);
loflf8 = filtfilt(b,a,lf8);

dsf = 20;

slope_loflf8 = [0;diff(loflf8)];

loflf8_z = loflf8 - mean(loflf8);
loflf8_z = loflf8_z/std(loflf8_z);


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
    if loflf8_z(up_points(i)) < 0
        trans_points(i) = find(loflf8_z(up_points(i):end)>0,1,'first')+up_points(i);
    else
        %if up point is above mean
        trans_points(i) = find(loflf8_z(1:up_points(i))<0,1,'last');
    end

end
