Fs = mean(CSC8_SampleFrequencies);
nyqf = Fs/2;
hif = 2/nyqf;
[b,a] = butter(4,hif,'low');
lofwcv = filtfilt(b,a,wcv_minus_spike);
hif = 1/nyqf;
[b,a] = butter(4,hif,'low');
lowlowwcv = filtfilt(b,a,wcv_minus_spike);
dsf = 20;

slope_lofwcv = [0;diff(lofwcv)];

lofwcv_z = lofwcv - mean(lofwcv);
lofwcv_z = lofwcv_z/std(lofwcv_z);

low_z = downsample(lowlowwcv,100);
low_z = low_z - mean(low_z);
low_z = low_z/std(low_z);
min_points = -findpeaks(-low_z,'minpeakheight',0);
ds_thresh = mean(min_points)+2*std(min_points);

slope_lofwcv_z = slope_lofwcv - mean(slope_lofwcv);
slope_lofwcv_z = slope_lofwcv_z/std(slope_lofwcv_z);

down_slope = downsample(slope_lofwcv_z,dsf);
[z,ds_up_points] = findpeaks(down_slope,'minpeakheight',1.5);

up_points = dsf*ds_up_points;

flat_points = find(abs(slope_lofwcv_z) < 0.1);
diff_flat_points = [0;diff(flat_points)];

flat_points(diff_flat_points == 1) = [];
down_points = flat_points;
down_points(lofwcv_z(down_points)>ds_thresh) = [];



%get rid of down_points before first up point
while min(down_points) < min(up_points)
    down_points(1) = [];
end


%make there is an up state in between every down state
for i= 1:length(up_points)-1
    
    nextDown = down_points(find(down_points>up_points(i),1,'first'));
    nextUp = up_points(i+1);
    down_points(find((down_points>nextDown)&(down_points<nextUp))) = [];
    
end

%make sure there is a down state between every up state

for i = 1:length(down_points)-1
    
    nextUp = up_points(find(up_points>down_points(i),1,'first'));
    nextDown = down_points(i+1);
    if isempty(nextUp)
        break
    end
    up_points(find((up_points>nextUp)&(up_points<nextDown))) = [];
    
end
