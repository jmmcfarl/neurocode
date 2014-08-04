function [up_points,down_points,up_start_sh,up_end_sh,lofwcv_z] = get_transition_times(wcv_minus_spike,Fs)


clear up_points down_points up_start_sh up_end_sh
% Fs = mean(CSC8_SampleFrequencies);
nyqf = Fs/2;
hif = 2/nyqf;
[b,a] = butter(4,hif,'low');
lofwcv = filtfilt(b,a,wcv_minus_spike);
dsf = 20;
slope_lofwcv = [0;diff(lofwcv)];

hif = 5/nyqf;
[b,a] = butter(4,hif,'low');
medwcv = filtfilt(b,a,wcv_minus_spike);
slope_medwcv = [0;diff(medwcv)];

lofwcv_z = lofwcv - mean(lofwcv);
lofwcv_z = lofwcv_z/std(lofwcv_z);

ds_thresh = -1;

% hif = 1/nyqf;
% [b,a] = butter(4,hif,'low');
% lowlowwcv = filtfilt(b,a,wcv_minus_spike);
% low_z = downsample(lowlowwcv,100);
% low_z = low_z - mean(low_z);
% low_z = low_z/std(low_z);
% min_points = -findpeaks(-low_z,'minpeakheight',0);
% ds_thresh = mean(min_points)+2*std(min_points);

slope_lofwcv_z = slope_lofwcv - mean(slope_lofwcv);
slope_lofwcv_z = slope_lofwcv_z/std(slope_lofwcv_z);

ds_slope = downsample(slope_lofwcv_z,dsf);
[z,ds_up_points] = findpeaks(ds_slope,'minpeakheight',1);

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
% for i= 1:length(up_points)-1
%     
%     nextDown = down_points(find(down_points>up_points(i),1,'first'));
%     nextUp = up_points(i+1);
%     down_points(find((down_points>nextDown)&(down_points<nextUp))) = [];
%     
% end


for i= 1:length(up_points)-1
    
    nextDown = down_points(find(down_points>up_points(i),1,'first'));
    nextUp = up_points(i+1);
    if isempty(nextDown)
        break
    end
    possible_downs_ind = find((down_points>=nextDown)&(down_points<nextUp));
    [downVal,minDown] = min(lofwcv_z(down_points(possible_downs_ind)));
    badDown = [1:length(possible_downs_ind)];
    badDown(minDown) = [];
    down_points(possible_downs_ind(badDown)) = [];
    clear possible_downs_ind minDown badDown
    
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

%make sure you end on a down state
while max(up_points) > max(down_points)
    up_points(end) = [];
end

%get rid of states where the up doesn't get above mean
badTrans = [];
for i = 1:length(up_points)-1
    
   upBeg = up_points(i);
   upEnd = down_points(i);
   maxVal = max(lofwcv_z(upBeg:upEnd));
   if maxVal < 0
badTrans = [badTrans i];
   end
    
end

up_points(badTrans) = [];
down_points(badTrans) = [];



%% NOW CALCULATE SHOULDER POINTS

%up state start shoulder is first point derivative becomes negative after
%up transition
minUpDur = round(Fs*0.1);
for i = 1:length(up_points)-1
    up_start_sh(i) = up_points(i) + find((slope_medwcv(up_points(i):end)<0)&(lofwcv_z(up_points(i):end)>-0.5),1,'first');
    up_end_sh(i) = find((slope_medwcv(1:down_points(i))>0)&(lofwcv_z(1:down_points(i))>-0.5),1,'last');
    if (up_end_sh(i) - up_start_sh(i)) < minUpDur
       up_end_sh(i) = down_points(i); 
    end
end
