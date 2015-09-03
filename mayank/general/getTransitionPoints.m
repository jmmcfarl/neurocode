function [up_trans_points,down_trans_points] = getTransitionPoints(Fs,lf8)


nyqf = Fs/2;
hif = 1/nyqf;
% lif = 0.1/nyqf;
[b,a] = butter(2,hif,'low');

lf8_z = lf8 - mean(lf8);
lf8_z = lf8_z/std(lf8_z);

loflf8_z = filtfilt(b,a,lf8_z);

% hif = 50/nyqf;
% lif = 0.1/nyqf;
% [b,a] = butter(2,[lif hif]);
% medlf8 = filtfilt(b,a,lf8);

dsf = 16;
Fsd = Fs/16;

slope_loflf8 = [0;diff(loflf8_z)];

% loflf8_z = loflf8 - mean(loflf8);
% loflf8_z = loflf8_z/std(loflf8_z);
% 
% medlf8 = medlf8 - mean(medlf8);
% medlf8 = medlf8/std(medlf8);

slope_loflf8_z = slope_loflf8- mean(slope_loflf8);
slope_loflf8_z = slope_loflf8_z/std(slope_loflf8_z);

down_slope = downsample(slope_loflf8_z,dsf);

[z,ds_up_points] = findpeaks(down_slope,'minpeakheight',1,'minpeakdistance',round(0.75*Fsd));
% [z,ds_down_points] = findpeaks(-down_slope,'minpeakheight',1,'minpeakdistance',round(0.75*Fsd));

up_points = dsf*ds_up_points;
% down_points = dsf*ds_down_points;

%set first transition to up
% while up_points(1)>down_points(1)
%     down_points(1) = [];
% end
% while up_points(end)>down_points(end)
%     up_points(end) = [];
% end

% %make sure transitions alternate
% badUps = [];
% for i = length(up_points)-1
% 
%     nextDown = find(down_points>up_points(i),1,'first');
%     badUps = [badUps (i + find(up_points(i:end) < nextDown))];
%     
% end
% 
% up_points(badUps) = [];
% 
% badDowns = [];
% for i = 1:length(down_points)-1
%     
%     nextUp = find(up_points>down_points(i));
%     badDowns = [badDowns (i+find(down_points(i:end)< up_points(i+1)))];
% 
% end
% 
% down_points(badDowns) = [];

bad_pts = [];
for i = 1:length(up_points)
   
    temp_pt = up_points(i) + find(slope_loflf8_z(up_points(i):end)<-1,1,'first');
    next_down = find(loflf8_z(temp_pt:end)<0,1,'first');
    if isempty(next_down)
        bad_pts = [bad_pts i];
    else
        down_points(i) = temp_pt+next_down;
    end
    
    clear temp_pt next_down
    
end
down_points = unique(down_points);
up_points(bad_pts) = [];

bad_ups = [];
%take as transition point the nearest 0 crossing around the up_point
for i = 1:length(up_points)

    %if up_point is below mean
    if lf8_z(up_points(i)) < 0
        temp_point = find(lf8_z(up_points(i):end)>0,1,'first')+up_points(i);
        if ~isempty(temp_point)
            up_trans_points(i) = temp_point;
        end
    else
        
        %if up point is above mean
        temp_point = find(lf8_z(1:up_points(i))<0,1,'last');
        if ~isempty(temp_point)
        up_trans_points(i) = temp_point;
        else
            up_trans_points(i) = 0;
        end
    end
    
    
    next_Down = find(down_points>up_points(i),1,'first');
    bad_ups = [bad_ups (i+find(up_points(i+1:end)<down_points(next_Down)))];

end
bad_ups = unique(bad_ups);
up_trans_points(bad_ups) = [];
up_trans_points(up_trans_points == 0) = [];


bad_downs = [];
for i = 1:length(down_points)

    %if down_point is below mean
    if lf8_z(down_points(i)) < 0
        temp_point = find(lf8_z(1:down_points(i))>0,1,'last');
        if ~isempty(temp_point)
            down_trans_points(i) = temp_point;
        end
    else
        
        %if up point is above mean
        temp_point = find(lf8_z(down_points(i):end)<0,1,'first')+down_points(i);
        if ~isempty(temp_point)
        down_trans_points(i) = temp_point;
        else
            down_trans_points(i) = 0;
        end
    end
    
end

down_trans_points(down_trans_points == 0) = [];

for i = 1:length(down_trans_points)
    
    next_Up = find(up_trans_points>down_trans_points(i),1,'first');
    bad_downs = [bad_downs (i+find(down_trans_points(i+1:end)<up_trans_points(next_Up)))];

end
bad_downs = unique(bad_downs);
down_trans_points(bad_downs) = [];

while min(down_trans_points) < min(up_trans_points)
    down_trans_points(1) = [];
end

