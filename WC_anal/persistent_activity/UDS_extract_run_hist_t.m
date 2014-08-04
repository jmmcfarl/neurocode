
function [up_state_dur,down_state_dur,up_trans,down_trans,sig_f2,int_thresh] = UDS_extract_run_hist(sig,threshold)

%minimum slope
min_slope = 1;
min_slope_dist = 10;

%minimum upstate duration
min_up_dur = 0.3;
min_down_dur = 0.3;

%first filter signal
niqf = 2016/2;
hcf1 = 1/niqf;
hcf2 = 5/niqf;
lcf = 0.1/niqf;

[b1,a1] = butter(2,[lcf hcf1]);
[b2,a2] = butter(2,[lcf hcf2]);
% [b1,a1] = butter(2,hcf1,'low');
sig_f1 = filtfilt(b1,a1,sig);
sig_f2 = filtfilt(b2,a2,sig);

%now downsample
dsf = 8;
Fsd = 2016/dsf;

sig_f1 = downsample(sig_f1,dsf);
sig_f2 = downsample(sig_f2,dsf);

m_sig = mean(sig_f1);
s_sig = std(sig_f1);

sig_f1 = (sig_f1 - m_sig)/s_sig;
sig_f2 = (sig_f2 - m_sig)/s_sig;

%interpolate threshold function
int_thresh = spline(1:length(threshold),threshold,linspace(1,length(threshold),length(sig_f1)));
int_thresh = int_thresh';

if size(sig_f1,1) ~= size(int_thresh,1)
    int_thresh = int_thresh';
end
thresh_det = sig_f1>int_thresh;

dtd = [0;diff(thresh_det)];
up_cross = find(dtd==1);
down_cross = find(dtd == -1);


%throw away crossings in the first second
up_cross(up_cross <=Fsd) = [];
down_cross(down_cross <= Fsd) = [];

%get rid of overlapping transitions (sub 100 ms states)
up_cross(find(diff(up_cross)<Fsd*.1)) = [];
down_cross(find(diff(down_cross)<Fsd*.1)) = [];

up_cross(end) = [];
down_cross(end) = [];

up_cross_times = up_cross/Fsd;
down_cross_times = down_cross/Fsd;


%cycle through threshold crossings and find where the less filtered signal
%crosses the interpolated threshold value
for i = 1:length(up_cross)-1

    if sig_f2(up_cross(i)) > int_thresh(up_cross(i))
        temp = find(sig_f2(up_cross(i)-Fsd:up_cross(i))< int_thresh(up_cross(i)),1,'last');
        if isempty(temp)
            up_trans(i) = NaN;
        elseif up_cross(i) > Fsd+temp
            up_trans(i) = up_cross(i)-Fsd+temp;
        else
            up_trans(i) = up_cross(i);
        end
    else
        temp = find(sig_f2(up_cross(i):up_cross(i)+Fsd)> int_thresh(up_cross(i)),1,'first');
        if isempty(temp)
            up_trans(i) = NaN;
        elseif up_cross(i)+temp < length(sig_f2)
            up_trans(i) = up_cross(i)+temp;
        else
            up_trans(i) = up_cross(i);
        end
    end

end

for i = 1:length(down_cross)-1

    if sig_f2(down_cross(i)) < int_thresh(down_cross(i))
        temp = find(sig_f2(down_cross(i)-Fsd:down_cross(i))>int_thresh(down_cross(i)),1,'last');
        if isempty(temp)
            down_trans(i) = NaN;
        elseif down_cross(i) > Fsd+temp
            down_trans(i) = down_cross(i)-Fsd+temp;
        else
            down_trans(i) = down_cross(i);
        end
    else
        temp = find(sig_f2(down_cross(i):down_cross(i)+Fsd)<int_thresh(down_cross(i)),1,'first');
        if isempty(temp)
            down_trans(i) = NaN;
        elseif down_cross(i) +temp < length(sig_f2)
            down_trans(i) = down_cross(i)+temp;
        else
            down_trans(i) = down_cross(i);
        end
    end

end
up_trans(isnan(up_trans)) = [];
down_trans(isnan(down_trans)) = [];

% up_trans(i+1) = up_cross(i+1);
% down_trans(i+1) = down_cross(i+1);

%get rid of overlapping transitions (sub 100 ms states)
up_trans(find(diff(up_trans)<Fsd*.1)) = [];
down_trans(find(diff(down_trans)<Fsd*.1)) = [];

%make sure start on up transition and end on down transition
while min(up_trans) > min(down_trans)
    down_trans(1) = [];
end
while max(up_trans) > max(down_trans)
    up_trans(end) = [];
end

%find slope at each transition
sig_slope = [0;diff(sig_f2)];
sig_slope = zscore(sig_slope);

up_trans_slope = sig_slope(up_trans);
bad_up_trans = [];
bad_down_trans = [];
for i = 1:length(up_trans)

    if up_trans_slope(i) < min_slope

        next_good_slope = find(sig_slope(up_trans(i):end)>min_slope,1,'first');
        %         prev_good_slope = up_trans(i)-find(sig_slope(1:up_trans(i))>min_slope,1,'last');

        %         if prev_good_slope/Fsd < min_slope_dist
        %
        %             up_trans(i) = up_trans(i) - prev_good_slope;
        %
        %         elseif next_good_slope < min_slope_dist
        %
        %             up_trans(i) = up_trans(i) + next_good_slope;
        %
        %         end
        next_down_slope = find(down_trans>up_trans(i),1,'first');
        if next_good_slope/Fsd < min_slope_dist & (next_good_slope+up_trans(i))<down_trans(next_down_slope)
            up_trans(i) = up_trans(i)+next_good_slope;
        else
           bad_up_trans = [bad_up_trans i];
           bad_down_trans = [bad_down_trans next_down_slope];
        end

    end

end

up_trans(bad_up_trans) = [];
down_trans(bad_down_trans) = [];

%get rid of overlapping transitions
badup = find(diff(up_trans)<Fsd*.1);
baddown = find(diff(down_trans)<Fsd*.1);

badup = [badup (baddown)];
baddown = [baddown (badup)];

up_trans(badup) = [];
down_trans(baddown) = [];

cur_state = 1;
cur_time = up_trans(1);

state_vec = zeros(size(sig_f1));
while cur_time < length(sig_f1)
    if cur_state == 1
        next_down = find(down_trans>cur_time,1,'first');
        if isempty(next_down)
            break
        end
        cur_state = 0;
        state_vec(cur_time:down_trans(next_down)) = 1;
        cur_time = down_trans(next_down);
    else
        next_up = find(up_trans>cur_time,1,'first');
        if isempty(next_up)
            break
        end
        cur_state = 1;
        cur_time = up_trans(next_up);
    end
end

bad_ups = find(state_vec(up_trans-1) == 1);
bad_downs = find(state_vec(down_trans) == 0);
up_trans(bad_ups) = [];
down_trans(bad_downs) = [];
% end
if length(up_trans) > length(down_trans)
    up_trans(end) = [];
end
up_state_dur = (down_trans - up_trans)/Fsd;
down_state_dur = (up_trans(2:end) - down_trans(1:end-1))/Fsd;

bad_downs = find(down_state_dur < min_down_dur);


disp(sprintf('Eliminated %d of %d down states',length(bad_downs),length(down_trans)+length(bad_downs)))

down_trans(bad_downs) = [];
up_trans(bad_downs+1) = [];

up_state_dur = (down_trans - up_trans)/Fsd;

bad_ups = find(up_state_dur < min_up_dur);

disp(sprintf('Eliminated %d of %d ups',length(bad_ups),length(up_trans)+length(bad_ups)))

up_trans(bad_ups) = [];
down_trans(bad_ups) = [];


up_state_dur = (down_trans - up_trans)/Fsd;
down_state_dur = (up_trans(2:end) - down_trans(1:end-1))/Fsd;




