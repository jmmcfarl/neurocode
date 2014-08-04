
function [up_state_dur,down_state_dur,up_trans,down_trans,sig_f1,sig_f2,threshold,effect_size] = UDS_extract_9_5(sig)

%first downsample to 252Hz
dsf = 8;
Fsd = 2016/dsf;
sig_d = downsample(sig,dsf);

min_effect_size = 2;

%zscore
sig_z = zscore(sig_d);

%filter low and medium
niqf = Fsd/2;
hcf1 = 1/niqf;
hcf2 = 5/niqf;
% lcf = .005/niqf;
[b1,a1] = butter(2,hcf1,'low');
[b2,a2] = butter(2,hcf2,'low');

sig_f1 = filtfilt(b1,a1,sig_z);
sig_f2 = filtfilt(b2,a2,sig_z);

segTime = 30;
segLength = segTime*Fsd;
winSlide = segLength-round(Fsd*.1);
numSegs = floor((length(sig_z)-segLength)/winSlide);
seg_Taxis = [1:numSegs]*segTime-round(segTime)/2;
up_cross = [];
down_cross = [];
for i = 1:numSegs

    segBeg = (i-1)*winSlide+1;
    segEnd = segBeg+segLength;

    wseg_f1 = sig_f1(segBeg:segEnd);
    wseg_f2 = sig_f2(segBeg:segEnd);

    %get GMM for sig in segment

    wseg_d2 = downsample(wseg_f2,2);
    [u,sig,t,iter] = fit_mix_gaussian(wseg_d2,2);
    effect_size(i) = get_effect_size(u,sig);
    if effect_size(i) > min_effect_size
        [threshold(i)] = find_GMM_xp(u,sig,t);
        if abs(threshold(i)) > 2
            threshold(i) = (max(wseg_f1)+min(wseg_f1))/2;
        end
    else
        threshold(i) = (max(wseg_f1) + min(wseg_f1))/2;
    end

    %find threshold crossings
    thresh_det = wseg_f1;
    above_points = find(thresh_det>threshold(i));
    below_points = find(thresh_det<threshold(i));
    thresh_det(above_points) = 1;
    thresh_det(below_points) = 0;

    dtd = [0;diff(thresh_det)];
    up_cross = [up_cross; (find(dtd==1) + segBeg)];
    down_cross = [down_cross; (find(dtd == -1) + segBeg)];

%     plot_mix_gaussian( u,sig,t,wseg_d2 )
%     pause
%     clf
    
%         plot(wseg_f1)
%         hold on
%         plot(down_cross,wseg_f1(down_cross),'ro')
%         plot(up_cross,wseg_f1(up_cross),'go')
%         pause;clf


end

%throw away crossings in the first second
up_cross(up_cross <=Fsd) = [];
down_cross(down_cross <= Fsd) = [];

%get rid of overlapping transitions
up_cross(find(diff(up_cross)<Fsd*1)) = [];
down_cross(find(diff(down_cross)<Fsd*.1)) = [];

up_cross(end) = [];
down_cross(end) = [];

up_cross_times = up_cross/Fsd;
down_cross_times = down_cross/Fsd;

%get a local threshold by interpolating between threshold points
up_thresh = spline(seg_Taxis,threshold,up_cross_times);
down_thresh = spline(seg_Taxis,threshold,down_cross_times);
up_thresh(up_cross_times<seg_Taxis(1)) = threshold(1);
up_thresh(up_cross_times>seg_Taxis(end)) = threshold(end);
down_thresh(down_cross_times<seg_Taxis(1)) = threshold(1);
down_thresh(down_cross_times>seg_Taxis(end)) = threshold(end);

%cycle through threshold crossings and find where the less filtered signal
%crosses the interpolated threshold value
for i = 1:length(up_cross)

    if sig_f2(up_cross(i)) > up_thresh(i)
        temp = find(sig_f2(up_cross(i)-Fsd:up_cross(i))<up_thresh(i),1,'last');
        if isempty(temp)
            up_trans(i) = NaN;
        else
            up_trans(i) = up_cross(i)-Fsd+temp;
        end
    else
        temp = find(sig_f2(up_cross(i):up_cross(i)+Fsd)>up_thresh(i),1,'first');
        if isempty(temp)
            up_trans(i) = NaN;
        else
            up_trans(i) = up_cross(i)+temp;
        end
    end

end

for i = 1:length(down_cross)

    if sig_f2(down_cross(i)) < down_thresh(i)
        temp = find(sig_f2(down_cross(i)-Fsd:down_cross(i))>down_thresh(i),1,'last');
        if isempty(temp)
            down_trans(i) = NaN;
        else
            down_trans(i) = down_cross(i)-Fsd+temp;
        end
    else
        temp = find(sig_f2(down_cross(i):down_cross(i)+Fsd)<down_thresh(i),1,'first');
        if isempty(temp)
            down_trans(i) = NaN;
        else
            down_trans(i) = down_cross(i)+temp;
        end
    end

end
up_trans(isnan(up_trans)) = [];
down_trans(isnan(down_trans)) = [];


%make sure start on up transition and end on down transition
while min(up_trans) > min(down_trans)
    down_trans(1) = [];
end
while max(up_trans) > max(down_trans)
    up_trans(end) = [];
end


%get rid of overlapping transitions
badup = find(diff(up_trans)<Fsd*.1);
baddown = find(diff(down_trans)<Fsd*.1);

badup = [badup (baddown-1)];
baddown = [baddown (badup-1)];

up_trans(badup) = [];
down_trans(baddown) = [];

% if length(up_trans) ~= length(down_trans)
    %make sure state transitions alternate
    cur_state = 1;
    cur_time = up_trans(1);
    % bad_downs = [];
    % bad_ups = [];
    % while cur_time < length(sig_f1)
    %     if cur_state == 1
    %         next_down = find(down_trans>cur_time,1,'first');
    %         if isempty(next_down)
    %             break
    %         end
    %         bad_ups = [bad_ups find(up_trans>cur_time & up_trans < down_trans(next_down))];
    %         cur_state = 0;
    %         cur_time = down_trans(next_down);
    %     else
    %         next_up = find(up_trans>cur_time,1,'first');
    %         if isempty(next_up)
    %             break
    %         end
    %         bad_downs = [bad_downs find(down_trans>cur_time & down_trans < up_trans(next_up))];
    %         cur_state = 1;
    %         cur_time = up_trans(next_up);
    %     end
    % end

    % up_trans(bad_ups) = [];
    % down_trans(bad_downs) = [];

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
up_state_dur = (down_trans - up_trans)/Fsd;
down_state_dur = (up_trans(2:end) - down_trans(1:end-1))/Fsd;


