function [cur_down_start, cur_down_stop] = down_state_loc_bim(dataseg,threshold)

Fs = 5000;
dsf = 10;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,3/niqf,'low');
[b2,a2] = butter(2,1.5/niqf,'low');

maxnoise = 0.1;
mindowndur = 0.75; %min dur in sec
minupdur = 1.0;
minpeakheight = 0.2;

filtseg = filtfilt(b,a,dataseg);
filtseg2 = filtfilt(b2,a2,dataseg);
diff_filt = [0;diff(filtseg)]*Fs;

cur_down_start = [];
cur_down_stop = [];
cur_up_start = [];
cur_up_stop = [];


sur_data = dataseg;
sur_data(filtseg2 > threshold) = 1;
sur_data(filtseg2 < threshold) = 0;
dsur_data = [0;diff(sur_data)];
thresh_cross_up = find(dsur_data == 1);
thresh_cross_down = find(dsur_data == -1);


%cycle through threshold crossings and find low derivative times
for u = 1:length(thresh_cross_up)
    prev_down_end = find(abs(diff_filt(1:thresh_cross_up(u))) < maxnoise,1,'last');
    if ~isempty(prev_down_end)
        cur_down_stop = [cur_down_stop prev_down_end];
    end
end
for d = 1:length(thresh_cross_down)
    next_down_start = thresh_cross_down(d)+find(abs(diff_filt(thresh_cross_down(d):end)) < maxnoise,1,'first');
    if ~isempty(next_down_start)
        cur_down_start = [cur_down_start next_down_start];
    end
end

if ~isempty(cur_down_stop) & ~isempty(thresh_cross_down)
    if cur_down_stop(1) < thresh_cross_down(1)
        cur_down_start = [cur_down_start 1];
        cur_down_start = sort(cur_down_start);
    end
elseif ~isempty(cur_down_stop) & isempty(cur_down_start)    
    cur_down_start = [cur_down_start 1];  
    cur_down_start = sort(cur_down_start);
end

if ~isempty(cur_down_start) & ~isempty(thresh_cross_up)
    if cur_down_start(end) > thresh_cross_up(end)
        cur_down_stop = [cur_down_stop length(dataseg)];
        cur_down_stop = sort(cur_down_stop);
    end
elseif ~isempty(cur_down_start) & isempty(cur_down_stop)
    cur_down_stop = [cur_down_stop length(dataseg)];
    cur_down_stop = sort(cur_down_stop);
end

cur_down_start = sort(cur_down_start);
cur_down_stop = sort(cur_down_stop);

if ~isempty(cur_down_start) & ~isempty(cur_down_stop)
    if cur_down_stop(1) < cur_down_start(1)
        cur_down_stop(1) = [];
    end
    if cur_down_stop(end) < cur_down_start(end)
        cur_down_start(end) = [];
    end
end

if length(cur_down_stop) ~= length(cur_down_start) 
    cur_down_stop = [];
    cur_down_start = [];
else
    down_state_dur = (cur_down_stop-cur_down_start)/Fs;
    cur_down_start(down_state_dur < mindowndur) = [];
    cur_down_stop(down_state_dur < mindowndur) = [];
end

