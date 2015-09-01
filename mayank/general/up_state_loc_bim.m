function [cur_up_start, cur_up_stop] = up_state_loc_bim(dataseg,threshold)

Fs = 5000;
dsf = 10;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,3/niqf,'low');
[b2,a2] = butter(2,1.5/niqf,'low');

maxnoise = 0.1;
minupdur = 1.5;
minpeakheight = 0.2;

filtseg = filtfilt(b,a,dataseg);
filtseg2 = filtfilt(b2,a2,dataseg);
diff_filt = [0;diff(filtseg)]*Fs;

cur_up_start = [];
cur_up_stop = [];


sur_data = dataseg;
sur_data(filtseg2 > threshold) = 1;
sur_data(filtseg2 < threshold) = 0;
dsur_data = [0;diff(sur_data)];
thresh_cross_up = find(dsur_data == 1);
thresh_cross_down = find(dsur_data == -1);

cur_up_start = thresh_cross_up;
cur_up_stop = thresh_cross_down;

if ~isempty(thresh_cross_up) & ~isempty(thresh_cross_down)
    if thresh_cross_up(1) > thresh_cross_down(1)
        cur_up_start = [1; cur_up_start];
    end
    if thresh_cross_up(end) > thresh_cross_down(end)
        cur_up_stop = [cur_up_stop; length(dataseg)];
    end
elseif isempty(thresh_cross_up) & ~isempty(thresh_cross_down)
    cur_up_start = [1 cur_up_start];
elseif ~isempty(thresh_cross_up) & isempty(thresh_cross_down)
    cur_up_stop = [cur_up_stop length(dataseg)];
end


up_durs = (cur_up_stop-cur_up_start)/Fs;
bad_states = find(up_durs < minupdur);

cur_up_start(bad_states) = [];
cur_up_stop(bad_states) = [];

