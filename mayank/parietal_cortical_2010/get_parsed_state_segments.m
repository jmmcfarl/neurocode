function [up_segs,down_segs,up_nsegs,down_nsegs] = get_parsed_state_segments(state_seq,Fs,win,trans,desynch_times)

nwin = round(Fs*win);

up_trans = find(state_seq(1:end-1)==1 & state_seq(2:end) == 2);
down_trans = find(state_seq(1:end-1)==2 & state_seq(2:end) == 1);

% %get rid of transitions near the start and end
% up_trans(up_trans < Fs) = [];
% down_trans(down_trans < Fs) = [];
% up_trans(up_trans > length(state_seq)-Fs) = [];
% down_trans(down_trans > length(state_seq)-Fs) = [];

up_trans(up_trans > down_trans(end)) = [];
down_trans(down_trans < up_trans(1)) = [];

bad_trans = [];
for i = 1:size(desynch_times,1)
    bad_trans = [bad_trans find(up_trans/Fs > desynch_times(i,1) & up_trans/Fs < desynch_times(i,2))];
end
up_trans(bad_trans) = [];
down_trans(bad_trans) = [];

up_start = up_trans;
up_stop = down_trans;

down_start = down_trans(1:end-1);
down_stop = up_trans(2:end);

up_state_ind = [up_start(:) up_stop(:)];
down_state_ind = [down_start(:) down_stop(:)];

up_state_ind(:,1) = up_state_ind(:,1) + round(trans*Fs);
up_state_ind(:,2) = up_state_ind(:,2) - round(trans*Fs);
down_state_ind(:,1) = down_state_ind(:,1) + round(trans*Fs);
down_state_ind(:,2) = down_state_ind(:,2) - round(trans*Fs);

up_durs = (up_state_ind(:,2)-up_state_ind(:,1))/Fs;
down_durs = (down_state_ind(:,2)-down_state_ind(:,1))/Fs;

up_state_ind(up_durs < win,:) = [];
down_state_ind(down_durs < win,:) = [];

up_nsegs = zeros(size(up_state_ind,1),1);
down_nsegs = zeros(size(down_state_ind,1),1);

up_segs = [];
for i = 1:length(up_nsegs)
    %     cur_up_segs = up_state_ind(i,1):nwin:(up_state_ind(i,2)-nwin+1);
%     cur_up_segs = up_state_ind(i,1)+nwin:nwin:(up_state_ind(i,2)-nwin+1);
    cur_up_segs = up_state_ind(i,1);
    up_segs = [up_segs cur_up_segs];
    up_nsegs(i) = length(cur_up_segs);
end

down_segs = [];
for i = 1:length(down_nsegs)
    %     cur_down_segs = down_state_ind(i,1):nwin:(down_state_ind(i,2)-nwin+1);
%     cur_down_segs = down_state_ind(i,1)+nwin:nwin:(down_state_ind(i,2)-nwin+1);
    cur_down_segs = down_state_ind(i,1);
    down_segs = [down_segs cur_down_segs];
    down_nsegs(i) = length(cur_down_segs);
end
