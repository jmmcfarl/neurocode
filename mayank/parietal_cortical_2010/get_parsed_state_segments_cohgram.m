function [up_segs,down_segs,n_upsegs,n_downsegs] = get_parsed_state_segments_cohgram...
    (state_seq,Fs,movingwin,desynch_times,edge_buffer,state_edges,back_wins)

nwin = round(Fs*movingwin(1));
nwinslide = round(Fs*movingwin(2));

up_trans = find(state_seq(1:end-1)==1 & state_seq(2:end) == 2);
down_trans = find(state_seq(1:end-1)==2 & state_seq(2:end) == 1);

%get rid of transitions near the start and end
up_trans(up_trans < edge_buffer) = [];
down_trans(down_trans < edge_buffer) = [];
up_trans(up_trans > length(state_seq)-edge_buffer) = [];
down_trans(down_trans > length(state_seq)-edge_buffer) = [];

up_trans(up_trans > down_trans(end)) = [];
down_trans(down_trans < up_trans(1)) = [];

%get rid of state trans during desynchronized epochs
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

up_state_ind(:,1) = up_state_ind(:,1)+state_edges(1);
up_state_ind(:,2) = up_state_ind(:,2)-state_edges(2);

down_state_ind(:,1) = down_state_ind(:,1)+state_edges(1);
down_state_ind(:,2) = down_state_ind(:,2)-state_edges(2);

up_durs = (up_state_ind(:,2)-up_state_ind(:,1))/Fs;
down_durs = (down_state_ind(:,2)-down_state_ind(:,1))/Fs;

up_state_ind(up_durs < movingwin(1),:) = [];
down_state_ind(down_durs < movingwin(1),:) = [];

up_nstates = size(up_state_ind,1);
down_nstates = size(down_state_ind,1);

max_up = max(up_durs*Fs);
n_segvecs = 1 + back_wins + ceil((max_up-nwin)/nwinslide);
for n = 1:n_segvecs
    up_segs{n} = [];
end
n_upsegs = zeros(n_segvecs,1);
for i = 1:up_nstates
    cur_nsegs = 1 + back_wins + floor((up_state_ind(i,2) - up_state_ind(i,1) - nwin)/nwinslide);
    for n = 1:cur_nsegs
        up_segs{n} = [up_segs{n} (up_state_ind(i,1) - back_wins*nwinslide + (n-1)*nwinslide)];
        n_upsegs(n) = n_upsegs(n) + 1;
    end
end

max_down = round(Fs*max(down_durs));
n_segvecs = 1 + back_wins + ceil((max_down-nwin)/nwinslide);
for n = 1:n_segvecs
    down_segs{n} = [];
end
n_downsegs = zeros(n_segvecs,1);
for i = 1:down_nstates
    cur_nsegs = 1 + back_wins + floor((down_state_ind(i,2) - down_state_ind(i,1) - nwin)/nwinslide);
    for n = 1:cur_nsegs
        down_segs{n} = [down_segs{n} (down_state_ind(i,1) -back_wins*nwinslide + (n-1)*nwinslide)];
        n_downsegs(n) = n_downsegs(n) + 1;
    end
end
