function [up_segs,down_segs,n_upsegs,n_downsegs] = get_parsed_state_segments_cohgram_seg...
    (state_seq,orig_Fs,Fs,movingwin,UDS_segs,edge_buffer,state_edges,back_wins)

us_fac = Fs/orig_Fs;
UDS_segs = round(UDS_segs*us_fac);
nwin = round(Fs*movingwin(1));
nwinslide = round(Fs*movingwin(2));

[state_durations] = compute_state_durations_seg(state_seq,orig_Fs);
max_up = max(state_durations{2})*Fs;
max_down = max(state_durations{1})*Fs;
nup_segvecs = 1 + back_wins + ceil((max_up-nwin)/nwinslide);
for n = 1:nup_segvecs
    up_segs{n} = [];
end
n_upsegs = zeros(nup_segvecs,1);
ndown_segvecs = 1 + back_wins + ceil((max_down-nwin)/nwinslide);
for n = 1:ndown_segvecs
    down_segs{n} = [];
end
n_downsegs = zeros(ndown_segvecs,1);

for ns = 1:size(UDS_segs,1)
    
    up_trans = find(state_seq{ns}(1:end-1)==1 & state_seq{ns}(2:end) == 2);
    down_trans = find(state_seq{ns}(1:end-1)==2 & state_seq{ns}(2:end) == 1);
    
    up_trans = round(up_trans*us_fac);
    down_trans = round(down_trans*us_fac);
    
    %get rid of transitions near the start and end
    up_trans(up_trans < edge_buffer) = [];
    down_trans(down_trans < edge_buffer) = [];
    up_trans(up_trans > length(state_seq{ns})-edge_buffer) = [];
    down_trans(down_trans > length(state_seq{ns})-edge_buffer) = [];
    
    up_trans(up_trans > down_trans(end)) = [];
    down_trans(down_trans < up_trans(1)) = [];
        
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
    
    for i = 1:up_nstates
        prev_down_trans = down_state_ind(find(down_state_ind(:,1) < up_state_ind(i,1),1,'last'),1);
        cur_nsegs = 1 + back_wins + floor((up_state_ind(i,2) - up_state_ind(i,1) - nwin)/nwinslide);
        for n = 1:cur_nsegs
            cur_win = up_state_ind(i,1)-back_wins*nwinslide+(n-1)*nwinslide;
            if cur_win > prev_down_trans
                up_segs{n} = [up_segs{n} cur_win];
                n_upsegs(n) = n_upsegs(n) + 1;
            end
        end
    end
    
    for i = 1:down_nstates
        prev_up_trans = up_state_ind(find(up_state_ind(:,1) < down_state_ind(i,1),1,'last'),1);
        cur_nsegs = 1 + back_wins + floor((down_state_ind(i,2) - down_state_ind(i,1) - nwin)/nwinslide);
        for n = 1:cur_nsegs
            cur_win = down_state_ind(i,1)-back_wins*nwinslide+(n-1)*nwinslide;
            if cur_win > prev_up_trans
                down_segs{n} = [down_segs{n} cur_win];
                n_downsegs(n) = n_downsegs(n) + 1;
            end
        end
    end        
end
