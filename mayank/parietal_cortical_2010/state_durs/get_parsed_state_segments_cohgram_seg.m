function [up_segs,down_segs,n_upsegs,n_downsegs] = get_parsed_state_segments_cohgram_seg...
    (state_seq,Fs,orig_Fs,movingwin,UDS_segs,edge_buffer,state_edges,back_wins)

us_fac = Fs/orig_Fs;
UDS_segs = round(UDS_segs*us_fac);
nwin = round(Fs*movingwin(1));
nwinslide = round(Fs*movingwin(2));

up_segs = [];
down_segs = [];

for ns = 1:size(UDS_segs,1)
    
    up_trans = find(state_seq{ns}(1:end-1)==1 & state_seq{ns}(2:end) == 2);
    down_trans = find(state_seq{ns}(1:end-1)==2 & state_seq{ns}(2:end) == 1);
    
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
    
    max_up = max(up_durs*Fs);
    n_segvecs = 1 + back_wins + ceil((max_up-nwin)/nwinslide);
    for n = 1:n_segvecs
        cup_segs{n} = [];
    end
    n_upsegs{ns} = zeros(n_segvecs,1);
    for i = 1:up_nstates
        cur_nsegs = 1 + back_wins + floor((up_state_ind(i,2) - up_state_ind(i,1) - nwin)/nwinslide);
        for n = 1:cur_nsegs
            cup_segs{n} = [cup_segs{n} (up_state_ind(i,1) - back_wins*nwinslide + (n-1)*nwinslide)];
            n_upsegs{ns}(n) = n_upsegs{ns}(n) + 1;
        end
    end
    
    max_down = round(Fs*max(down_durs));
    n_segvecs = 1 + back_wins + ceil((max_down-nwin)/nwinslide);
    for n = 1:n_segvecs
        cdown_segs{n} = [];
    end
    n_downsegs{ns} = zeros(n_segvecs,1);
    for i = 1:down_nstates
        cur_nsegs = 1 + back_wins + floor((down_state_ind(i,2) - down_state_ind(i,1) - nwin)/nwinslide);
        for n = 1:cur_nsegs
            cdown_segs{n} = [cdown_segs{n} (down_state_ind(i,1) -back_wins*nwinslide + (n-1)*nwinslide)];
            n_downsegs{ns}(n) = n_downsegs{ns}(n) + 1;
        end
    end
    
    up_segs = [up_segs cup_segs];
    down_segs = [down_segs cdown_segs];
    
end
