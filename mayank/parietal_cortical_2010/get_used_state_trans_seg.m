function [up_state_ind,down_state_ind] = get_used_state_trans_seg(state_seq,Fs,orig_Fs,UDS_segs,edge_buffer)

us_fac = Fs/orig_Fs;
N_segs = size(UDS_segs,1);
UDS_segs = round(UDS_segs*us_fac);
up_state_ind = [];
down_state_ind = [];
for i = 1:N_segs
    up_trans = find(state_seq{i}(1:end-1)==1 & state_seq{i}(2:end) == 2);
    down_trans = find(state_seq{i}(1:end-1)==2 & state_seq{i}(2:end) == 1);
    
    %get rid of transitions near the start and end
    up_trans(up_trans < edge_buffer) = [];
    down_trans(down_trans < edge_buffer) = [];
    up_trans(up_trans > length(state_seq{i})-edge_buffer) = [];
    down_trans(down_trans > length(state_seq{i})-edge_buffer) = [];
    
    up_trans(up_trans > down_trans(end)) = [];
    down_trans(down_trans < up_trans(1)) = [];
        
    up_start = UDS_segs(i,1) + up_trans;
    up_stop = UDS_segs(i,1) + down_trans;
    
    down_start = UDS_segs(i,1) + down_trans(1:end-1);
    down_stop = UDS_segs(i,1) + up_trans(2:end);
    
    up_state_ind = [up_state_ind; up_start(:) up_stop(:)];
    down_state_ind = [down_state_ind; down_start(:) down_stop(:)];
end
