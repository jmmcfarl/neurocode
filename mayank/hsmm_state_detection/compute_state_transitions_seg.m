function [up_trans_inds,down_trans_inds,up_trans_segnum,down_trans_segnum] = ...
    compute_state_transitions_seg(seg_inds,state_seq)

up_trans_inds = [];
down_trans_inds = [];
up_trans_segnum = []; 
down_trans_segnum = [];
for i = 1:length(state_seq)
    up_trans = 1+find(state_seq{i}(2:end)==2 & state_seq{i}(1:end-1)==1);
    down_trans = 1+find(state_seq{i}(2:end)==1 & state_seq{i}(1:end-1)==2);
    
    if ~isempty(up_trans) && ~isempty(down_trans)
    %make sure you start on an up-transition and end on a down transition
    down_trans(down_trans < up_trans(1)) = [];
    up_trans(up_trans > down_trans(end)) = [];

    up_trans = up_trans(:);
    down_trans = down_trans(:);
    
    up_trans_inds = [up_trans_inds; (up_trans + seg_inds(i,1)-1)];
    down_trans_inds = [down_trans_inds; (down_trans + seg_inds(i,1)-1)];
    up_trans_segnum = [up_trans_segnum; ones(size(up_trans))*i];
    down_trans_segnum = [down_trans_segnum; ones(size(down_trans))*i];
    end
end