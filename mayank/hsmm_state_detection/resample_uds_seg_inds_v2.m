function [new_seg_inds] = resample_uds_seg_inds_v2(UDS_segs,Fs_orig,Fs_new,state_seq)

Nsegs = size(UDS_segs,1);

for i = 1:Nsegs
    UDS_new_start = round((UDS_segs(i,1)-1)*Fs_new/Fs_orig)+1;
    new_seg_inds(i,:) = [UDS_new_start UDS_new_start+length(state_seq{i})-1];
end
