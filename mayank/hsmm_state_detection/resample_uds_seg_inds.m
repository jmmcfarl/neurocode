function [new_seg_inds] = resample_uds_seg_inds(UDS_segs,Fs_orig,Fs_new,newT)

Nsegs = size(UDS_segs,1);

for i = 1:Nsegs
    UDS_new_start = round((UDS_segs(i,1)-1)*Fs_new/Fs_orig)+1;
    curT = UDS_segs(i,2)-UDS_segs(i,1)+1;
    enddiff = newT-UDS_new_start+1;
    curT_new = min(round(Fs_new/Fs_orig*(curT+1)),enddiff);
    new_seg_inds(i,:) = [UDS_new_start UDS_new_start+curT_new-1];
end
