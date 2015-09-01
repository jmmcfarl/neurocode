function [new_seg_inds] = resample_uds_seg_inds(UDS_segs,Fs_orig,Fs_new,newT)

%convert index matrix specifying beginning and end of each UDS segment to a
%different sample-frequency. (Two sample-frequencies must be related by an
%integer down-sampling factor in this case)

%%INPUTS
%    UDS_segs: original Nx2 UDS segment matrix
%    Fs_orig: old sample-frequency (Hz)
%    Fs_new: new sample-frequency (Hz)
%    newT: number of samples in the new data vector
%%OUTPUTS
%    new_seg_inds: Nx2 UDS segment matrix matching the new sample-frequency
   

Nsegs = size(UDS_segs,1);

for i = 1:Nsegs
    UDS_new_start = round((UDS_segs(i,1)-1)*Fs_new/Fs_orig)+1;
    curT = UDS_segs(i,2)-UDS_segs(i,1)+1;
    enddiff = newT-UDS_new_start+1;
    curT_new = min(round(Fs_new/Fs_orig*(curT+1)),enddiff);
    new_seg_inds(i,:) = [UDS_new_start UDS_new_start+curT_new-1];
end
