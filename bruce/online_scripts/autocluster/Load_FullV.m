function [V,Vtime,Fs] = Load_FullV(name, add_Vmean, filt_cutoff,use_chs)

if nargin < 2 || isempty(add_Vmean)
   add_Vmean = true; 
end
if nargin < 3 || isempty(filt_cutoff)
    filt_cutoff = [nan nan];
end
if nargin < 4 || isempty(use_chs) || any(isnan(use_chs))
    use_chs = 1;
end

load(name);

V = double(FullV.V(use_chs,:))';
lfp_int2V = FullV.intscale(1)/FullV.intscale(2);

if add_Vmean
    first_dot = find(name == '.',1);
    Vmean_fname = [name(1:(first_dot-1)) 'FullVmean.mat'];
    if ~exist(Vmean_fname,'file')
        error('Cant find FullVmean file!');
    end
    load(Vmean_fname);
    V = V + FullV.sumscale*sumv;
end

V = V*lfp_int2V;
Fs = 1/FullV.samper;

if any(~isnan(filt_cutoff))
    niqf = Fs/2;
   if all(~isnan(filt_cutoff)) %bandpass
       [bb,aa] = butter(2,filt_cutoff/niqf);
   elseif ~isnan(filt_cutoff(1))
       [bb,aa] = butter(2,filt_cutoff(1)/niqf,'high');
   elseif ~isnan(filt_cutoff(2))
       [bb,aa] = butter(2,filt_cutoff(2)/niqf,'low');
   end
   V = filtfilt(bb,aa,V);
end

if nargout > 1
    first = 1;
    Vtime = nan(size(V,1),1);
    for j = 1:length(FullV.blklen)
        last = first+FullV.blklen(j)-1;
        Vtime(first:last) = FullV.blkstart(j)+[1:FullV.blklen(j)].*FullV.samper;
        first = last+1;
    end
end

bad_pts = find(isnan(Vtime));
if ~isempty(bad_pts)
%     fprintf('Eliminating %d of %d bad V samples\n',length(bad_pts),length(Vtime));
    
    V(bad_pts,:) = [];
    Vtime(bad_pts) = [];
end