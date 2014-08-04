function peak_stats = locate_peaks_jmm(sig,peak_threshold,edge_threshold,min_width,max_width,usable)

if nargin < 5
    usable = ones(size(sig));
end

sig = sig(:);

thresh_up_crossings = find(sig(1:end-1) < edge_threshold & sig(2:end) >= edge_threshold);
thresh_down_crossings = find(sig(1:end-1) >= edge_threshold & sig(2:end) < edge_threshold);

if sig(1) >= edge_threshold
    thresh_up_crossings = [1; thresh_up_crossings];
end
if sig(end) >= edge_threshold
    thresh_down_crossings = [thresh_down_crossings; length(sig)];
end

peak_durs = thresh_down_crossings-thresh_up_crossings;

bad_peaks = find(peak_durs < min_width | peak_durs > max_width);

thresh_up_crossings(bad_peaks) = [];
thresh_down_crossings(bad_peaks) = [];
peak_durs(bad_peaks) = [];

cur_n_peaks = length(thresh_up_crossings);
peak_amps = nan(cur_n_peaks,1);
peak_locs = nan(cur_n_peaks,1);
usable_peaks = ones(cur_n_peaks,1);
for i = 1:cur_n_peaks
    cur_range = thresh_up_crossings(i):thresh_down_crossings(i);
    [peak_amps(i),cur_loc] = max(sig(cur_range));
    peak_locs(i) = cur_loc + cur_range(1)-1;
    
    if any(usable(cur_range)==0)
        usable_peaks(i) = 0;
    end
end
usable_peaks = logical(usable_peaks);
usable_peaks(peak_amps < peak_threshold) = 0;

n_peaks = sum(usable_peaks);

if n_peaks >= 1
    peak_stats = struct(...
        'start_loc',mat2cell(thresh_up_crossings(usable_peaks),ones(n_peaks,1)),...
        'stop_loc',mat2cell(thresh_down_crossings(usable_peaks),ones(n_peaks,1)),...
        'peak_amp',mat2cell(peak_amps(usable_peaks),ones(n_peaks,1)),...
        'peak_loc',mat2cell(peak_locs(usable_peaks),ones(n_peaks,1)),...
        'peak_dur',mat2cell(peak_durs(usable_peaks),ones(n_peaks,1)));
else
    peak_stats.start_loc = [];
    peak_stats.stop_loc = [];
    peak_stats.peak_amp = [];
    peak_stats.peak_loc = [];
    peak_stats.peak_dur = [];
end

