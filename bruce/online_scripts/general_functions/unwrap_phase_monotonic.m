function unwrapped_phase = unwrap_phase_monotonic(phase)

NT = length(phase);
unwrapped_phase = unwrap(phase);
rev_points = 1 + find(diff(unwrapped_phase(1:end-1)) > 0 & diff(unwrapped_phase(2:end)) < 0);
int_inds = [];
for i = 1:length(rev_points)
    next_use_point = rev_points(i) + find(unwrapped_phase(rev_points(i)+1:end) >= unwrapped_phase(rev_points(i)),1,'first');
    if ~isempty(next_use_point)
    cur_inds = (rev_points(i)+1):next_use_point;
    else
        cur_inds = (rev_points(i)+1):NT;
    end
    int_inds = [int_inds cur_inds];
end
all_inds = 1:NT;
use_inds = setdiff(all_inds,int_inds);
unwrapped_phase = interp1(all_inds(use_inds),unwrapped_phase(use_inds),all_inds,'pchip');