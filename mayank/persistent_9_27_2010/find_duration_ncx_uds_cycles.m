function [mp_updurs_lfpc,mp_downdurs_lfpc] = find_duration_ncx_uds_cycles(mp_uptrans,mp_downtrans,...
    lfp_uplags,lfp_downlags,lfp_cycle_vec)

mp_updurs_lfpc = nan(size(mp_uptrans));
mp_downdurs_lfpc = nan(size(mp_uptrans));
for i = 1:length(mp_uptrans)
    if ~isnan(lfp_uplags(i)) && mp_downtrans(i) > lfp_uplags(i) && mp_downtrans(i)-lfp_uplags(i) < length(lfp_cycle_vec)
        cur_up_end = lfp_cycle_vec(mp_downtrans(i)-lfp_uplags(i));
        cur_up_beg = lfp_cycle_vec(mp_uptrans(i) - lfp_uplags(i));
        mp_updurs_lfpc(i) = cur_up_end - cur_up_beg;
        if length(mp_uptrans) > i && ~isnan(lfp_downlags(i)) && mp_uptrans(i+1) - lfp_downlags(i) < length(lfp_cycle_vec)
            cur_down_end = lfp_cycle_vec(mp_uptrans(i+1) - lfp_downlags(i));
            cur_down_beg = lfp_cycle_vec(mp_downtrans(i) - lfp_downlags(i));
            mp_downdurs_lfpc(i) = cur_down_end - cur_down_beg;
        end
    end
end