function [bs_tslu,bs_ttnu,bs_tsld,bs_ttnd] = get_bs_reltime_samps(up_trans8,down_trans8,datalen)

num_samps = 1000;
Fsd = 252;

rand_inds = ceil(rand(num_samps,1)*datalen);

bs_tslu = nan(num_samps,1);
bs_ttnu = nan(num_samps,1);
bs_tsld = nan(num_samps,1);
bs_ttnd = nan(num_samps,1);

for i = 1:num_samps

    cur_trans = rand_inds(i);

    [dummy,nearest_lfp_up] = min(abs(up_trans8-cur_trans));
    if ~isempty(nearest_lfp_up)
        bs_ttnu(i) = (cur_trans-up_trans8(nearest_lfp_up))/Fsd;
    end

    previous_lfp_up = find(up_trans8 < cur_trans,1,'last');
    if ~isempty(previous_lfp_up)
        bs_tslu(i) = (cur_trans-up_trans8(previous_lfp_up))/Fsd;
    end

    [dummy,nearest_lfp_down] = min(abs(down_trans8-cur_trans));
    if ~isempty(nearest_lfp_down)
        bs_ttnd(i) = (cur_trans-down_trans8(nearest_lfp_down))/Fsd;
    end

    previous_lfp_down = find(down_trans8 < cur_trans,1,'last');
    if ~isempty(previous_lfp_down)
        bs_tsld(i) = (cur_trans - down_trans8(previous_lfp_down))/Fsd;
    end

end