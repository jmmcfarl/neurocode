function [mp_upskipped,mp_downskipped] = greedy_find_skipped_ncx_states(...
    corresp_lf8_upinds,corresp_lf8_downinds,lf8_updurs,lf8_downdurs,thresh_down_dur,thresh_up_dur)

%store the LFP transition indices, number of skipped LFP states, robust
%number of skipped states, and the minimum duration of skipped states for
%each MP state
mp_upskipped.inds = cell(length(corresp_lf8_upinds),1);
mp_downskipped.inds = cell(length(corresp_lf8_upinds),1);
mp_upskipped.min_dur = nan(length(corresp_lf8_upinds),1);
mp_downskipped.min_dur = nan(length(corresp_lf8_upinds),1);
mp_upskipped.num_skipped = nan(length(corresp_lf8_upinds),1);
mp_downskipped.num_skipped = nan(length(corresp_lf8_upinds),1);
mp_upskipped.rnum_skipped = nan(length(corresp_lf8_upinds),1);
mp_downskipped.rnum_skipped = nan(length(corresp_lf8_upinds),1);

for i = 1:length(corresp_lf8_upinds)
    if ~isnan(corresp_lf8_upinds(i)) %if the mp up state was not itself skipped
        %for up state
        cur_skipped_downs = corresp_lf8_upinds(i):(corresp_lf8_downinds(i)-1);
        if isnan(cur_skipped_downs) %if the next down state was skipped
            mp_upskipped.inds{i} = [];
            mp_upskipped.num_skipped(i) = 0;
            mp_upskipped.rnum_skipped(i) = 0;
        else
            mp_upskipped.inds{i} = cur_skipped_downs;
            mp_upskipped.num_skipped(i) = length(cur_skipped_downs);
            mp_upskipped.rnum_skipped(i) = length(cur_skipped_downs);
            if ~isempty(cur_skipped_downs)
                mp_upskipped.min_dur(i) = min(lf8_downdurs(cur_skipped_downs));
                mp_upskipped.rnum_skipped(i) = mp_upskipped.rnum_skipped(i) - sum(lf8_downdurs(cur_skipped_downs) < thresh_down_dur);
            end
        end
    end
    if i < length(corresp_lf8_upinds) && ~isnan(corresp_lf8_downinds(i)) %if the MP down state was not itself skipped
        %for down state
        cur_skipped_ups = (corresp_lf8_downinds(i)+1):(corresp_lf8_upinds(i+1)-1);
        if isnan(cur_skipped_ups) %if the next MP up state was skipped
            mp_downskipped.inds{i} = [];
            mp_downskipped.num_skipped(i) = 0;
            mp_downskipped.rnum_skipped(i) = 0;
        else
            mp_downskipped.inds{i} = cur_skipped_ups;
            mp_downskipped.num_skipped(i) = length(cur_skipped_ups);
            mp_downskipped.rnum_skipped(i) = length(cur_skipped_ups);
            if ~isempty(cur_skipped_ups)
                mp_downskipped.min_dur(i) = min(lf8_updurs(cur_skipped_ups));
                mp_downskipped.rnum_skipped(i) = mp_downskipped.rnum_skipped(i) - sum(lf8_updurs(cur_skipped_ups) < thresh_up_dur);
            end
        end
    end
end

