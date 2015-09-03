function [mp_upskipped,mp_downskipped] = find_skipped_ncx_states(mp_upinds,...
    mp_downinds,lf8_upinds,lf8_downinds,corresp_lf8_upinds,corresp_lf8_downinds,...
    lf8_uplags,lf8_downlags,lf8_updurs,lf8_downdurs,thresh_updur,thresh_downdur)

mp_upskipped.lf8ups = cell(length(mp_upinds),1);
mp_downskipped.lf8ups = cell(length(mp_upinds),1);
mp_upskipped.lf8downs = cell(length(mp_upinds),1);
mp_downskipped.lf8downs = cell(length(mp_upinds),1);
mp_upskipped.num_lf8ups = nan(length(mp_upinds),1);
mp_upskipped.num_lf8downs = nan(length(mp_upinds),1);
mp_downskipped.num_lf8ups = nan(length(mp_upinds),1);
mp_downskipped.num_lf8downs = nan(length(mp_upinds),1);
mp_upskipped.min_updur = nan(length(mp_upinds),1);
mp_upskipped.min_downdur = nan(length(mp_upinds),1);
mp_downskipped.min_updur = nan(length(mp_upinds),1);
mp_downskipped.min_downdur = nan(length(mp_upinds),1);

for i = 1:length(mp_upinds)
   
    %for up state
    cur_skipped_ups = find(lf8_downinds > lf8_upinds(corresp_lf8_upinds(i)) & ...
        lf8_downinds < mp_downinds(i) - lf8_uplags(i));
    cur_skipped_downs = find(lf8_upinds > lf8_upinds(corresp_lf8_upinds(i)) & ...
        lf8_upinds < mp_downinds(i) - lf8_uplags(i));
    if ~isempty(cur_skipped_downs)
        cur_skipped_downs = cur_skipped_downs-1;
    end
    mp_upskipped.lf8ups{i} = cur_skipped_ups;
    mp_upskipped.lf8downs{i} = cur_skipped_downs;
    mp_upskipped.num_lf8ups(i) = sum(lf8_updurs(mp_upskipped.lf8ups{i}) > thresh_updur);
    mp_upskipped.num_lf8downs(i) = sum(lf8_downdurs(mp_upskipped.lf8downs{i}) > thresh_downdur);
    if ~isempty(mp_upskipped.lf8ups{i})
        mp_upskipped.min_updur(i) = min(lf8_updurs(mp_upskipped.lf8ups{i}));
    end
    if ~isempty(mp_upskipped.lf8downs{i})
        mp_upskipped.min_downdur(i) = min(lf8_downdurs(mp_upskipped.lf8downs{i}));
    end
    
    if i < length(mp_upinds)
        %for down state
        cur_skipped_downs = find(lf8_upinds > lf8_downinds(corresp_lf8_downinds(i)) & ...
            lf8_upinds < mp_upinds(i+1) - lf8_downlags(i));
        cur_skipped_ups = find(lf8_downinds > lf8_downinds(corresp_lf8_downinds(i)) & ...
            lf8_downinds < mp_upinds(i+1) - lf8_downlags(i));
        if ~isempty(cur_skipped_downs)
            cur_skipped_downs = cur_skipped_downs - 1;
        end
        mp_downskipped.lf8ups{i} = cur_skipped_ups;
        mp_downskipped.lf8downs{i} = cur_skipped_downs;
        mp_downskipped.num_lf8ups(i) = sum(lf8_updurs(mp_downskipped.lf8ups{i}) > thresh_updur);
        mp_downskipped.num_lf8downs(i) = sum(lf8_downdurs(mp_downskipped.lf8downs{i}) > thresh_downdur);
        if ~isempty(mp_downskipped.lf8ups{i})
            mp_downskipped.min_updur(i) = min(lf8_updurs(mp_downskipped.lf8ups{i}));
        end
        if ~isempty(mp_downskipped.lf8downs{i})
            mp_downskipped.min_downdur(i) = min(lf8_downdurs(mp_downskipped.lf8downs{i}));
        end     
    end
end

