function [corresp_lfp_upinds,corresp_lfp_downinds] = greedy_find_corresponding_ncx_state_transitions(...
    mp_uptrans,mp_downtrans,lfp_uptrans,lfp_downtrans)

num_mp_ups = length(mp_uptrans);
num_lfp_ups = length(lfp_uptrans);
N = min(num_mp_ups,num_lfp_ups);

corresp_lfp_upinds = nan(num_mp_ups,1);
corresp_lfp_downinds = nan(num_mp_ups,1);
corresp_lfp_dist = nan(num_mp_ups,1);

%make a greedy assignment of all MP up transitions to one LFP up-transition
available_lfp_upinds = 1:length(lfp_uptrans);
remaining_mp_ups = 1:length(mp_uptrans);
for i = 1:N
    mp_updists = nan(size(remaining_mp_ups));
    corr_lfp_upinds = nan(size(remaining_mp_ups));
    for j = 1:length(remaining_mp_ups)
        %find nearest lfp up transition for each mp up transition
        [mp_updists(j),corr_lfp_upinds(j)] = min(abs(mp_uptrans(remaining_mp_ups(j)) - lfp_uptrans(available_lfp_upinds)));
        %         [~,corresp_lfp_upinds(i)] =  min(abs(mp_uptrans(i) - lfp_uptrans));
    end
    [~,winner_loc] = min(mp_updists);
    %    if winner_loc == 1
    %        1
    %    end
    corresp_lfp_upinds(remaining_mp_ups(winner_loc)) = available_lfp_upinds(corr_lfp_upinds(winner_loc));
    corresp_lfp_dist(remaining_mp_ups(winner_loc)) = mp_updists(winner_loc);
    remaining_mp_ups(winner_loc) = [];
    available_lfp_upinds(corr_lfp_upinds(winner_loc)) = [];
    %     corresp_lfp_upinds(1)
    %     pause
end

%now break crossed links
for i = 1:num_mp_ups-1
    temp = i+find(corresp_lfp_upinds(i+1:end) <= corresp_lfp_upinds(i));
    if ~isempty(temp)
        [min_new_dist,min_new_loc] = min(corresp_lfp_dist(temp));     
        if corresp_lfp_dist(i) < min_new_dist
            corresp_lfp_upinds(temp) = nan;
            corresp_lfp_dist(temp) = nan;
        else
            corresp_lfp_upinds(i) = nan;
            corresp_lfp_dist(i) = nan;
        end       
    end
end


all_lfp_trans = [lfp_uptrans ones(size(lfp_uptrans))];
all_lfp_trans = [all_lfp_trans; lfp_downtrans 2*ones(size(lfp_downtrans))];
[~,order] = sort(all_lfp_trans(:,1));
all_lfp_trans = all_lfp_trans(order,:);

corresp_lfp_dist = nan(num_mp_ups,1);
for i = 1:num_mp_ups
    %if the MP down transition occurs before the subsequent LFP down
    %transition, assign the next LFP down
    if ~isnan(corresp_lfp_upinds(i))
        if mp_downtrans(i) <= lfp_downtrans(corresp_lfp_upinds(i))
            corresp_lfp_downinds(i) = corresp_lfp_upinds(i);
        else
            %otherwise, check if the down transition happened during a cortical up state.
            prev_trans_ind = find(all_lfp_trans(:,1) < mp_downtrans(i),1,'last');
            
            %IF so
            if all_lfp_trans(prev_trans_ind,2) == 1
                %check whether there is another MP up-transition assigned
                %to this supsequent LFP up state
                or_loc = find(lfp_uptrans == all_lfp_trans(prev_trans_ind,1));
                temp = find(corresp_lfp_upinds == or_loc);
                %if there is then assign the current MP down-transition to
                %the previous LFP down transition
                if ~isempty(temp)
                    corresp_lfp_downinds(i) = find(lfp_downtrans < mp_downtrans(i),1,'last');
                    %otherwise, assign the current MP down-transition to the
                    %next LFP down transition
                else
                    corresp_lfp_downinds(i) = find(lfp_downtrans >= mp_downtrans(i),1,'first');
                end
            else
                %If not, assign the previous LFP down
                corresp_lfp_downinds(i) = find(lfp_downtrans < mp_downtrans(i),1,'last');
            end
        end
        corresp_lfp_dist(i) = abs(mp_downtrans(i) - lfp_downtrans(corresp_lfp_downinds(i)));
    end
end

%now break crossed links
for i = 1:num_mp_ups-1
    if corresp_lfp_downinds(i) >= corresp_lfp_downinds(i+1)
        if corresp_lfp_dist(i) > corresp_lfp_dist(i+1)
            corresp_lfp_downinds(i) = nan;
            corresp_lfp_dist(i) = nan;
        else
            corresp_lfp_downinds(i+1) = nan;
            corresp_lfp_upinds(i+1) = nan;
            corresp_lfp_dist(i+1) = nan;
        end
    end
end

%now scan and check for skipped MP down states
for i = 1:num_mp_ups-1
    %if the corresponding LFP down transition occurs after the next MP
    %up-transition
    if ~isnan(corresp_lfp_downinds(i))
        if lfp_downtrans(corresp_lfp_downinds(i)) >= mp_uptrans(i+1)
            [shortest_dist,shortest_loc] = min(abs(mp_downtrans-lfp_downtrans(corresp_lfp_downinds(i))));
            %if the current MP down transition is not best for the
            %corresponding LFP down transition, then consider reassignment
            if shortest_loc ~= i
                if isnan(corresp_lfp_downinds(shortest_loc))
                    corresp_lfp_downinds(shortest_loc) = corresp_lfp_downinds(i);
                    corresp_lfp_downinds(i) = nan;
                end
            end
        end
    end
end

% %if either the up or down transition has no correspnding LFP transition
% %then set them both to nans
% for i = 1:num_mp_ups
%     if isnan(corresp_lfp_upinds(i)) || isnan(corresp_lfp_downinds(i))
%         corresp_lfp_upinds(i) = nan;
%         corresp_lfp_downinds(i) = nan;
%     end
% end

%final check for any crossed links
crossed_ups = [];
crossed_downs = [];
for i = 1:num_mp_ups-1
    temp = i+find(corresp_lfp_upinds(i+1:end) <= corresp_lfp_upinds(i));
    crossed_ups = [crossed_ups; temp];
    temp = i+find(corresp_lfp_downinds(i+1:end) <= corresp_lfp_downinds(i));
    crossed_downs = [crossed_downs; temp];
end
fprintf('%d crossed UPs.  %d crossed DOWNs.\n',length(crossed_ups),length(crossed_downs));
