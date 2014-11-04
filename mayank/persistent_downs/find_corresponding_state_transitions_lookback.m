function [corresp_lfp_upinds,corresp_lfp_downinds] = find_corresponding_state_transitions_lookback(...
    mp_uptrans,mp_downtrans,lfp_uptrans,lfp_downtrans)

num_mp_ups = length(mp_uptrans);
num_lfp_ups = length(lfp_uptrans);
N = min(num_mp_ups,num_lfp_ups);

%initialize vectors for storing assigned state transitions
corresp_lfp_upinds = nan(num_mp_ups,1);
corresp_lfp_downinds = nan(num_mp_ups,1);
corresp_lfp_dist = nan(num_mp_ups,1);

%make a greedy assignment of all MP up transitions to one LFP up-transition
available_lfp_upinds = 1:length(lfp_uptrans); %vector of unassigned LFP up transition indices
remaining_mp_ups = 1:length(mp_uptrans); %vector of unassigned MP up transition indices
for i = 1:N
    %initialize vectors to store the minimized LFP up-transition distance and index value for each MP up transition 
    mp_updists = nan(size(remaining_mp_ups));
    corr_lfp_upinds = nan(size(remaining_mp_ups));
    for j = 1:length(remaining_mp_ups)
        %find nearest lfp up transition for each mp up transition
        [mp_updists(j),corr_lfp_upinds(j)] = min(abs(mp_uptrans(remaining_mp_ups(j)) - lfp_uptrans(available_lfp_upinds)));
    end
    [~,winner_loc] = min(mp_updists); %index of the overall closest MP-LFP pair
    corresp_lfp_upinds(remaining_mp_ups(winner_loc)) = available_lfp_upinds(corr_lfp_upinds(winner_loc));
    corresp_lfp_dist(remaining_mp_ups(winner_loc)) = mp_updists(winner_loc);
    %remove assigned LFP and MP up transitions
    remaining_mp_ups(winner_loc) = [];
    available_lfp_upinds(corr_lfp_upinds(winner_loc)) = [];
end

%now break crossed links
for i = 1:num_mp_ups-1 %cycle through each MP up transition (up to the second to last)
    %find any corresponding LFP up transitions for later MP up transitions
    %that get assigned before the current corresponding LFP up 
    temp = i+find(corresp_lfp_upinds(i+1:end) <= corresp_lfp_upinds(i));  
    if ~isempty(temp)
        %find the shortest pair distance of the overlapping set
        [min_new_dist,min_new_loc] = min(corresp_lfp_dist(temp));   
        %if the current pair distance is shorter eliminate all links of the
        %proposed set
        if corresp_lfp_dist(i) < min_new_dist
            corresp_lfp_upinds(temp) = nan;
            corresp_lfp_dist(temp) = nan;
        %otherwise get rid of the current linkage and continue
        else
            corresp_lfp_upinds(i) = nan;
            corresp_lfp_dist(i) = nan;
        end       
    end
end

%% For all MP down transitions
corresp_lfp_dist = nan(num_mp_ups,1);
assigned_mp_ups = find(~isnan(corresp_lfp_upinds)); %the set of MP up transitions which had corresponding LFP ups
for i = 1:length(assigned_mp_ups)
%     cur_mp_down_set = assigned_mp_ups(i):(assigned_mp_ups(i+1)-1);
    %the potential corresponding LFP downs are the ones that are within the
    %corresponding range of LFP up transitions
    if i < length(assigned_mp_ups)
        cur_lfp_down_set = find(lfp_downtrans > lfp_uptrans(corresp_lfp_upinds(assigned_mp_ups(i))) & ...
            lfp_downtrans < lfp_uptrans(corresp_lfp_upinds(assigned_mp_ups(i+1))));
    else
        cur_lfp_down_set = find(lfp_downtrans > lfp_uptrans(corresp_lfp_upinds(assigned_mp_ups(i))));
    end
    
    %these are the LFP down-transitions which occur after the corresponding
    %MP down-transition
    after_dtrans = find(lfp_downtrans(cur_lfp_down_set) > mp_downtrans(assigned_mp_ups(i)));
    
    if length(cur_lfp_down_set) > 1
        %if all LFP down trans are after the MP downtrans, use the first
        %one
        if length(after_dtrans) == length(cur_lfp_down_set)
            cur_lfp_down_set = cur_lfp_down_set(1);
        else %otherwise, use the nearest one occuring after the MP down trans
            cur_lfp_down_set(lfp_downtrans(cur_lfp_down_set) > mp_downtrans(assigned_mp_ups(i))) = [];
        end
    end
    if isempty(cur_lfp_down_set) %if there are no potential corresponding LFP downs
        corresp_lfp_downinds(assigned_mp_ups(i)) = nan;
    else
        %find nearest lfp up transition for each mp up transition
        [mp_downdists,corr_lfp_downinds] = min(abs(mp_downtrans(assigned_mp_ups(i)) - lfp_downtrans(cur_lfp_down_set)));
        corresp_lfp_downinds(assigned_mp_ups(i)) = cur_lfp_down_set(corr_lfp_downinds);
        corresp_lfp_dist(assigned_mp_ups(i)) = mp_downdists;
    end
end

%now break crossed links
for i = 1:num_mp_ups-1
    temp = i+find(corresp_lfp_downinds(i+1:end) <= corresp_lfp_downinds(i));
    if ~isempty(temp)
        [min_new_dist,min_new_loc] = min(corresp_lfp_dist(temp));     
        if corresp_lfp_dist(i) < min_new_dist
            corresp_lfp_downinds(temp) = nan;
            corresp_lfp_dist(temp) = nan;
        else
            corresp_lfp_downinds(i) = nan;
            corresp_lfp_dist(i) = nan;
        end       
    end
end

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
