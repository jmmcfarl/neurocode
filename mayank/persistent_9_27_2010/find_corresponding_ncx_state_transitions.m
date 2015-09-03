function [corresp_lfp_upinds,corresp_lfp_downinds] = find_corresponding_ncx_state_transitions(...
    mp_uptrans,mp_upsegnums,mp_downtrans,mp_downsegnums,lfp_uptrans,lfp_upsegnums,lfp_downtrans,lfp_downsegnums)

num_mp_ups = length(mp_uptrans);

corresp_lfp_upinds = nan(num_mp_ups,1);
corresp_lfp_downinds = nan(num_mp_ups,1);
for i = 1:num_mp_ups
    %find nearest lfp up transition
   [~,corresp_lfp_upinds(i)] =  min(abs(mp_uptrans(i) - lfp_uptrans));
   
%    prev_lfp_up = find(lfp_uptrans < mp_uptrans(i),1,'last');
%    if mp_uptrans(i) < lfp_downtrans(prev_lfp_up)
%       corresp_lfp_upinds(i) = prev_lfp_up; 
%    end
   
   %if the MP down transition occurs before the subsequent LFP down
   %transition, assign the next LFP down
   if mp_downtrans(i) <= lfp_downtrans(corresp_lfp_upinds(i))
       corresp_lfp_downinds(i) = corresp_lfp_upinds(i);
   else
       %otherwise, assign the previous LFP down
       corresp_lfp_downinds(i) = find(lfp_downtrans < mp_downtrans(i),1,'last');
   end  
end

if any(mp_upsegnums ~= lfp_upsegnums(corresp_lfp_upinds))
    disp('error: Up trans from different UDS segments aligned!')
end
if any(mp_downsegnums ~= lfp_downsegnums(corresp_lfp_downinds))
    disp('error: Down trans from different UDS segments aligned~')
end