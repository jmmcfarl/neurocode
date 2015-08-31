function [corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts,all_eye_vals,all_eye_ts,all_t_axis,all_blockvec,block_bar_oris,used_inds,lin_correction)
% [corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts,all_eye_vals,all_eye_ts,all_t_axis,all_blockvec,block_bar_oris,used_inds,lin_correction)
% do basic preprocessing on ET data and align to bar oriented coordinates
% INPUTS:
%     Expts: Expts struct
%     all_eye_vals: TxN matrix of ET samples
%     all_eye_ts: ET time samples
%     all_t_axis: expt time axis
%     all_blockvec: corresponding block #s of each expt time sample
%     block_bar_oris: bar orientation in each block
%     used_inds: set of time indices used for analysis
%     <lin_correction>: do we want to remove linear relationship between parallel and orth eye position dims? [default is true]
% OUTPUTS:
%     corrected_eye_vals: processed eye position data (at ET resolution)
%     corrected_eye_vals_interp: processed eye position data, interpolated onto expt time axis
% 
if nargin < 8 
    lin_correction = true;
end

fprintf('Computing corrected ET data\n');

interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_t_axis); %interpolate raw ET data onto expt time axis
interp_eye_blocks = round(interp1(all_t_axis,all_blockvec,all_eye_ts)); %find block numbers of each ET sample

corrected_eye_vals = all_eye_vals;
corrected_eye_vals_interp = interp_eye_vals;

%subtract off fixation point location
expt_fx = cellfun(@(x) x.Stimvals.fx,Expts,'UniformOutput',true);
expt_fy = cellfun(@(x) x.Stimvals.fy,Expts,'UniformOutput',true);
expt_fx = unique(expt_fx);
expt_fy = unique(expt_fy);
if length(expt_fx) > 1 || length(expt_fy) > 1
    warning('RF position not constant across blocks!\n');
    expt_fx = mode(expt_fx);
    expt_fy = mode(expt_fy);
end
corrected_eye_vals_interp(:,[1 3]) = corrected_eye_vals_interp(:,[1 3]) - expt_fx;
corrected_eye_vals_interp(:,[2 4]) = corrected_eye_vals_interp(:,[2 4]) - expt_fy;
corrected_eye_vals(:,[1 3]) = corrected_eye_vals(:,[1 3]) - expt_fx;
corrected_eye_vals(:,[2 4]) = corrected_eye_vals(:,[2 4]) - expt_fy;

block_bar_oris = block_bar_oris/180*pi; %block bar oris in rads
unique_bar_oris = unique(block_bar_oris);

%cycle over unique bar orientations
for bb = 1:length(unique_bar_oris)
   cur_blocks = find(block_bar_oris == unique_bar_oris(bb)); 
   cur_inds = used_inds(ismember(all_blockvec(used_inds),cur_blocks));
   cur_eye_inds = ismember(interp_eye_blocks,cur_blocks);
   
   %rotate data into bar-oriented coordinate system
   Rmat = [cos(unique_bar_oris(bb)) sin(unique_bar_oris(bb)); -sin(unique_bar_oris(bb)) cos(unique_bar_oris(bb))];
   corrected_eye_vals_interp(cur_inds,[1 2]) = corrected_eye_vals_interp(cur_inds,[1 2])*Rmat;
   corrected_eye_vals_interp(cur_inds,[3 4]) = corrected_eye_vals_interp(cur_inds,[3 4])*Rmat;
   corrected_eye_vals(cur_eye_inds,[1 2]) = corrected_eye_vals(cur_eye_inds,[1 2])*Rmat;
   corrected_eye_vals(cur_eye_inds,[3 4]) = corrected_eye_vals(cur_eye_inds,[3 4])*Rmat;
  
   %subtract off median in each dimension
   corrected_eye_vals_interp(cur_inds,:) = bsxfun(@minus,corrected_eye_vals_interp(cur_inds,:),nanmedian(corrected_eye_vals_interp(cur_inds,:)));
   corrected_eye_vals(cur_eye_inds,:) = bsxfun(@minus,corrected_eye_vals(cur_eye_inds,:),nanmedian(corrected_eye_vals_interp(cur_inds,:)));

   if lin_correction
       %remove linear relationship between parallel and orthoganol eye
       %positions
       b = robustfit(corrected_eye_vals_interp(cur_inds,1),corrected_eye_vals_interp(cur_inds,2));
       pred_eyevals = corrected_eye_vals(used_inds,1)*b(2) + b(1);
       corrected_eye_vals_interp(used_inds,2) = corrected_eye_vals_interp(used_inds,2) - pred_eyevals;
       pred_eyevals = corrected_eye_vals(cur_eye_inds,1)*b(2) + b(1);
       corrected_eye_vals(cur_eye_inds,2) = corrected_eye_vals(cur_eye_inds,2) - pred_eyevals;
       
       b = robustfit(corrected_eye_vals_interp(cur_inds,3),corrected_eye_vals_interp(cur_inds,4));
       pred_eyevals = corrected_eye_vals_interp(used_inds,3)*b(2) + b(1);
       corrected_eye_vals_interp(used_inds,4) = corrected_eye_vals_interp(used_inds,4) - pred_eyevals;
       pred_eyevals = corrected_eye_vals(cur_eye_inds,3)*b(2) + b(1);
       corrected_eye_vals(cur_eye_inds,4) = corrected_eye_vals(cur_eye_inds,4) - pred_eyevals;
       
       %subtract off median again
       corrected_eye_vals_interp(cur_inds,:) = bsxfun(@minus,corrected_eye_vals_interp(cur_inds,:),nanmedian(corrected_eye_vals_interp(cur_inds,:)));
       corrected_eye_vals(cur_eye_inds,:) = bsxfun(@minus,corrected_eye_vals(cur_eye_inds,:),nanmedian(corrected_eye_vals_interp(cur_inds,:)));
   end
end

