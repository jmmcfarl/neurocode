function new_X = tb_proc_stim(X,tent_space,flen)
%new_X = tb_proc_stim(X,tent_space,flen)
% applies SPATIAL tent-basis processing to the stimulus
% INPUTS:
%     X: time-embedded stimulus matrix [N x K] where K is total number of stim dims
%     tent_space: spacing (in pixels) of tent-bases
%     flen: number of stimulus time lags
% OUTPUTS:
%     new_X: stimulus projected onto the spatial tent-bases
    
if nargin < 3 || isempty(flen)
    flen = 1;
end

%create a tent-basis (triangle) filter
tent_filter = [(1:tent_space)/tent_space 1-(1:tent_space-1)/tent_space]/tent_space;

nPix = size(X,2)/flen;
X = reshape(X,size(X,1),flen,nPix);

%apply to the stimulus
new_X = zeros(size(X));
for i = 1:length(tent_filter)
    new_X = new_X + shift_matrix_Nd(X,i-tent_space,3)*tent_filter(i);
end
new_X = new_X(:,:,1:tent_space:end); %no apply spatial downsampling
new_X = reshape(new_X,size(new_X,1),[]);
