function stim_up = up_sample_stimulus(stim,usfac)
%stim_up = up_sample_stimulus(stim,usfac)
%takes a one-d stimulus (not time-embedded), and applies spatial up-sampling, assuming stimulus is
%constant within the interpolated regions.
if usfac > 1
    stim_up = zeros(size(stim,1),size(stim,2)*usfac);
    for ii = 1:size(stim,2)
        for jj = 1:usfac
            stim_up(:,usfac*(ii-1)+jj) = stim(:,ii);
        end
    end
elseif usfac == 1
    stim_up = stim;
else
   error('invalid usfac'); 
end
