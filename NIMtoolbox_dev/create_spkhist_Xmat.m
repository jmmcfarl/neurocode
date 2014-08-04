function Xspkhst = create_spkhist_Xmat(Robs,bin_edges)
%
% Xspkhst = create_spkhist_Xmat(Robs,bin_edges)
%
% Creates an X-matrix out of observed spike train Robs and specified 'bin edges'

%%
NT = length(Robs);
maxlag = max(bin_edges);
spkbns = convert_to_spikebins(Robs); %spike bins

%%
Tmat = zeros(NT + maxlag,length(bin_edges)-1);
for i = 1:length(spkbns)
    for j = 1:length(bin_edges)-1
        Tmat(spkbns(i)+(bin_edges(j):(bin_edges(j+1)-1)),j) = Tmat(spkbns(i)+(bin_edges(j):(bin_edges(j+1)-1)),j) + 1;
    end
end

Xspkhst =  Tmat(1:NT,:); %concatenate onto X_matrix

