function peak_locs = find_gram_peaks(gram,thresh_value,f_range)

temp = [zeros(size(gram,1),1) diff(gram,[],2)];
dgram = zeros(size(gram));
dgram(:,f_range) = temp(:,f_range);
peak_locs = find(dgram(:,1:end-1) > 0 & dgram(:,2:end) < 0);
peak_locs(gram(peak_locs) < thresh_value) = [];