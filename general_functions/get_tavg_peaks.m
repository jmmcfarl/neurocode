function [peak_vals,peak_locs] = get_tavg_peaks(tavg_data,lags,search_range)

poss_lagrange = find(lags >= search_range(1) & lags <= search_range(2));
[peak_vals,peak_locs] = deal(nan(size(tavg_data,1),1));
for ii = 1:size(tavg_data,1)
    if all(~isnan(tavg_data(ii,poss_lagrange)))
        [temp,temploc] = findpeaks(tavg_data(ii,poss_lagrange),'sortstr','descend');
        if ~isempty(temp)
            peak_vals(ii) = temp(1); 
            peak_locs(ii) = temploc(1);
        end
    end
end
peak_locs(~isnan(peak_locs)) = lags(poss_lagrange(peak_locs(~isnan(peak_locs))));
peak_vals(isinf(peak_vals)) = nan;