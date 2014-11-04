clear all
close all

load C:/WC_germany/final_pdown_analysis/compiled_data.mat

%%
for dd = 1:length(data)
    cd(data(dd).dir)
    pwd
    
    if exist('./mua_data3.mat','file')
        load('./mua_data3.mat');
        load ./sync_times.mat
        
        avg_rates(dd,:) = cellfun(@(x) length(x),mua_amps)/range(synct)*1e6;
        
        %find peaks and valleys of avg spike waveforms for each channel
        [wave_peaks(dd,:),peak_locs] = max(avg_waveform(:,4:20),[],2);
        peak_locs = peak_locs + 3;
        [~,min_loc] = min(avg_waveform(:,6:end),[],2);
        min_loc = min_loc + 5;
        wave_spkwidths(dd,:) = min_loc - peak_locs;

        %compute different measures of spike width
        avg_spkwidths(dd,:) = nan(1,8);
        for ii = 1:8
            avg_spkwidths(dd,ii) = mean(mua_widths{ii});
            dp = find(avg_waveform(ii,peak_locs(ii):end) < wave_peaks(dd,ii)/2,1,'first');
            sp = find(avg_waveform(ii,1:peak_locs(ii)) < wave_peaks(dd,ii)/2,1,'last');
            if ~isempty(dp) && ~isempty(sp)
                avg_spkwidth_fwhm(dd,ii) = peak_locs(ii)+dp-sp;
            else
                avg_spkwidth_fwhm(dd,ii) = nan;
            end
            dp = find(avg_waveform(ii,peak_locs(ii):end) < wave_peaks(dd,ii)/4,1,'first');
            sp = find(avg_waveform(ii,1:peak_locs(ii)) < wave_peaks(dd,ii)/4,1,'last');
            if ~isempty(dp) && ~isempty(sp)
                avg_spkwidth_fwqm(dd,ii) = peak_locs(ii)+dp-sp;
            else
                avg_spkwidth_fwqm(dd,ii) = nan;
            end
        end
        avg_waveforms(dd,:,:) = avg_waveform;
                
        %find max of rate vs depth profile
        temp = avg_rates(dd,:); temp(1) = 0;
        [peak_hpcmua_rate(dd),peak_hpcmua_loc(dd)] = max(temp(3:4));
        peak_hpcmua_loc(dd) = peak_hpcmua_loc(dd) + 2;
        
        %make sure the max is actually a peak
        if avg_rates(dd,peak_hpcmua_loc(dd)+1) >= peak_hpcmua_rate(dd)
            peak_hpcmua_loc(dd) = nan;
        elseif avg_rates(dd,peak_hpcmua_loc(dd)-1) >= peak_hpcmua_rate(dd)
            peak_hpcmua_loc(dd) = nan;
        end
        
    else
        
        peak_hpcmua_loc(dd) = nan;
        peak_hpcmua_rate(dd) = nan;
        avg_rates(dd,:) = nan(1,8);
        avg_spkwidths(dd,:) = nan(1,8);
        peak_mua_rate(dd) = nan;
        peak_mua_loc(dd) = nan;
        avg_waveforms(dd,:,:) = nan(8,32);
        avg_spkwidth_fwhm(dd,:) = nan(1,8);
        avg_spkwidth_fwqm(dd,:) = nan(1,8);
    end
end

%%
for ii = 67:size(avg_rates,1)
    plot(2:7,avg_rates(ii,2:7),'o-')
    ii
    pause
    clf
end
%%
bad_peaks = [25 26 27 29 30 67 68 78]; %these peaks in the rate vs depth profile are not clear
% bad_peaks = [25 26 27 29 30]; %these peaks in the rate vs depth profile are not clear
peak_hpcmua_loc(bad_peaks) = nan;

%usable MUA on all channels
usable_mua = (avg_spkwidths >= 10 & avg_spkwidth_fwhm <= 10 & avg_spkwidth_fwqm <= 11);
usable_mua(98,8) = 0;

%%
cd C:\WC_Germany\final_pdown_analysis\
save mua_classification2 usable_mua avg_rates peak_hpc* avg_rates avg_*