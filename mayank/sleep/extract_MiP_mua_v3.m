function [mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua_v3(Fname,amp_threshold)

ExtractHeader = 1;
ExtractMode = 1;
if nargin < 3
    skip_ch1 = 1;
end

minFsize = 1e6;
D = dir(Fname);
if D.bytes < minFsize
   mua_times = [];
   mua_amps = [];
   mua_widths = [];
   avg_waveform = nan;
   std_waveform = nan;
    
    
    return
end

chunk_size = 50e3;
FieldSelection = [1 0 0 0 0];
[TimeStamps, Header] = Nlx2MatSpike_v3(Fname, FieldSelection, ExtractHeader, ExtractMode);
tot_n_spks = length(TimeStamps);
n_chunks = ceil(tot_n_spks/chunk_size);

FieldSelection = [1 0 0 0 1];
avg_waveform  = [];
std_waveform = [];
mua_widths = [];
mua_times = [];
mua_amps = [];
for nn = 1:n_chunks
    cur_range = [(nn-1)*chunk_size+1 nn*chunk_size];
    cur_range(cur_range > tot_n_spks) = tot_n_spks;
    [TimeStamps, DataPoints, Header] = Nlx2MatSpike_v3(Fname, FieldSelection, ExtractHeader, 2,cur_range);
    
    %extract bit conversion factors from header
    if strcmp(Header{16}(1:12),'-ADBitVolts ')
        bitconv = Header{16};
    else
        error('Header format');
    end
    conv_factors = bitconv(13:end);
    conv_factors = str2num(conv_factors);
    [n_pts,n_ch,n_spks] = size(DataPoints);
    DataPoints = DataPoints.*repmat(conv_factors,[n_pts,1,n_spks])*1e6; %in uV
    [peaks,peak_locs] = max(DataPoints(4:20,:,:),[],1);
    peaks = squeeze(peaks); peak_locs = squeeze(peak_locs);
    
    cur_spikes = find(peaks > amp_threshold);
    cur_spike_wave = squeeze(DataPoints(:,1,cur_spikes));
    cur_avg_waveform = mean(cur_spike_wave,2);
    cur_std_waveform = std(cur_spike_wave,[],2);
    
    [~,min_loc] = min(cur_spike_wave(6:end,:));
    min_loc = min_loc + 5;
    
    mua_widths = [mua_widths; min_loc' - peak_locs(cur_spikes)];
    mua_times = [mua_times; TimeStamps(cur_spikes)'];
    mua_amps = [mua_amps; peaks(cur_spikes)];
    
    avg_waveform = cat(3,avg_waveform,cur_avg_waveform);
    std_waveform = cat(3,std_waveform,cur_std_waveform);
end
avg_waveform = squeeze(mean(avg_waveform,3));
std_waveform = squeeze(mean(std_waveform,3));
