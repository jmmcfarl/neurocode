function [mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua_v2(amp_threshold,max_overlap,skip_ch1)

ExtractHeader = 1;
ExtractMode = 1;
if nargin < 3
    skip_ch1 = 1;
end

% FieldSelection = [1 0 0 0 0];
% Filename = 'Sc1.ntt';
% [TimeStamps1] = Nlx2MatSpike(Filename, FieldSelection, 0, ExtractMode);
% Filename = 'Sc2.ntt';
% [TimeStamps2] = Nlx2MatSpike(Filename, FieldSelection, 0, ExtractMode);
% [bad_set,bad_inds1,bad_inds2] = intersect(TimeStamps1,TimeStamps2);

chunk_size = 50e3;
FieldSelection = [1 0 0 0 0];
Filename = 'Sc1.ntt';
[TimeStamps, Header] = Nlx2MatSpike_v3(Filename, FieldSelection, ExtractHeader, ExtractMode);
tot_n_spks = length(TimeStamps);
n_chunks = ceil(tot_n_spks/chunk_size);

FieldSelection = [1 0 0 0 1];
sc1_avg_waveform  = [];
sc1_std_waveform = [];
for i = 1:4
    mua_widths{i} = [];
    mua_times{i} = [];
    mua_amps{i} = [];
end
for nn = 1:n_chunks
    cur_range = [(nn-1)*chunk_size+1 nn*chunk_size];
    cur_range(cur_range > tot_n_spks) = tot_n_spks;
    [TimeStamps, DataPoints, Header] = Nlx2MatSpike_v3(Filename, FieldSelection, ExtractHeader, 2,cur_range);
    
    %locate and remove spike events where the spike waveforms are too similar
    %across channels
    if skip_ch1
        spk_corrs = zeros(3,length(TimeStamps));
        cnt = 1;
        for aa = 2:4
            cur_set_a = squeeze(DataPoints(:,aa,:));
            cur_set_a = bsxfun(@minus,cur_set_a,mean(cur_set_a));
            for bb = aa+1:4
                cur_set_b = squeeze(DataPoints(:,bb,:));
                cur_set_b = bsxfun(@minus,cur_set_b,mean(cur_set_b));
                
                shared_var = sum(cur_set_a.*cur_set_b);
                spk_corrs(cnt,:) = shared_var./sqrt(sum(cur_set_a.^2).*sum(cur_set_b.^2));
                
                cnt = cnt + 1;
            end
        end
    else
        spk_corrs = zeros(6,length(TimeStamps));
        cnt = 1;
        for aa = 1:4
            cur_set_a = squeeze(DataPoints(:,aa,:));
            cur_set_a = bsxfun(@minus,cur_set_a,mean(cur_set_a));
            for bb = aa+1:4
                cur_set_b = squeeze(DataPoints(:,bb,:));
                cur_set_b = bsxfun(@minus,cur_set_b,mean(cur_set_b));
                
                shared_var = sum(cur_set_a.*cur_set_b);
                spk_corrs(cnt,:) = shared_var./sqrt(sum(cur_set_a.^2).*sum(cur_set_b.^2));
                
                cnt = cnt + 1;
            end
        end
    end
    mean_spk_corr = mean(spk_corrs);
    bad_spikes = find(mean_spk_corr > max_overlap);
    TimeStamps(bad_spikes) = [];
    DataPoints(:,:,bad_spikes) = [];
    
    %extract bit conversion factors from header
    bitconv = Header{15};
    conv_factors = bitconv(13:end);
    conv_factors = str2num(conv_factors);
    [n_pts,n_ch,n_spks] = size(DataPoints);
    DataPoints = DataPoints.*repmat(conv_factors,[n_pts,1,n_spks])*1e6; %in uV
    [peaks,peak_locs] = max(DataPoints(4:20,:,:),[],1);
    peaks = squeeze(peaks); peak_locs = squeeze(peak_locs);
    
    for i = 1:4
        cur_spikes = find(peaks(i,:) > amp_threshold);
        cur_spike_wave = squeeze(DataPoints(:,i,cur_spikes));
        cur_avg_waveform(i,:) = mean(cur_spike_wave,2);
        cur_std_waveform(i,:) = std(cur_spike_wave,[],2);
        [~,min_loc] = min(cur_spike_wave(6:end,:));
        min_loc = min_loc + 5;
        mua_widths{i} = [mua_widths{i} min_loc - peak_locs(i,cur_spikes)];
        mua_times{i} = [mua_times{i} TimeStamps(cur_spikes)];
        mua_amps{i} = [mua_amps{i} peaks(i,cur_spikes)];
    end
    sc1_avg_waveform = cat(3,sc1_avg_waveform,cur_avg_waveform);
    sc1_std_waveform = cat(3,sc1_std_waveform,cur_std_waveform);
end
sc1_avg_waveform = squeeze(mean(sc1_avg_waveform,3));
sc1_std_waveform = squeeze(mean(sc1_std_waveform,3));

FieldSelection = [1 0 0 0 0];
Filename = 'Sc2.ntt';
[TimeStamps, Header] = Nlx2MatSpike_v3(Filename, FieldSelection, ExtractHeader, ExtractMode);
tot_n_spks = length(TimeStamps);
n_chunks = ceil(tot_n_spks/chunk_size);

FieldSelection = [1 0 0 0 1];
sc2_avg_waveform  = [];
sc2_std_waveform = [];
for i = 5:8
    mua_widths{i} = [];
    mua_times{i} = [];
    mua_amps{i} = [];
end
for nn = 1:n_chunks
    cur_range = [(nn-1)*chunk_size+1 nn*chunk_size];
    cur_range(cur_range > tot_n_spks) = tot_n_spks;
    [TimeStamps, DataPoints, Header] = Nlx2MatSpike_v3(Filename, FieldSelection, ExtractHeader, 2,cur_range);
    
    spk_corrs = nan(6,length(TimeStamps));
    n_chunks = ceil(length(TimeStamps)/chunk_size);
    for cc = 1:n_chunks
        cur_chunk = (cc-1)*chunk_size + (1:chunk_size);
        cur_chunk(cur_chunk > length(TimeStamps)) = [];
        cnt = 1;
        for aa = 1:3
            cur_set_a = squeeze(DataPoints(:,aa,cur_chunk));
            cur_set_a = bsxfun(@minus,cur_set_a,mean(cur_set_a));
            for bb = aa+1:3
                cur_set_b = squeeze(DataPoints(:,bb,cur_chunk));
                cur_set_b = bsxfun(@minus,cur_set_b,mean(cur_set_b));
                
                shared_var = sum(cur_set_a.*cur_set_b);
                spk_corrs(cnt,cur_chunk) = shared_var./sqrt(sum(cur_set_a.^2).*sum(cur_set_b.^2));
                
                cnt = cnt + 1;
            end
        end
    end
    mean_spk_corr = nanmean(spk_corrs);
    bad_spikes = find(mean_spk_corr > 0.5);
    % bad_spikes = find(ismember(TimeStamps,bad_set));
    TimeStamps(bad_spikes) = [];
    DataPoints(:,:,bad_spikes) = [];
    
    %extract bit conversion factors from header
    bitconv = Header{15};
    conv_factors = bitconv(13:end);
    conv_factors = str2num(conv_factors);
    [n_pts,n_ch,n_spks] = size(DataPoints);
    DataPoints = DataPoints.*repmat(conv_factors,[n_pts,1,n_spks])*1e6; %in uV
    [peaks,peak_locs] = max(DataPoints(4:20,:,:),[],1);
    peaks = squeeze(peaks); peak_locs = squeeze(peak_locs);
    
    for i = 1:4
        cur_spikes = find(peaks(i,:) > amp_threshold);
        cur_spike_wave = squeeze(DataPoints(:,i,cur_spikes));
        cur_avg_waveform(i,:) = mean(cur_spike_wave,2);
        cur_std_waveform(i,:) = std(cur_spike_wave,[],2);
        [~,min_loc] = min(cur_spike_wave(6:end,:));
        min_loc = min_loc + 5;
        mua_widths{i+4} = [mua_widths{i+4} min_loc - peak_locs(i,cur_spikes)];
        mua_times{i+4} = [mua_times{i+4} TimeStamps(cur_spikes)];
        mua_amps{i+4} = [mua_amps{i+4} peaks(i,cur_spikes)];
    end
    sc2_avg_waveform = cat(3,sc2_avg_waveform,cur_avg_waveform);
    sc2_std_waveform = cat(3,sc2_std_waveform,cur_std_waveform);
end
sc2_avg_waveform = squeeze(mean(sc2_avg_waveform,3));
sc2_std_waveform = squeeze(mean(sc2_std_waveform,3));

avg_waveform = [sc1_avg_waveform; sc2_avg_waveform];
std_waveform = [sc1_std_waveform; sc2_std_waveform];