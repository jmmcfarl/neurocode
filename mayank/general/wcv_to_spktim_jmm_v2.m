load all_eeg_data

[len,widwcv] = size(CSC1_Samples);
% [lfplen] = length(CSC8_Samples(1,:));

% Get the data.

wcv = reshape(CSC1_Samples,len*widwcv,1);
f2 = CSC8_SampleFrequencies(1);
f = CSC1_SampleFrequencies(1);

clear CSC1* CSC2* CSC3* CSC4* CSC5* CSC6* CSC7* CSC8_Ch* CSC8_Num*

    lof = 20;  %default 300
    if f > 16000
        hif = 7999;
    else
        hif = f/2.1;
    end
    nyqf = f/2;
    lof = lof/nyqf; hif = hif/nyqf;
    [b,a] = butter(2, [lof hif]);
    chunk_size = 300;

    wcv = -wcv; % this sign flip ensures that spikes go upwards.

    load sync_times synct synct1id


    wcvlen = length(wcv);

    time_dur = wcvlen/f;
    num_chunks = ceil(time_dur/chunk_size);
    spk_locs = [];

    for i = 1:num_chunks
        seg_beg = chunk_size*f*(i-1)+1;
        seg_end = seg_beg+chunk_size*f;
        if i < num_chunks
            wcv_piece = wcv(seg_beg:seg_end);
        else
            wcv_piece = wcv(seg_beg:end);
        end
        wcv_f = filtfilt(b,a,wcv_piece);
        if i < num_chunks
            cur_max_spike = max(wcv_f);
        end
        wcv_f = wcv_f/cur_max_spike;
        if ~isempty(find(wcv_f > 0.4))
            [temp_peaks,temp_spk_locs] = findpeaks(wcv_f,'minpeakheight',0.4);
            slope = abs([0;diff(wcv_piece)]);
            bad_spikes = [];
            for s = 1:length(temp_spk_locs)-1
                temp = sum(slope(temp_spk_locs(s):temp_spk_locs(s+1)));
                if temp <1
                    bad_spikes = [bad_spikes s+1];
                end
            end
            temp_spk_locs(bad_spikes) = [];
            spk_isis = diff(temp_spk_locs)/f*1000;
            temp_spk_locs(spk_isis < 1.5) = [];

            spk_locs = [spk_locs temp_spk_locs+seg_beg];
        end
    end
        i
    
    
    spk_locs(spk_locs > max(synct1id)) = [];
    spk_ids{d} = spk_locs;
    clear spk_locs
    save spike_time_jmm spk_ids

