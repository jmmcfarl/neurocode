function [trial_ids,seed_nums,frame_trialid,frame_nums,left_pix,right_pix,trial_ds,has_mtrS] = parse_random_bar_file_dispExpt(fname)
%
%% Same as parse_random_bar_file, but adjusted to handle new data from disparity expts. These files have extra
% In these files there's an extra L0R0 at the end of the lines. In Bruce's words:
% The extra L0R0 lines have been added to deal with cases that did not happen that day.  
% [ If the disparity is not a integer multipe of the bar width, then (a) the pattern will start with a partial bar, and (b) the size of this will be different in the two eyes.  The L/R numbers record this. But I was careful that the disparity was always an even multiple. 
% So, this should be handled properly at some point! (Bruce probably has code to do this

%open rls.rc file and read in text data
fid = fopen(fname);
data = textscan(fid,'%s','delimiter','\n','bufsize',2^13-1);
data = data{1};
n_lines = length(data);
fclose(fid);

if length(data) < 5 %if there are too few lines
    trial_ids =[];
    seed_nums = [];
    left_pix = [];
    right_pix = [];
    frame_trialid = [];
    frame_nums = [];
    trial_ds = [];
    has_mtrS = nan;
    fprintf('%s has too few data, ignoring...\n',fname);
    return
end


%% first scan through text and identify lines
ds_lines = find(strncmp('ds',data,2));
mtrS_lines = find(strncmp('mtrS',data,4));
Resp_lines = find(strncmp('Resp',data,4));
trial_id_lines = find(strncmp('id',data,2));
frame_lines = setdiff(1:n_lines,[ds_lines; mtrS_lines; trial_id_lines; Resp_lines]);
line_lengths = cellfun(@length,data);
colon_locs = cell2mat(strfind(data(frame_lines),':'));
L_frames = cellfun(@(x) ~isempty(strfind(x,'L')),data(frame_lines));
nan_frames = cellfun(@(x) ~isempty(strfind(x,'N')),data(frame_lines));

L_locs = line_lengths(frame_lines);
has_L = cellfun(@(x) ~isempty(strfind(x,'L')),data(frame_lines));
L_locs(has_L) = cell2mat(strfind(data(frame_lines(has_L)),'L'));
L_locs(nan_frames) = nan;
% L_locs(~nan_frames) = cell2mat(strfind(data(frame_lines(~nan_frames)),'L'));
n_pix = nan(length(frame_lines),1);
n_pix(~nan_frames) = L_locs(~nan_frames) - colon_locs(~nan_frames)-1;
max_npix = nanmax(n_pix);

trial_ids = nan(length(trial_id_lines),1);
seed_nums = nan(length(trial_id_lines),1);
for ii = 1:length(trial_id_lines)
    a = sscanf(data{trial_id_lines(ii)},'id%dse%d');
    trial_ids(ii) = a(1); seed_nums(ii) = a(2);
end


%%
% is_ds = false(n_lines,1);
% is_mtrS = false(n_lines,1);
% trial_ids = nan(n_lines,1);
% seed_nums = nan(n_lines,1);
% frame_nums = nan(n_lines,1);
% n_pix = nan(n_lines,1);
% cloc = nan(n_lines,1);
% for ii = 1:n_lines
%     if strcmp(data{ii}(1:2),'id')
%         seloc = find(data{ii} == 's');
%         trial_ids(ii) = str2num(data{ii}(3:seloc-1));
%         seed_nums(ii) = str2num(data{ii}(seloc+2:end));
%     elseif strcmp(data{ii}(1:2),'ds')
%         is_ds(ii) = true;
%     elseif strcmp(data{ii}(1:4),'mtrS')
%         is_mtrS(ii) = true;
%     else
%         cloc(ii) = find(data{ii} == ':');
%         frame_nums(ii) = hex2dec(data{ii}(1:cloc(ii)-1));
%         n_pix(ii) = length(data{ii}) - cloc(ii);
%     end
% end
% max_npix = nanmax(n_pix);

%% NOW READ IN FRAME DATA
fprintf('Max pixel size %d\n',max_npix);
n_frames = length(frame_lines);
frame_nums = nan(n_frames,1);
buffers = floor((max_npix - n_pix)/2);
pix_data = char(zeros(n_frames,max_npix));
for ii = 1:n_frames
    frame_nums(ii) = hex2dec(data{frame_lines(ii)}(1:(colon_locs(ii)-1)));
    if ~nan_frames(ii)
    pix_data(ii,(1:n_pix(ii)) + buffers(ii)) = ...
        data{frame_lines(ii)}(colon_locs(ii)+1:L_locs(ii)-1);
    end
end

%COMPUTE LEFT AND RIGHT PIXEL DATA
left_pix = -999*ones(n_frames,max_npix,'int8');
right_pix = -999*ones(n_frames,max_npix,'int8');

%Hex   Left  Right
%0   =    B      B
%1   =    G      B
%2   =    W      B
%4   =    B      G
%5   =    G      G
%6   =    W      G
%8   =    B      W
%9   =    G      W
%a   =    W      W

left_pix(pix_data == '0') = -1;
right_pix(pix_data == '0') = -1;

left_pix(pix_data == '1') = 0;
right_pix(pix_data == '1') = -1;

left_pix(pix_data == '2') = 1;
right_pix(pix_data == '2') = -1;

left_pix(pix_data == '4') = -1;
right_pix(pix_data == '4') = 0;

left_pix(pix_data == '5') = 0;
right_pix(pix_data == '5') = 0;

left_pix(pix_data == '6') = 1;
right_pix(pix_data == '6') = 0;

left_pix(pix_data == '8') = -1;
right_pix(pix_data == '8') = 1;

left_pix(pix_data == '9') = 0;
right_pix(pix_data == '9') = 1;

left_pix(pix_data == 'a') = 1;
right_pix(pix_data == 'a') = 1;

is_binoc = any(left_pix(:) ~= right_pix(:));

%%
% trial_ids = trial_ids(~isnan(trial_ids));
% seed_nums = seed_nums(~isnan(seed_nums));
n_trials = length(trial_ids);
frame_trialid = zeros(n_frames,1);

frame_nums = frame_nums(~isnan(frame_nums));
trial_start_inds = find(frame_nums == 0);
for ii = 1:n_trials-1
    frame_trialid(trial_start_inds(ii):trial_start_inds(ii+1)) = ii;
end
frame_trialid(trial_start_inds(end):end) = n_trials;

trial_ds = false(n_trials,1);
for ii = 1:n_trials-1
    temp = find(ds_lines > trial_id_lines(ii) & ds_lines < trial_id_lines(ii+1));
    if ~isempty(temp)
        trial_ds(ii) = true;
    end
end
temp = find(ds_lines > trial_id_lines(n_trials));
if ~isempty(temp)
    trial_ds(end) = true;
end


has_ds = ~isempty(ds_lines);
has_mtrS = ~isempty(mtrS_lines)
fprintf('%d trials,  binoc = %d, Guidedsac = %d, hasmtrS = %d \n',n_trials,is_binoc,has_ds,has_mtrS);


