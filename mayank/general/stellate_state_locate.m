load('C:\WC_Germany\JMM_analysis_ste\stellate_heka_dir.mat')
% load('C:\WC_Germany\JMM_analysis_pyr\pyr_heka_dir.mat')

d = 2;
load(f_loc{d})
dat_name = f_loc{d}(25:end);

eval(['data = ' dat_name '_MP;'])
eval(['time = ' dat_name '_sampletimes;'])
Fs = 5000;
dsf = 10;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,3/niqf,'low');
[b2,a2] = butter(2,[7/niqf 20/niqf]);

sweep_t = find(diff(time) < 0);
low = min(data);
high = max(data)-0.3;
figure
dataseg = data(1:sweep_t(1));
filtseg = filtfilt(b,a,dataseg);

minpeakheight = 0.1;
maxnoise = 0.04;
mindowndur = 0.75; %min dur in sec

down_start = [];
down_stop = [];

for i = 1:length(sweep_t)-1
    dataseg = data(sweep_t(i)+1:sweep_t(i+1));
    filtseg = filtfilt(b,a,dataseg);
    % filtseg2 = filtfilt(b2,a2,dataseg);
    down_filt = downsample(filtseg,dsf);
    diff_filt = [0;diff(down_filt)]*Fsd;
    
    [dummy,curpeaks] = findpeaks(diff_filt,'minpeakheight',minpeakheight);
    cur_down_start = [];
    cur_down_stop = [];
    
    %cycle through all threshold crossings
    for p = 1:length(curpeaks)
        
        %find most recent instance where slope was below noise level
        prev_quiet = find(diff_filt(1:curpeaks(p))<maxnoise,1,'last');
        if ~isempty(prev_quiet)
            
            %find previous instance where magnitude of slope was above
            %noise level
            prev_noise = find(abs(diff_filt(1:prev_quiet))>maxnoise,1,'last');            
            if ~isempty(prev_noise)
                
                %if its a negative slope crossing
                if diff_filt(prev_noise) < 0
                    
                    %if down state is long enough
                    if (prev_quiet-prev_noise)/Fsd > mindowndur
                        cur_down_start = [cur_down_start prev_noise];
                        cur_down_stop = [cur_down_stop prev_quiet];
                    end
                end
                
            elseif prev_quiet/Fsd > mindowndur
               cur_down_start = [cur_down_start 1];
               cur_down_stop = [cur_down_stop prev_quiet];
            end
            
        end
        
    end
    
    subplot(2,1,1)
       plot(time(sweep_t(i)+1:sweep_t(i+1)),data(sweep_t(i)+1:sweep_t(i+1)))
       hold on
          plot(time(sweep_t(i)+1:sweep_t(i+1)),filtseg,'r','linewidth',2)
    %       plot(time(sweep_t(i)+1:sweep_t(i+1)),filtseg2*3+mean(dataseg),'g','linewidth',2)
       ylim([low high])
          subplot(2,1,2)
    plot(diff_filt,'k')
    hold on
    plot(curpeaks,diff_filt(curpeaks),'ro')
    plot(cur_down_start,diff_filt(cur_down_start),'go')
    plot(cur_down_stop,diff_filt(cur_down_stop),'mo')
    line([0 length(diff_filt)],[minpeakheight minpeakheight],'Color','r')
    line([0 length(diff_filt)],[maxnoise maxnoise],'Color','b')
    line([0 length(diff_filt)],[-maxnoise -maxnoise],'Color','b')

   grid
   pause
   clf
   
end