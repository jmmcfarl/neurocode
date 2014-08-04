function [results] = despikeLFPbyChunks(y,S,Bs,g,opts,overlap)
    %function [results] = despikeLFPByChunks(y,S,Bs,g,opts,overlap)
    %
    %Despikes an LFP by splitting the wideband signal into chunks of length
    %length(g), calling despikeLFP for each chunk, and sewing the results
    %back together. 
    %
    %The arguments y,S,Bs,g,opts are as in despikeLFP except that
    %length(y) > length(g) since the chunking algorithm takes cares of
    %setting the lengths properly. opts is simply forwarded to despikeLFP
    %except for opts.displaylevel which is used to set the amount of
    %feedback from both despikeLFP and despikeLFPbyChunks. It can take
    %values 0 (no feedback), 1 (some feedback, default), 2 (full debug
    %info)
    %
    %overlap specifies the extent of overlap between chunks. If unspecified 
    %it is set to 1/32. It has to be smaller than 1/6. There should be no
    %reason to want to change this parameter.
    %
    %Chunking reconstructs the signal z by summing partial zs estimated in 
    %chunks. If size(zs) == [length(y),nchunks], the length of a chunk, then 
    %z is given by:
    %
    % z = sum(recwindows.*zs)
    %
    %And the reconstruction windows (basis functions) look like this:
    %
    % First:
    % ------\
    %        \
    %         \---------------------------------------------------------
    % Second:
    %         /-------\
    %        /         \
    % ------/           \-----------------------------------------------
    % Third:
    %                  /-------\
    %                 /         \
    % ---------------/           \--------------------------------------
    % (etc.)
    % Last:
    %                                                            /------
    %                                                           /
    % ---------------------------------------------------------/
    %
    %Since the basis functions have compact support, we can estimate one of the zs
    %by cutting off a slightly larger portion of the signal than the size of
    %the basis function and feeding that into despikeLFP to obtain the zs
    %for that portion of the signal. Thus the signal is split into
    %overlapping chunks, with the support of each chunk being larger by
    %4*overlap*length(g) than the support of its corresponding basis. 
    %The basis functions overlap by overlap*length(g). If the numbers don't
    %work out, the procedure adjusts the overlap upwards to make everything
    %fit properly.
    %
    %To avoid edge artifacts, within each chunk we add a tapered mirror image 
    %of the content at an edge to the other edge. The taper has width
    %overlap*length(g). The first and last chunks are dealt with similarly
    %but the taper is longer and we discard more of the signal.
    %
    %The function returns a struct called results which contains the
    %following keys:
    %
    %z: the despiked wideband signal, sewn back together from the chunks
    %phi,gamma,mu,sigma: as in despikeLFP but now as arrays.
    %opts: The options given to the function
    %
    %Because the function preassigns every variable, if you don't run out
    %of memory by the second or third chunk you should be fine.
    %
    % History:
    %
    % 07/07/2010: changed which spikes were tagged in each chunk from those
    % for which win = 1 to those for which recwin >0
    % 23/06/2010: v1.0
    %
    % Created by Patrick Mineault
    %
    %See also despikeLFP
    
    %Sanitize input
    if nargin < 5
        opts = struct();
    end
    
    if ~isfield(opts,'displaylevel')
        opts.displaylevel = 1;
    end
    if nargin < 6
        overlap = 1/32;
    end
    
    %Count number of chunks
    n = length(g);
    N = length(y);
    nchunks = ceil((N - 2*n*(1-overlap*5.5))/(n-5*overlap))+2;
    
    %equispace chunks as much as possible
    overlap = (nchunks*n-N)/(11*n+5*n*(nchunks-2));
    w = round(overlap*n);
    knot1 = floor( (N - (nchunks-2)*(n-5*w) - w)/2 ) - 2*w;
    
    %Preassign variables to keep out of memory errors at the beginning
    results.z = zeros(size(y));
    results.phi = zeros([nchunks,size(Bs,2),size(S,2)]);
    results.gamma = zeros(nchunks,1);
    results.sigma = zeros(nchunks,1);
    results.mu = zeros(nchunks,1);
    results.zs = zeros(n,nchunks);
    
    %Deal with first chunk first
    win = ones(n,1);
    win(n -4*w + (1:w*2)) = 1-(1:w*2)/(w*2);
    win((n-2*w+1):end) = 0;
    
    recwin = ones(n,1);
    recwin(knot1 + 2*w + (1:w)) = 1-(1:w)/(w);
    recwin((knot1+3*w+1):end) = 0;
    
    if opts.displaylevel > 0
        fprintf('Despiking chunk #%d/%d\n\n',1,nchunks);
    end
    
    ys = y(1:n).*win + y(n:-1:1).*(1-win);
    Ss = bsxfun(@times,S(1:n,:),recwin>0);
    result = despikeLFP(ys,Ss,Bs,g,opts);
    
    ii = 1;
    results.zs(:,ii) = result.z;
    results.phi(ii,:,:) = result.phi;
    results.mu(ii) = result.mu;
    results.gamma(ii) = result.gamma;
    results.sigma(ii) = result.sigma;
    
    %opts.gamma = result.gamma;
    %opts.sigma = result.sigma;
    

    results.z(1:n) = (result.z+result.mu).*recwin;
    
    %Create the windows for the second to second-to-last chunks
    win = ones(n,1);
    win((1:w)) = 1/2 + 1/2*(1:w)/(w);
    win(n-w + (1:w)) = 1-1/2*(1:w)/(w);
    
    recwin = ones(n,1);
    recwin(1:2*w) = 0;
    recwin(2*w + (1:w)) = (1:w)/(w);
    recwin(n-3*w + (1:w)) = 1-(1:w)/(w);
    recwin(n - 2*w + (1:2*w)) = 0;

    for ii = 2:nchunks-1
        
        if opts.displaylevel > 0
            fprintf('Despiking chunk #%d/%d\n\n',ii,nchunks);
        end
        
        rg = knot1 + (ii-2)*(n-5*w) + (1:n);
        ys = y(rg).*win + y(rg(end:-1:1)).*(1-win);
        Ss = bsxfun(@times,S(rg,:),recwin>0);
        result = despikeLFP(ys,Ss,Bs,g,opts);
        results.zs(:,ii) = result.z;
        results.z(rg) = results.z(rg) + (result.z+result.mu).*recwin;

        results.zs(:,ii) = result.z;
        results.phi(ii,:,:) = result.phi;
        results.mu(ii) = result.mu;
        results.gamma(ii) = result.gamma;
        results.sigma(ii) = result.sigma;

        %opts.gamma = result.gamma;
        %opts.sigma = result.sigma;
    end
    
    win = ones(n,1);
    win(n-4*w + (1:w*2)) = 1-(1:w*2)/(w*2);
    win((n-2*w+1):end) = 0;
    
    win = win(end:-1:1);
    
    if opts.displaylevel > 0
        fprintf('Despiking chunk #%d/%d\n\n',nchunks,nchunks);
    end

    cutoff = n - (N - (knot1 + 2*w+ (n-5*w)*(nchunks-2)));
    
    recwin = zeros(n,1);
    recwin(cutoff + (1:w)) = (1:w)/(w);
    recwin((cutoff + w + 1):end) = 1;
    
    %Now for the last chunk
    rg = N-n + (1:n);
    ys = y(rg).*win + y(rg(end:-1:1)).*(1-win);
    Ss = bsxfun(@times,S(rg,:),recwin>0);
    result = despikeLFP(ys,Ss,Bs,g,opts);
    
    results.z(rg) = results.z(rg) + (result.z+result.mu).*recwin;
    
    results.z = results.z - mean(results.z);
    
    ii = nchunks;
    results.zs(:,ii) = result.z;
    results.phi(ii,:,:) = result.phi;
    results.mu(ii) = result.mu;
    results.gamma(ii) = result.gamma;
    results.sigma(ii) = result.sigma;

    results.phi = squeeze(results.phi);
    results.opts = opts;

end