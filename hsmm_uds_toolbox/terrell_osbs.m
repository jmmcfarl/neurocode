function sm_bandwidth = terrell_osbs(data)

%computes Terrell's oversmoothing bandwidth selector from:
%Terrell, G. (1990). "The maximal smoothing principle in density estimation." J. Amer. Statist. Assoc. 85: 470-477.
	
data = data(:);
N = length(data);
C = 3*35^(-1/5)/(2*sqrt(pi))^(1/5);
sm_bandwidth = C *std(data)/N^(1/5);
