function nim_out = NIMfit_alt( nim, Robs, Xstim, XLin, targets, silent, LLtolerance, desired_optim_params )
%
% Usage: nim_out = NIMfit_alt( nim, Robs, Xstim, <XLin>, <targets>, <silent>, <LLtolerance>, <desired_optim_params> )
% 
% will alternate between fitting filters and upstreamNLs until better improvements 
% between successive iterations is less than LLtolerance (default = 0.01)

MAXiter = 10;

if nargin < 4
	XLin = [];
end
if nargin < 5
	targets = [];
end
if nargin < 6
	silent = 1;
end
if nargin < 7
	LLtolerance = 0.01;
end
if nargin < 8
	desired_optim_params = [];
end

LLpast = -1e6;
LL = nim.LL_seq(end);
nim_out = nim;
iter = 0;

if ~silent
	fprintf( 'Beginning LL = %f\n', LL )
end

while (((LL-LLpast) > LLtolerance) && (iter < MAXiter))
	
	nim_out = NIMfit_filters( nim_out, Robs, Xstim, XLin, targets, 1, desired_optim_params );
	nim_out = NIMfit_upstreamNLs( nim_out, Robs, Xstim, XLin, targets, [], 1, desired_optim_params );
	LLpast = LL;
	LL = nim_out.LL_seq(end);
	
	iter = iter + 1;

	if ~silent
		fprintf( '  Iter %2d: LL = %f\n', iter, LL )
	end

end
