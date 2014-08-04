function nim_out = NMMfit_alt( nim, Robs, Xstims, Gmults, targets, silent, LLtolerance, desired_optim_params, regmat_custom )
%
% Usage: nim_out = NMMfit_alt( nim, Robs, Xstims, <Gmults>, <targets>, <silent>, <LLtolerance>, <desired_optim_params>, <regmat_custom> )
% 
% will alternate between fitting filters and upstreamNLs until better improvements 
% between successive iterations is less than LLtolerance (default = 0.01)

MAXiter = 10;

if nargin < 4
	Gmults = [];
end
if nargin < 5
	targets = [];
end
if (nargin < 6) || isempty(silent)
	silent = 1;
end
if (nargin < 7) || isempty(LLtolerance)
	LLtolerance = 0.01;
end
if nargin < 8
	desired_optim_params = [];
end
if nargin < 9
	regmat_custom = [];
end

LLpast = -1e6;
LL = nim.LL_seq(end);
nim_out = nim;
iter = 0;

if ~silent
	fprintf( 'Beginning LL = %f\n', LL )
end

while (((LL-LLpast) > LLtolerance) && (iter < MAXiter))
	
	nim_out = NMMfit_filters( nim_out, Robs, Xstims, Gmults, targets, 1, desired_optim_params, regmat_custom );
	nim_out = NMMfit_upstreamNLs( nim_out, Robs, Xstims, Gmults, targets, [], 1, desired_optim_params, regmat_custom );
	LLpast = LL;
	LL = nim_out.LL_seq(end);
	
	iter = iter + 1;

	if ~silent
		fprintf( '  Iter %2d: LL = %f\n', iter, LL )
	end

end
