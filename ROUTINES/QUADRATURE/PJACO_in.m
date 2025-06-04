function [v, intv2] = PJACO(n, alf,bt, x, varargin)
%PJACO returns the Jacobi polynomial
%
%	USAGE:
%		[v, intv2] = PJACO(n, alf,bt, x);
%	INPUTS:
%		n, alf,bt, x
%	OUTPUTS:
%		v, intv2
    if length(varargin)==2
        a = varargin{1};
        b = varargin{2};
    else
        a = 0;
        b = 1;
    end
    x = 2*(x-a)/(b-a)-1;
    
    v = jacobiP(n, alf,bt, x);
    intv2 = ((b-a)^(alf+bt+1)./(2*n+bt+alf+1).*...
             gamma(n+bt+1).*gamma(n+alf+1)./...
             (gamma(n+bt+alf+1).*factorial(n)));
end
