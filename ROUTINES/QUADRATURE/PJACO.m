function [v, intv2] = PJACO ( n, alf, bt, x, varargin )
%*****************************************************************************80
%
%% j_polynomial() evaluates the Jacobi polynomials J(N,A,B,X).
%
%  Differential equation:
%
%    (1-X*X) Y'' + (ALPHA-BETA-(BETA+ALPHA+2) X) Y' + N (N+BETA+ALPHA+1) Y = 0
%
%  Recursion:
%
%    P(0,BETA,ALPHA,X) = 1,
%
%    P(1,BETA,ALPHA,X) = ( (2+BETA+ALPHA)*X + (BETA-ALPHA) ) / 2
%
%    P(N,BETA,ALPHA,X)  = 
%      ( 
%        (2*N+BETA+ALPHA-1) 
%        * ((BETA^2-ALPHA^2)+(2*N+BETA+ALPHA)*(2*N+BETA+ALPHA-2)*X) 
%        * P(N-1,BETA,ALPHA,X)
%        -2*(N-1+BETA)*(N-1+ALPHA)*(2*N+BETA+ALPHA) * P(N-2,BETA,ALPHA,X)
%      ) / 2*N*(N+BETA+ALPHA)*(2*N-2+BETA+ALPHA)
%
%  Restrictions:
%
%    -1 < BETA
%    -1 < ALPHA
%
%  Norm:
%
%    Integral ( -1 <= X <= 1 ) ( 1 - X )^BETA * ( 1 + X )^ALPHA 
%      * P(N,BETA,ALPHA,X)^2 dX 
%    = 2^(BETA+ALPHA+1) * Gamma ( N + BETA + 1 ) * Gamma ( N + ALPHA + 1 ) /
%      ( 2 * N + BETA + ALPHA ) * N! * Gamma ( N + BETA + ALPHA + 1 )
%
%  Special values:
%
%    P(N,BETA,ALPHA)(1) = (N+BETA)!/(N!*BETA!) for integer BETA.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    18 March 2012
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Milton Abramowitz, Irene Stegun,
%    Handbook of Mathematical Functions,
%    US Department of Commerce, 1964.
%
%  Input:
%
%    integer M, the number of evaluation points.
%
%    integer N, the highest order polynomial to compute.  Note
%    that polynomials 0 through N will be computed.
%
%    real BETA, ALPHA, the parameters.
%    -1 < BETA, ALPHA.
%
%    real X(M,1), the evaluation points.
%
%  Output:
%
%    real V(M,1:N+1), the values of the first N+1 Jacobi
%    polynomials at the point X.
%
    if ( bt <= -1.0 )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'J_POLYNOMIAL - Fatal error!\n' );
        fprintf ( 1, '  Illegal input value of BETA = %f\n', bt );
        fprintf ( 1, '  But BETA must be greater than -1.\n' );
        error ( 'J_POLYNOMIAL - Fatal error!' );
    end
    
    if ( alf <= -1.0 )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'J_POLYNOMIAL - Fatal error!\n' );
        fprintf ( 1, '  Illegal input value of ALPHA = %f\n', alf );
        fprintf ( 1, '  But ALPHA must be greater than -1.\n' );
        error ( 'J_POLYNOMIAL - Fatal error!' );
    end
    if length(varargin)==2
        a = varargin{1};
        b = varargin{2};
    else
        a = 0;
        b = 1;
    end
    xor = x(:);

    x = (x-a)/(b-a)*2-1;
    m = length(x);

    % For consistency with MATLAB's internal notation
    alf = alf-1;
    bt = bt-1;
    
    if ( n < 0 )
        v = [];
        return
    end
    v = ones(m, n+1);
    ns = (0:n);    
    L2 = sqrt((b-a)^(bt+alf+1)./(2*ns+bt+alf+1).*...
              gamma(ns+bt+1).*gamma(ns+alf+1)./...
              (gamma(ns+bt+alf+1).*factorial(ns)))/...
              1;
    % sqrt(beta(bt+1,alf+1));
    if alf==bt && alf==-0.5
        L2(1) = sqrt(pi);
    end

    if ( n == 0 )
        v = v;
        intv2 = L2^2;
        return
    end

    v(:, 2) = ( 1.0 + 0.5 * ( bt + alf ) ) * x  + 0.5 * ( bt - alf );
    
    for i = 2 : n

        c1 = 2 * i * ( i + bt + alf ) * ( 2 * i - 2 + bt + alf );
        c2 = ( 2 * i - 1 + bt + alf ) * ( 2 * i + bt + alf ) ...
             * ( 2 * i - 2 + bt + alf );
        c3 = ( 2 * i - 1 + bt + alf ) * ( bt + alf ) * ( bt - alf );
        c4 = - 2 * ( i - 1 + bt ) * ( i - 1 + alf )  * ( 2 * i + bt + alf );

        
        v(:, i+1) = ( (c3+c2*x ) .* v(:, i) + c4 * v(:, i-1) ) / c1;
    end
    
    v = v(:, end);  % This is dumb but it's okay for now
    intv2 = L2(end)^2;
end
