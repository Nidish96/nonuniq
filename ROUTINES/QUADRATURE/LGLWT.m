function [x,w,P]=LGLWT(N, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lglnodes.m
%
% Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde 
% matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
% integration and spectral methods. 
%
% Reference on LGL nodes and weights: 
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 04/17/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin==1
        a=-1;
        b=1;
    else
        a=varargin{1};
        b=varargin{2};
    end
    
    % Truncation + 1
    N1=N+1;
    % Use the Chebyshev-Gauss-Lobatto nodes as the first guess
    x=cos(pi*(N:-1:0)/N)';
    % The Legendre Vandermonde Matrix
    P=zeros(N1,N1);
    % Compute P_(N) using the recursion relation
    % Compute its first and second derivatives and 
    % update x using the Newton-Raphson method.
    xold=2;
    while max(abs(x-xold))>eps
        xold=x;

        P(:,1)=1;    P(:,2)=x;

        for k=2:N
            P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
        end

        x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );

    end
    w=2./(N*N1*P(:,N1).^2);
    
	% Linear map from[-1,1] to [a,b]
    x=(a*(1-x)+b*(1+x))/2;      
    w=(b-a)/2*w;
end