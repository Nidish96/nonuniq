function [Fnl, Jnl] = JENK2D(UtUn, h, Nt, fpars, varargin)
%JENK2D Returns the nonlinear forces corresponding to a 2D Jenkins Element
%given the tangential and normal displacements (everything in frequency
%domain).
%   UtUn harmonics are ordered as, [A0; A1; B1; A2; B2; ...]
%       where Ai = [Ut_Ai; Un_Ai];
%       and Fourier Expansion is
%   utun(t) = A0 + sum_i Ai*cos(i*tau) + Bi*sin(i*tau);
%
%   USAGE:
%       [Fnl, Jnl] = JENK2D(UtUn, h, Nt, pars);
%   INPUTS:
%       UtUn    : (2*Nhc, 1) Harmonics of tangential and normal DOFs
%       h       : (Nh, 1)
%       Nt      : (1, 1)
%       fpars   : (4,1) array containing [kt;kn;mu;m] such that,
%           kt, kn, mu  : tangential stiffness, normal stiffness,
%                         coefficient of friction
%           m          : initial multiplier (for iterations)
%   OUTPUTS:
%       Fnl     : (2*Nhc, 1)
%       Jnl     : (2*Nhc, 2*Nhc)

    if isempty(varargin)
        tol = 1e-10;
    else
        tol = varargin{1};
    end

    Nhc = sum((h==0)+2*(h~=0));
    ftfn = zeros(Nt,2);

    % Parameters 
    kt = fpars(1);
    kn = fpars(2);
    mu = fpars(3);
    m  = fpars(4);
    
    % Frequency-Time Transformation
    utun = TIMESERIES_DERIV(Nt, h, reshape(UtUn,2,Nhc)', 0);
    cst  = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);

    % Normal Contact Force
    ftfn(:,2) = max(kn*utun(:,2), 0);
    dfndUn = kn*(utun(:,2)>0).*cst;

    % Tangential Contact Force
    ftfn(Nt,1) = m*mu*ftfn(Nt,2);  %% Setting Scaling
    dftdUt = zeros(Nt, Nhc);
    dftdUn = zeros(Nt, Nhc);
    dftdUn(Nt,:) = m*mu*dfndUn(Nt,:);
    its = 0;
    while its==0 || abs(fprev-ftfn(end,1))/mu>tol
        fprev = ftfn(end,1);
        for ti=1:Nt
            tim1 = mod(ti-1-1,Nt)+1;
    
            % Stick Prediction
            fsp = kt*(utun(ti,1)-utun(tim1,1)) + ftfn(tim1,1);
            if abs(fsp)<mu*ftfn(ti,2)  % stick
                ftfn(ti, 1) = fsp;
                dftdUt(ti, :) = kt*(cst(ti,:)-cst(tim1,:)) + dftdUt(tim1,:);
                dftdUn(ti, :) = dftdUn(tim1,:);
            else  % slip
                ftfn(ti, 1) = mu*ftfn(ti, 2)*sign(fsp);
                dftdUt(ti, :) = 0;
                dftdUn(ti, :) = mu*dfndUn(ti, :)*sign(fsp);
            end
        end
        its = its + 1;
        if its>10
            % warning('Non Convergence');
            % keyboard;
            break;
        end
    end

    % Time-Frequency Conversion
    Fnl = reshape(GETFOURIERCOEFF(h, ftfn)', 2*Nhc, 1);
    Jnl = zeros(2*Nhc);
    Jnl(1:2:end, 1:2:end) = GETFOURIERCOEFF(h, dftdUt);
    Jnl(1:2:end, 2:2:end) = GETFOURIERCOEFF(h, dftdUn);
    Jnl(2:2:end, 2:2:end) = GETFOURIERCOEFF(h, dfndUn);
end
