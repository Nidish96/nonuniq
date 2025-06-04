function [R, dRdUxw, dRdl] = EPMCRESFUN(uxwl, Fstat, Lb, h, Nt, fpars, mpars, varargin)
%EPMCRESFUN returns the residue vector and appropriate derivatives for
%EPMC
%
%   USAGE:
%       [R, dRdUxw, dRdl] = EPMCRESFUN(uxwl, Fstat, Fl, h, Nt, fpars, mpars);
%   INPUTS:
%       uxwl    : (Nd*Nhc+2, 1) 
%       Fstat   : (Nd*Nhc, 1)
%       Lb      : (Nd*Nhc+2, Nd*Nhc+1) [Lb*uxwl(1:end-1);uxwl(end)] -> [U; xi; w; log(A)]
%       h       : (Nh, 1)
%       Nt      : (1, 1)
%       fpars   : (5,1) array containing [alpha;kt;kn;mu;f0] such that,
%           alpha       : half angle of triangular wedge
%           kt, kn, mu  : tangential stiffness, normal stiffness,
%                         coefficient of friction
%           f0          : initial force (for iterations)
%       mpars   : Model parameters struct containing (at least)
%           M,C,K       : Linear mass, damping, stiffness matrices
%           Lrels       : (2*Nhc, Nd*Nhc, 2) with pages containing
%                           selectors for the two contacts
%   OUTPUTS:
%       R       : (Nd*Nhc+2, 1)
%       dRdUxw  : (Nd*Nhc+2, Nd*Nhc+2)
%       dRdl    : (Nd*Nhc, 1)

    Uxwl = [Lb*uxwl(1:end-1); uxwl(end)];

    Nhc = sum((h==0)+2*(h~=0));
    Nd = (size(Uxwl,1)-3)/Nhc;

    xi = Uxwl(end-2);
    w = Uxwl(end-1);

    lA = Uxwl(end);
    A = 10^lA;
    dAdl = A*log(10);
    % if ~isempty(varargin)
    %     A = varargin{1}*A;
    %     dAdl = varargin{1}*dAdl;
    % end

    % Uh = Asc.*U;
    [~, ~, zinds, rinds, iinds] = HINDS(Nd, h);
    Asc = zeros(Nd*Nhc,1);  dAscdl = zeros(Nd*Nhc,1);
    Asc(zinds) = 1.0;
    Asc([rinds iinds]) = A;
    dAscdl([rinds iinds]) = dAdl;

    [E, dEdw] = HARMONICSTIFFNESS(mpars.M, mpars.C-xi*mpars.M, mpars.K, w, h);
    dEdxi = HARMONICSTIFFNESS(zeros(Nd), -mpars.M, zeros(Nd), w, h);

    FNL = zeros(Nd*Nhc, 1);
    JNL = zeros(Nd*Nhc, 1);
    Lrels = mpars.Lrels(fpars(1));
    for i=1:size(Lrels,3)  % Iterate over Nonlinearities
        UtUn = Lrels(:,:,i)*reshape(Asc.*Uxwl(1:end-3), Nd, Nhc);
        [Fnl, Jnl] = JENK2D(UtUn, h, Nt, fpars(2:end));

        FNL = FNL + kron(eye(Nhc), Lrels(:,:,i))'*Fnl;
        JNL = JNL + kron(eye(Nhc), Lrels(:,:,i))'*Jnl*kron(eye(Nhc), Lrels(:,:,i));
    end
    % JNL is wrt Uh, not U.

    % Modal Amplitude
    mAmp = Uxwl(rinds(1:Nd))'*mpars.M*Uxwl(rinds(1:Nd)) + Uxwl(iinds(1:Nd))'*mpars.M*Uxwl(iinds(1:Nd));
    dmAmpdU = zeros(1, Nhc*Nd);
    dmAmpdU(rinds(1:Nd)) = 2*Uxwl(rinds(1:Nd))'*mpars.M;
    dmAmpdU(iinds(1:Nd)) = 2*Uxwl(iinds(1:Nd))'*mpars.M;

    % Assemble final residue
%     R      = [E*(Asc.*Uxwl(1:end-3)) + FNL - Fstat;
%               A^2*(mAmp-1.0);
%               Fl'*Uxwl(1:end-3)];
%     dRdUxw = [(E+JNL).*Asc' dEdxi*(Asc.*Uxwl(1:end-3)) dEdw*(Asc.*Uxwl(1:end-3));
%               A^2*(dmAmpdU) zeros(1,2);
%               Fl' zeros(1,2)];
%     dRdl   = [(E+JNL)*(dAscdl.*Uxwl(1:end-3)); 
%         2*A*dAdl*(mAmp-1.0); 0];
    
    % R      = [E*(Asc.*Uxwl(1:end-3)) + FNL - Fstat;
    %           A^2*(mAmp-1.0)];
    % dRdUxw = [(E+JNL).*Asc' dEdxi*(Asc.*Uxwl(1:end-3)) dEdw*(Asc.*Uxwl(1:end-3));
    %           A^2*(dmAmpdU) zeros(1,2)]*Lb;
    % dRdl   = [(E+JNL)*(dAscdl.*Uxwl(1:end-3)); 
    %           2*A*dAdl*(mAmp-1.0)];

    R      = [E*(Asc.*Uxwl(1:end-3)) + FNL - Fstat;
              mAmp-1.0];
    dRdUxw = [(E+JNL).*Asc' dEdxi*(Asc.*Uxwl(1:end-3)) dEdw*(Asc.*Uxwl(1:end-3));
              dmAmpdU zeros(1,2)]*Lb;
    dRdl   = [(E+JNL)*(dAscdl.*Uxwl(1:end-3)); 
              0];

    % if isfield(mpars, 'dK')  
    %     dE = kron(speye(Nhc), mpars.dK);

    %     % full eqn: (E+dE)*(Asc.*Uxwl(1:end-3)) - dE*(Asc.*Uxwl(1:end-3)) + FNL - Fstat;
    %     % Displacement coords: (Asc.*Uxwl(1:end-3)) + (E+dE)\(-dE*(Asc.*Uxwl(1:end-3)) + FNL - Fstat);
        
    %     rv = -dE*(Asc.*Uxwl(1:end-3))+FNL-Fstat;
    %     drvdU = (-dE+JNL).*Asc';

    %     R = [Asc.*Uxwl(1:end-3) + (E+dE)\rv;
    %          (mAmp-1.0)];
    %     dRdUxw = [diag(Asc)+(E+dE)\drvdU, -(E+dE)\[dEdxi*((E+dE)\rv), dEdw*((E+dE)\rv)];
    %               dmAmpdU 0 0]*Lb;
    %     dRdl = [dAscdl.*Uxwl(1:end-3) + (E+dE)\((-dE+JNL)*(dAscdl.*Uxwl(1:end-3)));
    %             0];

    %     % R(1:Nd) = E(1:Nd, :)*(Asc.*(Uxwl(1:end-3))) + FNL(1:Nd) - Fstat(1:Nd);
    %     % dRdUxw(1:Nd, :) = [(E(1:Nd, :)+JNL(1:Nd, :)).*Asc' dEdxi(1:Nd,:)*(Asc.*Uxwl(1:end-3)) dEdw(1:Nd,:)*(Asc.*Uxwl(1:end-3))]*Lb;
    % end
end
