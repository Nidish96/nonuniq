function [R, dRdU, dRdw] = HBRESFUN(Uw, Fl, h, Nt, fpars, mpars)
%HBRESFUN returns the residue vector and appropriate derivatives for
%Harmonic Balance
%
%   USAGE:
%       [R, dRdU, dRdw] = HBRESFUN(Uw, Fl, h, Nt, fpars, mpars);
%   INPUTS:
%       Uw      : (Nd*Nhc+1, 1)
%       Fl      : (Nd*Nhc, 1)
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
%       R       : (Nd*Nhc, 1)
%       dRdU    : (Nd*Nhc, Nd*Nhc)
%       dRdw    : (Nd*Nhc, 1)

    Nhc = sum((h==0)+2*(h~=0));
    Nd = (size(Uw,1)-1)/Nhc;

    [E, dEdw] = HARMONICSTIFFNESS(mpars.M, mpars.C, mpars.K, Uw(end), h);

    FNL = zeros(Nd*Nhc, 1);
    JNL = zeros(Nd*Nhc, 1);
    Lrels = mpars.Lrels(fpars(1));
    for i=1:size(Lrels,3)  % Iterate over Nonlinearities
        UtUn = Lrels(:,:,i)*reshape(Uw(1:end-1), Nd, Nhc);
        [Fnl, Jnl] = JENK2D(UtUn, h, Nt, fpars(2:end));

        FNL = FNL + kron(eye(Nhc), Lrels(:,:,i))'*Fnl;
        JNL = JNL + kron(eye(Nhc), Lrels(:,:,i))'*Jnl*kron(eye(Nhc), Lrels(:,:,i));
    end

    % Assemble final residue
    R    = E*Uw(1:end-1) + FNL - Fl;
    dRdU = E + JNL;
    dRdw = dEdw*Uw(1:end-1); 
end