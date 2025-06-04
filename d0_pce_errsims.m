% clc
clear all
addpath('./ROUTINES/HARMONIC/')
addpath('./ROUTINES/SOLVERS/')
addpath('./ROUTINES/QUADRATURE/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

savfig = true;
%% Model Parameters (Modified, from Erhan)
pars = struct('m', 10, 'mD', 0.1, 'c', 20, ...
              'k', 1e5, 'kxy', 3e5, 'kD', 1e6);

M = diag([repmat(pars.m, 1,4) repmat(pars.mD, 1,2)]);
K = diag([repmat(pars.k, 1,4) zeros(1,2)]);

K(1:2,1:2) = K(1:2,1:2) + 3*pars.kxy*[1 -1;-1 1];
K(3:4,3:4) = K(3:4,3:4) + 1*pars.kxy*[1 -1;-1 1];
K(5:6,5:6) = K(5:6,5:6) + pars.kD*[1 -1;-1 1];
C = diag([repmat(pars.c, 1,4) zeros(1,2)]);
F0 = [0; 0; 0; 0; 0; 120];
Fv = [0; -20; 5; 0; 0; 0];

Lrel1 = @(al) kron([-1 0 1], [cos(al) sin(al);-sin(al) cos(al)]);
Lrel2 = @(al) kron([0 -1 1], [cos(al) -sin(al);sin(al) cos(al)]);

ang_alpha = 40;  % Angle in degrees.
mpars = struct('M', M, 'C', C, 'K', K, ...
               'Lrels', @(al) cat(3, Lrel1(al), Lrel2(al)));

mi = 3; % Mode of interest

%% Initialize 
sims = {'Sim1_0011', 'Sim2_1001', 'Sim3_0101', 'Sim4_1101', 'H1Sim1_1111'};
symnam = sims{5};
Nqpce = 5;  % Number of quadrature points (per dimension)
uqpars = symnam(end-3:end);
uqparsi = find(uqpars=='1');
Nnq = sum(uqpars=='1');
pttls = {'k_t', 'k_n', '\mu', 'm'};

savdir = sprintf('./DATS/%s_Nq=%d', symnam, Nqpce);
load([savdir '/sim-data.mat'], 'uqpars', 'XIs', 'WIs', 'ktdist', 'kndist', ...
    'mudist', 'mdist', 'erris', 'h', 'Nt', 'Nhc', 'Fstat', 'Lb');
[inds0, indsh, zinds, rinds, iinds] = HINDS(6, h(:));
polfuns = {ktdist.pol,kndist.pol,mudist.pol,mdist.pol};
pdffuns = {@(x) pdf(ktdist.pdf,x,0,1), @(x) pdf(kndist.pdf,x,0,1), @(x) pdf(mudist.pdf,x,1), @(x) pdf(mdist.pdf,x,-1,1)};

germgens = @(xi, i) eval(ktdist.genfun).*(i==1) + ...
    eval(kndist.genfun).*(i==2) + ...
    eval(mudist.genfun).*(i==3) + ...
    eval(mdist.genfun).*(i==4);

XIgs = zeros(Nqpce^2, Nnq);
WIgs = zeros(Nqpce^2, 1);
for isim=1:Nqpce^2
    qinds = arrayfun(@(a) base2dec(a, Nqpce), dec2base(isim-1, Nqpce, Nnq))+1;

    XIgs(isim,:) = XIs(sub2ind([Nqpce Nnq], qinds, 1:Nnq));
    WIgs(isim) = prod(WIs(sub2ind([Nqpce Nnq], qinds, 1:Nnq)));
end

%% Harmonics Constraint (use smaller number of harmonics)
hn = (0:1);
Nhcn = sum((hn==0) + 2*(hn~=0));
Lhn = zeros(Nhcn,Nhc);
[~, ni] = find(hn(hn~=0)'==h(h~=0));

[~,~,zi,ri,ii] = HINDS(1, h(:));
[~,~,zin,rin,iin] = HINDS(1, hn(:));
Lhn(zi,zin) = 1;
Lhn(rin, ri(ni)) = eye(sum(hn~=0));
Lhn(iin, ii(ni)) = eye(sum(hn~=0));

Lhn = kron(Lhn, speye(6));

Lh0 = [Lhn zeros(Nhcn*6, 2);
       zeros(2, Nhc*6) eye(2)];
Lh0(:, rinds(5)) = [];

[~,~,zindsn,rindsn,iindsn] = HINDS(6, hn(:));
Lbn = speye(Nhcn*6+2);
Lbn(:, rindsn(5)) = [];

%%
Ast_h = 1e-1;
Aen_h = 1e2;
DSS = [0.25 0.1 0.2 0.5 0.05 0.02 0.025];  % List of 'ds' with priority to achieve convergence

errinds = find(erris);
nerris = zeros(length(errinds),1);
parsave = @(i, uxwL, xis, wis, length_scale) save(sprintf('%s/sim-%d.mat', savdir, i), 'uxwL', 'xis', 'wis', 'length_scale');
parfor ii=1:length(errinds)
    isim = errinds(ii);

    qinds = arrayfun(@(a) base2dec(a, Nqpce), dec2base(isim-1, Nqpce, Nnq))+1;
    
    xisim = XIs(sub2ind([Nqpce Nnq], qinds, 1:Nnq));
    wisim = prod(WIs(sub2ind([Nqpce Nnq], qinds, 1:Nnq)));
    %% Contact Parameters from Germs
    fpars = [deg2rad(ang_alpha); ktdist.mu; kndist.mu; mudist.mu; mdist.mu];  % Initialize all to means
    for i=1:length(xisim)
        fpars(1+uqparsi(i)) = germgens(xisim(i), uqparsi(i));
    end

    %% Set Initial Guess using linear modal analysis
    K0 = mpars.K + Lrel1(fpars(1))'*diag([fpars(2) fpars(3)])*Lrel1(fpars(1)) ...
         + Lrel2(fpars(1))'*diag([fpars(2) fpars(3)])*Lrel2(fpars(1));
    
    [V, Wst] = eig(K0, M);
    [Wst, si] = sort(sqrt(diag(Wst)));
    V = V(:,si);
    V = V./sqrt(diag(V'*M*V)');
    
    Uxw0 = zeros(6*Nhc+2,1);
    Uxw0(zinds) = K0\F0;
    %     Uxw0([rinds(1:6) iinds(1:6)]) = [V(:, mi); V(:, mi)]/sqrt(2);
    Uxw0(iinds(1:6)) = V(:,mi);
    Uxw0(end-1) = V(:,mi)'*C*V(:,mi);
    Uxw0(end) = Wst(mi);

    %% Obtain length scaling
    utunstat = [Lrel1(fpars(1))*(K0\F0) Lrel2(fpars(1))*(K0\F0)];
    utunms = [Lrel1(fpars(1))*V(:,mi) Lrel2(fpars(1))*V(:,mi)];

    fn0 = fpars(3)*utunstat(2,:);
    fnm = fpars(3)*utunms(2,:);

    ft0 = fpars(2)*utunstat(1,:);
    ftm = fpars(2)*utunms(1,:);
    ftm = ftm.*sign(ft0);
    ft0 = ft0.*sign(ft0);

    slims = fpars(4)*fpars(3)*utunstat(2,:);  % mu*kn*un0 -> fslip
    aslips = (fpars(4)*fn0)./(abs(ftm)-fpars(4)*abs(fnm));

    length_scale = min(abs(aslips));

    Ast = log10(Ast_h*length_scale);
    Aen = log10(Aen_h*length_scale);    
    
    try
        %% Conduct Simulation
        Sopt = struct('stepmax', 10000, 'jac', 'xl', 'dsmin', 0.001, ...
                      'dynamicDscale', 1, 'Display', false);
        Sopt.Dscale = [ones(6*Nhcn-1,1)*1e-4; (V(:,mi)'*C*V(:,mi)); Wst(mi); 1];        

        %         Sopt.Display = true;
        flg = false;
        Amaxs = zeros(length(DSS),1);
        for ik=1:length(DSS)
            ds = DSS(ik);
            UxwL = solve_and_continue(Lbn'*Lh0*Lb'*Uxw0, ...
                                      @(Uxwl) EPMCRESFUN(Uxwl, Lhn*Fstat, Lbn, hn, Nt, ...
                                                         fpars, mpars), ...
                                      Ast, Aen, ds, Sopt);
            UxwL = [Lh0'*Lbn*UxwL(1:end-1,:); UxwL(end,:)];

            UxwL = solve_and_continue(Lbn'*Lh0*Lb'*Uxw0, @(Uxwl) EPMCRESFUN(Uxwl, Lhn*Fstat, Lbn, hn, Nt, fpars, mpars), Ast, Aen, ds, Sopt);
            
            Amaxs(ik) = UxwL(end);
            if UxwL(end)>=Aen
                flg = true;
                break;
            end
        end
        if ~flg
            % Try in reverse
            Amins = zeros(length(DSS),1);
            for ik=1:length(DSS)
                ds = DSS(ik);
                UxwL = solve_and_continue(Lbn'*Lh0*Lb'*Uxw0, ...
                                          @(Uxwl) EPMCRESFUN(Uxwl, Lhn*Fstat, Lbn, hn, Nt, ...
                                                             fpars, mpars), ...
                                          Aen, Ast, ds, Sopt);
                UxwL = [Lh0'*Lbn*UxwL(1:end-1,:); UxwL(end,:)];
                
                Amins(ik) = UxwL(end);
                if UxwL(end)<=Ast
                    UxwL = UxwL(:, end:-1:1);
                    flg = true;
                    break;
                end
            end
        end
        if UxwL(end)<=Aen
            disp(table(fpars))
            error('isim %d didn''t converge: (%f, %f), %f\n', ...
                  isim, fpars(1+uqparsi(1)), fpars(1+uqparsi(2)), UxwL(end));
        end
        UxwL = [Lb*UxwL(1:end-1,:); UxwL(end,:)];
        nerris(ii) = 0;
        
        %% Save Results
        parsave(isim, UxwL, xisim, wisim, length_scale);
    catch me
        nerris(ii) = 1;
        disp(me.message)
    end
    
    %% 
    fprintf('Done %d/%d\n', isim, Nqpce^Nnq)    
end
erris(errinds) = nerris;

save([savdir '/sim-data.mat'], 'erris', '-append');
