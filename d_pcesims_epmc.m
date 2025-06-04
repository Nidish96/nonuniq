%% * Preamble
distcomp.feature( 'LocalUseMpiexec', false )
% clc
clear all
addpath('./ROUTINES/HARMONIC/')
addpath('./ROUTINES/SOLVERS/')
addpath('./ROUTINES/QUADRATURE/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

showfigs = false;

%% * Model Parameters (Modified, from Erhan)
pars = struct('m', 10, 'mD', 0.1, 'c', 20, ...
              'k', 1e5, 'kxy', 3e5, 'kD', 1e6);

M = diag([repmat(pars.m, 1,4) repmat(pars.mD, 1,2)]); 
K = diag([repmat(pars.k, 1,4) zeros(1,2)]);

K(1:2,1:2) = K(1:2,1:2) + 3*pars.kxy*[1 -1;-1 1];
K(3:4,3:4) = K(3:4,3:4) + 1*pars.kxy*[1 -1;-1 1];
K(5:6,5:6) = K(5:6,5:6) + pars.kD*[1 -1;-1 1];
C = diag([repmat(pars.c, 1,4) zeros(1,2)]);
F0 = [0; 0; 0; 0; 0; 120]*2;
Fv = [0; -20; 5; 0; 0; 0];

Lrel1 = @(al) kron([-1 0 1], [cos(al) sin(al);-sin(al) cos(al)]);
Lrel2 = @(al) kron([0 -1 1], [cos(al) -sin(al);sin(al) cos(al)]);

ang_alpha = 40;  % Angle in degrees.
mpars = struct('M', M, 'C', C, 'K', K, ...
               'Lrels', @(al) cat(3, Lrel1(al), Lrel2(al)));

mi = 3; % Mode of interest

%% * Setup HB
h = 0:5;
Nhc = sum((h==0)+2*(h~=0));
[inds0, indsh, zinds, rinds, iinds] = HINDS(6, h(:));
Nt = 256;

Fstat = zeros(6*Nhc,1);
Fstat(zinds(1:6)) = F0;

% Construct matrix for explicitly enforcing phase constraint
Lb = speye(Nhc*6+2);  Lb(:, rinds(5)) = [];

Ast_h = 1e-2;
Aen_h = 1e2;

DSS = [0.25 0.1 0.2 0.5 0.05 0.02 0.06 0.025];  % List of 'ds' with priority to achieve convergence

%% * Setup PCE Germs & Weights
% (see ./ROUTINES/QUADRATURE/quadcheck.m for quadrature rules for weighted integration)
Nqpce = 5;  % Number of quadrature points (per dimension)
logsck = true;

% 1. Tangential Stiffness
kmu = 3e5;
ksg = kmu*.3;
if logsck
    lkmu = log(kmu^2/sqrt(kmu^2+ksg^2));
    lksg = sqrt(log(1+ksg^2/kmu^2));

    kmu = lkmu;
    ksg = lksg;
end

ktdist  = struct('mu', kmu, 'sig', ksg, 'quad', 'he_prob', ...
                 'genfun', 'ktdist.mu+ktdist.sig*xi', ...
                 'lims', [], ...
                 'pol', @(n,x) PHERM(n,x), 'pdf', 'norm');                       % Gaussian (Probabilist's Hermite Polynomial)
                                                                                 % 2. Normal Stiffness
kndist  = struct('mu', kmu, 'sig', ksg, 'quad', 'he_prob', ...
                 'genfun', 'kndist.mu+kndist.sig*xi', ...
                 'lims', [], ...
                 'pol', @(n,x) PHERM(n,x), 'pdf', 'norm');                       % Gaussian (Probabilist's Hermite Polynomial)
                                                                                 % 3. Coefficient of Friction
                                                                                 % mudist  = struct('mu', 0.2, 'quad', 'lague', ...
                                                                                 %                  'genfun', 'mudist.mu*xi', ...
                                                                                 %                  'lims', [], ...
                                                                                 %                  'pol', @(n,x) PLAGU(n,x), 'pdf', 'exp');                        % Exponential (Laguerre Polynomial)
mudist = struct('ab', [0 1.25], 'mu', [], 'quad', 'legen', ...
                'genfun', 'mudist.ab(1)+(xi+1)*diff(mudist.ab)/2', ...
                'lims', [0 1.25], ...
                'pol', @(n,x) PLEGE(n,x), 'pdf', 'unif');
mudist.mu = mean(mudist.ab);
% mudist = struct('ab', [0 1], 'mu', 0.5, 'quad', 'legen', ...
%                 'genfun', 'mudist.ab(1)+(xi+1)*diff(mudist.ab)/2', ...
%                 'pol', @(n,x) PLEGE(n,x), 'pdf', 'unif');                     % Uniform (Legendre (preferred) or Lobatto Polynomial)
% 4. Contact Multiplier
% mdist   = struct('ab', [-1 1], 'mu', 0, 'quad', 'legen', ...
%                  'genfun', 'mdist.ab(1)+(xi+1)*diff(mdist.ab)/2', ...
%                  'pol', @(n,x) PLEGE(n,x), 'pdf', 'unif');                       % Uniform (Legendre (preferred) or Lobatto Polynomial)
mdist   = struct('alpha', .5, 'beta', .5, 'quad', 'jacobi', ...
                 'mu', [], 'lims', [-1 1], ...
                 'genfun', '2*xi-1', ...
                 'pol', [], 'pdf', 'beta');                       % Beta (Jacobi Polynomial)

% beta pars: [.5, .5]. [3.5, 3.5]. [2, 5]. [5, 2]. [.2, .2]. [.5, .2]. [.2, .5]. [1, 1]
mdist.mu = mdist.alpha/(mdist.alpha+mdist.beta);
mdist.mu = 2*mdist.mu-1;
mdist.pol = @(n,x) PJACO(n, mdist.alpha,mdist.beta, x);

uqpars = '1111';  % kt, kn, mu, m
% mudist.mu = 0.5;
uqparsi = find(uqpars=='1');
Nnq = sum(uqpars=='1');

germgens = @(xi, i) eval(ktdist.genfun).*(i==1) + ...
    eval(kndist.genfun).*(i==2) + ...
    eval(mudist.genfun).*(i==3) + ...
    eval(mdist.genfun).*(i==4);

% Generate Germs & Weights (points on standard domains)
[xi1, wi1] = GAUSS_QUAD(Nqpce, ktdist.quad);
[xi2, wi2] = GAUSS_QUAD(Nqpce, kndist.quad);
if ~strcmp(mudist.quad, 'jacobi')
    [xi3, wi3] = GAUSS_QUAD(Nqpce, mudist.quad);
else
    [xi3, wi3] = GAUSS_QUAD(Nqpce, mudist.quad, mudist.alpha,mudist.beta);
end
if ~strcmp(mdist.quad, 'jacobi')
    [xi4, wi4] = GAUSS_QUAD(Nqpce, mdist.quad);
else
    [xi4, wi4] = GAUSS_QUAD(Nqpce, mdist.quad, mdist.alpha, mdist.beta);
end

% Germs for [kt kn muN m]
XIs = [xi1(:) xi2(:) xi3(:) xi4(:)];
% Corresponding quadrature weights
WIs = [wi1(:) wi2(:) wi3(:) wi4(:)];

% Use only those required for the uq propagation
XIs = XIs(:, uqparsi);
WIs = WIs(:, uqparsi);

% Create Directory for data
[~,sysnam]=system('hostname');
sysnam = deblank(sysnam);
savdir = ['./DATS/' sysnam '_' uqpars '_' datestr(datetime, 'dd-mmm_at_HH:MM') '_Nq=' num2str(Nqpce)];
mkdir(savdir);
disp(['Saving to ' savdir])

save([savdir '/sim-data.mat'], 'uqpars', 'XIs', 'WIs', 'ktdist', 'kndist', ...
     'mudist', 'mdist', 'h', 'Nt', 'Nhc', 'Fstat', 'Lb', 'logsck');
parsave = @(i, uxwL, xis, wis, length_scale) ...
          save(sprintf('%s/sim-%d.mat', savdir, i), 'uxwL', 'xis', 'wis', ...
               'length_scale');

%% * Mean Model Simulation
xisim = [eval(['fzero(@(xi)' ktdist.genfun '-' num2str(ktdist.mu) ',0)']),
         eval(['fzero(@(xi)' kndist.genfun '-' num2str(kndist.mu) ',0)']),
         eval(['fzero(@(xi)' mudist.genfun '-' num2str(mudist.mu) ',0)']),
         eval(['fzero(@(xi)' mdist.genfun '-' num2str(mdist.mu) ',0)'])];
disp('hey');
wisim = [nan nan nan nan];
fpars = [deg2rad(ang_alpha); ktdist.mu; kndist.mu; mudist.mu; mdist.mu];
if logsck
    fpars(2:3) = exp(fpars(2:3));
end
% Set Initial Guess using linear modal analysis
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

% Obtain length scaling
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
try % Conduct Simulation
    Sopt = struct('stepmax', 10000, 'jac', 'xl', 'dsmin', 0.001, ...
                  'dynamicDscale', 1, 'Display', true);
    % -> Dscale can also be used to weight the unknowns. Set Asc = 0.1 to get better resolution in the NMA.
    Sopt.Dscale = [ones(6*Nhc-1,1)*1e-4; (V(:,mi)'*C*V(:,mi)); Wst(mi); 1];

    flg = false;
    % Try forward continuation (preferred)
    Amaxs = zeros(length(DSS),1);
    for ik=1:length(DSS)
        ds = DSS(ik);
        UxwL = solve_and_continue(Lb'*Uxw0, ...
                                  @(Uxwl) EPMCRESFUN(Uxwl, Fstat, Lb, h, Nt, ...
                                                     fpars, mpars), ...
                                  Ast, Aen, ds, Sopt);
        Amaxs(ik) = UxwL(end);
        if UxwL(end)>=Aen
            flg = true;
            break;
        end
    end
    % if ~flg % Try in reverse
    %     Amins = zeros(length(DSS),1);
    %     for ik=1:length(DSS)
    %         ds = DSS(ik);
    %         UxwL = solve_and_continue(Lb'*Uxw0, ...
    %                                   @(Uxwl) EPMCRESFUN(Uxwl, Fstat, Lb, h, Nt, ...
    %                                                      fpars, mpars), ...
    %                                   Aen, Ast, ds, Sopt);
    %         Amins(ik) = UxwL(end);
    %         if UxwL(end)<=Ast
    %             UxwL = UxwL(:, end:-1:1);
    %             flg = true;
    %             break;
    %         end  
    %     end
    % end
    if ~flg  % Continuation broke
        fprintf('isim %d didn''t converge:\n', isim);
        disp(table(fpars))
        error('isim %d didn''t converge:\n', isim);
        % error('isim %d didn''t converge: (%f, %f), %f\n', ...
        %       isim, fpars(1+uqparsi(1)), fpars(1+uqparsi(2)), UxwL(end));
    end
    UxwL = [Lb*UxwL(1:end-1,:); UxwL(end,:)];

    % Save Results
    parsave(0, UxwL, xisim, wisim, length_scale);
catch me
    % erris(isim) = 1;
    parsave(0, [], xisim, wisim, length_scale);
    disp(me.message)
end

%% * Initiate Simulation seeded by isim 
% isim represents the cartesian index corresponding to the quadrature point
% in 4D space.

parsave = @(i, uxwL, xis, wis, length_scale) save(sprintf('%s/sim-%d.mat', savdir, i), 'uxwL', 'xis', 'wis', 'length_scale');
erris = zeros(Nqpce^Nnq,1);

w0s = zeros(1, Nqpce^Nnq);
dw0s = zeros(1, Nqpce^Nnq);
parfor isim=1:Nqpce^Nnq
    qinds = arrayfun(@(a) base2dec(a, Nqpce), dec2base(isim-1, Nqpce, Nnq))+1;
    
    xisim = XIs(sub2ind([Nqpce Nnq], qinds, 1:Nnq));
    wisim = prod(WIs(sub2ind([Nqpce Nnq], qinds, 1:Nnq)));
    % Contact Parameters from Germs
    fpars = [deg2rad(ang_alpha); ktdist.mu; kndist.mu; mudist.mu; mdist.mu];  % Initialize all to means
    for i=1:length(xisim)
        fpars(1+uqparsi(i)) = germgens(xisim(i), uqparsi(i));
    end
    if logsck
        fpars(2:3) = exp(fpars(2:3));
    end

    % Set Initial Guess using linear modal analysis
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

    % Initial Guess from mean model
    prv=load([savdir '/sim-0.mat'], 'uxwL');
    Uxw0 = prv.uxwL(1:end-1,1);

    % Obtain length scaling
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
    
    try % Conduct Simulation
        Sopt = struct('stepmax', 10000, 'jac', 'xl', 'dsmin', 0.001, ...
                      'dynamicDscale', 1, 'Display', false);
        % -> Dscale can also be used to weight the unknowns. Set Asc = 0.1 to get better resolution in the NMA.
        Sopt.Dscale = [ones(6*Nhc-1,1)*1e-4; (V(:,mi)'*C*V(:,mi)); Wst(mi); 1];

        flg = false;
        % Try forward continuation (preferred)
        Amaxs = zeros(length(DSS),1);
        for ik=1:length(DSS)
            ds = DSS(ik);
            UxwL = solve_and_continue(Lb'*Uxw0, ...
                                      @(Uxwl) EPMCRESFUN(Uxwl, Fstat, Lb, h, Nt, ...
                                                         fpars, mpars), ...
                                      Ast, Aen, ds, Sopt);
            Amaxs(ik) = UxwL(end);
            if UxwL(end)>=Aen
                flg = true;
                break;
            end
        end
        if ~flg % Try in reverse
            Uxw0 = prv.uxwL(1:end-1,end);
            Amins = zeros(length(DSS),1);
            for ik=1:length(DSS)
                ds = DSS(ik);
                UxwL = solve_and_continue(Lb'*Uxw0, ...
                                          @(Uxwl) EPMCRESFUN(Uxwl, Fstat, Lb, h, Nt, ...
                                                             fpars, mpars), ...
                                          Aen, Ast, ds, Sopt);
                Amins(ik) = UxwL(end);
                if UxwL(end)<=Ast
                    UxwL = UxwL(:, end:-1:1);
                    flg = true;
                    break;
                end  
            end
        end
        if ~flg  % Continuation broke
            fprintf('isim %d didn''t converge:\n', isim);
            disp(table(fpars))
            error('isim %d didn''t converge:\n', isim);
            % error('isim %d didn''t converge: (%f, %f), %f\n', ...
            %       isim, fpars(1+uqparsi(1)), fpars(1+uqparsi(2)), UxwL(end));
        end
        UxwL = [Lb*UxwL(1:end-1,:); UxwL(end,:)];

        % Save Results
        dw0s(isim) = diff(UxwL(end-1, 1:2));
        w0s(isim) = UxwL(end-1,1);
        parsave(isim, UxwL, xisim, wisim, length_scale);
    catch me
        erris(isim) = 1;
        parsave(isim, [], xisim, wisim, length_scale);
        disp(me.message)
    end
    fprintf('Done %d/%d\n', isim, Nqpce^Nnq)
end
iols = find(~isoutlier(dw0s, 'mean'));
% Do the second step only if necessary
% iols = [setdiff(1:length(dw0s), iols) iols(find(isoutlier(dw0s(iols), 'mean')))];
iols = setdiff(1:length(dw0s), iols);

iols2 = find(isoutlier(w0s));

% iols = unique([iols iols2]);
% iols = iols2;
iols = intersect(iols, iols2);

erris(iols) = 1;  % Spurious run
save([savdir '/sim-data.mat'], 'erris', '-append');

%% * Impute Missing Data (assuming there's nothing contiguous that is missing)
parsave = @(i, uxwL, xis, wis, length_scale) save(sprintf('%s/sim-%d.mat', savdir, i), 'uxwL', 'xis', 'wis', 'length_scale');

load([savdir '/sim-data.mat'], 'erris');

errinds = find(erris);
imputedis = [];
for isim=errinds(:)'
    fprintf('Imputing %d\n', isim);
    
    % Get Current Point Features
    qinds = arrayfun(@(a) base2dec(a, Nqpce), dec2base(isim-1, Nqpce, Nnq))+1;
    load(sprintf('%s/sim-%d.mat', savdir, isim), 'uxwL', 'xis', 'wis', 'length_scale');
    save(sprintf('%s/sim-%d_o.mat', savdir, isim), 'uxwL', 'xis', 'wis', 'length_scale');
    
    % Get (Valid) Neighboring indices
    qinds_ngb = [qinds+eye(Nnq); qinds-eye(Nnq)];
    qinds_ngb = qinds_ngb(all(qinds_ngb~=0 & qinds_ngb<=Nqpce, 2), :);
    isim_ngb = base2dec(int2str((qinds_ngb-1)*(10.^(Nnq-1:-1:0)')), Nqpce)+1;
    isim_ngb = setdiff(isim_ngb, errinds);
    
    % Load data from neighbors
    ngbs = arrayfun(@(i) load(sprintf('%s/sim-%d.mat', savdir, i), ...
                              'uxwL', 'xis', 'wis', 'length_scale'), isim_ngb);
    Napts = fix(mean(arrayfun(@(n) size(n.uxwL,2), ngbs))*1.1);
    ass = logspace(log10(Ast_h), log10(Aen_h), Napts);
    uxwL = zeros(size(ngbs(1).uxwL,1)-1, Napts);
    for i=1:length(isim_ngb)
        uxwL = uxwL + pchip(ngbs(i).uxwL(end,:)-log10(ngbs(i).length_scale), ...
                            ngbs(i).uxwL(1:end-1,:), log10(ass));
    end
    uxwL = uxwL/length(isim_ngb);
    uxwL = [uxwL; log10(ass*length_scale)];
    
    % Save average of data
    save(sprintf('%s/sim-%d.mat', savdir, isim), 'uxwL', 'xis', 'wis', 'length_scale');
    imputedis = [imputedis; isim];

    errinds = errinds(2:end);
end
save([savdir '/sim-data.mat'], 'imputedis', '-append');

%% * Plotting to check
if showfigs
    load([savdir '/sim-data.mat'], 'erris');
    
    figure(20)
    clf()
    set(gcf, 'Color', 'white')
    % for isim=find(~erris)'
    for isim=1:Nqpce^Nnq
	load(sprintf('%s/sim-%d.mat', savdir, isim), 'uxwL', 'xis', 'wis', 'length_scale');
	UxwL = uxwL;
        w0s(isim) = UxwL(end-1, 1);
        dw0s(isim) = diff(UxwL(end-1, 1:2));
	
	subplot(2,1,1)
	semilogx(10.^UxwL(end,:)/length_scale, UxwL(end-1,:)/2/pi, '-'); hold on
        % plot(UxwL(end-1,:)/2/pi,'.-'); hold on
	
	subplot(2,1,2)
	semilogx(10.^UxwL(end,:)/length_scale, ...
                 UxwL(end-2,:)./(2*UxwL(end-1,:))*100, '-'); hold on
	% drawnow
    end
    subplot(2,1,1)
    ylabel('Natural Frequency (Hz)')
    subplot(2,1,2)
    ylabel('Damping Factor (\%)')
    xlabel('Modal Amplitude')
    
    % answer = questdlg('Delete Data?');
    % if strcmp(answer, 'Yes')
    %     rmdir(savdir, 's');
    % end
end
