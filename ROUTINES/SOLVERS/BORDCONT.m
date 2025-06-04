function [Us, efls] = BORDCONT(func, U0, Wst, Wen, ds, varargin)
%BORDCONT
%
%	USAGE:
%		[Us] = BORDCONT(func, U0, Wst, Wen, ds, Copt);
%	INPUTS:
%		func		: @(Uw)->[R, dRdU, dRdw]
%    	        U0		: (Nun, 1) initial guess
%   		Wst, Wen	: Start and end Parameters
%    		ds		: Initial Arclength
%		Copt		: Options structure with fields,
%			stepmax		: [100] maximum steps
%			dsmax		: [inf] max ds
%			dsmin		: [ds/5] min ds
%			startdir	: [sign(Wen-Wst)] starting dxn
%			itopt		: [6] Optimal number of iterations
%			stepadapt	: [true] adaptive step size
%			predictor	: ['t'] 't','tangent' or 's','secant'
%			Display		: [true]
%			dynsc		: [true] Dynamic unknown scaling
%			rsc		: [true] Dynamic residue scaling
%			brconds		: [[]] @(Uw) function handle or cell array of such
%						Return true to break.    
%	OUTPUTS:
%		Us		: (Nun+1, Npts) Unknowns
%		efls		: (1, Npts) Error flags    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVENTION                                                                            %
% Commenting Style: #.		Numbers denote sequence                           	%
%                   [ch]:	Denotes a choice made by me (possibly debatable)        %
%    		    [Q]: 	Denotes questions to be regarded                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% 0. Default options
    Copt = struct('stepmax', 100, 'dsmax', inf, 'dsmin', ds/5, ...
                  'startdir', sign(Wen-Wst), ...
                  'itopt', 6, ...
                  'stepadapt', true, ...
                  'predictor', 'tangent', ...
                  'Display', true, 'dynsc', true, ...
                  'rsc', true, ...
                  'brconds', []);
    if nargin>=6
        nflds = fieldnames(varargin{1});
        for i=1:length(nflds)
            Copt.(nflds{i}) = varargin{1}.(nflds{i});
        end
    end
    if ~isempty(Copt.brconds) && ~iscell(Copt.brconds)
        Copt.brconds = {Copt.brconds};
    end
    if nargin==7
        parm = lower(varargin{2});
        switch parm
          case {'orthogonal', 'orth', 'normal', 'riks', 'o', 'n', 'r'}
            parm = 'orthogonal';
          case {'arclength', 'arc', 'arcl', 'a' }
            parm = 'arclength';
          otherwise
            error('Unkown parameterization: %s', parm);
        end
    else
        parm = 'orthogonal';
    end
    Copt.predictor = lower(Copt.predictor);
    switch Copt.predictor
      case {'tangent', 'tgt', 't'}
        Copt.predictor = 'tangent';
      case {'secant', 'sec', 's'}
        Copt.predictor = 'secant';
      otherwise
        error('Unknown predictor: %s', Copt.predictor);
    end
    opt = optimoptions('fsolve', 'specifyobjectivegradient', true, 'Display', 'iter', ...
                       'MaxIter', 2*Copt.itopt);

    % Setup quantities to save/return
    % [ch]: Initializing as empty and growing with iterations
    Nun = size(U0,1);
    Us = [];
    alphs = [];
    efls = [];
    NITs = [];
    sel = ones(Nun+1,1);

    %% 1. Initial Point Correction
    % [ch]: No scaling used here.
    [U0, ~, efl, out] = fsolve(@(U) func([U;Wst]), U0, opt);
    Us = [Us [U0;Wst]];
    efls = [efls efl];
    NITs = [NITs out.iterations];
    opt.Display = 'off';
    
    % 1a. Tangent at Initial Point
    % [ch]: Tangent computed using unscaled Jacobians
    [~, dRdU, dRdw] = func(Us(:,end));
    z = -dRdU\dRdw;
    % Normalization for physical tangent
    alph = sqrt(1+z'*z)*Copt.startdir;
    alphs = [alphs alph];
    
    % 1b. Right preconditioning: Scaling unknowns
    if Copt.dynsc
        % Determine a selector based on absolute value of tangent
        % [ch]: Use the tangent in physical (unscaled) space for this
        sel = zeros(Nun+1,1);
        [~, si] = sort(abs(z), 'descend');
        sid = si;
        sel([si(1); end]) = 1.0;

        % Scale selected unknowns based on solution
        % Don't scale unselected unknowns
        dsc = sel.*abs(Us(:,end)) + (~sel);
    else
        dsc = ones(Nun+1,1);
    end
    % 1c. Left preconditioning: Scaling residue
    if Copt.rsc
        % Choose residue scaling to minimize the condition number of Jacobian @ soln
        % [Q]: Is there a better strategy? Just 1/diag?
        rsc = 1./mean(abs([dRdU dRdw].*dsc'),2);
        rsc(rsc==0) = 1.0;
    else
        rsc = ones(size(dRdU,1),1);
    end

    % 1d. Express tangents in scaled space
    zhat = dsc(end).*(z./dsc(1:end-1));
    % Use only selected tangent elements for normalization
    alphhat = sqrt(1+zhat'*(sel(1:end-1).*zhat))*sign(alph);

    % 1e. Precompute "optimal" quantities for step length adaptation
    if Copt.stepadapt
        % [ch]: It is assumed that the user wants the spacing between
        % all points on solution curve resembling that at the starting point.
        % [Q]: This may be too restrictive. Can we update the optimal quantities?

        % Start with tangent initial guess
        Uguess = Us(:,end)./dsc + ds/alphhat*[zhat;1];
        [R0, dRdUw0] = BORDRES(func, parm, Uguess, Us(:,end)./dsc, ...
                               zhat, alphhat, ds, dsc, rsc);
        U1 = Uguess - dRdUw0\R0;
        R1 = BORDRES(func, parm, U1, Us(:,end)./dsc, ...
                     zhat, alphhat, ds, dsc, rsc);
        
        kap0 = norm(dRdUw0\R1)./norm(dRdUw0\R0);  % contraction rate
        del0 = norm(dRdUw0\R0);  % norm of step
    end
    kap = kap0;
    while true
        %% 2. Predict & Correct Point
        % 2a. Make prediction in scaled domain
        switch Copt.predictor
          case 'tangent'
            Uguess = Us(:,end)./dsc+ds/alphhat*[zhat;1];
          case 'secant'
            if size(Us,2)>1
                zc = (Us(1:end-1,end)-Us(1:end-1,end-1))/(Us(end,end)-Us(end,end-1));
                zhatc = dsc(end).*(zc./dsc(1:end-1));
                alphhatc = sqrt(1+zhatc'*zhatc)*sign(alphhat);
            else
                zhatc = zhat;
                alphhatc = alphhat;
            end
            Uguess = Us(:,end)./dsc+ds/alphhatc*[zhatc;1];
          otherwise
            error('Unknown predictor %s', Copt.predictor);
        end

        % 2b. Iteratively Correct using Bordered Residue
        [U0, ~, efl, out] = fsolve(@(Uw) BORDRES(func, parm, Uw, Us(:,end)./dsc, ...
                                                 zhat, alphhat, ds, dsc, rsc), ...
                                   Uguess, opt);
        % 2b.i. Nonconvergence Handling
        if efl<=0
            % [ch]: If unconverged, first try to reduce step size.
            % [ch]: If still unconverged after step size==dsmin,
            % 		then increase number of iterations.
            opt.MaxIter = 20;
            if ds~=Copt.dsmin
                ds = max(ds/2, Copt.dsmin);
                continue;
            else
                opt.MaxIter = opt.MaxIter*2;
                for i=1:2
                    % [ch]: Try solving for more iterations (everything else same)
                    [U0, ~, efl, out] = fsolve(@(Uw) BORDRES(func, parm, Uw, ...
                                                             Us(:,end)./dsc, zhat, ...
                                                             alphhat, ds, dsc, rsc), ...
                                               Uguess, opt);
                    if efl>0
                        break;
                    end
                    % [ch]: Try solving without left preconditioning
                    [U0, ~, efl, out] = fsolve(@(Uw) BORDRES(func, parm, Uw, ...
                                                             Us(:,end)./dsc, zhat, ...
                                                             alphhat, ds, dsc, 1.0), ...
                                               Uguess, opt);
                    if efl>0
                        break;
                    end
                    
                    opt.MaxIter = 400;
                end
                if efl<=0
                    warning('No Convergence (%d). Exiting.', efl)
                    break;
                end
                opt.MaxIter = 2*Copt.itopt;
            end
        end
        % 2b.ii. [ch]: Conduct one Newton iteration if fsolve refuses to do any
        if out.iterations==0
            [R, dRdUw] = BORDRES(func, parm, Uguess, Us(:,end)./dsc, ...
                                 zhat, alphhat, ds, dsc, rsc);
            U0 = Uguess - dRdUw\R;
            out.iterations = out.iterations+1;
        end

        % 2c. Accept solution, append to global list.
        Us = [Us dsc.*U0];
        efls = [efls efl];
        NITs = [NITs out.iterations];

        % 2d. Check for breaking conditions
        if ~isempty(Copt.brconds)
            bflag = (sum(cellfun(@(c) c(Us(:,end)), Copt.brconds))~=0);
            if bflag
                break;
            end
        end

        % 2d. Estimate new tangent
        % [ch]: Do this in physical (unscaled) space
        [~, dRdU, dRdw] = func(Us(:,end));
        zn = -dRdU\dRdw;
        % Normalize tangent
        alph = sqrt(1+zn'*zn)*sign((1+z'*zn)/alph);
        alphs = [alphs alph];

        % 2d.i. Estimate new secant
        zd = Us(:,end)-Us(:,end-1);
        if zd(end)==0
            zd = zd(1:end-1);
        else
            zd = zd(1:end-1)/zd(end);
        end

        % 2e. Scaling the unknowns (Right preconditioning)
        if Copt.dynsc
            dscp = dsc; % Save previous scaling

            % Determine selector vector
            sel = zeros(Nun+1,1);
            sip = si;
            sidp = sid;
            [~, si] = sort(abs(zn), 'descend');
            [~, sid] = sort(abs(zd), 'descend');
            % [ch]: Use top value indices from current and previous tangents & secants
            % [ch]: Let the parameter always be selected.
            sel(unique([si(1); sid(1); sip(1); sidp(1); end])) = 1.0;

            % Determine Scaling.
            % [ch]: Scaling is 1 for all the unselected ones
            dsc = sel.*abs(Us(:,end)) + (~sel);
        else
            dscp = dsc;
            dsc = ones(Nun+1,1);
        end
        % 2f. Left preconditioning: Scaling residue
        if Copt.rsc
            % Choose residue scaling to minimize the condition number of Jacobian @ soln
            % [Q]: Is there a better strategy? Just 1/diag?
            rsc = 1./mean(abs([dRdU dRdw].*dsc'),2);
            rsc(rsc==0) = 1.0;
        else
            rsc = ones(size(dRdU,1),1);
        end

        % 2g. Express tangents & secants in scaled space
        zhat  = dsc(end).*(zn./dsc(1:end-1));  % new tangent with new scaling
        zhatd = dsc(end).*(zd./dsc(1:end-1));  % new secant with new scaling
        zhatp = dsc(end).*(z./dsc(1:end-1));   % old tangent with new scaling

        % [ch]: Normalize scaled tangent with direction fixed from tangents scaled
        % previous scaling. I found this to be able to robustly detect correct
        % turning points better.
        zhh  = dscp(end).*(zn./dscp(1:end-1)); % new tangent with old scaling
        zhhp = dscp(end).*(z./dscp(1:end-1));  % old tangent with old scaling
        
        alphhatp = alphhat;
        alphhat = sqrt(1+zhat'*(sel(1:end-1).*zhat))*sign((1+zhh'*zhhp)/alphhatp);
        % Update previous tangent vector 
        z = zn;
        
        % 3. Check if whole path has been traced & break if necessary
        if Us(end)*sign(Wen-Wst)>Wen*sign(Wen-Wst) || size(Us, 2)>Copt.stepmax
            break;
        end

        % 3a. Step size adaptation
        if Copt.stepadapt
            % 3a.i. Determine metrics we would expect if we did Newton with current
            % step size ds. (See Sec. 6.1 in Allgower & George 2003 for definitions)
            Uguess = Us(:,end)./dsc + ds/alphhat*[zhat;1]; % Tangent initial guess
            [R0, dRdUw0] = BORDRES(func, parm, Uguess, Us(:,end)./dsc, ...
                                   zhat, alphhat, ds, dsc, rsc);
            U1 = Uguess - dRdUw0\R0;  % Newton update
            R1 = BORDRES(func, parm, U1, Us(:,end)./dsc, ...
                         zhat, alphhat, ds, dsc, rsc);
            
            kap = norm(dRdUw0\R1)./norm(dRdUw0\R0); % Contraction rate
            del = norm(dRdUw0\R0);  % norm of update
            % [Q]: Is it appropriate to compute angle this way?
            % [ch]: I'm using the scaled angle
            ang = acos((zhatp'*(sel(1:end-1).*zhat)+1)/...
                       (norm(sel.*[zhatp;1])*norm(sel.*[zhat;1])));  % angle b/w tgts


            % [Q] Optimal quantities: how to best set these???
            kapopt = kap0; % optimal contraction rate
            delopt = del0; % optimal norm of update
            angopt = 1e-1; % optimal angle

            xi1 = sqrt(kapopt/kap);  		% 1. contraction rate criterion
            xi2 = Copt.itopt/out.iterations;	% 2. convergence iterations criterion
            xi3 = sqrt(delopt/del);  		% 3. update norm criterion
            xi4 = angopt/ang;  			% 4. tangent angle criterion

            % 3a.ii. Synthesize [xi1, xi2, xi3, xi4] to obtain xi.
            % [Q]: How to synthesize [x1 x2 x3 x4] to give recommended xi?
            % The textbook recommends minimum of all.
            % But this makes the code too conservative.
            xis = [xi1 xi2 xi3 xi4];
            % [ch]: After trial & error I found using the geometric mean of
            % criterions 1 & 2 (contract rate & iterations) perform well.
            xi = sqrt(prod(xis(1:2)));

            % 3a.iii. Update step size
            ds = min(max(ds*xi, Copt.dsmin), Copt.dsmax);
        end

        % 4. Display Progress if necessary
        if Copt.Display
            fprintf('(%d) %f : %d (%d, %f), %f: %f\n', size(Us,2), Us(end), ...
                    out.iterations, efl, kap, ds, alphhat);
        end            
    end
end

%% Bordered Residue
function [R, dRdUw] = BORDRES(func, parm, Uw, Uw0, z, alph, ds, dsc, varargin)
%BORDCONT Return the bordered residue function.
%
%	USAGE:
%		[R, dRdUw] = BORDCONT(func, parm, Uw, Uw0, z, alph, ds, dsc, varargin);
%	INPUTS:
%		func, parm, Uw, Uw0, z, alph, ds, dsc, varargin
%    			Optional argument allows you to scale the residue.
%	OUTPUTS:
%		R, dRdUw    
    % Residue of new point
    [R, dRdU, dRdw] = func(dsc.*Uw);

    switch parm
      case 'orthogonal'
        % Orthogonal Parameterization
        R = [R; ...
             [z;1]'*(Uw-Uw0)-ds*alph];
        dRdUw = [[dRdU dRdw].*dsc';
                 [z;1]'];
      case 'arclength'
        % Arc-Length Parameterization
        R = [R;
             (Uw-Uw0)'*(Uw-Uw0)-ds^2];
        dRdUw = [[dRdU dRdw].*dsc';
                 2*(Uw-Uw0)'];
      otherwise
        error('Unkown parameterization: %s', parm);
    end
    if ~isempty(varargin)
        try
        R(1:end-1) = varargin{1}.*R(1:end-1);
        dRdUw(1:end-1,:) = varargin{1}.*dRdUw(1:end-1,:);
        catch
            keyboard
        end
    end
end

%% Fsolve Exitflag meaning
%  1  FSOLVE converged to a root.
%  2  Change in X too small.
%  3  Change in residual norm too small.
%  4  Computed search direction too small.
%  0  Too many function evaluations or iterations.
% -1  Stopped by output/plot function.
% -2  Converged to a point that is not a root.
% -3  Trust region radius too small (Trust-region-dogleg).
