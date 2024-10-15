classdef EnGMF < datools.filter.ensemble.EnF
    % Ensemble Gaussian Mixture Filter
    % BRUF ---DOI: https://doi.org/10.23919/FUSION52260.2023.10224134
    % IEKF ---DOI: https://doi.org/10.1109/9.250476
    % UKF  ---Book: Simo Sarkka (2013). Bayesian Filtering and Smoothing. Cambridge University Press.
    
    properties
        BandwidthScale
        BRUFSteps
        GMMUpdateType
        UseRobustSampling
        SamplingType
        GMMConstruction
        UseAdaptiveEnsembleSize = false
        UseAKDE

        Name = 'Ensemble Gaussian Mixture Filter'
    end

    properties (Access = private)
        QMCStream
    end

    methods
        function obj = EnGMF(varargin)

            p = inputParser;
            p.KeepUnmatched = true;

            % this is recommended according to the Reich book

            addParameter(p, 'BandwidthScale', 1);
            addParameter(p, 'BRUFSteps', 1);
            addParameter(p, 'GMMUpdateType', 'EKF');
            addParameter(p, 'UseRobustSampling', false);
            addParameter(p, 'SamplingType', "Gaussian");
            addParameter(p, 'GMMConstruction', 'KDE');
            addParameter(p, 'UseAKDE', false);

            parse(p, varargin{2:end});

            s = p.Results;

            kept = p.Unmatched;

            obj@datools.filter.ensemble.EnF(varargin{1}, kept);

            obj.BandwidthScale = s.BandwidthScale;
            obj.BRUFSteps = s.BRUFSteps;
            obj.GMMUpdateType = s.GMMUpdateType;
            obj.UseRobustSampling = s.UseRobustSampling;
            obj.SamplingType = s.SamplingType;
            obj.GMMConstruction = s.GMMConstruction;
            obj.UseAKDE = s.UseAKDE;

        end
    end

    methods

        function analysis(obj, obs)

            obsgmm = obs.Uncertainty.asGMM();

            ys = obsgmm.Mu;
            Rs = obsgmm.Sigma;
            wobs = obsgmm.W;

            ensM = size(ys, 2);

            Xb = obj.Ensemble;

            % inflation
            Xbm = mean(Xb, 2);
            Xb = obj.Inflation*(Xb - Xbm) + Xbm;

            ensN = obj.NumEnsemble;
            n = size(Xb, 1);

            if obj.SamplingType == "QMC" && isempty(obj.QMCStream)
                obj.QMCStream = qrandstream('sobol', n + 1, 'Skip', 100);
            end

            beta = obj.BandwidthScale*((4/(n + 2))^(2/(n + 4)))*((ensN)^(-2/(n + 4)));

            if iscell(obj.GMMConstruction)
                gmmconst = obj.GMMConstruction{1};
                % MATLAB does not support cell array slicing for some
                % strange reason, thus we have to do it ourselves
                gmmconstparams = cell(1, (numel(obj.GMMConstruction) - 1));
                for i = 1:(numel(obj.GMMConstruction) - 1)
                    gmmconstparams{i} = obj.GMMConstruction{i + 1};
                end
            else
                gmmconst = obj.GMMConstruction;
                gmmconstparams = {};
            end


            switch gmmconst
                case 'ELocalization'

                    if isempty(gmmconstparams)
                        radiusscale = 1;
                        alpha = 0;
                        epsilon1 = -inf;
                        epsilon2 = -inf;
                    else

                        nparams = numel(gmmconstparams);
                        if nparams < 1
                            radiusscale = 1;
                        else
                            radiusscale = gmmconstparams{1};
                        end
                        if nparams < 2
                            alpha = 0;
                        else
                            alpha = gmmconstparams{2};
                        end
                        if nparams < 3
                            epsilon1 = -inf;
                        else
                            epsilon1 = gmmconstparams{3};
                        end
                        if nparams < 4
                            epsilon2 = -inf;
                        else
                            epsilon2 = gmmconstparams{4};
                        end
                    end
                    % calculate the distances
                    dist = zeros(ensN, ensN);
                    wD = zeros(ensN, ensN);
                    rs = zeros(1, 1, ensN);

                    for i = 1:ensN
                        for j = 1:ensN
                            dist(i, j) = norm(Xb(:, i) - Xb(:, j));
                        end
                        [~, I] = sort(dist(i, :));
                        r = radiusscale*dist(i, I(round(sqrt(ensN))));
                        rs(1, 1, i) = r;

                        wD(i, :) = exp(-0.5*(dist(i, :)/r).^2);
                    end
                    wD = wD./sum(wD, 2);

                    wD = (1 - alpha)*wD + alpha*ones(ensN)/ensN;

                    Btilde = zeros(n, n, ensN);
                    Xtilde = zeros(n, ensN);

                    for i = 1:ensN
                        wDi = wD(i, :).';
                        Xbm = Xb*wDi;
                        c = 1/(1 - sum(wDi.^2));
                        %c
                        Pbi = c*((Xb - Xbm)*(diag(wDi))*(Xb - Xbm).');

                        Si = (rs(:, :, i).^2)*eye(n);

                        A = (Si - Pbi);

                        [U, S] = eig(A);
                        U = real(U);
                        S = real(S);

                        s = diag(S);
                        %s
                        s(s < epsilon2) = epsilon2;
                        %A = U*diag(s)*U.';
                        %Sfrak = Pbi*(A\Si);

                        Sfrak = Pbi*(U*(s.\(U.'*Si)) );

                        Sfrak = (Sfrak + Sfrak.')/2;

                        [U, S] = eig(Sfrak);
                        U = real(U);
                        S = real(S);
                        s = diag(S);
                        %s
                        s(s < epsilon1) = epsilon1;
                        Sfrak = U*diag(s)*U.';

                        Btilde(:, :, i) = beta*Sfrak;

                        Xtilde(:, i) = Xb(:, i);

                    end
                case 'FitGMM'
                    gmm = [];
                    bestBIC = inf;

                    if isempty(gmmconstparams)
                        fullncomp = 1:(ensN - 1);
                    else
                        fullncomp = gmmconstparams{1};
                    end

                    ncomp = fullncomp(end);
                    for tmpcomp = fullncomp
                        tmpgmm = fitgmdist(Xb.', tmpcomp,...
                            "RegularizationValue", 1e-5, ...
                            'Replicates', 10, ...
                            'ProbabilityTolerance', 0);
                        if tmpgmm.BIC < bestBIC
                            gmm = tmpgmm;
                            ncomp = tmpcomp;
                            bestBIC = tmpgmm.BIC;
                        end
                    end

                    Btilde = gmm.Sigma;
                    obj.Weights = gmm.ComponentProportion;

                    Xtilde = gmm.mu.';

                    obj.UseAdaptiveEnsembleSize = true;
                    ensNa = ensN;

                    ensN = ncomp;

                case 'KDE'
                    Xbm = mean(Xb, 2);
                    Ab = Xb - repmat(Xbm, 1, ensN);
                    Ab = sqrt(beta)*Ab/sqrt(ensN - 1);
                    B = Ab*Ab.';
                    Btilde = repmat(B, 1, 1, ensN);
                    Xtilde = Xb;
            end

            as = zeros(ensN, 1);
            w = obj.Weights.';

            %% Perform AKDE

            if obj.UseAKDE
                % calculate geometrix means of the samples at the currently
                % estimated distribution

                %naive for now
                logfxi = zeros(1, ensN);
                for i = 1:ensN
                    as = zeros(1, ensN);
                    for j = 1:ensN
                        Bj = Btilde(:, :, j);
                        logdetBtilde = sum(log(diag(chol(Bj))));
                        z = (Xb(:, i) - Xb(:, j));
                        as(j) = log(w(j)) - 0.5*z.'*(Bj\z) ...
                            - 0.5*logdetBtilde - (n/2)*log(2*pi);
                    end
                    m = max(as);
                    logfxi(i) = m + log(sum(exp(as-m)));
                end
                logg = mean(logfxi);

                gamma = -1/n;
                lambdas = exp(logfxi - logg).^gamma;

                lambdaScale = reshape(lambdas.^2, 1, 1, ensN);

                Btilde = lambdaScale.*Btilde;

            end

            if ~isempty(obj.Localization)
                H = eye(n);
                rho = obj.Localization(mean(Xtilde, 2), H);
                for i = 1:ensN
                    Btilde(:, :, i) = rho.*Btilde(:, :, i);
                end
            end

            %% Perform the Gaussian Mixture update
            gmmupdatetype = obj.GMMUpdateType;
            gmmupdateparam = [];
            if iscell(gmmupdatetype)
                gmmupdateparam = gmmupdatetype{2};
                gmmupdatetype = gmmupdatetype{1};
            end

            Xtilde = repmat(Xtilde, 1, ensM);
            Btilde = repmat(Btilde, 1, 1, ensM);
            w = repmat(w, 1, ensM);

            switch gmmupdatetype
                case 'none'

                case "EKF"

                    m = size(Rs, 1);

                    H = obs.linearization(Xtilde);
                    S = pagemtimes(pagemtimes(H, Btilde), 'none', H, 'transpose') ...
                        + reshape(repmat(Rs, 1, ensN), m, m, []);
                    S = (S + pagectranspose(S))/2;

                    ysfull = reshape(repmat(ys, ensN, 1), m, []);
                    
                    HXtilde = obs.observeWithoutError(Xtilde);
                    d = reshape((HXtilde - ysfull), m, 1, []);

                    Sid = pagemldivide(S, d);

                    inov = pagemtimes(H, 'transpose', Sid, 'none');
                    inov = reshape(pagemtimes(Btilde, inov), n, []);

                    Xtilde = Xtilde - inov;

                    Ba = pagemtimes(H, Btilde);
                    Ba = pagemldivide(S, Ba);
                    Ba = pagemtimes(H, 'transpose', Ba, 'none');
                    Btilde = Btilde - pagemtimes(Btilde, Ba);

                    Btilde = (Btilde + pagectranspose(Btilde))/2;

                    wobsfull = reshape(repmat(wobs.', ensN, 1), 1, []);

                    L = pageeig(S);

                    aTa = reshape(pagemtimes(d, 'transpose', Sid, 'none'), 1, []);
                    as = -0.5*aTa + log(w) + log(wobsfull) ...
                        - reshape(0.5*sum(log(L), 1), 1, []);
                    as = as.';

                    m = max(as);
                    w = exp(as-(m + log(sum(exp(as-m))))).';
                    w = w/sum(w);
                case "UKF"
                    alpha = 1;
                    kappa = 3 - n;
                    beta = 2;
                    lambda = (alpha^2)*(n + kappa) - n;

                    as = zeros(ensN*ensM, 1);
                    for j = 1:ensM
                        y = ys(:, j);
                        R = Rs(:, :, j);
                        for i = 1:ensN
                            ind = sub2ind([ensN, ensM], i, j);

                            % set the local mean and covariance
                            Xmu = Xtilde(:, ind);
                            Btildei =  Btilde(:, :, ind);
                            % sample sigma points
                            sqBt = sqrt(n + lambda)*sqrtm(Btildei);
                            Xi = [Xmu, Xmu + sqBt, Xmu - sqBt];
                            HXi = obs.observeWithoutError(Xi);

                            Wm = [lambda/(lambda + n), (1/(2*(n + lambda)))*ones(1, 2*n)];
                            Wc = [lambda/(lambda + n) + (1 - alpha^2 + beta), ...
                                (1/(2*(n + lambda)))*ones(1, 2*n)];

                            HXimu = HXi*(Wm.');
                            A = (Xi - Xmu);
                            Z = (HXi - HXimu);

                            S = ((Wc.*Z)*Z.' + R);
                            S = (S + S.')/2;

                            RS = chol(S);
                            d = (HXimu - y);
                            BHt = (Wc.*A)*(Z.');
                            Xtilde(:, ind) = Xtilde(:, ind) - BHt*(S\d);

                            % recompute Btilde
                            Btildei = (Btildei - BHt*(S\(BHt.')));
                            Btilde(:, :, ind) = (Btildei + Btildei.')/2;

                            a = (RS.'\d);
                            as(ind) = -0.5*(a.'*a) - sum(log(diag(RS))) ...
                                + log(w(ind)) + log(wobs(j));
                        end
                    end
                    m = max(as);
                    w = exp(as-(m + log(sum(exp(as-m))))).';
                    w = w/sum(w);
                case "IEKF"

                    if isempty(gmmupdateparam)
                        % set the tolerance
                        gmmupdateparam = 1e-6;
                    end

                    h = @(x) obs.observeWithoutError(x);
                    

                    as = zeros(ensN*ensM, 1);
                    for j = 1:ensM
                        y = ys(:, j);
                        R = Rs(:, :, j);

                        V = @(x, xb, P) 0.5*(h(x) - y).'*(R\(h(x) - y)) ...
                            + 0.5*(x - xb).'*(P\(x - xb));

                        for i = 1:ensN
                            ind = sub2ind([ensN, ensM], i, j);
                            Xbi = Xb(:, ind);
                            Xtildei = Xtilde(:, ind);
                            Btildei =  Btilde(:, :, ind);

                            Vcur = V(Xtildei, Xbi, Btildei);

                            gnorm = inf;
                            linesearchdone = false;
                            while gnorm > gmmupdateparam
                                H = obs.linearization(Xtildei);

                                S = (H*Btildei*H.' + R);
                                S = (S + S.')/2;

                                deltax = Xbi - Xtildei + Btildei*(H.'*(S\(y - h(Xtildei) ...
                                    - H*(Xbi - Xtildei))));
                                gnorm = norm(deltax);

                                % perform naive line search
                                alpha = 1;

                                while true
                                    Xtildeinew = Xtildei + alpha*deltax;
                                    Vnew = V(Xtildeinew, Xbi, Btildei);

                                    if alpha < eps
                                        linesearchdone = true;
                                        break
                                    end

                                    if Vnew < Vcur
                                        Xtildei = Xtildeinew;
                                        Vcur = Vnew;
                                        break;
                                    else
                                        alpha = alpha/(1.5);
                                    end
                                end

                                if linesearchdone
                                    break;
                                end
                            end

                            % update x
                            Xtilde(:, ind) = Xtildei;

                            % update the covariance
                            H = obs.linearization(Xtildei);

                            S = (H*Btildei*H.' + R);
                            S = (S + S.')/2;

                            Btildei = (Btildei - Btildei*H.'*(S\(H*Btildei)));
                            Btilde(:, :, ind) = (Btildei + Btildei.')/2;

                            % compute the log weights

                            RS = chol(S);
                            d = h(Xb(:, ind)) - y;

                            a = (RS.'\d);
                            as(ind) = -0.5*(a.'*a) - sum(log(diag(RS))) ...
                                + log(w(ind))  + log(wobs(j));
                        end
                    end

                    % update the weights
                    m = max(as);
                    w = exp(as-(m + log(sum(exp(as-m))))).';
                    w = w/sum(w);

                case "BRUF"
                    K = gmmupdateparam;
                    for k = 1:K
                        wold = w;
                        HXtilde = obs.observeWithoutError(Xtilde);
                        for i = 1:ensN
                            H = obs.linearization(Xtilde(:, i));

                            Btildei =  Btilde(:, :, i);

                            S = (H*Btildei*H.' + K*R);
                            S = (S + S.')/2;
                            RS = chol(S);
                            d = (HXtilde(:, i) - y);
                            Xtilde(:, i) = Xtilde(:, i) - Btildei*(H.'*(S\d));

                            % recompute Btilde
                            Btildei = (Btildei - Btildei*(H.'*(S\(H*Btildei))));
                            Btilde(:, :, i) = (Btildei + Btildei.')/2;

                            a = (RS.'\d);
                            as(i) = -0.5*(a.'*a) - sum(log(diag(RS))) + log(wold(i));
                        end
                        m = max(as);
                        w = exp(as-(m + log(sum(exp(as-m))))).';
                        w = w/sum(w);
                    end

            end

            % resample
            Xa = Xb;

            if obj.UseAdaptiveEnsembleSize

            else
                ensNa = ensN;
            end

            switch obj.SamplingType
                case "QMC"
                    for  k = 1:ensNa

                        samples = obj.QMCStream.rand(1, n + 1).';
                        ind = find(samples(1) <= cumsum(w), 1, 'first');

                        sqBa = sqrtm(Btilde(:, :, ind));

                        if obj.UseRobustSampling
                            Xa(:, k) = Xtilde(:, ind) + sqBa*robustnorminv(samples(2:end));
                        else
                            Xa(:, k) = Xtilde(:, ind) + sqBa*norminv(samples(2:end));
                        end

                    end

                case "Gaussian"
                    for  k = 1:ensNa

                        ind = find(rand <= cumsum(w), 1, 'first');

                        sqBa = sqrtm(Btilde(:, :, ind));

                        if obj.UseRobustSampling
                            Xa(:, k) = Xtilde(:, ind) + sqBa*robustnorminv(rand(n, 1));
                        else
                            Xa(:, k) = Xtilde(:, ind) + sqBa*randn(n, 1);
                        end
                    end
                case "Rejection"
                    error('Not supported yet');

                case "None"
                    obj.Ensemble = Xtilde;

                    obj.Weights = w.';
                    return
            end

            if any(~isreal(Xa), "all")
                keyboard
            end

            obj.Ensemble = Xa;

            obj.Weights = ones(ensNa, 1) / ensNa;

        end

    end

end


function [c, g] = cost(L, Pbi, Si)

W = Si + L*L.';

d = Pbi - Si + Si*(W\Si);

c = sum(d.^2, 'all');

if nargout > 1
    g = -4*(W\(Si*(Pbi - Si + Si*(W\Si))*Si*(W\L)));
end

end

function x = robustnorminv(x)

x(x < 0) = 0;
x(x > 1) = 1;

I1 = x <= 1/6;
I2 = and(x > 1/6, x <= 5/6);
I3 = x > 5/6;

x(I1) = (-3)+2.*6.^(1/3).*x(I1).^(1/3);
x(I2) = 2.*3.^(1/2).*sin((1/3).*atan(((-2)+4.*x(I2)).*((-1)+(-16).*((-1)+x(I2)).* ...
    x(I2)).^(-1/2)));
x(I3) = 3+(-2).*(6+(-6).*x(I3)).^(1/3);

x = real(x);

end
