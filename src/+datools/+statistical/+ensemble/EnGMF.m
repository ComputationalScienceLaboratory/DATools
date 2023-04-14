classdef EnGMF < datools.statistical.ensemble.EnF

    properties
        BandwidthScale
        BRUFSteps
        UseRobustSampling
        SamplingType
        UseELocalization
        UseAdaptiveEnsembleSize = false
        RadiusScale
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
            addParameter(p, 'UseRobustSampling', false);
            addParameter(p, 'SamplingType', "Gaussian");
            addParameter(p, 'UseELocalization', false);
            addParameter(p, 'RadiusScale', 1);

            parse(p, varargin{2:end});

            s = p.Results;

            kept = p.Unmatched;

            obj@datools.statistical.ensemble.EnF(varargin{1}, kept);

            obj.BandwidthScale = s.BandwidthScale;
            obj.BRUFSteps = s.BRUFSteps;
            obj.UseRobustSampling = s.UseRobustSampling;
            obj.SamplingType = s.SamplingType;
            obj.UseELocalization = s.UseELocalization;
            obj.RadiusScale = s.RadiusScale;

        end
    end

    methods

        function analysis(obj, R, y)

            tc = obj.Model.TimeSpan(1);

            Xb = obj.Ensemble;
            ensN = obj.NumEnsemble;
            n = size(Xb, 1);

            if obj.SamplingType == "QMC" && isempty(obj.QMCStream)
                obj.QMCStream = qrandstream('sobol', n + 1, 'Skip', 100);
            end

            beta = obj.BandwidthScale*((4/(n + 2))^(2/(n + 4)))*((ensN)^(-2/(n + 4)));

            if obj.UseELocalization
                Xbm = mean(Xb, 2);
                Ab = Xb - repmat(Xbm, 1, ensN);
                Ab = Ab/sqrt(ensN - 1);

                radiusscale = obj.RadiusScale;
                epsilon = 1e-4;
                % calculate the distances
                dist = zeros(ensN, ensN);
                wD = zeros(ensN, ensN);
                for i = 1:ensN
                    for j = 1:ensN
                        dist(i, j) = norm(Xb(:, i) - Xb(:, j));
                    end
                    [~, I] = sort(dist(i, :));
                    r = radiusscale*dist(i, I(round(sqrt(ensN))));

                    wD(i, :) = exp(-0.5*(dist(i, :)/r).^2) + epsilon;
                end
                wD = wD./sum(wD, 2);

                Pb = Ab*Ab.';
                trPb = trace(Pb);

                Btilde = zeros(n, n, ensN);
                trPbis = zeros(ensN, 1);

                for i = 1:ensN
                    wDi = wD(i, :).';
                    Xbm = Xb*wDi;
                    c = 1/(1 - sum(wDi.^2));
                    Pbi = c*((Xb - Xbm)*(diag(wDi))*(Xb - Xbm).');
                    trPbis(i) = trace(Pbi);

                    Btilde(:, :, i) = beta*Pbi;
                end

                Btilde = (trPb/mean(trPbis))*Btilde;

            else
                Xbm = mean(Xb, 2);
                Ab = Xb - repmat(Xbm, 1, ensN);
                Ab = sqrt(beta)*Ab/sqrt(ensN - 1);
                B = Ab*Ab.';
                Btilde = repmat(B, 1, 1, ensN);
            end

            Xtilde = Xb;
            as = zeros(ensN, 1);
            w = obj.Weights.';

            K = obj.BRUFSteps;
            for k = 1:K
                wold = w;
                HXtilde = obj.Observation.observeWithoutError(tc, Xtilde);
                for i = 1:ensN
                    H = obj.Observation.linearization(tc, Xtilde(:, i));

                    Btildei =  Btilde(:, :, i);

                    S = (H*Btildei*H.' + K*R);
                    S = (S + S.')/2;
                    RS = chol(S);
                    d = (HXtilde(:, i) - y);
                    Xtilde(:, i) = Xtilde(:, i) - Btildei*H.'*(S\d);

                    % recompute Ptilde
                    Btildei = (Btildei - Btildei*H.'*(S\(H*Btildei)));
                    Btilde(:, :, i) = (Btildei + Btildei.')/2;

                    a = (RS.'\d);
                    as(i) = -0.5*(a.'*a) - sum(log(diag(RS)));
                end
                m = max(as);
                w = exp(as-(m + log(sum(exp(as-m))))).';
                w = w/sum(w);
                w = w.*wold;
                w = w/sum(w);
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
                case "None"
                    obj.Ensemble = Xtilde;

                    obj.Weights = w.';
                    return
            end

            obj.Ensemble = Xa;

            obj.Weights = ones(ensNa, 1) / ensNa;

        end

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
