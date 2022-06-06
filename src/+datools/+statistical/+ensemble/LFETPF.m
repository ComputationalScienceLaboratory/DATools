classdef LFETPF < datools.statistical.ensemble.EnF

    properties
        B
        Bsqrt
        SurrogateEnsN
        Laplace
        LocalizationEnsembleDistance
    end

    methods
        function obj = LFETPF(varargin)

            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'B', 1);
            addParameter(p, 'SurrogateEnsembleSize', 2);
            addParameter(p, 'Laplace', false);
            addParameter(p, 'LocalizationEnsembleDistance', ...
                @(~, ~, inds, k) ...
                spdiags([zeros(k-1, 1); 1; zeros(numel(inds)-k+1, 1)], 1, numel(inds), numel(inds)));

            parse(p, varargin{2:end});

            s = p.Results;

            kept = p.Unmatched;

            obj@datools.statistical.ensemble.EnF(varargin{1}, kept);

            obj.B = s.B;
            if iscell(obj.B)
                Bs = cell(numel(obj.B), 1);
                for i = 1:numel(obj.B)
                    Bs{i} = sqrtm(obj.B{i});
                end
                obj.Bsqrt = Bs;
            else
                obj.Bsqrt = sqrtm(s.B);
            end
            obj.SurrogateEnsN = s.SurrogateEnsembleSize;
            obj.Laplace = s.Laplace;

            obj.LocalizationEnsembleDistance = s.LocalizationEnsembleDistance;

        end
    end

    methods

        function analysis(obj, R, y)

            inflation = obj.Inflation;
            tau = obj.Rejuvenation;

            ensM = obj.SurrogateEnsN;

            tc = obj.Model.TimeSpan(1);

            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;

            AeqT = [kron(speye(ensN), ones(1, ensN+ensM)); kron(ones(ensN, 1), speye(ensN+ensM)).'];
            lbT = zeros((ensN + ensM)*ensN, 1);
            optsT = optimoptions('linprog', 'Display', 'off');


            xfm = mean(xf, 2);
            Af = xf - xfm;

            n = size(xf, 1);

            Covsqrt = (1 / sqrt(ensN-1) * Af);
            Bsqrtf_all = obj.Bsqrt;

            NN = ensN - 1;

            if iscell(Bsqrtf_all)

                Bi = 0;
                Uhopt = inf;

                for i = 1:numel(Bsqrtf_all)

                    Bsqrtfi = Bsqrtf_all{i};
                    s = svd(Bsqrtfi \ Covsqrt);
                    trC = sum(s.^2);
                    tr2C = trC * trC;
                    trC2 = sum(s.^4);

                    % calculate sphericity
                    Uh = (n * trC2 / tr2C - 1) / (n - 1);

                    if Uh < Uhopt
                        Uhopt = Uh;
                        Bi = i;
                    end

                end
                Bsqrtf = Bsqrtf_all{Bi};

            else

                Bsqrtf = Bsqrtf_all;
                s = svd(Bsqrtf \ Covsqrt);
                trC = sum(s.^2);
                tr2C = trC * trC;
                trC2 = sum(s.^4);

                % calculate sphericity
                Uhopt = (n * trC2 / tr2C - 1) / (n - 1);

            end

            gamma = min((NN - 2)/(NN * (NN + 2))+((n + 1) * NN - 2)/(Uhopt * NN * (NN + 2) * (n - 1)), 1);

            if ensM == 0
                gamma = 0;
            end
            mu = trC / n;

            %mu =  sum(s1.^2)/sum(s2.^2);
            Asynth = sqrt(mu) * Bsqrtf * randn(n, ensM);

            laplace = obj.Laplace;
            if laplace
                Z = exprnd(1, 1, ensM);
                Asynth = sqrt(Z) .* Asynth;
            end

            Asynth = Asynth - mean(Asynth, 2);

            Afrak = [Af, inflation * Asynth];
            chiF = xfm + Afrak;

            Hchif = obj.Observation.observeWithoutError(tc, chiF);

            invR = spdiags(1./diag(R), 0, size(R, 1), size(R, 2));
            Hi = obj.Observation.Indices;

            xa = xf;

            for k = 1:n

                if isempty(obj.Localization)
                    C = speye(size(invR));
                else
                    C = obj.Localization(tc, mean(xf, 2), Hi, k);
                end

                if isempty(obj.LocalizationEnsembleDistance)
                    rho = ones(n, 1);
                else
                    rho = obj.LocalizationEnsembleDistance(tc, mean(chiF, 2), 1:n, k);
                    rho = full(diag(rho));
                end

                xdist = zeros(ensN+ensM, ensN+ensM);

                for i = 1:(ensN + ensM)
                    xtemp = chiF - chiF(:, i);
                    xdist(i, :) = sum(rho.*(xtemp.^2), 1);
                end

                t0 = Hchif - y;

                % more efficient way of calculating weights
                as = (-0.5 * sum(t0.*((C * invR) * t0), 1)).';
                m = max(as);
                w = exp(as-(m + log(sum(exp(as-m)))));

                xdist = xdist(:, 1:ensN);

                w(1:ensN) = w(1:ensN) * (1 - gamma) / (ensN);
                w((ensN + 1):end) = w((ensN + 1):end) * gamma / (ensM);

                w = w / sum(w);

                beqT = [ones(ensN, 1) / ensN; w];
                f = xdist(:);
                Tx = linprog(f, [], [], AeqT, beqT, lbT, [], optsT);
                Tx = ensN * reshape(Tx, ensN+ensM, ensN);

                xa(k, :) = chiF(k, :) * Tx;

            end

            obj.Ensemble = xa;
            obj.Weights = ones(ensN, 1) / ensN;

            % rejuvenation of this sort is not performed by default,
            % however it is included in this algorithm purely for
            % completeness
            obj.rejuvenate(tau, xf);

            obj.Model.update(0, obj.BestEstimate);

        end

    end

end
