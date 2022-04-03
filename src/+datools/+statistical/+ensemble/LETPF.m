classdef LETPF < datools.statistical.ensemble.EnF

    properties
        LocalizationEnsembleDistance
    end


    methods
        function obj = LETPF(varargin)

            p = inputParser;
            p.KeepUnmatched = true;

            % this is recommended according to the Reich book

            addParameter(p, 'LocalizationEnsembleDistance', ...
                @(~, ~, inds, k) ...
                spdiags([zeros(k-1, 1); 1; zeros(numel(inds)-k+1, 1)], 1, numel(inds), numel(inds)));

            parse(p, varargin{2:end});

            s = p.Results;

            kept = p.Unmatched;

            obj@datools.statistical.ensemble.EnF(varargin{1}, kept);

            obj.LocalizationEnsembleDistance = s.LocalizationEnsembleDistance;

        end
    end

    methods

        function analysis(obj, R, y)

            tau = obj.Rejuvenation;

            tc = obj.Model.TimeSpan(1);

            xf = obj.Ensemble;
            n = size(xf, 1);
            ensN = obj.NumEnsemble;

            AeqT = [kron(speye(ensN), ones(1, ensN)); kron(ones(ensN, 1), speye(ensN)).'];
            lbT = zeros((ensN)*ensN, 1);
            optsT = optimoptions('linprog', 'Display', 'off');

            Hxf = obj.Observation.observeWithoutError(tc, xf);

            xdist = zeros(ensN, ensN);


            t0 = Hxf - y;

            xa = xf;

            invR = spdiags(1./diag(R), 0, size(R, 1), size(R, 2));
            Hi = obj.Observation.Indices;

            for k = 1:n

                if isempty(obj.Localization)
                    C = speye(size(invR));
                else
                    C = obj.Localization(tc, mean(xf, 2), Hi, k);
                end

                if isempty(obj.LocalizationEnsembleDistance)
                    rho = ones(n, 1);
                else
                    rho = obj.LocalizationEnsembleDistance(tc, mean(xf, 2), 1:n, k);
                    rho = full(diag(rho));
                end

                for i = 1:ensN
                    xtemp = xf - repmat(xf(:, i), 1, ensN);
                    xdist(i, :) = sum(rho.*(xtemp.^2), 1);
                end

                % more efficient way of calculating weights
                as = (-0.5 * sum(t0.*((C * invR) * t0), 1)).';
                m = max(as);
                w = exp(as-(m + log(sum(exp(as-m)))));

                beqT = [ones(ensN, 1) / ensN; w];
                f = xdist(:);
                Tx = linprog(f, [], [], AeqT, beqT, lbT, [], optsT);
                Tx = ensN * reshape(Tx, ensN, ensN);

                xa(k, :) = xf(k, :) * Tx;

            end

            obj.Ensemble = xa;
            obj.Weights = ones(ensN, 1) / ensN;
            obj.rejuvenate(tau, xf);

            obj.Model.update(0, obj.BestEstimate);

        end

    end

end
