classdef ETPF < datools.statistical.ensemble.EnF

    properties
        SinkhornKnoppLambda = 0;
        SinkhornKnoppIterations = 10;
        SecondOrderCorrection = false;
        LocalizationEnsembleDistance
        Name = "Ensemble Transform Particle Filter"
    end

    methods

        function obj = ETPF(varargin)

            p = inputParser;
            p.KeepUnmatched = true;
            p.PartialMatching = false;

            addParameter(p, 'SinkhornKnoppLambda', 0);
            addParameter(p, 'SinkhornKnoppIterations', 10);
            addParameter(p, 'SecondOrderCorrection', false);

            % this is recommended according to the Reich book

            %addParameter(p, 'LocalizationEnsembleDistance', ...
            %    @(~, inds, k) ...
            %    spdiags([zeros(k-1, 1); 1; zeros(numel(inds)-k+1, 1)], 0, numel(inds), numel(inds)));

            addParameter(p, 'LocalizationEnsembleDistance', ...
                @(~, inds, k) ...
                [zeros(k-1, 1); 1; zeros(numel(inds)-k, 1)]);


            parse(p, varargin{2:end});

            s = p.Results;

            kept = p.Unmatched;

            obj@datools.statistical.ensemble.EnF(varargin{1}, kept);

            obj.LocalizationEnsembleDistance = s.LocalizationEnsembleDistance;
            obj.SinkhornKnoppLambda = s.SinkhornKnoppLambda;
            obj.SinkhornKnoppIterations = s.SinkhornKnoppIterations;
            obj.SecondOrderCorrection = s.SecondOrderCorrection;

        end

        function analysis(obj, obs)

            %y = obs.Mean;
            %R = obs.Covariance;

            tau = obj.Rejuvenation;

            %dR = decomposition(R, 'chol');

            X = obj.Ensemble;
            ensN = obj.NumEnsemble;

            HX = obj.Observation.observeWithoutError(X);

            xdist = zeros(ensN, ensN);

            if isempty(obj.Localization)

                for i = 1:ensN
                    xtemp = X - repmat(X(:, i), 1, ensN);
                    xdist(i, :) = vecnorm(xtemp).^2;
                end

                % more efficient way of calculating weights
                as = obs.log(HX).';
                m = max(as);
                w = exp(as-(m + log(sum(exp(as-m)))));

                lambda = obj.SinkhornKnoppLambda;

                if lambda == 0
                    AeqT = [kron(speye(ensN), ones(1, ensN)); kron(ones(ensN, 1), speye(ensN)).'];
                    lbT = zeros((ensN)*ensN, 1);
                    optsT = optimoptions('linprog', 'Display', 'off');
                    beqT = [ones(ensN, 1) / ensN; w];
                    f = xdist(:);
                    Tx = linprog(f, [], [], AeqT, beqT, lbT, [], optsT);
                    Tx = ensN * reshape(Tx, ensN, ensN);
                else
                    its = obj.SinkhornKnoppIterations;
                    Tx = datools.utils.sinkhornknopp(xdist, ones(ensN, 1), ensN*w, lambda, its);
                end

                D = zeros(ensN, ensN);
                if obj.SecondOrderCorrection
                    % Delta
                    W = diag(w);
                    Bric = Tx - w * ones(1, ensN);
                    Aric = ensN * (W - w * w.') - Bric * Bric.';

                    rs = @(D) reshape(D, ensN, ensN);
                    urs = @(D) reshape(D, ensN*ensN, 1);
                    f = @(t, D) urs(-Bric*rs(D)-rs(D)*(Bric.')+Aric-rs(D)*rs(D));

                    D = rs(datools.utils.dpETPF2(f, D(:)));
                end


                xa = X * (Tx + D);

            else
                xa = X;
                R = obs.Covariance;
                y = obs.Mean;
                invR = spdiags(1./diag(R), 0, size(R, 1), size(R, 2));
                Hi = obj.Observation.Indices;

                t0 = HX - y;

                n = size(X, 1);
                            
                for k = 1:n

                    C = obj.Localization(mean(X, 2), Hi, k);

                    if isempty(obj.LocalizationEnsembleDistance)
                        rho = ones(n, 1);
                    else
                        rho = obj.LocalizationEnsembleDistance(mean(X, 2), 1:n, k);
                        %rho = full(diag(rho));
                    end

                    for i = 1:ensN
                        xtemp = X - repmat(X(:, i), 1, ensN);
                        xdist(i, :) = sum(rho.*(xtemp.^2), 1);
                    end

                    % more efficient way of calculating weights
                    as = (-0.5 * sum(t0.*((C * invR) * t0), 1)).';
                    m = max(as);
                    w = exp(as-(m + log(sum(exp(as-m)))));

                    lambda = obj.SinkhornKnoppLambda;

                    if lambda == 0
                        AeqT = [kron(speye(ensN), ones(1, ensN)); kron(ones(ensN, 1), speye(ensN)).'];
                        lbT = zeros((ensN)*ensN, 1);
                        optsT = optimoptions('linprog', 'Display', 'off');
                        beqT = [ones(ensN, 1) / ensN; w];
                        f = xdist(:);
                        Tx = linprog(f, [], [], AeqT, beqT, lbT, [], optsT);
                        Tx = ensN * reshape(Tx, ensN, ensN);
                    else
                        its = obj.SinkhornKnoppIterations;
                        Tx = datools.utils.sinkhornknopp(xdist, ones(ensN, 1), ensN*w, lambda, its);
                    end

                    D = zeros(ensN, ensN);
                    if obj.SecondOrderCorrection
                        % Delta
                        W = diag(w);
                        Bric = Tx - w * ones(1, ensN);
                        Aric = ensN * (W - w * w.') - Bric * Bric.';

                        rs = @(D) reshape(D, ensN, ensN);
                        urs = @(D) reshape(D, ensN*ensN, 1);
                        f = @(t, D) urs(-Bric*rs(D)-rs(D)*(Bric.')+Aric-rs(D)*rs(D));

                        D = rs(datools.utils.dpETPF2(f, D(:)));
                    end

                    xa(k, :) = X(k, :) * (Tx + D);

                end

            end

            obj.Ensemble = xa;
            obj.Weights = ones(ensN, 1) / ensN;
            obj.rejuvenate(tau, X);

            obj.Model.update(0, obj.MeanEstimate);


        end

    end

end
