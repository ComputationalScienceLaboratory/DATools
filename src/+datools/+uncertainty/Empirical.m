classdef Empirical < datools.uncertainty.Uncertainty

    properties
        States
    end

    properties (Dependent)
        Covariance
        Mean
    end

    methods

        function obj = Empirical(varargin)
            p = inputParser;
            addOptional(p, 'States', 0);
            parse(p, varargin{:});

            s = p.Results;

            obj.States = s.States;
        end

        function x = sampleWithReplacement(obj, N)
            Xs = obj.States;
            NS = size(Xs, 2);
            x = Xs(:, randi(NS, 1, N));
        end

        function x = sampleWithoutReplacement(obj, N)
            Xs = obj.States;
            NS = size(Xs, 2);
            if N > NS
                error('Cannot sample without replacement more states than exist.')
            end
            x = Xs(:, randperm(NS, N));
        end

        function x = sample(obj, N)
            x = obj.sampleWithoutReplacement(N);
        end

        function xp = addError(obj, x)
            xp = x + obj.sample(size(x, 2));
        end

        function lp = log(obj, x)
            mu = obj.Mean;
            sigma = obj.Covariance;

            n = size(x, 1);

            sigmaj = sigma;
            R = chol(sigmaj);
            as = sum(-0.5*(x - mu).*(sigmaj\(x - mu)), 1) ...
                - sum(log(diag(R))) - (n/2)*log(2*pi);
            lp = as;

        end

        function px = pdf(obj, x)
            px = exp(obj.log(x));
        end

        function set.Covariance(obj, Cnew)
            C = obj.Covariance;
            X = obj.States;
            A = X - obj.Mean;
            A = sqrtm(Cnew)*(sqrtm(C)\A);
            obj.States = obj.Mean + A;
        end

        function C = get.Covariance(obj)
            X = obj.States;
            N = size(X, 2);
            A = X - mean(X, 2);
            C = (A*A.')/(N - 1);
        end

        function set.Mean(obj, munew)
            X = obj.States;
            A = X - obj.Mean;
            obj.States = munew + A;
        end

        function mu = get.Mean(obj)
            X = obj.States;
            mu = mean(X, 2);
        end

        function g = asGaussian(obj)
            g = datools.uncertainty.Gaussian(obj.Mean, obj.Covariance);
        end

        function gmm = asGMM(obj, ~, sbeta)
            if nargin < 3 || isempty(sbeta)
                sbeta = 1;
            end
            % use Silverman's rule of thumb
            X = obj.States;
            [n, N] = size(X);
            beta2 = sbeta*(((4/(n + 2))^(2/(n + 4)))*((N)^(-2/(n + 4))));
            P = beta2*(obj.Covariance);
            Sigmas = repmat(P, 1, 1, N);
            ws = ones(1, N)/N;

            gmm = datools.uncertainty.GMM(X, Sigmas, ws);
        end

    end

end
