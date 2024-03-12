classdef Gaussian < datools.uncertainty.Uncertainty

    properties
        Mean
    end

    properties (GetAccess=private)
        InternalCovariance
        InternalCovarianceSqrt
    end

    properties (Dependent)
        Covariance
    end

    methods

        function obj = Gaussian(varargin)
            p = inputParser;
            addOptional(p, 'Covariance', 1);
            addOptional(p, 'Mean', 0);
            parse(p, varargin{:});

            s = p.Results;

            obj.InternalCovariance = s.Covariance;
            obj.InternalCovarianceSqrt = sqrtm(s.Covariance);
            obj.Mean = s.Mean;
        end

        function x = sample(obj, N)
            A = obj.InternalCovarianceSqrt*randn(size(obj.InternalCovarianceSqrt, 1), N);
            x = A + obj.Mean;
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

        function set.Covariance(obj, C)
            obj.InternalCovariance = C;
            obj.InternalCovarianceSqrt = sqrtm(C);
        end

        function C = get.Covariance(obj)
            C = obj.InternalCovariance;
        end

        function obj = asGaussian(obj)

        end

        function gmm = asGMM(obj)
            gmm = datools.uncertainty.GMM(obj.Mean, obj.Covariance, 1);
        end

    end

end
