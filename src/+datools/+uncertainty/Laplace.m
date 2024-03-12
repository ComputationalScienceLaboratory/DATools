classdef Laplace < datools.uncertainty.Uncertainty

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

        function obj = Laplace(varargin)
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
            Z = exprnd(1, 1, N);
            A = sqrt(Z) .* A;
            x = A + obj.Mean;
        end

        function xp = addError(obj, x)
            xp = x + obj.sample(size(x, 2));
        end

        function lp = log(obj, x)
            mu = obj.Mean;
            sigma = obj.Covariance;

            n = size(x, 1);
            v = (2 - n)/2;

            sigmaj = sigma;
            R = chol(sigmaj);
            b = sum((x - mu).*(sigmaj\(x - mu)), 1);
            c = (v/2)*log(b/2);
            z = sqrt(2*b);
            d = log(besselk(v, z, 1)) - z;
            as = log(2) + c + d ...
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

        function g = asGaussian(obj)
            g = datools.uncertainty.Gaussian("Mean", obj.Mean, "Covariance", obj.Covariance);
        end

        function gmm = asGMM(obj, N)
            if nargin < 2 || isempty(N)
                N = 4;
            end

            pkm1 = datools.utils.Polynomial(1);
            pk = datools.utils.Polynomial([-1 1]);

            for k = 1:N
                tmp = datools.utils.Polynomial([-1, 2*k+1]);
                pkp1 = (1/(k + 1))*(tmp*pk - k*pkm1);

                pkm1 = pk;
                pk = pkp1;
            end

            zs = roots(pkm1.p);

            w = zs./(((N + 1).^2)*(polyval(pk.p, zs).^2));
            w = w./sum(w);

            zs = reshape(zs, 1, 1, []);
            
            Sigma = zs.*(obj.Covariance);
            means = repmat(obj.Mean, 1, N);

            gmm = datools.uncertainty.GMM(means, Sigma, w);

        end

    end

end
