classdef GMM < datools.uncertainty.Uncertainty

    properties
        Mu
        W
    end

    properties(GetAccess=private)
        InternalSigma
        InternalSigmaSqrt
    end

    properties(Dependent)
        Sigma
        SigmaSqrt
    end

    properties(Dependent)
        Mean
        Covariance
    end

    methods
        function obj = GMM(mu, sigma, w)
            obj.Mu = mu;
            obj.Sigma = sigma;
            obj.W = w;
        end

        function X = sample(obj, N)
            mu = obj.Mu;
            sigmasqrt = obj.SigmaSqrt;
            w = obj.W;
            n = size(mu, 1);
            w = w(:);
            
            cw = cumsum(w);
            inds = sum(rand(1, N) <= cw, 1);

            sqBas =  sigmasqrt(:, :, inds);

            X = mu(:, inds) + reshape(pagemtimes(sqBas, randn(n, 1, N)), n, N);
        end

        function m = get.Mean(obj)
            w = obj.W;
            mu = obj.Mu;
            w = reshape(w, 1, []);

            m = sum(w.*mu, 2);
        end

        function c = get.Covariance(obj)
            w = obj.W;
            mu = obj.Mu;
            sigma = obj.Sigma;
            w = reshape(w, 1, []);

            m = obj.Mean;
            A = sqrt(w).*(mu - m);
            wc = reshape(w, 1, 1, []);
            c = sum(wc.*sigma, 3) + A*A.';
        end

        function lp = log(obj, x)
            mu = obj.Mu;
            sigma = obj.Sigma;
            w = obj.W;

            n = size(x, 1);
            N = size(x, 2);
            M = size(mu, 2);

            as = zeros(N, M);
            for j = 1:M
                sigmaj = sigma(:, :, j);
                R = chol(sigmaj);
                as(:, j) = (-0.5*sum((x - mu(:, j)).*(sigmaj\(x - mu(:, j))), 1) ...
                    - sum(log(diag(R))) - (n/2)*log(2*pi) + log(w(j)))';
            end

            ma = max(as, [], 2);
            lp = ma + log(sum(exp(as - ma), 2));
            lp = lp.';

        end

        function px = pdf(obj, x)
            px = exp(obj.log(x));
        end

        function set.Sigma(obj, s)
            obj.InternalSigma = s;
            for i = 1:size(s, 3)
                obj.InternalSigmaSqrt(:, :, i) = sqrtm(s(:, :, i));
            end
        end

        function s = get.Sigma(obj)
            s = obj.InternalSigma;
        end

        function set.SigmaSqrt(obj, s)
            obj.InternalSigmaSqrt = s;
            for i = 1:size(s, 3)
                obj.InternalSigma(:, :, i) = s(:, :, i)*s(:, :, i).';
            end
        end

        function s = get.SigmaSqrt(obj)
            s = obj.InternalSigmaSqrt;
        end

    end

end
