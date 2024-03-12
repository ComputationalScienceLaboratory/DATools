classdef Uniform < datools.uncertainty.Uncertainty

    properties
        Lower
        Upper
    end

    properties (Dependent)
        Mean
        Covariance
    end

    methods

        function obj = Uniform(varargin)
            p = inputParser;
            addOptional(p, 'Lower', 0);
            addOptional(p, 'Upper', 1);
            parse(p, varargin{:});

            s = p.Results;

            obj.Lower = reshape(s.Lower, [], 1);
            obj.Upper = reshape(s.Upper, [], 1);
        end

        function x = sample(obj, N)
            n = size(obj.Lower, 1);
            d = obj.Upper - obj.Lower;
            x = (d.*rand(n, N)) + obj.Lower;
        end

        function xp = addError(obj, x)
            xp = x + obj.sample(size(x, 2));
        end

        function lp = log(obj, x)
            N = size(x, 2);
            bin = all(x <= obj.Upper, 1) & all(x >= obj.Lower, 1);
            
            volume = prod(obj.Upper - obj.Lower);
            
            lp = -inf*ones(1, N);
            lp(bin) = -log(volume);
        end

        function px = pdf(obj, x)
            px = exp(obj.log(x));
            px(isnan(px)) = 0;
        end

        function set.Covariance(obj, C)
            mu = obj.Mean;
            % we can only set the variance
            c = diag(C);
            dnew = sqrt(12*c);
            obj.Lower = mu - 0.5*dnew;
            obj.Upper = mu + 0.5*dnew;
        end

        function C = get.Covariance(obj)
            C = diag((1/12)*((obj.Upper - obj.Lower).^2));
        end

        function set.Mean(obj, munew)
            mu = 0.5*(obj.Lower + obj.Upper);
            obj.Lower = obj.Lower - mu + munew;
            obj.Upper = obj.Upper - mu + munew;
        end

        function mu = get.Mean(obj)
            mu = 0.5*(obj.Lower + obj.Upper);
        end

        function g = asGaussian(obj)
            g = datools.uncertainty.Gaussian('Mean', obj.Mean, ...
                'Covariance', obj.Covariance);
        end

        function gmm = asGMM(obj, N)
            if nargin < 2
                N = 25;
            end

            n = numel(obj.Lower);

            Ngrid = ceil(N^(1/n));

            ls = linspace(0, 1, Ngrid);

            xs = (1 - ls).*(obj.Lower) + ls.*(obj.Upper);

            grid = creategrid(xs, n);

            em = datools.uncertainty.Empirical('State', grid);
            
            gmm = em.asGMM();
        end

    end

end

function grid = creategrid(xs, n)
if n == 1
    grid = xs(1, :);
else
    grid = creategrid(xs, n - 1);
    Nnew = size(grid, 2);
    grid = repmat(grid, 1, size(xs, 2));

    xnew = reshape(repmat(xs(n, :), Nnew, 1), 1, []);
    grid = [grid; xnew];
end
end
