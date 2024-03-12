classdef Nonlinear < datools.observation.Observation

    properties
        F
        Jacobian
        Indices
    end

    methods

        function obj = Nonlinear(nvars, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'F', @(x) x);
            addParameter(p, 'Jacobian', @(x) speye(numel(x)));
            addParameter(p, 'Indices', 1);
            parse(p, varargin{:});

            s = p.Results;

            obj@datools.observation.Observation(nvars, p.Unmatched);

            obj.F = s.F;
            obj.Jacobian = s.Jacobian;
            obj.Indices = s.Indices;

        end

        function y = observeWithoutError(obj, x)
            yfull = obj.F(x);
            y = yfull(obj.Indices, :);
        end

        function H = linearization(obj, x, v)
            if nargin < 4
                H = obj.Jacobian(x);
                if size(H, 1) == size(H, 2)
                    if size(H, 3) == 1
                        H = H(obj.Indices, :);
                    else
                        H = H(obj.Indices, :, :);
                    end
                end
            else
                H = obj.Jacobian(x, v);
                H = H(obj.Indices, :);
            end
        end

    end

end