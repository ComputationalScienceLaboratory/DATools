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
            addParameter(p, 'F', @(~, x) x);
            addParameter(p, 'Jacobian', @(~, x) speye(numel(x)));
            addParameter(p, 'Indices', 1);
            parse(p, varargin{:});

            s = p.Results;

            obj@datools.observation.Observation(nvars, p.Unmatched);

            obj.F = s.F;
            obj.Jacobian = s.Jacobian;
            obj.Indices = s.Indices;

        end

        function y = observeWithoutError(obj, t, x)
            yfull = obj.F(t, x);
            y = yfull(obj.Indices, :);
        end

        function H = linearization(obj, t, x)
            H = obj.Jacobian(t, x);
            if size(H, 1) == size(H, 2)
                H = H(obj.Indices, :);
            end
        end

    end

end