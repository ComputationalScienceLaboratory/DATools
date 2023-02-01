classdef Nonlinear < datools.observation.operator.Observation
%NONLINEAR Defines the non linear observation operator

    properties
        F
        Jacobian    % Jacobian of non-linear Observation operator
        Indicies    % Indices of observed state
    end

    methods

        function obj = Nonlinear(nvars, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'F', @(~, x) x);
            addParameter(p, 'Jacobian', @(~, x) speye(numel(x)));
            addParameter(p, 'Indicies', 1);
            parse(p, varargin{:});

            s = p.Results;

            obj@datools.observation.operator.Observation(nvars, p.Unmatched);

            obj.F = s.F;
            obj.Jacobian = s.Jacobian;
            obj.Indicies = s.Indicies;

        end

        function y = observeWithoutError(obj, t, x)
            yfull = obj.F(t, x);
            y = yfull(obj.Indicies, :);
        end

        function H = linearization(obj, t, x)
            H = obj.Jacobian(t, x);
            if size(H, 1) == size(H, 2)
                H = H(obj.Indices, :);
            end
        end

    end

end