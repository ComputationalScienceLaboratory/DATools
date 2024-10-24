classdef Indexed < datools.observation.Observation
    % Derived from the "observation" parent class. This can
    % be used (for ex.) to observe the state space via indexing.
    properties
        Indices
    end

    methods

        function obj = Indexed(nvars, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'Indices', 1);
            parse(p, varargin{:});

            s = p.Results;

            obj@datools.observation.Observation(nvars, p.Unmatched);

            obj.Indices = s.Indices;

        end

        function y = observeWithoutError(obj, x)
            y = x(obj.Indices, :);
        end

        function H = linearization(obj, x)
            %I = speye(obj.NumVars);
            %H = I(obj.Indices, :);
            N = size(x, 2);
            I = eye(obj.NumVars);
            H = repmat(I(obj.Indices, :), 1, 1, N);
        end

    end

end