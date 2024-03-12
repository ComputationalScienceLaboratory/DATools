classdef IndexedObservation < datools.observation.operator.Observation
%INDEXED Defines the indexed observation operator

    properties
        Indices    % Indices of observed state
    end

    methods

        function obj = IndexedObservation(nvars, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'Indices', 1);
            parse(p, varargin{:});

            s = p.Results;

            obj@datools.observation.operator.Observation(nvars, p.Unmatched);

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