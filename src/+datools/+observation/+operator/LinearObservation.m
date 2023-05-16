classdef LinearObservation < datools.observation.operator.Observation
%LINEAR Defines the linear observation operator
%   H is already linearized

    properties
        H    % Observation operator
    end

    methods

        function obj = LinearObservation(nvars, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'H', speye(nvars));
            parse(p, varargin{:});

            s = p.Results;

            obj@datools.observation.operator.Observation(nvars, p.Unmatched);

            obj.H = s.H;

        end

        function y = observeWithoutError(obj, ~, x)
            y = obj.H * x;
        end

        function H = linearization(obj, ~, ~)
            H = obj.H;
        end

    end

end