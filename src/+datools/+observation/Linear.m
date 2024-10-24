classdef Linear < datools.observation.Observation
    % Derived from the "observation" parent class. This can be chosen when
    % the observation operation is linear.
    properties
        H
    end

    methods

        function obj = Linear(nvars, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'H', speye(nvars));
            parse(p, varargin{:});

            s = p.Results;

            obj@datools.observation.Observation(nvars, p.Unmatched);

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