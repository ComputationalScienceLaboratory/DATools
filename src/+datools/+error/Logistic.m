classdef Logistic < datools.error.Uncertainty

    properties
        State
        Scale
        R
    end

    methods

        function obj = Logistic(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'NumVars', 1);
            addParameter(p, 'Scale', 2*sqrt(3));
            addParameter(p, 'R', 3.99999);
            parse(p, varargin{:});

            s = p.Results;

            obj.Scale = s.Scale;
            obj.R = s.R;

            kept = p.Unmatched;

            p = inputParser;
            addParameter(p, 'InitialState', rand(s.NumVars, 1));
            parse(p, kept);

            s = p.Results;

            obj.State = s.InitialState;

        end

        function xp = adderr(obj, ~, x)
            r = obj.R;

            z = obj.State;

            z = r * z .* (1 - z);

            obj.State = z;

            err = obj.Scale * (z - 0.5);

            xp = x + err;
        end

    end

end
