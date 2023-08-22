classdef Tent < datools.error.Error

    properties
        State
        Scale
        Mu
    end

    methods

        function obj = Tent(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'NumVars', 1);
            addParameter(p, 'Scale', 2*sqrt(3));
            addParameter(p, 'Mu', 1.99999);
            parse(p, varargin{:});

            s = p.Results;

            obj.Scale = s.Scale;
            obj.Mu = s.Mu;

            kept = p.Unmatched;

            p = inputParser;
            addParameter(p, 'InitialState', rand(s.NumVars, 1));
            parse(p, kept);

            s = p.Results;

            obj.State = s.InitialState;
        end

        function xp = adderr(obj, ~, x)
            mu = obj.Mu;

            z = obj.State;

            z = mu * min(z, 1-z);

            obj.State = z;

            err = obj.Scale * (z - 0.5);

            xp = x + err;
        end

    end

end
