classdef Triangle < datools.error.Uncertainty

    properties
        Lower
        Middle
        Upper
    end

    methods

        function obj = Triangle(varargin)
            p = inputParser;
            addOptional(p, 'Lower', -1);
            addOptional(p, 'Middle', 0);
            addOptional(p, 'Upper', 1);

            parse(p, varargin{:});

            s = p.Results;

            obj.Lower = s.Lower;
            obj.Middle = s.Middle;
            obj.Upper = s.Upper;
        end

        function xp = adderr(obj, ~, x)
            td = makedist('Triangular', obj.Lower, obj.Middle, obj.Upper);

            A = td.random([size(x, 1), size(x, 2)]);
            xp = x + A;
        end

    end

end
