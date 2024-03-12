classdef Observation < handle

    properties
        NumVars
        Uncertainty
    end

    methods

        function obj = Observation(varargin)
            p = inputParser;
            addRequired(p, 'NumVars');
            addParameter(p, 'Uncertainty', datools.uncertainty.NoUncertainty);
            parse(p, varargin{:});

            s = p.Results;

            obj.NumVars = s.NumVars;
            obj.Uncertainty = s.Uncertainty;
        end

        function y = observeWithError(obj, x)
            y = obj.Uncertainty.addError(obj.observeWithoutError(x));
        end

        function y = observeWithoutError(~, x)
            y = x;
        end

        function H = linearization(obj, ~)
            H = speye(obj.NumVars);
        end

    end

end
