classdef Observation < handle
    %OBSERVATION This is the base class for all observation 
    %   Derive from this class and implement methods/functions as required
    %   Deriving from handle base class allows an object of this class to be
    %   passed by reference.

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