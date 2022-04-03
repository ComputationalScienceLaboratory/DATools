classdef Observation < handle

    properties
        NumVars
        ErrorModel
    end

    methods

        function obj = Observation(varargin)
            p = inputParser;
            addRequired(p, 'NumVars');
            addParameter(p, 'ErrorModel', datools.error.Error);
            parse(p, varargin{:});

            s = p.Results;

            obj.NumVars = s.NumVars;
            obj.ErrorModel = s.ErrorModel;
        end

        function y = observeWithError(obj, t, x)
            y = obj.ErrorModel.adderr(t, obj.observeWithoutError(t, x));
        end

        function y = observeWithoutError(~, ~, x)
            y = x;
        end

        function H = linearization(obj, ~)
            H = speye(obj.NumVars);
        end

    end

end
