classdef Observation < handle
%OBSERVATION Base class for every observation operator
%   other details 

    properties
        Y             % Current Observation of States (maybe sparse)
        NumVars      % number of state variables
        ErrorModel   % observation error object
        NumObs        % Number of Observations(Make sure > 0)
    end

    methods

        function obj = Observation(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'Y', []);
            addRequired(p, 'NumVars');
            addParameter(p, 'ErrorModel', []);
            validationFunction = @(x) (x>0);
            addParameter(p, 'NumObs', 1, validationFunction);

            parse(p, varargin{:});

            s = p.Results;

            obj.NumVars = s.NumVars;
            obj.ErrorModel = s.ErrorModel;
            obj.Y = s.Y;
            obj.NumObs = s.NumObs;
        end

        function y = observeWithError(obj, t, x)
            y = obj.ErrorModel.addError(t, obj.observeWithoutError(t, x));
        end

        function updateY(obj, y)
            obj.Y = y;
        end

        function y = observeWithoutError(~, ~, x)
            y = x;
        end

        function H = linearization(obj, ~)
            H = speye(obj.NumVars);
        end

    end

end
