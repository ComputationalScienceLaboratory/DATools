classdef Variational < datools.DABase
    % This is the base class for variational methods. Derive from this class
    % and implements methods as required. Abstract methods and properties
    % should be defined in the child class

    properties
        Model % type of ODE solver (ode45/Runge Kutta) and the model (eg: Lorenz63)
        State % current values for all the states
        B % background covariance
        OptAlg % lbfgs or newton (we can try Poblano too!!)
        CovarianceEstimate % needed bc of class abstraction
        Type = "Variational"
    end

    properties (Access = protected)
        BDecomposition % Cholesky decomposition of B
    end

    properties (Abstract)
        Name % name of the filter
    end

    methods (Abstract)
        analysis(obj, obs) % assimilation algorithm
    end

    methods
        function obj = Variational(varargin)
            % The constructor initializes the properties/attributes
            %
            %   OBJ = Variational(VARARGIN) accepts variable length argument list
            %   VARARGIN and updates the properties/attributes of the object
            %   (OBJ) of this class or a derived class
            optAlgValFcn = @(x) (strcmp(x, 'lbfgs') || strcmp(x, 'newton'));

            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p, 'Model');
            addParameter(p, 'InitialState', []);
            addParameter(p, 'BackgroundCovariance', []);
            addParameter(p, 'OptimizationType', 'lbfgs', optAlgValFcn);

            addParameter(p, 'CovarianceEstimate', []);
            parse(p, varargin{:});

            s = p.Results;

            obj.Model = s.Model;
            obj.State = s.InitialState;
            obj.B = s.BackgroundCovariance;
            obj.BDecomposition = decomposition(obj.B, 'chol');
            obj.OptAlg = s.OptimizationType;
        end

        function forecast(obj)
            %FORECAST Method to propagate the model forward in time
            %
            %   FORECAST(OBJ) propoagates the model one step in time
            %   using a user defined time integration method

            [time , yend] = obj.Model.solve([], obj.State);
            obj.State = obj.Model.Uncertainty.addError(yend);

            obj.Model.update(time, obj.State);
        end
    end

end