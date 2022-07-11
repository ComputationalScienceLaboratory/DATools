classdef ThreeDVar < handle
    % 

    properties
        Model % type of ODE solver (ode45/Runge Kutta) and the model (eg: Lorenz63)
        ModelError % type err
        Observation % type of obervation
        State % current state values for all the states
        B % background covariance
    end

    properties (Dependent)
        BestEstimate % Current estimate of the particles/ensembles
    end

    properties (Access=protected)
        BDecomposition
    end

    methods
        function obj = ThreeDVar(varargin)
            %ENF   The constructor initializes the properties/attributes
            %
            %   OBJ = ENF(VARARGIN) accepts variable length argument list
            %   VARARGIN and updates the properties/attributes of the object
            %   (OBJ) of this class or a derived class

            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p, 'Model', @(x) isa(x, 'datools.Model'));
            addParameter(p, 'ModelError', datools.error.Error);
            addParameter(p, 'InitialState', []);
            addParameter(p, 'BackgroundCovariance', []);
            parse(p, varargin{:});

            s = p.Results;

            obj.Model = s.Model;
            obj.ModelError = s.ModelError;
            obj.State = s.InitialState;
            obj.B = s.BackgroundCovariance;
            obj.BDecomposition = decomposition(obj.B, 'chol');

            kept = p.Unmatched;

            p = inputParser;
            addParameter(p, 'Observation', datools.observation.Observation(s.Model.NumVars));
            parse(p, kept);

            s = p.Results;

            obj.Observation = s.Observation;
        end

        function forecast(obj)
            %FORECAST Method to propagate the model forward in time
            %
            %   FORECAST(OBJ) propoagates the model one step in time
            %   using a user defined time integration method

            [~ , yend] = obj.Model.solve([], obj.State);

            obj.State = obj.ModelError.adderr(obj.Model.TimeSpan(end), yend);


        end

        function analysis(obj, R, y)
            dB = obj.BDecomposition;
            dR = decomposition(R, 'chol');
            
            xb = obj.State;
            obs = obj.Observation;

            t = obj.Model.TimeSpan(end);

            H = @(x) obs.observeWithoutError(t, x);
            Hadjoint = @(x) obs.linearization(t, x).';
            
            J = @(x) cost(x, xb, dB, y, dR, H, Hadjoint);

            opts = optimoptions('fmincon','Display','none','SpecifyObjectiveGradient',true, ...
                'HessianApproximation', {'lbfgs', 50});
            xa = fmincon(J, xb, [],[],[],[],[],[],[], opts);

            obj.State = xa;

            function [c, g] = cost(x, xb, dB, y, dR, H, Hadjoint)
                
                Hx = H(x);

                Bix = dB\(x - xb);

                Riy = dR\(Hx - y);
                
                c = 0.5*(x - xb).'*Bix + 0.5*(Hx - y).'*Riy;

                g = Bix + Hadjoint(x)*Riy;
                
            end


        end

        function x = get.BestEstimate(obj)
            %GET.BESTESTIMATE Method to estimate of ensemble values
            %
            %   X = GET.BESTESTIMATE(OBJ) uses an in-built getter method,
            %   derived from handle class, to return the best estimate of
            %   the information from the current ensembles of states and
            %   its corresponding weights

            x = obj.State;

        end

    end

end
