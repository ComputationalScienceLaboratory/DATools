classdef ThreeDVar < datools.filter.variational.Variational
    %

    properties
        Name = 'ThreeDVar'
    end

    properties (Dependent)
        MeanEstimate
    end

    methods
        function analysis(obj, obs)

            % A constrained problem like double pendulum will need a
            % constraint coming in from the object.

            dB = obj.BDecomposition;

            y = obs.Uncertainty.Mean;
            R = obs.Uncertainty.Covariance;
            dR = decomposition(R, 'chol');

            % if strcmp(class(R), "decomposition")
            %     dR = R;
            % else
            %     dR = decomposition(R, 'chol');
            % end

            xb = obj.State;

            t = obj.Model.ODEModel.TimeSpan(end);

            H = @(x) obs.observeWithoutError(x);  % make it robust
            Hadjoint = @(x) obs.linearization(x)';  % make it robust

            J = @(x) cost(x, xb, dB, y, dR, H, Hadjoint);

            if strcmp(obj.OptAlg, 'lbfgs')
                opts = optimoptions('fmincon','Display','none', ...
                    'SpecifyObjectiveGradient',true, ...
                    'HessianApproximation', {'lbfgs', 50});
            else

                hessv = @(Hax, v) dB\v + Hax*(dR\(Hax'*v));
                opts = optimoptions('fmincon','Display','none', ...
                    'Algorithm','trust-region-reflective', ...
                    'SpecifyObjectiveGradient',true, ...
                    'HessianMultiplyFcn', hessv, ...
                    'SubproblemAlgorithm', 'cg');

            end

            xa = fmincon(J, xb, [], [], [], [], [], [], [], opts);

            obj.State = xa;

            function [c, g, hessinfo] = cost(x, xb, dB, y, dR, H, Hadjoint)

                dx = x - xb;
                dy = H(x) - y;
                Hax = Hadjoint(x);

                Bix = dB\dx;
                Riy = dR\dy;

                c = 0.5*((dx'*Bix) + (dy'*Riy));
                g = Bix + Hax*Riy;
                hessinfo = Hax;

            end


        end


        function x = get.MeanEstimate(obj)
            x = obj.State;
        end

    end

end
