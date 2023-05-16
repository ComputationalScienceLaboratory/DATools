classdef SIR < datools.statistical.ensemble.EnF

    methods

        function analysis(obj, observation)
            %ANALYSIS   Method to overload the analysis function
            %
            %   ANALYSIS(OBJ) assimilates the current observation with the
            %   background/prior information to get a better estimate
            %   (analysis/posterior)
            
            tau = obj.Rejuvenation;
            tc = obj.Model.TimeSpan(1);

            xf = obj.Ensemble;
            xa = xf;

            ensN = obj.NumEnsemble;
            wf = obj.Weights;
            
            R = observation.ErrorModel.Covariance;

            Hxf = observation.observeWithoutError(tc, xf);
            t0 = Hxf - observation.Y;

            dR = decomposition(R, 'chol');
            as = exp(-0.5*sum(t0.*(dR \ t0), 1)).';

            w = wf .* as;
            w = w / sum(w);

            ensEff = 1 / sum(w.^2);

            if ensEff < ensN * obj.ResamplingThreshold

                what = cumsum(w);
                a = rand(1, ensN);

                for i = 1:ensN
                    ind = find(a(i) < what, 1);
                    xa(:, i) = xf(:, ind);
                end

                w = ones(ensN, 1) / ensN;

            end


            obj.Ensemble = xa;
            obj.Weights = w;
            obj.rejuvenate(tau, xf);

            obj.Model.update(0, obj.BestEstimate);


        end

    end

end
