classdef SIR < datools.statistical.ensemble.EnF

    methods

        function analysis(obj, obs)

            tau = obj.Rejuvenation;

            xf = obj.Ensemble;
            xa = xf;

            ensN = obj.NumEnsemble;
            wf = obj.Weights;

            Hxf = obj.Observation.observeWithoutError(xf);
            %t0 = Hxf - y;

            %dR = decomposition(R, 'chol');
            %as = exp(-0.5*sum(t0.*(dR \ t0), 1)).';

            as = obs.log(Hxf).';
            m = max(as);
            as = exp(as-(m + log(sum(exp(as-m)))));

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

            obj.Model.update(0, obj.MeanEstimate);


        end

    end

end
