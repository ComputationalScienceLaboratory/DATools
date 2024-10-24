classdef BPF < datools.filter.ensemble.EnF
    % Bootstrap Particle Filter
    % citations/references
    properties
        Name = "Bootstrap Particle Filter"
    end

    methods

        function analysis(obj, obs)
            %ANALYSIS   Method to overload the analysis function
            %
            %   ANALYSIS(OBJ) assimilates the current observation with the
            %   background/prior information to get a better estimate
            %   (analysis/posterior)

            tau = obj.Rejuvenation;

            xf = obj.Ensemble;
            xa = xf;

            ensN = obj.NumEnsemble;
            wf = obj.Weights;
            
            R = obs.Uncertainty.Covariance;

            Hxf = obs.observeWithoutError(xf);

            as = obs.Uncertainty.log(Hxf).';
            m = max(as);
            as = exp(as-(m + log(sum(exp(as-m)))));

            w = wf.*as;
            w = w/sum(w);
            
            if any(isnan(w))
                w = ones(ensN, 1) / ensN;
            end

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

            obj.Model.update(0, obj.MeanEstimate); % t=0 here because we do analysis


        end

    end

end
