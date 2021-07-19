classdef Hybrid < datools.statistical.ensemble.EnF
    
    properties
        Alphas
        Filters
    end
    
    methods
        function obj = Hybrid(varargin)
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'Alphas', 1);
            addParameter(p, 'Filters', {});
            
            parse(p, varargin{2:end});
            
            s = p.Results;
            
            kept = p.Unmatched;
            
            obj@datools.statistical.ensemble.EnF(varargin{1}, kept);
            
            obj.Alphas = s.Alphas;
            obj.Filters = s.Filters;
            
        end
    end
    
    
    methods
        
        function analysis(obj, R, y)
            
            tau = obj.Rejuvenation;
            
            X = obj.Ensemble;
            alphas = obj.Alphas;
            
            Xf = X;
            
            for i = 1:numel(obj.Filters)
                alpha = alphas(i);
                if alpha > 0
                    obj.Filters{i}.Weights = obj.Weights;
                    obj.Filters{i}.Ensemble = X;
                    obj.Filters{i}.analysis(R/alpha, y);
                    obj.Weights = obj.Filters{i}.Weights;
                    X = obj.Filters{i}.Ensemble;
                end
            end
            
            obj.Ensemble = X;
            obj.rejuvenate(tau);
            
            obj.Model.update(0, obj.BestEstimate);

        end
        
    end
    
end
