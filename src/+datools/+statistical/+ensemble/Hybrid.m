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
                obj.Filters{i}.Weights = obj.Weights;
                obj.Filters{i}.Ensemble = X;
                obj.Filters{i}.analysis(R/alpha, y);
                obj.Weights = obj.Filters{i}.Weights;
                X = obj.Filters{i}.Ensemble;
            end
            
            [n, ensN]  = size(X);
            
            if n < ensN + 2
                Af = (Xf - mean(Xf))/sqrt(ensN -1);
                                
                vs = sqrt(sum(Af.^2, 2));
                
                Xi = sqrt(tau)*vs.*rand(n, ensN);
                Xi = Xi - mean(Xi, 2);

                X = X + Xi;
            else
                P = sqrt(tau/(ensN - 1))*(eye(ensN) - ones(ensN)/ensN)*randn(ensN)*(eye(ensN) - ones(ensN)/ensN);
                X = X + Xf*P;
            end
            
            obj.Ensemble = X;

        end
        
    end
    
end
