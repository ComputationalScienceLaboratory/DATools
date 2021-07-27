function RH(model, xt)
    %function to plot the rank histogram
        
    flag = isa(model, 'datools.statistical.ensemble.EnF');
    if ~flag
        frpintf('Model should be ensemble type');
        exit;
    end
    if nargin < 2
        fprintf('need to pass the truth at this time step');
        exit;
    end
    var = length(model.RankHistogram);
    % define variables
    %fprintf('Ens = %.3f\n', model.NumEnsemble);
    
    xf = model.Ensemble;
    [stateN,ensN] = size(xf);
    
    %rankVal = zeros(length(var), ensN+2);
   
   for i = 1:length(var)
       tempVar = var(i);
       model.RankValue(i, ensN+2) = tempVar;
       ensFor = sort(xf(tempVar, :));
       
       if xt(tempVar) < ensFor(1)
           model.RankValue(i, 1) = model.RankValue(i, 1) + 1;
       elseif xt(tempVar) > ensFor(length(ensFor))
           model.RankValue(i, ensN+1) = model.RankValue(i, ensN+1) + 1;
       else
           index = find(ensFor > xt(tempVar), 1);
           model.RankValue(i, index) = model.RankValue(i, index) + 1;
       end
       
   end
    
end



