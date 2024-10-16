function RH(model, xt, y, Hxa, observationvariable)
%function to plot the rank histogram

flag = isa(model, 'datools.filter.ensemble.EnF');

if ~flag
    frpintf('Model should be ensemble type');
    exit;
end
if nargin < 2
    fprintf('need to pass the truth at this time step');
    exit;
end
var = length(model.RankHistogram);

xa = model.Ensemble;
[~, ensN] = size(xa);

% remove this
observationvariable = 'Truth';

switch observationvariable
    case 'Truth'
        for i = 1:var
            tempVar = model.RankHistogram(i);
            model.RankValue(i, ensN+2) = tempVar;
            ensFor = sort(xa(tempVar, :));
            
            if xt(tempVar) < ensFor(1)
                model.RankValue(i, 1) = model.RankValue(i, 1) + 1;
            elseif xt(tempVar) > ensFor(length(ensFor))
                model.RankValue(i, ensN+1) = model.RankValue(i, ensN+1) + 1;
            else
                index = find(ensFor > xt(tempVar), 1);
                model.RankValue(i, index) = model.RankValue(i, index) + 1;
            end
            
        end
    case 'Observation'
        for i = 1:numel(varobservtion)
            tempVar = model.RankHistogram(varobservtion(i));
            model.RankValue(i, ensN+2) = tempVar;
            ensFor = sort(Hxa(i, :));
            
            if y(i) < ensFor(1)
                model.RankValue(i, 1) = model.RankValue(i, 1) + 1;
            elseif y(i) > ensFor(length(ensFor))
                model.RankValue(i, ensN+1) = model.RankValue(i, ensN+1) + 1;
            else
                index = find(ensFor > y(i), 1);
                model.RankValue(i, index) = model.RankValue(i, index) + 1;
            end
            
        end
end


end
