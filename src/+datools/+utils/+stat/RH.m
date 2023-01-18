function RH(filter, xt, y, Hxa, observationvariable)
%function to plot the rank histogram

flag = isa(filter, 'datools.statistical.ensemble.EnF');
if ~flag
    frpintf('Model should be ensemble type');
    exit;
end
if nargin < 2
    fprintf('need to pass the truth at this time step');
    exit;
end
var = length(filter.RankHistogram);
varobservtion = filter.Observation.Indices;
xa = filter.Ensemble;
[~, ensN] = size(xa);


switch observationvariable
    case 'Truth'
        for i = 1:var
            tempVar = filter.RankHistogram(i);
            filter.RankValue(i, ensN+2) = tempVar;
            ensFor = sort(xa(tempVar, :));
            
            if xt(tempVar) < ensFor(1)
                filter.RankValue(i, 1) = filter.RankValue(i, 1) + 1;
            elseif xt(tempVar) > ensFor(length(ensFor))
                filter.RankValue(i, ensN+1) = filter.RankValue(i, ensN+1) + 1;
            else
                index = find(ensFor > xt(tempVar), 1);
                filter.RankValue(i, index) = filter.RankValue(i, index) + 1;
            end
            
        end
    case 'Observation'
        for i = 1:numel(varobservtion)
            tempVar = filter.RankHistogram(varobservtion(i));
            filter.RankValue(i, ensN+2) = tempVar;
            ensFor = sort(Hxa(i, :));
            
            if y(i) < ensFor(1)
                filter.RankValue(i, 1) = filter.RankValue(i, 1) + 1;
            elseif y(i) > ensFor(length(ensFor))
                filter.RankValue(i, ensN+1) = filter.RankValue(i, ensN+1) + 1;
            else
                index = find(ensFor > y(i), 1);
                filter.RankValue(i, index) = filter.RankValue(i, index) + 1;
            end
            
        end
end


end
