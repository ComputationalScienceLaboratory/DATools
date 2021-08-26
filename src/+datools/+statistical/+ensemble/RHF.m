classdef RHF < datools.statistical.ensemble.EnF
    
    % Define the update method
    methods
        function analysis(obj, R, y)
            inflation = obj.Inflation;
            tc = obj.Model.TimeSpan(1);
            
            % define length of gaussian tail (make it a member of base class, if required)
            gtl = 20;
            discrete_pts = 1e5;
            
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            xfm = mean(xf, 2);
            xfv = var(xf, 0, 2);
            
            xas = zeros(size(xf));
            
            % Inflate the forecast
            Af = xf - repmat(xfm, 1, ensN);
            Af = inflation*Af;
            xf = repmat(xfm, 1, ensN) + Af;
            
            % Observable
            Hxf = obj.Observation.observeWithoutError(tc, xf);
            Hxfm = mean(Hxf, 2);
            
            % do for each state variable
            for i = 1:size(xf,1)
                [xfs, sortin] = sort(xf(i,:));
                
                prior_ht = cell(1, length(xfs)+1);
                % find the height of uniform part
                for ii = 2:length(xfs)
                    prior_ht{ii} = 1/ ((ensN + 1)*abs((xfs(ii) - xfs(ii-1))));
                end
                
                Var = xfv(i);
                xleft = linspace(xfs(1) - gtl, xfs(1), discrete_pts);
                xright = linspace(xfs(length(xfs)), xfs(length(xfs)) + gtl, discrete_pts);
                
                if strcmp(obj.Tail, 'Gaussian') == 1
                    options = optimoptions(@fminunc,'FunctionTolerance',1e-12, 'Display', 'none');
                    
                    pdfl = @(mu) (1/sqrt(Var*2*pi)) * exp(-0.5*(xleft-mu).^2/(Var));
                    minf = @(mu) (1/(ensN + 1) - riemannsum(xleft, pdfl(mu))).^2;
                    [muleft, ~] = fminunc(minf, xleft(end), options);
                    prior_ht{1} = pdfl(muleft);
                    
                    pdfr = @(mu) (1/sqrt(Var*2*pi)) * exp(-0.5*(xright-mu).^2/(Var));
                    minf = @(mu) (1/(ensN + 1) - riemannsum(xright, pdfr(mu))).^2;
                    [muright, ~] = fminunc(minf, xright(1), options);
                    prior_ht{length(xfs) + 1} = pdfr(muright);
                elseif strcmp(obj.Tail, 'Flat') == 1
                    prior_ht{1} = ones(1,discrete_pts)*1/((ensN+1) * gtl);
                    prior_ht{length(xfs) + 1} = ones(1,discrete_pts)*1/((ensN+1) * gtl);
                end
                                
                inn = y(i) - sort(Hxf(i,:));
                ntemp = length(y(i));
                likelihood = (diag((1/sqrt((2*pi)^ntemp*det(R(i,i)))) * ...
                    exp(-0.5*inn.'*(R(i,i)\inn)))).';
                l_ht = cell(1, ensN + 1);
                l_ht{1} = likelihood(1) * ones(1,length(prior_ht{1}));
                l_ht{ensN+1} = likelihood(ensN) *...
                    ones(1, length(prior_ht{ensN+1}));
                for ii = 2:ensN
                    l_ht{ii} = (likelihood(ii-1) + likelihood(ii))/2;
                end
                
                % find the posterior ht
                post_ht = cell(1, ensN + 1);
                for jj = 1:ensN+1
                    post_ht{jj} = prior_ht{jj}.*l_ht{jj};
                end
                
                % find total area to normalize
                area = zeros(1, ensN+1);
                area(1) = riemannsum(xleft, post_ht{1});
                area(ensN+1) = riemannsum(xright, post_ht{ensN+1});
                
                for kk = 2:ensN
                    area(kk) = post_ht{kk} * abs(xfs(kk) - xfs(kk-1));
                end
                for ll = 1:ensN+1
                    post_ht{ll} = post_ht{ll}/sum(area);
                end
                area = area/sum(area);
                % find the updated posterior points/particles
                tempa = findpos(post_ht, ensN, area, xfs, xleft, xright);
                %[val, in2] = sort(sortin);
                xas(i,:) = tempa(sortin);
            end
            % end of function
            
            obj.Ensemble = xas;
            obj.Model.update(0, obj.BestEstimate);
        end
        
        %end of method
    end
    %end of class
end

function pts = findpos(post_ht, ensN, area, xfs, xleft, xright)
pts = zeros(1,ensN);
for i =1:ensN
    temp = area;
    A = i/(ensN + 1);
    for j =1:ensN+1
        if temp(j) < A
            A = A - temp(j);
            temp(j) = 0;
        else
            if j == 1
                pts(i) = riemannsearch(xleft, post_ht{1}, A);
            elseif j == ensN+1
                pts(i) = riemannsearch(xright, post_ht{ensN+1}, A);
            else
                pts(i) = xfs(j-1) + A/post_ht{j};
                %temp(j) = temp(j) - A;
            end
            break;
        end
    end
end
end

function area = riemannsum(x, pdf)
dx = abs(x(2)-x(1));
temp = dx*pdf;
area = sum(temp(1:end-1));
end

function xp = riemannsearch(x, y, area)
dx = abs(x(2) - x(1));
darea = dx*y;
cumarea = cumsum(darea);
xp = x(find(cumarea - area>0, 1));
end