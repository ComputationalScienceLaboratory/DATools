classdef RHF < datools.statistical.ensemble.EnF
    
    properties
        Tail
        Truncate
        DiscretePoint
    end
    
    % Define the update method
    methods
        function obj = RHF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParameter(p, 'Tail', 'Flat');
            addParameter(p, 'Truncate', 10);
            addParameter(p, 'DiscretePoint', 1e4);
            parse(p, varargin{2:end});
            
            s = p.Results;
            
            kept = p.Unmatched;
            
            obj@datools.statistical.ensemble.EnF(varargin{1}, kept);
            obj.Tail = s.Tail;
            obj.Truncate = s.Truncate;
            obj.DiscretePoint = s.DiscretePoint;
        end
        
        function analysis(obj, R, y)
            inflation = obj.Inflation;
            tc = obj.Model.TimeSpan(1);
            
            % define length of gaussian tail
            gtl = obj.Truncate;
            discrete_pts = obj.DiscretePoint;
            
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
            
            ind1 = 1:ensN-1; ind2 = 2:ensN;
            
            % do for each state variable
            for i = 1:size(xf,1)
                [xfs, sortin] = sort(xf(i,:));
                
                % try inflating the background
                xfs = mean(xfs) + inflation*(xfs - mean(xfs));
                
                prior_ht = zeros(1, ensN + 1);
                % find the height of uniform part
                for ii = 2:ensN
                    prior_ht(ii) = 1/ ((ensN + 1)*abs((xfs(ii) - xfs(ii-1))));
                end
                
                Var = xfv(i);
                
                if strcmp(obj.Tail, 'Gaussian') == 1
                    muleft = xfs(1) - sqrt(2*Var)*erfinv(2/(ensN + 1) -1);
                    muright = xfs(end) - sqrt(2*Var)*erfinv(1 - 2/(ensN + 1));
                    mu = [muleft muright];
                else
                    fprintf('Error! Gaussian tail not provided!\n');
                    return;
                end
                                
                inn = y(i) - sort(Hxf(i,:));
                ntemp = length(y(i));
                likelihood = (diag((1/sqrt((2*pi)^ntemp*det(R(i,i)))) * ...
                    exp(-0.5*inn.'*(R(i,i)\inn)))).';
                likelihood = likelihood/min(likelihood); % for scaling
                l_ht = [likelihood(1), 0.5*(likelihood(ind1) + likelihood(ind2)), likelihood(end)];
                
                % find the posterior ht
                post_ht = zeros(1, ensN + 1);
                for jj = 2:ensN
                    post_ht(jj) = prior_ht(jj).*l_ht(jj);
                end
                
                % find total area to normalize
                area = [l_ht(1) * 0.5 * (1 + erf((xfs(1) - muleft)/(sqrt(2*Var)))),...
                    (xfs(ind2) - xfs(ind1)).*(post_ht(2:end-1)),...
                    l_ht(end).*(1 - 0.5*(1 + erf((xfs(end) - muright)/(sqrt(2*Var)))))];
                ta = sum(area);
                area = area/ta;
                
                post_ht = post_ht/ta;
                
                
                % find scaling actors for the tails
                tailscale = [l_ht(1) l_ht(end)];
                tailscale = tailscale/ta;
                % find the updated posterior points/particles
                tempa = findpos(post_ht, ensN, area, xfs, tailscale, Var, mu);
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

function pts = findpos(post_ht, ensN, area, xfs, tailscale, V, mu)
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
                pts(i) = sqrt(2*V)*erfinv(2*(A/tailscale(1)) - 1) + mu(1);
            elseif j == ensN+1
                pts(i) = sqrt(2*V)*erfinv(1 - 2*(A/tailscale(2))) + mu(2);
            else
                pts(i) = xfs(j-1) + A/post_ht(j);
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