classdef RHF < datools.statistical.ensemble.EnF

    properties
        Tail
        Truncate
        DiscretePoint
    end

    methods
        function obj = RHF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;

            addParameter(p, 'Tail', 'Flat');
            addParameter(p, 'Truncate', 10);
            addParameter(p, 'DiscretePoint', 1e4);
            parse(p, varargin{7:end});

            s = p.Results;

            kept = p.Unmatched;

            obj@datools.statistical.ensemble.EnF(varargin{1:6}, kept);
            obj.Tail = s.Tail;
            obj.Truncate = s.Truncate;
            obj.DiscretePoint = s.DiscretePoint;
        end

        function analysis(obj, observation)
            %ANALYSIS   Method to overload the analysis function
            %
            %   ANALYSIS(OBJ) assimilates the current observation with the
            %   background/prior information to get a better estimate
            %   (analysis/posterior)
            
            inflation = obj.Inflation;
            tc = obj.Model.TimeSpan(1);

            % get the current ensemble forecast
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;
            
            R = observation.ErrorModel.Covariance;
            
            y = observation.Y;

            % calculate mean and variance of the forecast
            xfm = mean(xf, 2);
            xfv = var(xf, 0, 2);

            xas = zeros(size(xf));

            % Inflate the forecast
            Af = xf - repmat(xfm, 1, ensN);
            Af = inflation * Af;
            xf = repmat(xfm, 1, ensN) + Af;

            % Observable
            Hxf = observation.observeWithoutError(tc, xf);
            % mean of observable
            Hxfm = mean(Hxf, 2);

            % will be used for slicing
            ind1 = 1:ensN - 1;
            ind2 = 2:ensN;

            % assimilate each state variable
            for i = 1:size(xf, 1)
                % sort the forecast ensemble
                %xfs = sorted forecast ensemble
                [xfs, sortin] = sort(xf(i, :));

                priorht = zeros(1, ensN+1);
                % find the height of uniform part
                for ii = 2:ensN
                    priorht(ii) = 1 / ((ensN + 1) * abs((xfs(ii) - xfs(ii-1))));
                end

                Var = xfv(i);

                if strcmp(obj.Tail, 'Gaussian') == 1
                    % find the mean of the left gaussian tail
                    muleft = xfs(1) - sqrt(2*Var) * erfinv(2/(ensN + 1)-1);
                    % find the mean of the right gaussian tail
                    muright = xfs(end) - sqrt(2*Var) * erfinv(1-2/(ensN + 1));
                    mu = [muleft, muright];
                else
                    error('Only Gaussian Tail is supported as of now');
                end

                % calculate the likelihood for this state
                d = y(i) - sort(Hxf(i, :));
                ntemp = length(y(i));
                likelihood = (diag((1 / sqrt((2 * pi)^ntemp*R(i, i)))* ...
                    exp(-0.5*d.'*(R(i, i) \ d)))).';
                likelihood = likelihood / min(likelihood); % for scaling
                likelihood = [likelihood(1), 0.5 * (likelihood(ind1) + likelihood(ind2)), likelihood(end)];

                % find the posterior ht
                postht = zeros(1, ensN+1);
                for jj = 2:ensN
                    postht(jj) = priorht(jj) .* likelihood(jj);
                end

                % find total area to normalize everything
                area = [likelihood(1) * 0.5 * (1 + erf((xfs(1) - muleft)/(sqrt(2*Var)))), ...
                    (xfs(ind2) - xfs(ind1)) .* (postht(2:end-1)), ...
                    likelihood(end) .* (1 - 0.5 * (1 + erf((xfs(end) - muright)/(sqrt(2*Var)))))];
                % find total sum of area
                ta = sum(area);
                % Normalize the area
                area = area / ta;

                % normalize the posterior height
                postht = postht / ta;

                % find scaling factors for the tails
                tailscale = [likelihood(1), likelihood(end)];
                tailscale = tailscale / ta;
                % find the updated posterior points/particles
                tempa = findpos(postht, ensN, area, xfs, tailscale, Var, mu);
                xas(i, :) = tempa(sortin);
            end
            % end of function

            obj.Ensemble = xas;
            obj.Model.update(0, obj.BestEstimate);
        end

        %end of method
    end
    %end of class
end

function pts = findpos(postht, ensN, area, xfs, tailscale, V, mu)
pts = zeros(1, ensN);
for i = 1:ensN
    temp = area;
    A = i / (ensN + 1);
    for j = 1:ensN + 1
        if temp(j) < A
            A = A - temp(j);
            temp(j) = 0;
        else
            if j == 1
                pts(i) = sqrt(2*V) * erfinv(2*(A / tailscale(1))-1) + mu(1);
            elseif j == ensN + 1
                pts(i) = sqrt(2*V) * erfinv(1-2*(A / tailscale(2))) + mu(2);
            else
                pts(i) = xfs(j-1) + A / postht(j);
            end
            break;
        end
    end
end
end


% extra functions not useful at this time.
function area = riemannsum(x, pdf)
dx = abs(x(2)-x(1));
temp = dx * pdf;
area = sum(temp(1:end-1));
end

function xp = riemannsearch(x, y, area)
dx = abs(x(2)-x(1));
darea = dx * y;
cumarea = cumsum(darea);
xp = x(find(cumarea-area > 0, 1));
end
