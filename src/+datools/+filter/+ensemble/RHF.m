classdef RHF < datools.filter.ensemble.EnF
    % Rank Histogram Filter (Anderson)
    % citation/reference
    properties
        Tail % the kind of tail it is
        Truncate % equal truncation length on each flat tail
        LeftBound % lower bound on the tail
        RightBound % upper bound on the tail
        Name = "Rank Histogram Filter"
    end

    methods
        function obj = RHF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;

            addParameter(p, 'Tail', 'Flat');
            addParameter(p, 'Truncate', 10);
            addParameter(p, 'LeftBound', -inf);
            addParameter(p, 'RightBound', inf);
            parse(p, varargin{8:end});

            s = p.Results;

            kept = p.Unmatched;

            obj@datools.filter.ensemble.EnF(varargin{1:7}, kept);
            obj.Tail = s.Tail;
            obj.Truncate = s.Truncate;
            obj.LeftBound = s.LeftBound;
            obj.RightBound = s.RightBound;
        end

        function analysis(obj, obs)
            %ANALYSIS   Method to overload the analysis function
            %
            %   ANALYSIS(OBJ) assimilates the current observation with the
            %   background/prior information to get a better estimate
            %   (analysis/posterior)

            inflation = obj.Inflation;

            % get the current ensemble forecast
            xf = obj.Ensemble;
            ensN = obj.NumEnsemble;

            R = obs.Uncertainty.Covariance;

            y = obs.Uncertainty.Mean;

            % calculate mean and variance of the forecast
            xfm = mean(xf, 2);
            xfv = var(xf, 0, 2);

            xas = zeros(size(xf));

            % Inflate the forecast
            Af = xf - repmat(xfm, 1, ensN);
            Af = inflation * Af;
            xf = repmat(xfm, 1, ensN) + Af;

            % Observable
            Hxf = obs.observeWithoutError(xf);
            % mean of observable
            Hxfm = mean(Hxf, 2);

            % will be used for slicing
            ind1 = 1:ensN - 1;
            ind2 = 2:ensN;

            % assimilate each state variable according to the tail
            if strcmp(obj.Tail, 'Gaussian') == 1
                for i = 1:size(xf, 1) % refactor to make only observed variable
                    % sort the forecast ensemble
                    [xfs, sortIndex] = sort(xf(i, :));
                    [~, index2] = sort(sortIndex);

                    priorht = zeros(1, ensN+1);
                    % find the height of uniform part
                    for ii = 2:ensN
                        priorht(ii) = 1 / ((ensN + 1) * abs((xfs(ii) - xfs(ii-1))));
                    end

                    Var = xfv(i);

                    % find the mean of the left gaussian tail
                    muleft = xfs(1) - sqrt(2*Var) * erfinv(2/(ensN + 1)-1);
                    % find the mean of the right gaussian tail
                    muright = xfs(end) - sqrt(2*Var) * erfinv(1-2/(ensN + 1));
                    mu = [muleft, muright];

                    % calculate the likelihood for this state
                    % we assume that observations for each state are considered
                    % independent
                    d = y(i) - sort(Hxf(i, :));
                    likelihood = 1 / sqrt(2*R(i, i)*pi) * exp(-0.5*(d.^2)/R(i, i));
                    % likelihood = likelihood / min(likelihood); % for scaling
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
                    tempa = findPosGauss(postht, ensN, area, xfs, tailscale, Var, mu);
                    xas(i, :) = tempa(index2);
                end
            elseif strcmp(obj.Tail, 'Flat') == 1
                for i = 1:size(xf, 1) % refactor to make only observed variable
                    % sort the forecast ensemble
                    [xfs, sortIndex] = sort(xf(i, :));
                    [~, index2] = sort(sortIndex);

                    % tail length
                    tailLengthLeft = 1 * (xfs(end) - xfs(1));
                    tailLengthRight = 1 * (xfs(end) - xfs(1));
                    if xfs(1) - tailLengthLeft < obj.LeftBound
                        tailLengthLeft = xfs(1) - obj.LeftBound;
                    end

                    if xfs(end)+tailLengthRight>obj.RightBound
                        tailLengthRight = obj.RightBound - xfs(end);
                    end

                    % length of the entire domain
                    lens = [tailLengthLeft, xfs(ind2) - xfs(ind1), tailLengthRight];

                    % calculate the height of the prior rectangles
                    priorht = 1 ./ (lens * (ensN + 1));

                    % calculate the likelihood for this state
                    % we assume that observations for each state are considered
                    % independent
                    d = y(i) - sort(Hxf(i, :));
                    likelihood = 1 / sqrt(2*R(i, i)*pi) * exp(-0.5*(d.^2)/R(i, i));
                    % likelihood = likelihood / min(likelihood); % for scaling
                    likelihood = [likelihood(1), 0.5 * (likelihood(ind1) + likelihood(ind2)), likelihood(end)];

                    % find the posterior height of the uniform parts
                    postht = priorht .* likelihood;
                    area = postht .* lens;

                    % normalizethe area
                    ta = sum(area);
                    area = area / ta;
                    postht = postht / ta;

                    % find the updated posterior points/particles
                    tempa = findPosFlat(postht, ensN, area, xfs, tailLengthLeft);
                    xas(i, :) = tempa(index2);
                end

            end

            % end of function

            obj.Ensemble = xas;
            obj.Model.update(0, obj.MeanEstimate);
        end

        %end of method
    end
    %end of class
end

function pts = findPosGauss(postht, N, area, xfs, tailscale, V, mu)
pts = zeros(1, N);
for i = 1:N
    temp = area;
    A = i / (N + 1);
    for j = 1:N + 1
        if temp(j) < A
            A = A - temp(j);
            temp(j) = 0;
        else
            if j == 1
                pts(i) = sqrt(2*V) * erfinv(2*(A / tailscale(1))-1) + mu(1);
            elseif j == N + 1
                pts(i) = sqrt(2*V) * erfinv(1-2*(A / tailscale(2))) + mu(2);
            else
                pts(i) = xfs(j-1) + A / postht(j);
            end
            break;
        end
    end
end
end

function pts = findPosFlat(postht, N, area, xfs, tailLength)
pts = zeros(1, N);
for i = 1:N
    temp = area;
    A = i / (N + 1);
    for j = 1:N + 1
        if temp(j) < A
            A = A - temp(j);
            temp(j) = 0;
        else
            if j == 1
                pts(i) = xfs(1) - tailLength + A / postht(j);
            else
                pts(i) = xfs(j-1) + A / postht(j);
            end
            break;
        end
    end
end
end
