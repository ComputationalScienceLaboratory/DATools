function LF = gasparicohnsqrtadaptive(d, m, theta)

if nargin < 2
   m = []; 
end

if nargin < 3
    theta = 1;
end

if isempty(theta)
    theta = 1.7386;
end

LF = @(obj, k, r) csl.dataassimilation.statistical.ensemble.localisation.functions.gasparicohnCTilde( ...
    obj.Problem.NumVars, r, d, obj.Problem.TimeSpan(1), obj.Problem.Y0, m, k, theta);

end
