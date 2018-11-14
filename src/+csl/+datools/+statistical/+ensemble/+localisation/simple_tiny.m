function LF = simple_tiny(r, d, m, theta)

if nargin < 3
   m = @(ri, rj) (ri + rj)/2; 
end

if nargin < 5
    theta = 1;
end

if isempty(theta)
    theta = 1.7386;
end

LF = @(obj, H) csl.dataassimilation.statistical.ensemble.localisation.functions.simplerho_tiny( ...
    obj.Problem.NumVars, r, d, obj.Problem.TimeSpan(1), obj.Problem.Y0, m, theta, H);

end


