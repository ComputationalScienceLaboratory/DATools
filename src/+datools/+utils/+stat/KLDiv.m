%function to find the KL divergence of two given distribution
% For more information look up: 
% https://en.wikipedia.org/wiki/Kullback-Leibler_divergence

%pdf1 = sample pdf
%pdf2 = reference pdf(Uniform)

%pval = value of the polynomial
%KLval = kl div. value
function [xs, pval, KLVal] = KLDiv(pdf1, pdf2)

if numel(pdf1) ~= numel(pdf2)
    fprintf('Both vectors must be of same length');
    return;
end

% check if the distributions are normalized
if sum(pdf1) ~= 1
    pdf1 = pdf1/sum(pdf1);
end
if sum(pdf2) ~= 1
    pdf2 = pdf2/sum(pdf2);
end

p = polyfit(1:1:length(pdf1), pdf1, 2);

%p = polyfit(linspace(0, 1, numel(pdf1)), pdf1, 2);

%a = p(1);

%KLVal =  sum(pdf1.*(log(pdf1./pdf2))) * sign(p(1));

xs = linspace(0, 1, numel(pdf1));
ys = numel(pdf1)*pdf1;

J = @(a) sum( ((a * xs.^2 - a * xs + 0.25*a + 1) - ys).^2 );
options = optimoptions(@fminunc, 'Display', 'none');
a = fminunc(J, 0, options);
%KLVal = a;
KLVal =  sum(pdf1.*(log(pdf1./pdf2))) * sign(a);
pval = a * xs.^2 - a * xs + 0.25*a + 1;

end
