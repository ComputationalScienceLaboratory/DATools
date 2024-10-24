%function to find the KL divergence of two given distribution
% For more information look up:
% https://en.wikipedia.org/wiki/Kullback-Leibler_divergence

%pdf1 = sample pdf
%pdf2 = reference pdf(Uniform)

%pval = value of the polynomial
%KLval = kl div. value

function KLVal = KLDivergence(pdf1, pdf2)

if numel(pdf1) ~= numel(pdf2)
    error('Both vectors must be of same length');
end

% check if the distributions are normalized
if sum(pdf1) ~= 1
    pdf1 = pdf1 / sum(pdf1);
end
if sum(pdf2) ~= 1
    pdf2 = pdf2 / sum(pdf2);
end

KLVal = sum(pdf1.*(log(pdf1./pdf2)));

end