function [xz, xzPlot] = normalNormalsq(N, mu, sigma, plotSampleNumber)
if nargin<2
    mu = [0,0];
end
if nargin<3
    sigma = [1,0.9;0.9,1];
end
if nargin<4
    plotSampleNumber = 1e7;
end

xz = generateSample(N, mu, sigma);
xzPlot = generateSample(plotSampleNumber, mu, sigma);
end


function x = generateSample(N, mu, sigma)
samples = mvnrnd(mu, sigma, N);
samples = samples.';
margMu1 = 1/N * sum(samples(1,:));
margMu2 = 1/N * sum(samples(2,:));
margStd1 = std(samples(1,:));
margStd2 = std(samples(2,:));
samples2 = [cdf('Normal', samples(1,:), margMu1, margStd1); cdf('Normal', samples(2,:), margMu2, margStd2)];
x = [icdf('Normal', samples2(1,:), 0,1); (icdf('Normal', samples2(2,:), 0, 1)).^2];
end