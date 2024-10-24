function [xz, xzPlot] = bimodalGaussian(N, mus, sigmas, ratio)
if nargin<2
    mus = [3,3; 10,10];
    sigmas(:,:,1) = [5 1.5; 1.5 3];
    sigmas(:,:,2) = [3 -1.5; -1.5 3];
    ratio = [0.5,0.5];
end
if nargin<3
    sigmas(:,:,1) = [5 1.5; 1.5 3];
    sigmas(:,:,2) = [3 -1.5; -1.5 3];
    ratio = [0.5, 0.5];
end
if nargin<4
    ratio = [0.5, 0.5];
end

mu1 = mus(1,:);
mu2 = mus(2,:);
sigma1 = sigmas(:,:,1);
sigma2 = sigmas(:,:,2);
ratio1 = ratio(1);
ratio2 = ratio(2);

xz = [mvnrnd(mu1,sigma1, ratio1*N); mvnrnd(mu2, sigma2, ratio2*N)].';
xzPlot = xz;
end