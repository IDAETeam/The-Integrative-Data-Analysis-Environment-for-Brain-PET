    function [MSD,SampleNo]=getMSD(x);

% To return mean and SDs.
%
%   usage:  MSD=getMSD(x);
%	    [MSD,SampleNo]=getMSD(x);
%
% X may be a matrix.
% MSD=[mean,SD,coefficients of variance]

[SampleNo,n]=size(x); MSD=zeros(3,n);
MSD=[mean(x,1);std(x,0,1);std(x,0,1)./sqrt(n);std(x,0,1).*100./mean(x,1)];


