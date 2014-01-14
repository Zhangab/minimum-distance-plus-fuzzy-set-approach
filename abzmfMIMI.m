function [a,b] = abzmfMIMI(targetspecies, nearestneighbor,IntraVariation095,InterVariation005)
% abzmftdr is to estimate two parameters [a,b] of zmf based on a DNA sequences matrix of 
% a target species, and the DNA sequences matrix of its nearest neighbor
% via TDR method (by calling tdrfun2.m)
% [a, b] - parameters of zmf
% targetspecies - DNA matrix of a target species 
% nearestneighbor - DNA matrix of the nearest neighbor of the target
% species
%   

%addpath C:\MATLAB\R2009a\work2009\TDR3;


[m,n] = size(targetspecies);
[m2,n2] = size(nearestneighbor);


[meanD01,SD1,meanD02,SD2,meanD012,SD3,miniD012,maxD012] = twospeciesdist3(targetspecies,nearestneighbor)
D0 = dn_k2p(targetspecies);
if (miniD012>max(max(D0,[],1)))
b = miniD012;
a = max(max(D0,[],1));

else 
  %  a = meanD01;
  %  b = meanD012;
  %IntraVariation095,InterVariation005
   a = IntraVariation095;
   b = InterVariation005;
  %a=0.0114;
   %b=0.0212;
   
end



%a =100 - mean(TDRintra);
%b =100 - mean(TDRinter);












%rmpath C:\MATLAB\R2009a\work2009\TDR3;

end

