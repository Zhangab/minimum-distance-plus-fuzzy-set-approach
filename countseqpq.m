function [P,Q]=countseqpq(S,gapratio)
%COUNTSEQPQ - Counts transition (P) and transversion (Q) for the given sequence pair S

% Molecular Biology & Evolution Toolbox, (C) 2006
% Author: James J. Cai
% Email: jamescai@hku.hk
% Website: http://bioinformatics.org/mbetoolbox/
% Last revision: 5/28/2005

if (isstruct(S)), S=S.seq; end
[n,m] = size(S);

TS = [0, 0, 1, 0;
      0, 0, 0, 1;
      1, 0, 0, 0;
      0, 1, 0, 0];

TV = [0, 1, 0, 1;
      1, 0, 1, 0;
      0, 1, 0, 1;
      1, 0, 1, 0];


P = zeros(n);
Q = zeros(n);


for i=1:n-1
for j=i+1:n
	[X,gap] = countntchange(S(i,:), S(j,:));
    
    if (gap./m<=gapratio)
	P(i,j) = sum(sum(TS.*X));% 数组乘法表示对应元素相乘。
	Q(i,j) = sum(sum(TV.*X));
    else
     P(i,j) = NaN ; % 
	 Q(i,j) = NaN ;
    end
P(j,i) = P(i,j);
Q(j,i) = Q(i,j);
end
end
