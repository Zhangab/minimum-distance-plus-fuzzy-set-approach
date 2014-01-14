function [meanD01,SD1,meanD02,SD2,meanD012,SD3,miniD012,maxD012] = twospeciesdist3(seqmatrix1,seqmatrix2)

% twospeciesdist is to produce intraspecific distances and interspecific distance and their std for two species
% by taking two sequence matrices of two species as input.

% meanD01-average intraspecific distance of spe1,SD1-std;
% meanD02-average intraspecific distance of spe2,SD2-std;
% meanD012-average interspecific distance of spe1 and spe2,SD3-std.
% miniD012-mini  interspecific distance of spe1 and spe2 via pair-wise
% comparision
% maxD012 - maximum interspecific distance of spe1 and spe2 via pair-wise
% comparision


% Note: the default genetic distance is k2p!.

% the following functions will be called: dn_k2p();table2row2().

[n1,m1]=size(seqmatrix1);
[n2,m2]=size(seqmatrix2);

seqmat=cat(1,seqmatrix1,seqmatrix2);

%since function dn_k2p() is located at the following folder!
%addpath d:\MATLAB6p5\work\gdsi;

D0=dn_k2p(seqmat);
%D0=dn_ntdiff(seqmat)

% substract intraspecific distance matrix for spe1 from D0.

% D01=D0([1:n1],[1:n2]); modified to the following line at 2009-7-26 9:33

D01=D0([1:n1],[1:n1]); % modified  at 2009-7-26 9:33


% substract intraspecific distance matrix for spe2 from D0.

D02=D0([end-n2+1:end],[end-n2+1:end]);

% substract interspecific distance matrix between spe1 and spe2 from D0.

D012=D0([end-n2+1:end],[1:n1]);

% summarize the results for output.

[non_diag01, diagele] = table2row2(D01);
[non_diag02, diagele] = table2row2(D02);

meanD01=mean(non_diag01);
SD1=std(non_diag01);

meanD02=mean(non_diag02);
SD2=std(non_diag02);

D12=reshape(D012,1,n1*n2);

meanD012=mean(D12);
SD3=std(D12);

[miniD012,I2]=min(D12);
[maxD012,I3]=max(D12);
