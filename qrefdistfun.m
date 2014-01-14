 function [minid,meanminid,maxd] = qrefdistfun(Refseqmatrix,que)

% % qrefdistfun is to generate three different distance of a que sequence to the reference matrix.
% note we use D0=dn_k2p(seqmat) here!!!!

% Refseqmatrix - A DNA sequences matrix of reference. 
% que - A single DNA query sequence.

% 
%addpath D:\Matlab6p5\work\GeneticDist;
%addpath C:\Matlab6p5\work\GeneticDist;
%temp;



%Refseqmatrix=Ref;
%que = Queryseqs.seq(1,:);



seqmat=cat(1,Refseqmatrix,que);
D0=dn_k2p(seqmat);
[mD0,nD0]=size(D0);

Dqr=D0(end,:);
Dqr1=Dqr(1:nD0-1);
[C,I]=max(Dqr1);

[C2,I2]=min(Dqr1);

meanminid=mean(Dqr1);

minid=C2;

maxd=C;



