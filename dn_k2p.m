function [D,VarD]=dn_k2p(aln,alpha)
%DN_K2P - Kimura 80 Distance
%
% Syntax: [D,VarD]=dn_k2p(aln,alpha)
%
% Inputs:
%    aln          - Alignment structure
%    alpha        - (optional) Shape parameter of gamma distribution
%
% Outputs:
%    D      - Distance matrix
%    VarD   - Variance of distance
%
% See also: 

% Molecular Biology & Evolution Toolbox, (C) 2006
% Author: James J. Cai
% Email: jamescai@hku.hk
% Website: http://bioinformatics.org/mbetoolbox/
% Last revision: 5/28/2005


% p1=(F(1,2)+F(2,1)+F(3,4)+F(4,3) )/N;
% p2= sum(sum(F(1:2,3:4)+F(3:4,1:2)))/N;

if (isstruct(aln))
    sq=aln.seq;
else
    sq=aln;
end
gapratio = 0.90;% gapratio=gaps/lengthofseq

[n,m] = size(sq);
[P,Q]=countseqpq(aln,gapratio);
P=P./m; Q=Q./m;
W1 = 1-2*P-Q;
W2 = 1-2*Q;

% where p1 and p2 are the frequencies of sites with transitional and
% transversional differences, respectively.

if(nargin==1)
	D=(-1/2)*i_safelog(W1)-(1/4)*i_safelog(W2);
elseif(nargin==2)
	D=(alpha/2)*((W1.^(-1/alpha) + (1/2)*W2.^(-1/alpha)-3/2));
end

if (nargout==2)
        W1=1./W1; W2=1./W2;
	W3=(W1+W2)./2;
	VarD=((W1.^2).*P + (W3.^2).*Q - (W1.*P+W3.*Q).^2)./m;
end

