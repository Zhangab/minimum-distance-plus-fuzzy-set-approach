function [NT,AA] = seqcode()
%SEQCODE - Return vector for mapping sequence letters to integers
%
% Syntax: [NT,AA] = seqcode
%
% See also:

% Molecular Biology & Evolution Toolbox, (C) 2006
% Author: James J. Cai
% Email: jamescai@hku.hk
% Website: http://bioinformatics.org/mbetoolbox/
% Last revision: 5/28/2005

NT = 'ACGT-';
if (nargout>1)
AA = 'ARNDCQEGHILKMFPSTWYV*-';
end

% AANames = {'ala' 'arg' 'asn' 'asp' 'cys' 'gln' 'glu' 'gly' 'his' 'ile' 'leu' 'lys' 'met' ...
%            'phe' 'pro' 'ser' 'thr' 'trp' 'tyr' 'val'};