%function [FMF,FMF3] = FMFvalue(tmpdatabase,singleQue)

% Function FMF is based on simufuzzyDNAbarcoding3.m
% 


% simufuzzyDNAbarcoding is to generate fuzzymembershipfunction values for
% each species in the reference database. It's only used as simulation
% study, not for final users.


% updated on 2010-09-22 9:12 
% changes since last version "simufuzzyDNAbarcoding.m"
% use smallest distance of nearest neighbor to the MRCA as interspecific variation but with
% mean distance of potential species to the MRCA
% 


%clear all;
% t0 = clock;   
% addpath C:\MATLAB\R2009a\work2009\TDR3;

% The following functions or script will be called sequecially:

% 1. produce input.fas, neighbor.fas, and outgroup.fas for each
%======= query
% 2.  calculate FMF (FuzzyMembershipFunction) values for each query 
%=======  and store them in Matrix FMF.
% 3. output the species identified by minimum genetic distance method and 
%  their FMF values for each query.
%
%
%
%
%
%
%
%



% 1. find out identified species by minimum genetic distances rules
% the indices of with minimum distance to the query is stored in vector I.





%speidentifyfun3;
%[I,I2ndmin] = speidentifyfunf(tmpdatabase,singleQue)

% 2. 




%for i=1:length(I)

%[Aa] = spegroupdividefun2(Ref,cleanspenamelist{I(i)});

%[tdr] = tdrfun2(Aa,Que.seq(i,:));

%TDR(i)=tdr;

%%%%% ####

%=========================================================================
%=======  1. produce input.fas, neighbor.fas, and outgroup.fas for each
%======= query
%======= 
%=========================================================================

% addpath C:\MATLAB\R2009a\work2009\TDR3;
%[References] = readfastaref('Ref.fas',1,0)   % for non-coding DNA/RNA
%[Referenceseqs] = readfastaref('epaf.fas',1,0)   % for non-coding DNA/RNA
%Ref2 = tmpdatabase.seq;
%[n2,m2]=size(Ref2);
%if (m2<2)
%error('Ref is a Matrix!');
%end

%[cleanspenamelist3]= removeReduntant(tmpdatabase.spenames);
%I=1;%?????????????
% 1. nearest neighbor searching for a specified species in the reference
% database.
%[Indexofspe] = nearestneighborsearchingf(tmpdatabase, cleanspenamelist3{I(i)});
%[I,I2ndmin] = speidentifyfunf(tmpdatabase,singleQue);
[I,I2ndmin] = speidentifyfunf(tmpdatabase,singleQue);

identifiedspe(iglobe) = I;


[cleanspenamelist3]= removeReduntant(tmpdatabase.spenames);

% 1. nearest neighbor searching for a specified species in the reference
% database.

[Indexofspe] = nearestneighborsearchingf(tmpdatabase, cleanspenamelist3{I});
% 2. estimation of parameters of a & b of zmf.


[onespe1] = select1species2(tmpdatabase, cleanspenamelist3{I});
[onespe2] = select1species2(tmpdatabase, cleanspenamelist3(Indexofspe));

%[a,b] = abzmftdr(onespe1.seq, onespe2.seq);% before 2010-10-23 8:55

%[a,b] = abzmfMIMI(onespe1.seq, onespe2.seq); %2010-10-23 8:55
[a,b] = abzmfMIMI(onespe1.seq, onespe2.seq,IntraVariation095,InterVariation005);
%3. estimate the zmf value of a specified query.

[nIndividual,tmp]=size(onespe1.seq);
[nIndividual2,tmp]=size(onespe2.seq);
if ( nIndividual>1 && nIndividual2>1)

%[TDRque] = tdrfun2(onespe1.seq,singleQue.seq);
[meanD011,SD11,meanD021,SD21,meanD0121,SD31,miniD0121,maxD0121] = twospeciesdist3(onespe1.seq,singleQue.seq)
queDist = meanD0121;

else
    [meanD0xx,SDxx,meanD0xx1,SDxx,queDist,SD31xx,miniD0121xx,maxD0121xx] = twospeciesdist3(onespe1.seq,singleQue.seq)
    
  %queDist=0; 
end

FMF(iglobe)=zmf(queDist,[a,b]);


