function [indexofspe2,I2ndmin] = speidentifyfunf(Ref,Que)

% speidentifyfun is to generate a list of species name identified by a sort
% of genetic distance specified by string 'distfun'

% Ref - a struct of reference sequences with both DNA matrix and species
% information, such as species names;
% Que - a struct of query sequences with both DNA matrix and species
% information, such as species names;

% speid - indexofspe 
%I-index of species with minimus distance to the query.

% it will call the following functions during runing


% removeReduntant() 
% qrefdistfun()


%========================================================================
%clear all;

%addpath C:\MATLAB\R2009a\work2009\TDR3;
%temp;
% 1. read reference sequences from a set of species.
%[Ref] = readfastaref('ref.fas',1,0)   % for non-coding DNA/RNA

%Ref = Referenceseqs.seq; % keep this!


%[Que] = readfasta('que.fas',1,0)   % for non-coding DNA/RNA

%que = Queryseqs.seq(1,:); % keep this!
% 
%========================================================================

%distfun='minid';

% #1  making a list of species names without reduntance in the reference database 


[cleanspenamelist2]=removeReduntant(Ref.spenames);

nspecies = length(cleanspenamelist2); % the number of species in the references sequences database.

% #2 
[m,n]=size(Ref.seq);
[m2,n2]=size(Que.seq);

for j=1:nspecies
    cleantemp1=cleanspenamelist2{j};
         k2=1;
    for i=1:m  
    if strcmp(Ref.spenames{i},cleantemp1)
      DNAmatrix1(k2,:)=Ref.seq(i,:);
      k2=k2+1;
  end
    end
    % substract a submatrix from DNAmatrix1

sub1=1:k2-1;
subDNAmatrix1=DNAmatrix1(sub1,:);
    
    for k=1:m2
        
        [minid,meanminid,maxd] = qrefdistfun(subDNAmatrix1,Que.seq(k,:));
        
    out.dist(j,k)=minid;
    out2.dist(j,k)=meanminid;
    out3.dist(j,k)=maxd;
    
end
end

[Y,I] =min(out.dist,[],1) ;
indexofspe2=I;

%[Ymean,Imean] =min(out2.dist,[],1) ;
%[Ymax,Imax] =max(out3.dist,[],1) ;

[Y3,I2ndmin] = secondmin(out.dist);

% speid=cleanspenamelist2{I}

%rmpath C:\MATLAB\R2009a\work2009\TDR3;
