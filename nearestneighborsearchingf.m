 function [Indexofspe] = nearestneighborsearchingf(References, speciesname)

% nearestneighborsearchingf.m is to find out the nearest neighbor for a
% target species (here we found out the nearest neighbor of all species in the reference database firstly, then select a targeted one.)
% by pointing out the index of nearest neighbor in the database of reference.
% note: nearest neighbor is defined by minimum average genetic distance between two species.

% References - a struct.
% inputs of the program include 
% 

% this matlab script is to take in a modified .fasta file (Ref.fas) as input file,then
% output two sorts of output files.


%%%% #def.
%%%% population size PS specified for scanning as references sequences or not.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%1.Take an modified file (Ref.fas) as input. the name of file could be changed!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% functions list:
% readfastaref();
% removeReduntant();
% twospeciesdist(); you need to choose different evolutionary model in this function!



%clear all;
%addpath C:\MATLAB6p5\work\GeneticDist;
%addpath D:\MATLAB6p5\work\GeneticDist;
%addpath D:\MATLAB6p5\work\GDSI;
% 1. read reference sequences from a set of species.
% #######
%[References] = readfastaref('30ne2Rep6gene0ref.fas',1,0)   % for non-coding DNA/RNA
%[References] = readfastaref('rep44.fas',1,0)   % for non-coding DNA/RNA
%[References] = readfastaref('epaftest.fas',1,0)   % for non-coding DNA/RNA
%[References] = readfastaref('epaftest2.fas',1,0)   % for non-coding DNA/RNA
%[References] = readfastaref('epaftest4.fas',1,0)   % for non-coding DNA/RNA
%[References] = readfastaref('epaftest5.fas',1,0)   % for non-coding DNA/RNA

%[References] = readfastaref('modifiedinput.fasta',1,0)   % for non-coding DNA/RNA


% 2. make species name list by removing reduntant names of species.
%    it calls function - removeReduntant

[cleanspenamelist]= removeReduntant(References.spenames);


Ref = References.seq;
[n2,m2]=size(Ref);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%2. find out nearest and farest species for each species (taken from multispeciesdist.m)
%%%%%%%%% and save them as C,C2, I,and I2 (the last two are the indices of max and min species to the current species checking!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3. start pairwise comparision among different species;

nspecies=length(cleanspenamelist);
% 4. based on the "cleanspenamelist" , start two loops (i & j)

k4=1;

for i=1:nspecies

% save population size of each species into matrxi/vector "popsize" respectively!


  for j=(i+1):nspecies
 
% search for DNA sequence matrix for species i, and for species j, and put into two 
% DNA matrices: DNAmatrix1 and DNAmatrix2, %than call twospeciesdist().% the output was firstly put into
% an rector rowdist, then formated into an matrix spedist.

% initializing DNAmatrix1 and DNAmatrix2!;
%DNAmatrix1=0;
%DNAmatrix2=0;

% 4.1 compare species i to References.spenames{i}


cleantemp1=cleanspenamelist{i};
cleantemp2=cleanspenamelist{j};
k2=1;
k3=1;
       for k=1:n2
    if strcmp(References.spenames{k},cleantemp1)
      DNAmatrix1(k2,:)=References.seq(k,:);
      k2=k2+1;
    end
       end
   % k2
    
       for k=1:n2
    if strcmp(References.spenames{k},cleantemp2)
      DNAmatrix2(k3,:)=References.seq(k,:);
      k3=k3+1;
    end
       end
    
      % k3
 % check if the sizes of   DNAmatrix1 and DNAmatrix2 are same as those of last loop since the program does not initializing these two matrices at the beginning of each loop!!!
 
%[tmp1,tmp2]=size(DNAmatrix1);
%[tmp3,tmp4]=size(DNAmatrix1);
 
 % if (tmp1==k2-1 & tmp3==k3-1)
  
 % [meanD01,SD1,meanD02,SD2,meanD012,SD3] = twospeciesdist(DNAmatrix1,DNAmatrix2);
    
 %   rowdist(k4)=meanD012;
  %  spedist(i,j)=meanD012;
 %   k4=k4+1;

 % else

% substract a submatrix from DNAmatrix1 and DNAmatrix2 

sub1=1:k2-1;
sub2=1:k3-1;
subDNAmatrix1=DNAmatrix1(sub1,:);
subDNAmatrix2=DNAmatrix2(sub2,:);

  [meanD01,SD1,meanD02,SD2,meanD012,SD3] = twospeciesdist(subDNAmatrix1,subDNAmatrix2);
    
    rowdist(k4)=meanD012;
    spedist(i,j)=meanD012;
    k4=k4+1;

    % end
  
    
    
   end

% save population size of each species into matrxi/vector "popsize" respectively!  
   % popsize(i)=k2-1;
    
    %popsize(i)=k3-1;
end

 %clear tmp1,tmp2,tmp3,tmp4; 
 
spedist2=cat(1,spedist,zeros(1,length(spedist)));

fullspedist=spedist2 + spedist2';

[C,I]=max(fullspedist);



temp1=reshape(fullspedist,1,length(fullspedist)*length(fullspedist));
tempmax=max(temp1);


temp2=ones(1,length(fullspedist))*tempmax;

tempfullspedist=fullspedist + diag(temp2);

[C2,I2]=min(tempfullspedist);


clear temp1;
clear temp2;
clear tempfullspedist;

%======speciesname
%cleantemp1=cleanspenamelist{i};
for j=1:nspecies
 if strcmp(cleanspenamelist{j},speciesname)
     Indexofspe= I2(j);
 end  
     end


 end

