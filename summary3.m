% this matlab script is to take in an modified .fasta file as input file,then
% output two sorts of output files, one is going to be input files for SAP (refsap.fasta and quesap.fasta),
% another is going to be input files for my own program TDR (reftdr.fasta, quetdr.fasta)



%clear all;

%diary('log.txt');



%t0 = clock;        

%%%% #def.
%%%% population size PS specified for scanning as references sequences or not.
%PS=3;
%Ratioquevsref=1./3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%1.Take an modified file as input. the name of file could be changed!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% functions list:
% readfastaref();
% removeReduntant();
% twospeciesdist(); you need to choose different evolutionary model in this function!



%clear all;
%addpath C:\MATLAB6p5\work\GeneticDist;
%addpath C:\MATLAB\R2009a\work2009\fuzzyset2\leaveOneOutSimuTDR;
%addpath D:\MATLAB6p5\work\GDSI;
% 1. read reference sequences from a set of species.
% #######
%[Referenceseqs] = readfastaref('30ne2Rep6gene0ref.fas',1,0)   % for non-coding DNA/RNA
%[Referenceseqs] = readfastaref('rep44.fas',1,0)   % for non-coding DNA/RNA
%[Referenceseqs] = readfastaref('epaftest.fas',1,0)   % for non-coding DNA/RNA
%[Referenceseqs] = readfastaref('epaftest2.fas',1,0)   % for non-coding DNA/RNA
%[Referenceseqs] = readfastaref('epaftest4.fas',1,0)   % for non-coding DNA/RNA
%[Referenceseqs] = readfastaref('epaftest5.fas',1,0)   % for non-coding DNA/RNA

%[Referenceseqs] = readfastaref('modifiedinput.fasta',1,0)   % for non-coding DNA/RNA

%[Referenceseqs] = readfastaref('test4.fas',1,0)   % for non-coding DNA/RNA
%[Referenceseqs] = readfastaref('epaf.fas',1,0)   % for non-coding DNA/RNA
%[Referenceseqs] = readfastaref('modifiedinput.fas',1,0)   % for non-coding DNA/RNA


%Ref = Referenceseqs.seq;
%[n2,m2]=size(Ref);
%if (m2<2)
%error('Ref is a Matrix!');
%end

% 2. make species name list by removing reduntant names of species.
%    it calls function - removeReduntant

[cleanspenamelist]= removeReduntant(mydatabase.spenames);




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

% 4.1 compare species i to Referenceseqs.spenames{i}


cleantemp1=cleanspenamelist{i};
cleantemp2=cleanspenamelist{j};
k2=1;
k3=1;
       for k=1:n2
    if strcmp(mydatabase.spenames{k},cleantemp1)
      DNAmatrix1(k2,:)=mydatabase.seq(k,:);
      k2=k2+1;
    end
       end
    %k2;
    
       for k=1:n2
    if strcmp(mydatabase.spenames{k},cleantemp2)
      DNAmatrix2(k3,:)=mydatabase.seq(k,:);
      k3=k3+1;
    end
       end
    
       %k3
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



    % end
  
    
    
   end

% save population size of each species into matrxi/vector "popsize" respectively!  
   % popsize(i)=k2-1;
    
    %popsize(i)=k3-1;
end

 %clear tmp1,tmp2,tmp3,tmp4; 
 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% # 3. start to output to SAP & TDR
%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% computate the population size for each species in the database
%% and put the sequences from the same species into a tmp seq matrix of a
%% struct Tmref.

    Tmpref.seqtype=1;
    Tmpref.geneticcode=0;
    
    Complispe.seqtype=1;
    Complispe.geneticcode=0;

for i=1:nspecies   %%%%#####
 k2=1;% k33=1;
    for k=1:n2
 
        if strcmp(mydatabase.spenames{k},cleanspenamelist{i})
 %  store the sequences from the same species into a tmp seq matrix of a struct
         
            k2=k2+1;
        else 
        % put all other sequences (which are not in cleanspenamelist{i} )into another struct "complispe" (comlimentary sequences other than species in the database ) 
        
        % =================
      
        
         % =================
         end
        
        
    end
    
%%%%%%%%% # 3.1 output the sequences of species which population size
%%%%%%%%% popsize<=PS to references for SAP refsap.fasta
    
    
    
  
   
%   Complispe.seq=Complispe.seqtmp(tmpI2,:)
   
   
    
    % output the complimentary sequences
    
    % writefastatdr(Complispe,'refcomplispe.fasta')
    %writefastatdr3(Complispe,'refcomplispe.fasta');%/2009-8-7 14:36 2009-8-7 11:31 
    
    
    
 %%%%%%%%% # 3.2 output the sequences of species which population size
%%%%%%%%% popsize>PS to references for SAP refsap.fasta   
    

%%%%%% # 3.2.1 started to randomly seperate the tmp struct into two parts:
%%%%%% Tmpref2. and Tmpque2.

%%%%%% generate a random vector myrand = randperm (n) in the range of 1 to
%%%%%% n which is equal to the popsize of the species.
%%%%%% 


    
    

    
      popsize(i)=k2-1;
      % k2=1;
         
end             %%%%#####
  







%popsize


%Timeused =etime(clock,t0)./60

%save logsaptdr;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%. # 4. output popsize information.
%%%%%%%%% 
%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%diary ('popsize.txt')
%Num=length(cleanspenamelist);
%for i=1:Num
%disp(cleanspenamelist{i})
%disp (':')
%disp(popsize(i))
%end


%diary off


%[n4,m5]=size(popsize);
%fid=fopen('output.txt','w');
%for k=1:m5
%outname=cell2mat(cleanspenamelist(k));
%end

% output:::::

%for j=1:m5
%outname=cell2mat(cleanspenamelist(j));


%out2=popsize(:,j);
%fprintf(fid,'%6.2f\t',out2);
%fprintf(fid,'%s\n',outname);

%fprintf(fid,'\r');
%end
%fclose(fid);


%rmpath C:\MATLAB\R2009a\work2009\fuzzyset2\leaveOneOutSimuTDR;