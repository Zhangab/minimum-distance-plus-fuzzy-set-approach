



%leaveOneOutSimu.m is to read a modified fasta file into a struct named
%"mydatabase",and 
% select random query to test " mini-distance + fuzzy set theory" method
% with a leave-one-out method for nRepp times.


%2010-10-23 8:55  note: 
% This method will use maximum intraspecific genetic distance of potential species 
% of a query, and the minimum interspecific genetic distance of the
% potential species and its nearest neighbor (NN) as the [a,b] of the curve
% zmf.


%%% Important Note1<:====================before 2010-10-23 8:55
% if you want to generate an output file named 'randomsingletons_first.txt' for
% the first time, please activate all lines marked with

% "%%to-generate-randomsingletons.txt##":
%%to-generate-randomsingletons.txt## fid2=fopen('randomsingletons_first.txt','w');
%%to-generate-randomsingletons.txt## fprintf(fid2,'%6.2f\t',r1(iglobe));
%%to-generate-randomsingletons.txt## fclose(fid2);
%%% Important Note1>:====================


%%% Important Note2<:====================
% if one just wants to generate a 'randomsingletons_first.txt' for later
% use or for other program to use, one'd better run this script without
% calling 'FMFvalueTDR.m'arround line 196;
% just use another scirpt named 'leaveOneOutSimuTDR2-randlist.m'
%%% Important Note2>:====================
%diary on;

%mydelete;
clear all;
%tStart=tic;
%fid5=fopen('log.txt','w');
%diary('log.txt');
%diary on;
% for test;
%global tmp_alpha;
%global tmp_j22;


%addpath C:\MATLAB\R2009a\work2009\fuzzyset2;
%addpath C:\MATLAB\R2009a\work2009\TDR3;
nRepp = 2; % the number of replications for simulation.
nRepp0 = nRepp;


% ### 1. read a modified fasta format file and save it to a construct called
% "mydatabase"

[mydatabase] = readfastaref('modifiedinput.fas',1,0)   % for non-coding DNA/RNA
% important notes: all sequences should contains only letters and / or
% numbers, nothing else! others, it will break.

%[mydatabase,nseq] = readfasta8('modifiedinput.fas',1,0)


%n=nseq; % the number of reference sequences
%n2 = 15; % the number of query sequenes

%[header,s]=sequence(n);
%d=[];


%for i=1:n         
 %   d1=dvcurve(s{i});
  %  d=[d;d1]; % reference data matrix.
% end



[IntraVariation095,InterVariation005] = IntraInterVariation(mydatabase,0.95,0.05);



Ref = mydatabase.seq;
[n2,m2]=size(Ref);
if (m2<2)
error('Ref is a Matrix!');
end

[nseq,nsite]=size(mydatabase.seq);

% ### 2. to find out singletons by calling "summary3.m", and save the
% indices
% of these singletons into a vector "indexSingletons", and those of nonsingletons
% into "indexNonsingletons"

summary3;

% ### 2.1 to find out the names of species for singletons, and save the
% species names into a cell named "listSingletons"

k=1;
for i=1:length(popsize)
    
    if (popsize(i)==1)
        
            listSingletons{k} = cleanspenamelist{i};
        k=k+1;
       
    else
        listSingletons = cell(0);
        
    end
    
end

% ### 2.1 compare the species names of sequences in the database with the
% "listSingletons", and to find out the indices of singletons in the whole
% database, and save these indices into a vector named "indexSingletons"



mydatabase.statusSingletons(length(mydatabase.seqnames))=0;
indexAll=[1:length(mydatabase.seqnames)];

if (isempty(listSingletons))
    %length_listSingletons=0;
    indexNonsingletons=indexAll;
else
    

length_listSingletons=length(listSingletons);

for j = 1:length_listSingletons
for    i = 1:length(mydatabase.seqnames)
    if strcmp(mydatabase.spenames{i},listSingletons{j})
        
        %if strcmp(mydatabase.spenames{k},cleanspenamelist{i})
        
    mydatabase.statusSingletons(i)=1;
    end
end
end


% ### 2.2 save into  "indexSingletons"
k = 1;

for i = 1:1:length(mydatabase.seqnames)
    
    if (mydatabase.statusSingletons(i)==1)
       indexSingletons(k) = i;
        k = k + 1;
    end
    
end

% save into "indexNonsingletons"
lengthIS=length(indexSingletons);

lengthIN=length(indexAll)-lengthIS;

for i=1:lengthIS
indexAll(indexSingletons(i))=0;
end

tmpAll=sort(indexAll,2,'descend');

indexNonsingletons=sort(tmpAll([1:lengthIN]));


r1 = randi(length(indexSingletons),nRepp,1);% for singletons


end  %endofIF-ELSE loop

% ### 3. to generate random replications for singletons (r1) and non-singletons
%(r2)

r2 = randi(length(indexNonsingletons),nRepp,1);% for non-singletons
r3 = randi(nseq,nRepp,1);% for both singletons and non-singletons



%  to check if file 'randomsingletons.txt'exit, if not, use r1, otherwise
%  use r1_fromfile to replace r1;
%reading data <:======================
fid4=fopen('randomNumberSpecified.txt','rt');%randomsingletons.txt
if fid4~=-1
 %   disp(msg);
  r3_fromfile=fscanf(fid4,'%f\t');
  r3=r3_fromfile;
  if max(r3_fromfile>nseq)
      disp('errors: Your input from randomsingletons.txt does not match that of modifiedinput.fas! ')
      return;
  end
end


% !!! to save computation time, only calculate FMF value for each unique query.

r3old=r3;
r3=unique(r3);
nRepp = length(r3);



% read data
%fclose(fid4);

% test reading data >:======================

%output r1 <:======================
%%to-generate-randomsingletons.txt## 
%fid22=fopen('randomNumber.txt','w');
%fid33=fopen('randomlist0.txt','w');

%fid2=fopen('uniquerandomNumber.txt','w');
%fid3=fopen('uniquerandomlist0.txt','w');


%fprintf(fid3,'%i\r\n',nRepp);
% do some cleanning
clc
myclear;
%tmpdatabase = mydatabase;
% ### 4. calcualate the FMF values 
% ### 4.1 for each singleton randomly selected (nRepp times)

%nRepp = 1;
% ### 4.2 output the title
fid=fopen('output.txt','w');
 
fprintf(fid,'%s','There were ');  

fprintf(fid,'%i',nRepp0);  

fprintf(fid,'%s',' replciations initially, but only ');  

fprintf(fid,'%i',nRepp);  

fprintf(fid,'%s',' non-redundant queries were performed to save computation time.');  
 fprintf(fid,'\r\n');
fprintf(fid,'%s',' The results are listed below: ');  


 fprintf(fid,'\r\n');
 
 
 fprintf(fid,'%s\t','No.');  
 fprintf(fid,'%s\t','Query');  
 fprintf(fid,'%s\t','Species identified by minimum genetic distance');  
 fprintf(fid,'%s\t','Fuzzymembershipvalue_MIMI');  
 %fprintf(fid,'%s\t','Fuzzymembershipvalue2');  
 %fprintf(fid,'%s\t','Fuzzymembershipvalueold');  
 fprintf(fid,'\r\n');


for iglobe = 1:nRepp
    tmpdatabase = mydatabase;
    tmpdatabase.seqnames(r3(iglobe))=[];
    tmpdatabase.spenames(r3(iglobe))=[];
    tmpdatabase.seq(r3(iglobe),:)=[];
    tmpdatabase.statusSingletons(r3(iglobe))=[];
%simufuzzyDNAbarcoding3


singleQue.seqtype = 1;
singleQue.geneticcode = 0;
singleQue.seqnames = mydatabase.spenames(r3(iglobe));
singleQue.seq = mydatabase.seq(r3(iglobe),:);

%simufuzzyDNAbarcoding3
%[indexofspe2] = speidentifyfunf(tmpdatabase,singleQue)

%[I,I2ndmin] = speidentifyfunf(tmpdatabase,singleQue)
FMFvalueMIMI; % Maximum Intra and Minimum Inter.

%end 

% ### 5. output

%fid=fopen('output.txt','a');
% fprintf(fid,'%s\t','No.');  
% fprintf(fid,'%s\t','Query');  
% fprintf(fid,'%s\t','Species identified by minimum genetic distance');  
% fprintf(fid,'%s\t','Fuzzymembershipvalue_mini_b&mean_a');  
% fprintf(fid,'%s\t','Fuzzymembershipvalue3_mini_b&max_a');  
% fprintf(fid,'\r');
 
 
 
 
%for i=1:nRepp
    
 %(0) output the number of query   
fprintf(fid,'%d\t',iglobe);  
 
%(1) output the names of query   
 fprintf(fid,'%s',mydatabase.seqnames{r3(iglobe)});
 if (mydatabase.statusSingletons(r3(iglobe))==1)
     { fprintf(fid,'%s\t','*');}
 else{ fprintf(fid,'\t');}
 end
 fprintf(fid,'%s\t',mydatabase.spenames{r3(iglobe)});
%(2) output the names of species indentified
fprintf(fid,'%s\t',cleanspenamelist3{identifiedspe(iglobe)});  
%fprintf(fid,'%s\t',cleanspenamelist3{I}');  
%(3)  output the FMF values  

fprintf(fid,'%6.2f\t',FMF(iglobe));

%fprintf(fid,'%6.2f\t',FMF3(iglobe));

%fprintf(fid,'%6.2f\t',FMFold(iglobe));

% ### 6. 

% ### 7. 

% ### 8. 
fprintf(fid,'\r\n');


%output r1 <:======================
%fid2=fopen('randomsingletons.txt','w');
%%to-generate-randomsingletons.txt## 
%fprintf(fid2,'%6.2f\t',r3(iglobe));

%fprintf(fid3,'%s\r\n',mydatabase.seqnames{r3(iglobe)});
%fprintf(fid3,'\r');
%fprintf(fid2,'\r');
%output r1 >:======================
end


% iglobe2-loop is going to output of list of sequences randomly selected!



%for iglobe2 = 1:length(r3old)
  %fprintf(fid22,'%d\t',iglobe2);   
  %fprintf(fid22,'%s\t',mydatabase.seqnames{r3old(iglobe2)});
 % fprintf(fid22,'%6.2f\t',r3old(iglobe2));
  %fprintf(fid33,'%d\t',iglobe2); 
  %fprintf(fid33,'%s\r\n',mydatabase.seqnames{r3old(iglobe2)});
%end





%tElapsed=toc(tStart)./60;
%fprintf(fid5,'%6.2f\t',tElapsed);
%fprintf(fid5,'%s','mins');
fclose(fid);
%%to-generate-randomsingletons.txt## 
%fclose(fid2);
%fclose(fid3);
%fclose(fid5);
%fclose(fid22);
%fclose(fid33);

fid4=fopen('modifiedinput.fas','rt');

fclose(fid4);

%diary off;

% test %reading data <:======================
%fid4=fopen('randomsingletons.txt','rt');
%if fid4==-1
 %   disp(msg);
 % return;
%end
% read data
%r1_fromfile=fscanf(fid4,'%f\t');
%fclose(fid4);

% test reading data >:======================



%rmpath C:\MATLAB\R2009a\work2009\fuzzyset2;
%rmpath C:\MATLAB\R2009a\work2009\TDR3;
%save output2010-10-05
%SendEMail('programrunningstatus@gmail.com', 'onetwothree123', 'smtp.gmail.com', 'zhangab2008@gmail.com', 'RunningstatusofPC01', 'The running on PC01 is over');
%important note:
%web('file:///C:\MATLAB\R2009a\work2009\fuzzyset2\leaveOneOutSimuTDR\leaveOneOutSimuTDR2_help.html', '-browser');

 
