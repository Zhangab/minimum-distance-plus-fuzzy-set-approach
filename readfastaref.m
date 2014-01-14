function [aln] = readfastaref(filename,seqtype,geneticcode)
%READFASTA - Reads data from a FASTA formatted file into a MATLAB structure
%A file with a FASTA format begins with a right angle bracket (>) and a single 
%line description. adopted from  Cai (2006).
%
% Syntax:  [aln] = readfasta('filename',seqtype,geneticcode)
%
% Inputs:
%    filename    - FASTA formatted file (ASCII text file).
%    seqtype     - Type of sequences. 1 - non-coding nucleotide; 2 - coding
%                  nucleotide; 3 - protein
%    geneticcode - Genetic code. 1 - Standard; 2 - Vertebrate Mitochondrial; 
%                  3 - Yeast Mitochondrial; 4 - Mold, Protozoan, Coelenterate 
%                  Mitochondrial, and Mycoplasma/Spiroplasma; 5 - Invertebrate 
%                  Mitochondrial; 6 - Ciliate, Dasycladacean, and Hexamita 
%                  Nuclear; 0 - non-coding
%
% Outputs:
%    aln    - alignment structure
%
% Examples:
%
% >> [aln] = readfasta('input.fas',1,0)   % for non-coding DNA/RNA
% >> [aln] = readfasta('input.fas',2,1)   % for protein-coding DNA
% >> [aln] = readfasta('input.fas',3,1)   % for protein           
%
% See also: 

% ++++++++++++++++++++++++++++++
% Author: Ai-bing ZHANG
% Email: zhangab2008@gmail.com
% Website: http://202.204.209.200/shizishow.php?idh=493
% Last revision: April/28/2008


MAXNAME = 200;

if nargin < 1
    [filename, pathname] = uigetfile( ...
       {'*.fasta;*.fas', 'FASTA Format Files (*.fasta, *.fas)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a FASTA file');
	if ~(filename), aln=[]; return; end
	filename=[pathname,filename];
end

%if nargin < 2
%	[seqtype, geneticcode]=selectSeqTypeAndGeneticCode;
%	if (isempty(seqtype)|isempty(geneticcode)), aln=[]; return; end
%end


% check input is char
% in a future version we may accept also cells
%if ~ischar(filename)
 %   error('mbetoolbox:InvalidInput','Input must be a character array')
%end

%if ~(exist(filename,'file') || exist(fullfile(cd,filename),'file')),
%  is a valid filename ?
  %  error('mbetoolbox:InvalidInput','Input must be a valid file')
%end



file = fopen(filename, 'r');
disp(['Reading ',filename]);

% Now we are looking for the maximum length of the sequence
n=0;    % the number of sequences
m=0;    % the maximum length
cm = 0; % current sequence length

while 1
	[x,nr] = fscanf(file,'%c',1);
   if nr == 0 break; end;
   if x =='>'  % new sequence started
		if cm > m m=cm; end;
		cm = 0;
		fgets(file);
		n=n+1;
	else
		if isletter(x) | x=='-'
			cm=cm+1;
		end;
	end;
end

if cm > m m=cm; end;

% go throught the file
Ss = char(m); S = [];
str = zeros(1,MAXNAME);
sizes = zeros(1,n);
frewind(file);
% names=[];
spenames={};
names={};
i=0;j=1;i1=0;
while 1
   [x,nr] = fscanf(file,'%c',1);
   if nr == 0 break; end;
	if x =='>'  % new sequence started
		if i~= 0 % save the sequence
			[x, sizes(i)]=size(Ss);
			S=strvcat(S,Ss);
			Ss = []; Ss = char(m);
		end;
		str=fgetl(file); % read the name, we remove the '>' symbol
		% names=strvcat(names,str);
		pos=find(str==' ');
		if ~(isempty(pos))
			str=str(1:pos(1,1));
		end
		i=i+1;
		names{i}=str;
		%disp(str);
		
		%
		
		str=fgetl(file); % read the name, we remove the '>' symbol
		% names=strvcat(names,str);
		pos=find(str==' ');
		if ~(isempty(pos))
			str=str(1:pos(1,1));
		end
		i1=i1+1;
		spenames{i1}=str;
		%disp(str);
		
		
		%
		
		
		
		
		j=1;
	else
		if isletter(x) | x=='-'
		% processing the sequence symbol
	   	Ss(j) = upper(x);
			j=j+1;
		end;
	end;
end

S=strvcat(S,Ss);
[x, sizes(i)]=size(Ss);

aln.seqtype = seqtype;
aln.geneticcode = geneticcode;
aln.seqnames = names;
aln.spenames = spenames;
aln.seq=S;	
aln=encodealn(aln);
fclose(file);