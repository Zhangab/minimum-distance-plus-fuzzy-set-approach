 function [onespe2] = select1species2(References, speciesname)
% function select1species() is to output a fasta file by taking the index
% of a species name in the cleanspenamelist (I2(1))
%  References - a struct
%  I2 - the index of a species name in the cleanspenamelist 
% output a fasta file.

%  make species name list by removing reduntant names of species.
%    it calls function - removeReduntant





[cleanspenamelist]= removeReduntant(References.spenames);


%===============================


 onespe.seqtype = 1;
    onespe.geneticcode = 0;



[m,n] = size(References.seq);
 % cleantemp1=cleanspenamelist{I2};
 %cleantemp1 = 'Chrotopterus auritus';
  cleantemp1 = speciesname;
 
  %cleantemp1=cleanspenamelist{2};
  
  
  k2=1;
for i=1:m
  
     if strcmp(References.spenames{i},cleantemp1)
      onespe.seq(k2,:)=References.seq(i,:);
      onespe.seqnames(k2)=References.seqnames(i);
      k2=k2+1;
    end
   
    %k2
    
    
        
end

%myoutput='onespe2.fas';
% writefasta(onespe,myoutput);  % 2009-10-11 10:54
 onespe2 = onespe;
 
 
 % writefasta(onespe,'onespe2.fas');  % 2009-10-11 10:54

end

