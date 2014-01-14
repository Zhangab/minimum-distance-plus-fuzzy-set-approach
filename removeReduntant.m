function [cleanlist]= removeReduntant(reduntantlist)

%sort the elements of a cell (reduntantlist) and remove the reduntant elements of cell, and store them in cell cleanlist.
%Note:reduntantlist and cleanlist are cells,not structs.



%sort(reduntantlist)

%sort the reduntantlist and remove the reduntant elements of it, and store them in cell cleanlist.


%diary off       

%diary('log.txt')

Marker=0; k3=1;

cleanlist{1}=reduntantlist{1};




for k1=2:length(reduntantlist)

   Marker=0; 
   for k2=1:k3


     if strcmp(reduntantlist{k1},cleanlist{k2})
   
   
  
       Marker=0;
       else Marker=Marker+1;
      end

    
     if Marker==k3
     cleanlist{k3+1}=reduntantlist{k1};
     k3=k3+1;  
     end
  
  end
 

end

