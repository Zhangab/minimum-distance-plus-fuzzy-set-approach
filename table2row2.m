function [non_diag, diagele] = table2row2(table)

%table2row2 is to transfer the elements of matrix (triu,seperating diag and non-diag elements) into two rows (rector).

[m,n]=size(table);

if (m==n & m>1)

k=1;
for j=1:n
    for i = 1:(j-1)
    
    arow(k) = table(i,j);
    k = k + 1;
    
   % i
   % j
    end
end


non_diag=arow';


diagele=diag(table);

else 
  
    non_diag=table(1,1);


diagele=diag(table);
    
    disp('m!=n')
    
end