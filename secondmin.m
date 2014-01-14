function [Y3,I2ndmin] = secondmin(A)
% secondmin() is to generate the second minimum values of each column
% vector and their indices.
%   Y3 - the second minmum values for each column vector.
%   I2ndmin - the indices of each second minmum value for each column
%   vector


%  The function will call the following functions
%  min
%  [] - delete an element of vector 


[m,n] = size(A);

for i=1:n
    oneColumn = A(:,i);
    TmpColumn = oneColumn;
    [Y1, I1] = min(oneColumn);
    
    TmpColumn(I1) = []; % delete the element with minmum value
    [Y2, I2] = min(TmpColumn);
    Y3(i) = Y2;
    if(I2>I1)
    I2ndmin(i) = I1;
    else 
         I2ndmin(i) = I2 + 1;
    end
    
end














end

