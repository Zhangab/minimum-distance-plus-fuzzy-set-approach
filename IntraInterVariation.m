function [IntraVariation095,InterVariation005] = IntraInterVariation(Ref, alpha_intra, beta_inter)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Ref - struct
% alpha-intra : 0.95 , 95% percentage for intraspecific variation
% beta-inter: 0.05, 5% percentage for interspecific variation

%alpha_intra=0.95;
%beta_inter=0.05;

%[Ref] = readfastaref('test.fas',1,0) 
D00=dn_k2p(Ref);
Vnan=diag(eye(size(D00)))*NaN;
[n,m]=size(D00);
MyOnes=ones(size(D00));
D000=diag(Vnan,0);
D0000=D00+D000;
IACM=MyOnes;%IntraCoeMat
IECM=MyOnes;%InterCoeMat
for i=1:n
    for j=1:m
     if strcmp(Ref.spenames{i},Ref.spenames{j})
   IECM(i,j)=NaN;
  
   else IACM(i,j)=NaN;
      end   
    end
end


IntraVariation=D0000.*IACM;
InterVariation=D0000.*IECM;

[n2,m2]=size(IntraVariation);
[n3,m3]=size(InterVariation);


IntraVariation2=sort(reshape(IntraVariation,1,n2*m2));
IntraVariation3=IntraVariation2(1:(n2*m2-nnz(isnan(IntraVariation2))));


InterVariation2=sort(reshape(InterVariation,1,n3*m3));
InterVariation3=InterVariation2(1:(n2*m2-nnz(isnan(InterVariation2))));

IntraVariation095=norminv(alpha_intra,mean(IntraVariation3),std(IntraVariation3));

InterVariation005=norminv(beta_inter,mean(InterVariation3),std(InterVariation3));

end

