function [FitnessValue, y]=Fitness(x)        
   
   global R;
   global nResource; %% N
   global TaskSize;
   global t_f;
   global t_c;
     
m=x;
j=1;
resource{1}=(1);
n=numel(x);
   for i=1:n
       if m(i)==1
           resource{j}=[resource{j} i+1];
           
       elseif m(i)==0
           j=j+1;
           resource{j}=[];
           resource{j}=[resource{j} i+1];
           
       end
   end
 y=[];
 for h=1:numel(resource)
     for k=1:numel(resource{h})
         y=[y h];
     end
 end

% numel(resource)
W=TaskSize;
CP=R(1,:); %% speed of resources
CB=R(3,:); %% band width between resources and scheduler
TC=zeros(1,nResource);
% EET=zeros(1,nResource);
Ly=numel(y);
for j=1:nResource
   for i=1:Ly
      if y(i)==j
%             
         TC(j)=TC(j)+(W(i)/CP(j))+(t_f(i)+t_c(i)/CB(j));
      end
   end

 end
%    FitnessValue=max(TC);
   FitnessValue=1/max(TC);
end