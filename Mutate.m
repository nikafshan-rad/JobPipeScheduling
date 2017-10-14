function y=Mutate(x,PT)

    global nResource;
    global PTnTask;
    
    nVar=numel(x);
    mu=rand(1);        
    m=ceil(mu*nVar);
    
    jj=randsample(nVar,m);
    
    y=x; 
    for i=1:numel(jj)
       if x(jj(i))==1 && sum(y)>PTnTask(PT)-nResource
           y(jj(i))=0;
       elseif x(jj(i))==0
           y(jj(i))=1;
       end          
    end
end
