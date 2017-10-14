function [ch1, ch2]=Crossover(p1,p2)

  x1=p1;
  x2=p2;
   
  nVar=numel(x1);
  c=ceil((nVar-1).*rand(1));
  
  y1=[x1(1:c) x2(c+1:end)];
  y2=[x2(1:c) x1(c+1:end)];
  
  ch1=y1;
  ch2=y2;

end