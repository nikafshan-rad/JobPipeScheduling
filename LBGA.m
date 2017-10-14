clc;
clear;
close all;

%% problem definition

global PTnTask;
global R;
global TaskSize; 
global MaxSize;
global nTournament;
global nPop;
global nResource;
global t_f;
global t_c;

nPT=3; %% number of paths
PTnTask=[60 70 50];    %% number of Tasks for every path
nResource=20;  %% number of Resources
MinSize=20;   %% minimum size of tasks
MaxSize=100;  %% maximum size of tasks
MinSpeedOfResource=2;     %% min speed of resources
MaxSpeedOfResource=10;    %% max speed of resources
% MinPriceOfResource=1;
% MaxPriceOfResource=5;
MinBandWidth=20;   %% (bps) min of band width between node and scheduler
MaxBandWidth=100;  %% (bps) max of band width between node and scheduler

R1=(MaxSpeedOfResource-MinSpeedOfResource).*rand(1,nResource)+MinSpeedOfResource;   %%speed of cpu
% R1=sort(R1);
% R2=(MaxPriceOfResource-MinPriceOfResource).*rand(1,nResource)+MinPriceOfResource;       %%cost of cpu
% R2=sort(R2);
R3=(MaxSize-MinSize).*rand(1,nResource)+MinSize;               %% initial tasks on resources
R4=(MaxBandWidth-MinBandWidth).*rand(1,nResource)+MinBandWidth;
R=[R1
%    R2
   R3
   R4];

%% GA parameters
MaxIt=40;
nPop=40;
nTournament=2;

pCrossover=0.8;
nCrossover=round(pCrossover*nPop/2)*2;

pMutation=0.1;
nMutation=round(pMutation*nPop);

%% Initialization

individual.Position=[];
individual.Cost=[];

pop=repmat(individual,nPop,1);


for PT=1:nPT
    %% generate size of both tasks and transfering files for path number pp
TaskSize=unifrnd(MinSize,MaxSize,[1 PTnTask(PT)]);  %% size of job workload (billion instructions)
t_f=TaskSize+ 0.01*TaskSize;  %% size(bit) of file for task i that is transfered to the node from the scheduler
t_c=TaskSize+ 0.02*TaskSize;  %% size(bit) of file for task i that is transfered to the scheduler from the node

 %% generate initial population for every path
    z=ceil((nResource-1).*rand(1))+1;
    jj=randsample(nResource,z);
    
for i=1:nPop
    pop(i).Position=ones(1,PTnTask(PT)-1);
    z=ceil((nResource-1).*rand(1));
    jj=randsample(PTnTask(PT)-1,z);
    pop(i).Position(jj)=0;
    pop(i).Cost=Fitness(pop(i).Position);
end

BestCost=zeros(1,MaxIt);
meanCost=zeros(1,MaxIt);

 %% GA main loop
for it=1:MaxIt
    
 %%crossover
   pop2=repmat(individual,nCrossover/2,2);
for k=1:nCrossover/2    
    
    L1=PTnTask(PT)-nResource-1;
    L2=L1;
    while L1<PTnTask(PT)-nResource || L2<PTnTask(PT)-nResource 
      i1=Tournament(pop);%ceil((nPop-1).*rand(1))+1;
      i2=Tournament(pop);%ceil((nPop-1).*rand(1))+1;
       
       p1=pop(i1).Position;
       p2=pop(i2).Position;
      
       [a, b]=Crossover(p1,p2);
       L1=sum(a);
       L2=sum(b);
    end
    
       ch1.Position=a ;
       ch2.Position=b; 
       ch1.Cost=Fitness(a);
       ch2.Cost=Fitness(b);
      
       pop2(k,1)=ch1;
       pop2(k,2)=ch2;       
end
       pop2=pop2(:);
       


       %%mutation 
%        mu=GeneVariance(pop);
%        r2=Variance(pop)/2;
%        pMutation=100*(1/(r2+1));
%        nMutation=ceil(pMutation*numel(pop));
              
%        G=GeneVariance(pop);
%        r1=Variance(pop);
%        r2=atan(r1);
%        fis=readfis('MuFuzzy');
%        mu=evalfis([G r2],fis);
       
%        fis=readfis('pMutateFuzzy');
%        pMutation=evalfis(r2,fis);
%        pMutation=mu;
%        nMutation=ceil(pMutation*numel(pop));
       
       pop3=repmat(individual,nMutation,1);
    for k=1:nMutation
        
        i3=Tournament(pop);%ceil((nPop-1)*rand(1)+1);
        q=Mutate(pop(i3).Position,PT);
        p.Position=q;
        p.Cost=Fitness(q);
              
        pop3(k)=p;      
    end
     
     %%merg population
      pop1=[pop 
            pop2  
            pop3];
    
     %%sort
    Costs=[pop1.Cost];
    [Costs , Sort_Order]=sort(Costs,'descend');   %% ascending order
    pop1=pop1(Sort_Order);
    pop1=pop1(1:nPop);                 %%delete extra individuals

%     
%     BestMakespan(it)=min(u1);
%     MeanMakespan(it)=mean(u1);
      pop=pop1; 
      
    BestCost(it)=pop(1).Cost;
    meanCost(it)=mean([pop.Cost]);

end
%% Display Results
[g1,g2]=Fitness(pop(1).Position);
disp(['path ' num2str(PT)]);
disp(['Best chromosome for path ' num2str(PT),' =' num2str(pop(1).Position)]);
disp(['Best shedule for path ' num2str(PT),' =' num2str(g2)]);
disp(['Best Fitness for path ' num2str(PT),' =' num2str(pop(1).Cost)]);

figure(PT);
plot(BestCost,'g-');
hold on
plot(meanCost,'b-');
legend('Best Fitness','Mean Fitness');
title(['Fitness of path ' num2str(PT), ' with ' num2str(PTnTask(PT)),' Tasks']);
xlabel('Iteration');
ylabel('Fitness');

%% choose 4 best solutions for every path in pipeline
bestsol{PT}={pop(1) pop(2) pop(3) pop(4)};  
end
bestsol
Greedy={};
for t=1:nPT 
    for k=1:4
        Greedy{k,t}=[bestsol{t}{k}.Position];
        
    end
end
