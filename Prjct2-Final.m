%Project # 2 for EE597
%authors: Mohit Aggarwal and Prathamesh Khedekar

size=10;                        %this variable plays controls the # of nodes
avgETX=zeros(size,1);           %variable to hold Avg ETX in case of Dijkstra
anyPathAvgETX=zeros(size,1);    %variable to hold Avg ETX for any path routing
EI=zeros(size,1);
anypathEI=zeros(size,1);
B=zeros(size,1);

for K=1:size
    N=K*K;   %# of wireless nodes for grid topology only
    try
    [topo , prM] = LinkLayerModel(1,N); % '1' == Grid Topology
    catch err
        disp(err.identifier);   % if link layer model throws an error the 
                                %error message will be displayed here
        break;
    end
    
    etx=zeros((N*N):1);         % preparing ETX vector
    x=zeros((N*N):1);
    y=zeros((N*N):1);
    k=1;
   
    for i=1:N
        for j=1:N
            etx(k)=1/prM(i,j);  % calculating ETX vector
            x(k)=i;
            y(k)=j;
            k=k+1;
        end
    end

    % creating sparse matrix
    DG=sparse(x, y, etx);
    cost=zeros(N,1);
    
    %calculating shortest path using dijkstra algotrihtm
    %iterating so that each node in the topology will get a chance to play
    %source node
    for b=1:N
        [dist,path,pred] = graphshortestpath(DG,b);
        cost(b)=mean(dist);
    end
    %calculating 95% confidence interval
    m=mean(cost);
    [H P C stat]=ttest(cost,m,0.05);
    EI(K)=stat.sd;
    B(K)=N;
    avgETX(K)=mean(cost);  %finally calculating avgETX value for given
                                %#of nodes
    
                                
                                
                                
    %----------------------------------------------------------------%
    %calculating shortest anypath route using Laufer et al algorithm
    %we are assuming node#1 to the destination node, and initializing the
    %algorithm
    
    %iterate over the number of nodes present in the topology such that
    %each node will get chance to act as destination
    anypathCost=zeros(N,1);
    
    %before calculating the anypath routes replace any value smaller than
    %1e-5 to '0' in the prM matrix
    for a=1:N
        for b=1:N
            if prM(a,b) < 1e-8
                prM(a,b)=0;
            end
            if a==b
                prM(a,b)=-1;
            end
        end
    end
    
    for b=1:N
        [Di Fi]=AnyPath(prM,b);
        %calculate the avg any path ETX for given number of nodes
        Di(Di==Inf)=0;
        anypathCost(b)=mean(Di);
    end
    
    %calculating 95% confidence interval
    m=mean(anypathCost);
    [H P C stat]=ttest(anypathCost,m,0.05);
    anypathEI(K)=stat.sd;
    anyPathAvgETX(K)=mean(anypathCost);     %finally calculating anyPathAvgETX 
                                            %value for given #of nodes
end
    
hold on;
p=errorbar(B,avgETX,EI);
q=errorbar(B,anyPathAvgETX,anypathEI);
grid on;
set(p,'Color','Blue','LineWidth',2);
set(q,'Color','Red','LineWidth',2);
title('Avg ETX vs # of nodes for Grid Topology');
xlabel('# of nodes');
ylabel('Avg ETX');
legend('Djikstra', 'Any Path');