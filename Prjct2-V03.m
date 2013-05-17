%Project # 2 for EE597
%authors: Mohit Aggarwal and Prathamesh Khedekar

size=8;                        %this variable plays controls the # of nodes
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
    for b=1:N
        Di=zeros(N,1);
        Fi=containers.Map('KeyType','int32','ValueType','any');
        Q=containers.Map('KeyType','int32','ValueType','int32');
        for i=1:N
            Di(i)=Inf;
            Q(i)=i;
            Fi(i)=[];
        end
   
        Di(b)=0;
        S=[];
        X=prM(:,b);     %fetching the vector of prrM corresponding to node#1
        e=1;
        while ~isempty(Q)
            %find the min distance node from Q
            if e==1
                p=b;
                e=e+1;
            else
                [j,p]=max(X);%finding the max value of prrM and its index in the 
            end              %vector. This index will be used to identify the
                            %node 'j'.  (do j <- EXTRACT-MIN(Q))
            X(p)=0;         %removing the max value from the vector 
            remove(Q,p);    %removing the max prrM value node from the Q vector
            S=[S p];        %Adding node 'p' in S vector (S <- S U {j})
        
            %find incoming edges (i,j) in E. For this purpose we can set the
            %threshold for deciding the neighboring the nodes to 1e-16
            Y = prM(:,p);   %fetching the vector of prrM corresponding to node 'p'.
        
            %use only incoming nodes
            if length(S)>1
                for u=1:length(S)
                    Y(S(u))=0;
                end
            end
        
            E=[];
            for l=1:length(Y)
                if Y(l)>1e-3 %if the prrM of a neighbor is greater than 1e-16 
                            % then include that node in 'E' vector (incoming
                            % edge (i, j) in E)
                E=[E l];
                end
            end
            J=[];
            for c=1:length(E)
                %creating forwarding set. (do J <- Fi U {j})
                J = [Fi(E(c)) p];
            
                if Di(E(c)) > Di(p) %if Di > Dj
                    PiJ=1;
                    DJ=0;
                    %calculating DJ
                    for u=1:length(J)
                        %calculate wij;
                        wij=1;
                        for o=u:-1:1
                            if o==u
                                wij=wij*prM(J(o),E(c));
                            else
                                wij=wij*(1-prM(J(o),E(c)));
                            end
                        end
                        DJ=DJ+wij*Di(J(u));
                        PiJ=PiJ*(1-prM(E(c),J(u)));  %calculating PiJ parameter
                    end
                    PiJ=1-PiJ;
                    diJ=1/PiJ;
                    DJ=DJ/PiJ;
               
                    Di(E(c)) = diJ+DJ;  %Di <- diJ + DJ
                    Fi(E(c))=J;         %Fi <- J
                end %end of IF block
            end %end of outer FOR loop
        
        end

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