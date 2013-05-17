%Project # 2 for EE597
%authors: Mohit Aggarwal and Prathamesh Khedekar

size=6;                        %this variable plays controls the # of nodes
avgETX=zeros(size,1);           %variable to hold Avg ETX in case of Dijkstra
anyPathAvgETX=zeros(size,1);    %variable to hold Avg ETX for any path routing
b=zeros(size,1);

for K=2:size
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
    %h = view(biograph(DG,[],'ShowWeights','on'));
    %calculating shortest path using dijkstra algotrihtm
    [dist,path,pred] = graphshortestpath(DG,1);

    count=length(dist);
    for i=1:count
            avgETX(K) = avgETX(K)+dist(i);
            b(K)=N;
    end
    avgETX(K)=avgETX(K)/count;  %finally calculating avgETX value for given
                                %#of nodes
    
    
    %calculating shortest anypath route using Laufer et al algorithm
    %we are assuming node#1 to the destination node, and initializing the
    %algorithm
    Di=zeros(N,1);
    Fi=containers.Map('KeyType','int32','ValueType','any');
    Q=containers.Map('KeyType','int32','ValueType','int32');
    for i=1:N
        Di(i)=Inf;
        Q(i)=i;
        Fi(i)=[];
    end
    
    Di(1)=0;
    S=[];
    X=prM(:,1);     %fetching the vector of prrM corresponding to node#1
    while ~isempty(Q)
        %find the min distance node from Q
        [j,p]=max(X);   %finding the max value of prrM and its index in the 
                        %vector. This index will be used to identify the
                        %node 'j'.  (do j <- EXTRACT-MIN(Q))
        X(p)=0;        %removing the max value from the vector 
        remove(Q,p);   %removing the max prrM value node from the Q vector
        S=[S p];        %Adding node 'p' in S vector (S <- S U {j})
        
        %find incoming edges (i,j) in E. For this purpose we can set the
        %threshold for deciding the neighboring the nodes to 1e-16
        Y = prM(:,p);   %fetching the vector of prrM corresponding to node 'p'.
        E=[];
        for l=1:length(Y)
           if Y(l)>1e-5 %if the prrM of a neighbor is greater than 1e-16 
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
                            wij=wij*prM(E(c),J(o));
                        else
                            wij=wij*(1-prM(E(c),J(o)));
                        end
                    end
                    DJ=DJ+wij*Di(J(u));
                    PiJ=PiJ*(1-prM(E(c),J(u)));  %calculating PiJ parameter
                end
                PiJ=1-PiJ;
                diJ=1/PiJ;
                DJ=DJ/PiJ;
                
                Di(E(c)) = diJ+DJ;  %Di <- diJ + DJ
                Fi(E(c))=[Fi(E(c)) J];         %Fi <- J
            end %end of IF block
        end %end of outer FOR loop
        
    end

    %calculate the avg any path ETX for given number of nodes
    for i=1:length(Di)
        anyPathAvgETX(K)=anyPathAvgETX(K)+Di(i);
        b(K)=N;
    end
    anyPathAvgETX(K)=anyPathAvgETX(K)/length(Di);
end
    
hold on;
p=plot(b,avgETX);
q=plot(b,anyPathAvgETX);
grid on;
set(p,'Color','Blue','LineWidth',2);
set(q,'Color','Red','LineWidth',2);
title('Avg ETX vs # of nodes for Grid Topology');
xlabel('# of nodes');
ylabel('Avg ETX');
legend('Djikstra', 'Any Path');