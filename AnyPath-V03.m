%anypath routing algorithm
function [Di] = AnyPath(prM1,b)

N=length(prM1);
Di=zeros(N,1);
Fi=containers.Map('KeyType','int32','ValueType','any');
Q=containers.Map('KeyType','int32','ValueType','int32');
    
%initialize the variables
for i=1:N
    Di(i)=Inf;
    Q(i)=i;
    Fi(i)=[];
end
   
Di(b)=0;
S=[];
while ~isempty(Q)
    %From Q find the node having min distance to destination node 'b'
    %finding the max value of prrM and its index in the vector. This index 
    %will be used to identify the node 'j'.  (do j <- EXTRACT-MIN(Q))
    min=Inf;
    p=1;
    k12=keys(Q);
    k1=cell2mat(k12);
    for i=1:length(k1)
        if Di(Q(k1(i)))<min
            min=Di(Q(k1(i)));
            p=k1(i);
        end
    end
    remove(Q,p);    %removing the max prrM value node from the Q vector
    S=[S p];        %Adding node 'p' in S vector (S <- S U {j})
    
    %find incoming edges (i,j) in E. For this purpose we can set the
    %threshold for deciding the neighboring the nodes to 1e-5
    Y = prM1(p,:);   %fetching the vector of prrM corresponding to node 'p'.
        
    %use only incoming nodes
    if length(S)>1
        for u=1:length(S)
            Y(S(u))=0;
        end
    end
        
    E=[];
    for l=1:length(Y)
        if Y(l)>0 %if the prrM of a neighbor is greater than 1e-5
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
                        wij=wij*prM1(J(o),E(c));
                    else
                        wij=wij*(1-prM1(J(o),E(c)));
                    end
                 end
                 DJ=DJ+wij*Di(J(u));
                 PiJ=PiJ*(1-prM1(J(u),E(c)));  %calculating PiJ parameter
             end
             PiJ=1-PiJ;
             diJ=1/PiJ;
             DJ=DJ/PiJ;
               
             Di(E(c)) = diJ+DJ;  %Di <- diJ + DJ
             Fi(E(c))=J;         %Fi <- J
        end %end of IF block
     end %end of outer FOR loop
    
end
