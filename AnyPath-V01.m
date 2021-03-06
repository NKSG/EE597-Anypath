%anypath routing algorithm
function [Di] = anypath(prM1,b)

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

X=prM1(:,b);     %fetching the vector of prrM corresponding to node#1
e=1;
while ~isempty(Q)
    %find the min distance node from Q
    if e==1
        p=b;
        e=e+1;
    else
       [j,p]=max(X);%finding the max value of prrM and its index in the 
    end             %vector. This index will be used to identify the
                    %node 'j'.  (do j <- EXTRACT-MIN(Q))
    X(p)=0;         %removing the max value from the vector 
    remove(Q,p);    %removing the max prrM value node from the Q vector
    S=[S p];        %Adding node 'p' in S vector (S <- S U {j})
    
    %find incoming edges (i,j) in E. For this purpose we can set the
    %threshold for deciding the neighboring the nodes to 1e-16
    Y = prM1(:,p);   %fetching the vector of prrM corresponding to node 'p'.
        
    %use only incoming nodes
    if length(S)>1
        for u=1:length(S)
            Y(S(u))=0;
        end
    end
        
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
                        wij=wij*prM1(J(o),E(c));
                    else
                        wij=wij*(1-prM1(J(o),E(c)));
                    end
                 end
                 DJ=DJ+wij*Di(J(u));
                 PiJ=PiJ*(1-prM1(E(c),J(u)));  %calculating PiJ parameter
             end
             PiJ=1-PiJ;
             diJ=1/PiJ;
             DJ=DJ/PiJ;
               
             Di(E(c)) = diJ+DJ;  %Di <- diJ + DJ
             Fi(E(c))=J;         %Fi <- J
        end %end of IF block
     end %end of outer FOR loop
        
end