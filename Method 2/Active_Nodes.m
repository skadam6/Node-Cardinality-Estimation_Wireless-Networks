%Active Nodes generation
%function [A,temp] = Active_Nodes(N,q,temp)
%Modified: Set of hash values for each type of user.
function temp = Active_Nodes(N,q,temp)
%A=[];
for i=1:N
    if (random('bino',1,q)) %Bern(q)
        %A = union(A,i);
        value = Hash(i,N)+1;
        temp(value) = temp(value)+1;
    end
end