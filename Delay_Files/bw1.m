function [Pset, Sset]= bw1(C)
% Simulates the behaviour in Broadcast Window 1. Returns Pset and Sset
% after simulating behaviour of primary user
P=[];
S=[];
M=length(C);
Temp= zeros(1,M);
for i = 1 : M
    Temp(i) = random('bino',1,C(i));%Bern(r)
    if Temp(i) == 1
        P = union(P,i);% Channel Busy.
    else
       S = union(S,i);% Channel Free.
    end
end
Pset= P;
Sset=S;
