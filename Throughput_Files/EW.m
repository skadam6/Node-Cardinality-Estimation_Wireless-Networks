function [B,a,b,c,d,e,f,g] =EW(Ea,Na,Pa,Hmax,Ne,Np,Nn)
% Simulates the behaviour of Estimation window
n1=length(Ea);
n2=length(Na);
n3=length(Pa);
temp=zeros(Hmax,3);
for i= 1:n1         %B((i),1)= Emergency Nodes
    value=Hash(Ea(i),Ne);
    temp((value+1),1) = temp((value+1),1)+1;
end
for i= 1:n2          %B((i),2)= Normal Nodes
    value=Hash(Na(i),Np);
    temp((value+1),2) = temp((value+1),2)+1;
end
for i= 1:n3          %B((i),3)= Periodic Nodes
    value=Hash(Pa(i),Nn);
    temp((value+1),3) = temp((value+1),3)+1;
end

count=1;
e_hat=0;
p_hat=0;
n_hat=0;
E2=zeros(1,4);
E3=zeros(1,4);

%You have the required set B saved in Temp set
for i=1:Hmax
    if temp(i,1)==1
        if (temp(i,2)==0 && temp(i,3)==0)  %alpha,alpha
            e_hat=e_hat+1;
        elseif (temp(i,2)>=1 && temp(i,3)==0) %alpha,collision
            e_hat=e_hat+1;
            n_hat=n_hat+1;
        elseif (temp(i,3)>=1 && temp(i,2)==0) %collision,alpha
            e_hat=e_hat+1;
            p_hat=p_hat+1;
        else                           %collision, collison, Phase 2 req      
            E2(count,1)=i;
            E2(count,2:4)=temp(i,1:3);
            count=count+1;   
        end
    elseif temp(i,1)==0
        if (temp(i,2)==1 && temp(i,3)==0)   %Empty,beta
            n_hat=n_hat+1;
        elseif (temp(i,2)>=1 && temp(i,3)==0) %Empty,collision
            %e_hat=e_hat+1;
            n_hat=n_hat+1;
        elseif (temp(i,2)==1 && temp(i,3)==1) %beta,beta
            p_hat=p_hat+1;
            n_hat=n_hat+1;
        elseif (temp(i,2)==1 && temp(i,3)>=1) %collision,beta
            p_hat=p_hat+1;
            n_hat=n_hat+1;
        elseif (temp(i,3)==1 && temp(i,2)==0) %beta,Empty
            p_hat=p_hat+1;
            %n_hat=n_hat+1;
        elseif (temp(i,3)>=1 && temp(i,2)==0) %collision,Empty
            p_hat=p_hat+1;
            %n_hat=n_hat+1;
        elseif (temp(i,2)>=1 && temp(i,3)==1) %beta,collision
            p_hat=p_hat+1;
            n_hat=n_hat+1;
        elseif (temp(i,2)==0 && temp(i,3)==0) %Empty, Empty
            
        else                  %collision, collison, Phase 2 req
            E2(count,1)=i;
            E2(count,2:4)=temp(i,1:3);
            count=count+1;   
        end
    else                      %collision, collison, Phase 2 req
        E2(count,1)=i;
        E2(count,2:4)=temp(i,1:3);
        count=count+1;   
    end
end

%-------------------------EW Phase 2------------------------------------%
count2=count-1;
count=1;
for i=1:count2
    if E2(i,2)==1               %1E, 1P, 1N; CC due to 1E 1P  1N
         e_hat=e_hat+1;
         p_hat=p_hat+1;
         n_hat=n_hat+1;
    elseif E2(i,2)==0               %2P, 2N; CC due to 0E 2P 2N
         p_hat=p_hat+1;
         n_hat=n_hat+1;
    else                      %Phase 3 req
        E3(count,1)=i;
        E3(count,2:4)=temp(i,1:3);
        count=count+1;   
    end
end      

%-------------------------EW Phase 3------------------------------------%
%count3=count-1;
for i=1:(count-1)
    if E3(i,4)==0               %0P
        if E3(i,3)==0            %1E, 0N; CC due to 2E 0P 0N
            e_hat=e_hat+1;
        else                    %1E, 1N; CC due to 2e 0P 2N
            e_hat=e_hat+1;
            n_hat=n_hat+1;
        end
    else                       %1P
        if E3(i,3)==0          %1E, 1P, 0N; CC due to 2E 1P 0N
             e_hat=e_hat+1;
             p_hat=p_hat+1;
        else                   %1E, 1P, 1N; CC due to 2E 1P 1N
            e_hat=e_hat+1;
            p_hat=p_hat+1;
            n_hat=n_hat+1;
        end
    end
end

%-------------------Actual calculation of estimates------------------
flag=0;
for j=1:Hmax         %Finds rightmost zero in bitmap for emergency nodes
    if(temp(j,1)==0)
        flag=1;
        break;
    end
end

if(flag==1)
    re=j-1;
else
    re=j;
end
e_hat=floor(1.2897*(2^re));     %Lof Estimate formula

flag=0;
for j=1:Hmax        %Finds rightmost zero in bitmap for normal nodes
    if(temp(j,2)==0)
        flag=1;
        break;
    end
end

if(flag==1)
    rn=j-1;
else
    rn=j;
end
n_hat=floor(1.2897*(2^rn));   %Lof Estimate formula

flag=0;
for j=1:Hmax       %Finds rightmost zero in bitmap for periodic nodes
    if(temp(j,3)==0)
        flag=1;
        break;
    end
end

if(flag==1)
    rp=j-1;
else
    rp=j;
end
p_hat=floor(1.2897*(2^rp));  %Lof Estimate formula

B=temp;
a=e_hat;
b=p_hat;
c=n_hat;
d=E2;
e=E3;
f=count2;
g=count-1;
