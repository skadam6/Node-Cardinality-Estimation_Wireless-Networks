clc;
%clear all;
%-------------------- Initialization ---------------------------------% 
C=[0.01 0.02 0.03 0.04 0.05 0.06 0.07]; %total channels
M=length(C); %total no of Channels
N_limit=16;
E=linspace(1,N_limit,N_limit); %Set of Emergency nodes
%P=linspace(1,100,100); %Set of Periodic nodes
%N=linspace(1,100,100); %Set of Normal nodes
P=E;
N=E;
Ne=length(E); %total emergency nodes
Nn=length(N); %total normal nodes
Np=length(P); %total periodic nodes
%EWP2=[];      %Variable used to plot K
%EWP3=[];      %Variable used to plot R

total_nodes=Ne+Np+Nn;  %total nodes

delta= 0.05;  %threshold probability
Pset = []; %Channels used by PUs
Sset = []; %SU usable channel list

Mf=0;  %Size of Sset
mf=0;   %Broadcasting Channel. i.e. first free channel 
%Pprev=0; % Periodic nodes in previous slots
lambda_e=1; %Packet arrival rate of emergency nodes
lambda_n=1; %Packet arrival rate of normal nodes
lambda_p=1; %Packet arrival rate of periodic nodes

Ea=[]; %Set of emergency nodes with active queue
Na=[]; %Set of Normal nodes with active queue
Pa=[]; %Set of periodic nodes with active queue
Hmax= ceil(log2(max(max(Ne,Nn),Np)))+1;% Maximum Hash Value. Plus one because log2(16)=4 but 16=10000
B=zeros(Hmax,3); %Set used in EW

Tframe=30;           %Total frame time
Tsense=0.1;           %Time spent for sensing the channel
%Tup=0.5;                %Time required for uploading 1 packet  
%Tack=0.2;             %Time req. for acknowlegdment/ downlink
Tspent=0;             %Remaining Time
%-------------------- Queues ---------------------------------% 
Qe=zeros(Ne,1);  % Queue for emergency nodes
Qn=zeros(Nn,1);  % Queue for normal nodes
Qp=zeros(Np,2);  % Queue for periodic nodes- structure description in the notebook
Stack=cell(M,1);  %Stack of every channel, needs to be 2D format=(node,frames reserved) 
Prev_set=[];        %Set of periodic nodes who have previosly reserved access 
F_limit=50;
%N_start=100;
TP_e=zeros(1,F_limit);
TP_n=zeros(1,F_limit);
TP_p=zeros(1,F_limit);

%h = waitbar(0,'Initializing..');
 for frames=1:F_limit
%-------------------- Packet Generation ---------------------------------% 
    for i=1:Ne    %Packet generation for emergency nodes
        k = poissrnd(lambda_e);
        if k ~= 0
            Ea= union(Ea,i);
            Qe(i)= Qe(i)+k;
        end
    end
    for i=1:Nn    %Packet generation for normal nodes
        k = poissrnd(lambda_n);
        if k ~= 0
            Na= union(Na,i);
             Qn(i)= Qn(i)+k;
        end
    end
   for i=1:Np     %Packet generation for periodic nodes
      k = poissrnd(lambda_p);
       if k ~= 0
          if (any(Prev_set==i))
              Qp(i,2)= Qp(i,2)+ k;
            else
                Pa= union(Pa,i);
                Qp(i,1)= k;
          end
       end
    end
      %Ea = randperm(N_limit,n);
      %Pa = randperm(N_limit,n);
      %Na = randperm(N_limit,n);
    
%-------------------- BW1 ---------------------------------%
    [Pset, Sset]= bw1(C);
    mf=Sset(1);
    Mf=length(Sset);
    Tbw1=Tsense*(mf + 1);  % Time spent in BW1


%-------------------- EW ---------------------------------%
     [B, e_hat, p_hat, n_hat, E2, E3, phase2, phase3]=EW(Ea,Na,Pa,Hmax,Ne,Np,Nn);
     Tewp1=ceil((ceil(log2(max(max(Ne,Nn),Np)))+1)/Mf);
     x=size(E2);
     Tewp2=ceil(x(1)/Mf);
     x=size(E3);
     Tewp3=ceil(x(1)/Mf);
     Tew=Tewp1+Tewp2+Tewp3;
     Tspent=Tbw1+Tew;
     Trem=Tframe-Tspent;
%-------------------- BW2 ---------------------------------%     
     [Se, Sp, Sn, Pe, Pp, Pn, CM, Stack]= bw2(e_hat, p_hat, n_hat, Ea, Pa, Na, Sset, M, Stack,Trem);
     
     Tbw2=10;       % Need to find out how much time is spent
    
     Tspent=Tspent+Tbw2;
     Trem=Trem-Tbw2;

%------------------------CDTW--------------------------------
    Pe_mat=ones(1,length(Se));         %Every channel gets its own Pe
    Pe_mat=Pe.*Pe_mat;
    Pn_mat=ones(1,length(Sn));          %Every channel gets its own Pn
    Pn_mat=Pn.*Pn_mat;
    Pp_mat=ones(1,length(Sp));         %Every channel gets its own Pp
    Pp_mat=Pp.*Pp_mat;
    [Qe, Qn, Qp, Stack, Prev_set,Ea,Na,Pa,temp_te,temp_tn,temp_tp]=CDTW(CM,Qe,Qn,Qp,Stack,Se,Sp,Sn,Pe_mat,Pp_mat,Pn_mat,Trem, Prev_set,Ea,Na,Pa);
    TP_e(frames)=temp_te/Ne; 
    TP_n(frames)=temp_tn/Nn;
    TP_p(frames)=temp_tp/Np;
 end
    
    %waitbar(n/ N_limit,h,sprintf('current n= %d%',n));






