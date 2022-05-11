clc;
%clear all;
%-------------------- Initialization ---------------------------------% 
%C=[0.01 0.02 0.03 0.04 0.05 0.06 0.07]; %total channels
C=0.15.*ones(1,30);
M=length(C); %total no of Channels
N_limit=50;
E=linspace(1,N_limit,N_limit); %Set of Emergency nodes
%P=linspace(1,100,100); %Set of Periodic nodes
%N=linspace(1,100,100); %Set of Normal nodes
P=E;
N=E;
Ne=length(E); %total emergency nodes
Nn=length(N); %total normal nodes
Np=length(P); %total periodic nodes
total_nodes=Ne+Np+Nn;  %total nodes
Pset = []; %Channels used by PUs
Sset = []; %SU usable channel list
alpha=0.2;
Mf=0;  %Size of Sset
mf=0;   %Broadcasting Channel. i.e. first free channel 
lambda_p=1; %Packet arrival rate of periodic nodes
lambda_e=lambda_p*(1-alpha); %Packet arrival rate of emergency nodes
lambda_n=lambda_p*(1+alpha); %Packet arrival rate of normal nodes

Ea=[]; %Set of emergency nodes with active queue
Na=[]; %Set of Normal nodes with active queue
Pa=[]; %Set of periodic nodes with active queue
Hmax= ceil(log2(max(max(Ne,Nn),Np)))+1;% Maximum Hash Value. Plus one because log2(16)=4 but 16=10000
B=zeros(Hmax,3); %Set used in EW

Tframe=50;           %Total frame time
Tsense=0.2;           %Time spent for sensing the channel
Tspent=0;             %Remaining Time
time=0;              %total time
%-------------------- Queues ---------------------------------% 
Qe=cell(Ne,1);  % Queue for emergency nodes
delay_e=cell(Ne,1); %Delay per packet for a node
Qn=cell(Nn,1);  % Queue for normal nodes
delay_n=cell(Nn,1);
Qp=cell(Np,2);  % Queue for periodic nodes- structure description in the notebook
delay_p=cell(Np,1);
Stack=cell(M,1);  %Stack of every channel,
Prev_set=[];        %Set of periodic nodes who have previosly reserved access 

%----------------------------Simulation Starts----------------------------
run_limit=50;
F_limit=50;
lambda_set=(0:0.3:3);
size_lambda=length(lambda_set);
mean_e=zeros(run_limit,size_lambda);
mean_n=zeros(run_limit,size_lambda);
mean_p=zeros(run_limit,size_lambda);
h = waitbar(0,'Initializing..');
h1 = waitbar(0,'Arrival rate');
for RUN=1:run_limit
for n=0:0.3:3
    lambda_p=n; %Packet arrival rate of periodic nodes
    %lambda_e=lambda_p*(1-alpha); %Packet arrival rate of emergency nodes
    %lambda_n=lambda_p*(1+alpha); %Packet arrival rate of normal nodes
    lambda_e=lambda_p; 
    lambda_n=lambda_p;
    Ea=[]; %Set of emergency nodes with active queue
    Na=[]; %Set of Normal nodes with active queue
    Pa=[]; %Set of periodic nodes with active queue
    Qe=cell(Ne,1);  % Queue for emergency nodes
    delay_e=cell(Ne,1);
    Qn=cell(Nn,1);  % Queue for normal nodes
    delay_n=cell(Nn,1);
    Qp=cell(Np,2);  % Queue for periodic nodes- structure description in the notebook
    delay_p=cell(Np,1);
    Stack=cell(M,1);  %Stack of every channel, needs to be 2D format=(node,frames reserved)
    Prev_set=[];        %Set of periodic nodes who have previosly reserved access
    time=0; 
    Tspent=0;

 for frames=1:F_limit
%-------------------- Packet Generation ---------------------------------% 
    for i=1:Ne    %Packet generation for emergency nodes
        k = poissrnd(lambda_e);
        if k ~= 0
            Ea= union(Ea,i);
            for ii=1:k
                Qe{i}= [Qe{i} [1,time]];
            end
        end
    end
   for i=1:Nn    %Packet generation for normal nodes
       k = poissrnd(lambda_n);
        if k ~= 0
            Na= union(Na,i);
            for ii=1:k
                Qn{i}= [Qn{i} [1,time]];
            end
        end
   end
   for i=1:Np     %Packet generation for periodic nodes
      k = poissrnd(lambda_p);
       if k ~= 0
          if (any(Prev_set==i))
              %Qp(i,2)= Qp(i,2)+ k;
              for ii=1:k
                  Qp{i,2}= [Qp{i,2} [1,time]];
              end
          else
                Pa= union(Pa,i);
                for ii=1:k
                    Qp{i,1}= [Qp{i,1} [1,time]];
                end
          end
       end
   end
          
%-------------------- BW1 ---------------------------------%
    [Pset, Sset]= bw1(C);
    mf=Sset(1);
    Mf=length(Sset);
    Tbw1=ceil(Tsense*(mf + 1));  % Time spent in BW1
    time=time+Tbw1;

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
     time=time+Tew;
%-------------------- BW2 ---------------------------------%     
     [Se, Sp, Sn, Pe, Pp, Pn, CM, Stack]= bw2(e_hat, p_hat, n_hat, Ea, Pa, Na, Sset, M, Stack,Trem);
     
     Tbw2=10;       % Need to find out how much time is spent
    
     Tspent=Tspent+Tbw2;
     Trem=Tframe-Tspent;
     time=time+Tbw2;

%------------------------CDTW--------------------------------
    Pe_mat=ones(1,length(Se));         %Every channel gets its own Pe
    Pe_mat=Pe.*Pe_mat;
    Pn_mat=ones(1,length(Sn));          %Every channel gets its own Pn
    Pn_mat=Pn.*Pn_mat;
    Pp_mat=ones(1,length(Sp));         %Every channel gets its own Pp
    Pp_mat=Pp.*Pp_mat;
    [Qe, Qn, Qp, Stack, Prev_set,Ea,Na,Pa,time,delay_e,delay_n,delay_p]=CDTW(CM,Qe,Qn,Qp,Stack,Se,Sp,Sn,Pe_mat,Pp_mat,Pn_mat,Trem, Prev_set,Ea,Na,Pa,time,delay_e,delay_n,delay_p);
 end
 index=find(abs(lambda_set-n)<1e-5);
     emergency_delay=zeros(1,Ne);
     for xy=1:Ne
         emergency_delay(xy)=mean(delay_e{xy});
     end
     indices=find(isnan(emergency_delay));
     for ii=1:length(indices)
            emergency_delay(indices(ii))=0;
     end
     mean_e(RUN,index)=mean(emergency_delay);

     normal_delay=zeros(1,Nn);
     for xy=1:Nn
         normal_delay(xy)=mean(delay_n{xy});
     end
     indices=find(isnan(normal_delay));
     for ii=1:length(indices)
            normal_delay(indices(ii))=0;
     end
     mean_n(RUN,index)=mean(normal_delay);

 waitbar(n/lambda_set(end),h1,sprintf('Arrival rate = %d%',n));
end

 waitbar(RUN/run_limit,h,sprintf('Current run = %d%',RUN));
end

ADE=mean(mean_e);
ADN=mean(mean_n);

figure;
x_axis=lambda_set;
plot(x_axis,ADE,'k-o',x_axis,ADN,'k--*');
axis tight;
%grid on;
title('Average Delay v/s \lambda');
xlabel('\lambda');
ylabel('Average Delay');
legend('Emergency Nodes','Normal Nodes');






