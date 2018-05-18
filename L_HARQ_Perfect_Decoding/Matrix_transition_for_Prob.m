function P=Matrix_transition_for_Prob(policy,M,R,SNR)
%% policy is vector of actions 
% M= number of discritization
% R rate
% SNR signal to noise ratio

snr=10^(SNR/10);
s=(M+1)*(M+3);

%% Initiate the matrix

%%%%%%%%%  G1-->Gi i=1,...6  %%%%%%
P11=zeros(M+1,M+1);
P12=zeros(M+1,(M+1)*(M+1));
P13=P11;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%  G2-->Gi i=1,...6  %%%%%%
P21=zeros((M+1)*(M+1),M+1);
P22=zeros((M+1)*(M+1),(M+1)*(M+1));
P23=P21;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%  G3-->Gi i=1,...6  %%%%%%
P31=P11;
P32=P12;
P33=P11;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% transition matrix

%%%%%%%%%  G1-->Gi i=1,...6  %%%%%%

for i=1:M+1
    for j=1:M+1
        P11(i,j)=transition_probability(i,j,1,1,policy(i),R,SNR,M);
    end
end


for i=1:M
    temp=1;
    for j=1:M+1
        for t=1:M+1
            P12(i,temp)=transition_probability(i,[j t],1,2,policy(i),R,SNR,M);
            temp=temp+1;
        end
    end
end

for i=1:M
    for j=i:M+1
        P13(i,j)=transition_probability(i,j,1,3,policy(i),R,SNR,M);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%  G2-->Gi i=1,...6  %%%%%%
temp=1;
for i=1:M+1
    for j=1:M+1
        for t=1:M+1
            P21(temp,t)=transition_probability([i j],t,2,1,policy(i*(M+1)+j),R,SNR,M);
        end
        temp=temp+1;
    end
end


temp1=1;
for i=1:M+1
    for j=1:M+1
        temp2=1;
        for t=1:M+1
            for k=1:M+1
                P22(temp1,temp2)=transition_probability([i j],[t k],2,2,policy(i*(M+1)+j),R,SNR,M);
                 temp2=temp2+1;
            end
        end
        temp1=temp1+1;
    end
end

temp=1;
for i=1:M+1
    for j=1:M+1
        for t=j:M+1
            P23(temp,t)=transition_probability([i j],t,2,3,policy(i*(M+1)+j),R,SNR,M);
        end
        temp=temp+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%  G3-->Gi i=1,...6  %%%%%%
for i=1:M+1
    for j=1:M+1
        P31(i,j)=transition_probability(i,j,3,1,policy((M+1)*(M+2)+i),R,SNR,M);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P=[P11 P12 P13;P21 P22 P23;P31 P32 P33];
