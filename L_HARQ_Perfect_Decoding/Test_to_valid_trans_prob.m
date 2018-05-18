clear
clc

%transition_probability(i,j,Gi,Gj,action,R,SNR,M)

M=3;
R=5;
SNR=15;

P=zeros(1,M+1);

for i=1:M+1

   P(i)=transition_probability(1,i,1,1,2,R,SNR,M);

end

P
sum(P,2)


P=zeros(M+1,M+1);


for i=1:M+1
    for j=1:M+1
        %P(i,j)=P(i,j)+transition_probability([2 3],[i j],2,3,[0 1],R,SNR,M);
        P(i,j)=transition_probability([3 2],[i j],5,2,0.4,R,SNR,M);
    end
end

%P
 sum(P,2)
 sum(sum(P,2))
 
