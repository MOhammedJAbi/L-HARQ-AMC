function [P_sta,P]=states_steady_probabilities(policy,M,R,SNR)

P=Matrix_transition_for_Prob(policy,M,R,SNR);
pas=1e-2;

[v ,d]=eig(P');
max(max(d))
min(min(d))
max(max(abs(d)))
min(min(abs(d)))
[a ,b]=find(d<=1+pas&d>1-pas);
P_sta=v(:,a)'/sum(v(:,a));