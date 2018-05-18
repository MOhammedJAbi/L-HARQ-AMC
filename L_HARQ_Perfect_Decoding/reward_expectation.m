function Rew=reward_expectation(i,alpha,R,SNR,M)

Rew=0;


%%%%%%%%%  G1  %%%%%%%%%
if i<=M+1
    Rew=Rew+R*transition_probability(i,M+1,1,1,alpha,R,SNR,M);
    for j=1:M
        Rew=Rew+R*transition_probability(i,[j M+1],1,2,alpha,R,SNR,M);
    end
    for j=1:M
        Rew=Rew+R*transition_probability(i,[M+1 j],1,2,alpha,R,SNR,M);
    end
    Rew=Rew+2*R*transition_probability(i,[M+1 M+1],1,2,alpha,R,SNR,M);
    Rew=Rew+R*transition_probability(i,M+1,1,3,alpha,R,SNR,M);
    return
end

%%%%%%%%%  G2  %%%%%%%%%
% if(i>(M+1)*(M+1)&&i<(M+1)*(M+2))
%     k=i-(M+1)*(M+1);
%     Rew=Rew+R*transition_probability([M+1 k],M+1,2,3,alpha,R,SNR,M);
%     for j=1:M
%         Rew=Rew+R*transition_probability([M+1 k],[j M+1],2,2,alpha,R,SNR,M);
%     end
%     for j=1:M
%         Rew=Rew+R*transition_probability([M+1 k],[M+1 j],2,2,alpha,R,SNR,M);
%     end
%     Rew=Rew+2*R*transition_probability([M+1 k],[M+1 M+1],2,2,alpha,R,SNR,M);   
% end
% if (i>(M+1)&&i<=(M+1)*(M+1))
%     k=ceil(i/(M+1));
%     t=i-(M+1)*(k-1);
%     %->G3
%     Rew=Rew+R*transition_probability([k-1 t],M+1,2,3,alpha,R,SNR,M);
% end

if(i<=(M+1)*(M+2)&&i>M+1)
    k=ceil(i/(M+1));
    t=i-(M+1)*(k-1);
    Rew=Rew+R*transition_probability([k-1 t],M+1,2,1,alpha,R,SNR,M);
    for j=1:M
        Rew=Rew+R*transition_probability([k-1 t],[j M+1],2,2,alpha,R,SNR,M);
    end
    for j=1:M
        Rew=Rew+R*transition_probability([k-1 t],[M+1 j],2,2,alpha,R,SNR,M);
    end
    Rew=Rew+2*R*transition_probability([k-1 t],[M+1 M+1],2,2,alpha,R,SNR,M);
    Rew=Rew+R*transition_probability([k-1 t],M+1,2,3,alpha,R,SNR,M);

    return
end

%%%%%%%%%  G3  %%%%%%%%%
if (i<=(M+1)*(M+3) && i>(M+1)*(M+2))
    k=i-(M+1)*(M+2);
    %->G1
    Rew=Rew+R*transition_probability(k,M+1,3,1,alpha,R,SNR,M);
    return
end





