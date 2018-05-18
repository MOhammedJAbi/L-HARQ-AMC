function Np=Np_expectation(i,alpha,R,M,SNR)

Np=0;

%%%%%%%%%  G1  %%%%%%%%%
if i<=M+1
    for j=1:M+1
        Np=Np+transition_probability(i,j,1,1,alpha,R,SNR,M);
    end
    for j=1:M+1
        for z=1:M+1
            Np=Np+transition_probability(i,[j z],1,2,alpha,R,SNR,M);
        end
    end
    return
end

%%%%%%%%%  G2  %%%%%%%%%
% if(i>(M+1)*(M+1)&&i<(M+1)*(M+2))
%     k=i-(M+1)*(M+1);
%     for j=1:M+1
%         for z=1:M+1
%             Np=Np+transition_probability([M+1 k],[j z],2,2,alpha,R,SNR,M);
%         end
%     end 
% end
% if(i<=(M+1)*(M+2)&&i>M+1)
%     k=ceil(i/(M+1));
%     t=i-(M+1)*(k-1);
%     for j=1:M+1
%         Np=Np+transition_probability([k-1 t],j,2,1,alpha,R,SNR,M);
%     end
%     return
% end

if(i<=(M+1)*(M+2)&&i>M+1)
    k=ceil(i/(M+1));
    t=i-(M+1)*(k-1);
    for j=1:M+1
        Np=Np+transition_probability([k-1 t],j,2,1,alpha,R,SNR,M);
    end
    for j=1:M+1
         for z=1:M+1
            Np=Np+transition_probability([k-1 t],[j z],2,2,alpha,R,SNR,M);
         end
    end
    return
end

%%%%%%%%%  G3  %%%%%%%%%

%->G1
if(i<=(M+1)*(M+3)&&i>(M+1)*(M+2))
    k=i-(M+1)*(M+2);
    for j=1:M+1
        Np=Np+transition_probability(k,j,3,1,alpha,R,SNR,M);
    end
    return
end
