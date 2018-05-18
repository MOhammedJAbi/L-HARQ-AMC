clear;
clc;

L_set=[10:10:100];
N_realization=1e4;
SNRdBmin=0;
path_loss=3;
SNR=@(d) d.^-path_loss;




for L_ind=1:length(L_set)
    L=L_set(L_ind);
    accumulated_rate=1e-6.*ones(1,L);
    accumulated_rate_r=0.*ones(1,L);
    agregate=zeros(1,L);
    cont=1
    for i=1:N_realization
        v_distance=sqrt(unifrnd(0,1,1,L));%sqrt(unifrnd(0,1,1,L));rand
        v_distance(v_distance<0.1)=0.1;
        SNRdB_vector=SNR(v_distance);
        SNR_lin_vector=exprnd(SNRdB_vector);
        rate_instan=log2(1+SNR_lin_vector);
        [Y,J]=max(rate_instan./accumulated_rate);
        agregate(J)=agregate(J)+1;
        accumulated_rate(J)=accumulated_rate(J)+rate_instan(J);
        accumulated_rate_r(cont)=accumulated_rate_r(cont)+rate_instan(cont);
        cont=cont+1;
        cont(cont>L)=1;
    end
    throughput(L_ind)=sum(accumulated_rate-1e-6.*ones(1,L))/N_realization
    throughput2(L_ind)=sum(accumulated_rate_r)/N_realization
end

