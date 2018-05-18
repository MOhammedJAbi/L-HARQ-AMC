function p=transition_probability(i,j,Gi,Gj,action,R,SNR,M)

%% M= number of discritization
% G=1 group of AMI with \zeta=(1,0)
% G=2 group of AMI with \zeta=(2,1)
% G=3 group of AMI with \zeta=(2,0)
% G=4 group of AMI with \zeta=(3,2)
% G=5 group of AMI with \zeta=(3,1)
% G=6 group of AMI with \zeta=(3,0)
% Gi=G initial state group
% Gj=G next state  group
% i number of the AMI in previous transmission (can be scalar or vector)
% i=1:M+1
% j number of the AMI in next transmission (can be scalar or vector)
% j=1:M+1
% action=[a1 a2], a1=0 TS, a1=1 1P
%%

p=0;

%%

if((Gi==1&&Gj==2&&action(1)~=2&&action(1)~=0&&i<M+1))
    Ii=(i-1)*R/M;
    
    if(i<=M&&j(1)<=M&&j(2)<=M)
        
        Ij=[(j(1)-1)*R/M (j(2)-1)*R/M];
        Ij_1=[j(1)*R/M j(2)*R/M];
        
        a=(Ij(1)-Ii)/action(1);
        b=(Ij_1(1)-Ii)/action(1);
        c=Ij(2)/(1-action(1));
        d=Ij_1(2)/(1-action(1));
        
        p=proba_joint_mutual_inform_1(SNR,a,b,c,d);
        return
    elseif(i<=M&&j(1)==M+1&&j(2)<=M)
        
        Ij=[(j(1)-1)*R/M (j(2)-1)*R/M];
        Ij_1=[Inf j(2)*R/M];
        
        
        a=(Ij(1)-Ii)/action(1);
        b=Inf;
        c=Ij(2)/(1-action(1));
        d=Ij_1(2)/(1-action(1));
        
        p=proba_joint_mutual_inform_1(SNR,a,b,c,d);
        return
    elseif(i<=M&&j(1)==M+1&&j(2)==M+1)
        
        Ij=[(j(1)-1)*R/M (j(2)-1)*R/M];
        
        a=(Ij(1)-Ii)/action(1);
        b=Inf;
        c=Ij(2)/(1-action(1));
        d=Inf;
        
        p=proba_joint_mutual_inform_1(SNR,a,b,c,d);
        return
    elseif(i<=M&&j(1)<=M&&j(2)==M+1)
        
        Ij=[(j(1)-1)*R/M (j(2)-1)*R/M];
        Ij_1=[j(1)*R/M Inf];
        
        
        a=(Ij(1)-Ii)/action(1);
        b=(Ij_1(1)-Ii)/action(1);
        c=Ij(2)/(1-action(1));
        d=Inf;
        
        p=proba_joint_mutual_inform_1(SNR,a,b,c,d);
        return
    end
end
%%
if((Gi==1&&Gj==3&&action(1)==2&&i<=M))
    Ii=(i-1)*R/M;
    if(j<=M)
        
        Ij=(j-1)*R/M;
        Ij_1=j*R/M;
        
        a=(Ij-Ii);
        b=(Ij_1-Ii);
        p=proba_joint_mutual_inform_1(SNR,a,b,a,b);
        return
    elseif(j==M+1)
        
        Ij=(j-1)*R/M;
        
        a=(Ij-Ii);
        p=proba_joint_mutual_inform_1(SNR,a,Inf,a,Inf);
        return
        
    end
    
end


%%
if((Gi==1&&Gj==1&&action(1)==0)||(Gi==2&&Gj==1&&action(1)==0)||(Gi==3&&Gj==1&&action(1)==0))
    if(j<=M)
        
        Ij=(j-1)*R/M;
        Ij_1=j*R/M;
        
        p=proba_joint_mutual_inform_1(SNR,Ij,Ij_1,Ij,Ij_1);
        return
        
    elseif(j==M+1)
        p=proba_joint_mutual_inform_1(SNR,R,Inf,R,Inf);
        return
    end
end

%%

if(Gi==2&&Gj==2&&i(2)<=M&&action(1)~=2&&action(1)~=0)
    Ii=(i-1)*R/M;
    
    if(j(1)<=M&&j(2)<=M)
        
        Ij=[(j(1)-1)*R/M (j(2)-1)*R/M];
        Ij_1=[j(1)*R/M j(2)*R/M];
        
        a=(Ij(1)-Ii(2))/action(1);
        b=(Ij_1(1)-Ii(2))/action(1);
        c=Ij(2)/(1-action(1));
        d=Ij_1(2)/(1-action(1));
        
        p=proba_joint_mutual_inform_1(SNR,a,b,c,d);
        return
    elseif(j(1)==M+1&&j(2)<=M)
        
        
        Ij=[(j(1)-1)*R/M (j(2)-1)*R/M];
        Ij_1=[Inf j(2)*R/M];
        
        a=(Ij(1)-Ii(2))/action(1);
        b=Inf;
        c=Ij(2)/(1-action(1));
        d=Ij_1(2)/(1-action(1));
        
        p=proba_joint_mutual_inform_1(SNR,a,b,c,d);
        return
    elseif(j(1)==M+1&&j(2)==M+1)
        
        Ij=[(j(1)-1)*R/M (j(2)-1)*R/M];
        
        a=(Ij(1)-Ii(2))/action(1);
        b=Inf;
        c=Ij(2)/(1-action(1));
        d=Inf;
        
        p=proba_joint_mutual_inform_1(SNR,a,b,c,d);
        return
    elseif(j(1)<=M&&j(2)==M+1)
        
        Ij=[(j(1)-1)*R/M (j(2)-1)*R/M];
        Ij_1=[j(1)*R/M Inf];
        
        a=(Ij(1)-Ii(2))/action(1);
        b=(Ij_1(1)-Ii(2))/action(1);
        c=Ij(2)/(1-action(1));
        d=Inf;
        
        p=proba_joint_mutual_inform_1(SNR,a,b,c,d);
        return
    end
end

%%


if(Gi==2&&Gj==3&&i(2)<=M&&action(1)==2)
    Ii=(i(2)-1)*R/M;
    if(j<=M)
        
        Ij=(j-1)*R/M;
        Ij_1=j*R/M;
        
        a=(Ij-Ii);
        b=(Ij_1-Ii);
        p=proba_joint_mutual_inform_1(SNR,a,b,a,b);
        return
    elseif(j==M+1)
        
        Ij=(j-1)*R/M;
        
        a=(Ij-Ii);
        p=proba_joint_mutual_inform_1(SNR,a,Inf,a,Inf);
        return
        
    end
    
end

%%








