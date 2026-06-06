function [XX,s]=RN_GMM_err(M0,y0,m,T,t_max,t_max1,MU)
%GMM单层/符号网络重构
%Copyrigt:马闯 chuang_m@126.com
%输入：M0{}(1:tt,:)XX=y0{i}(1:tt,:
%      m：表示解类型的个数（符号网络，3类型；单层网络：2类型）
%   t_max：最大迭代次数
%  t_max1：最大迭代次数
%      MU：迭代终止误差
%输出：
%    XX： 方程结果
%     s：类型值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load MY,
% m=3;   %值的个数
% tt=20; %方程组个数
% t_max=1000;
% t_max1=1000;
% MU=0.01;
%%%%%%%%%%%%%%%%%%%%%%
%%%% 参数设置%%%%%%%%
sig_0=0.01;
sig_1=0.01;
n=length(y0);
A_all=cell(1,n);
y_all=cell(1,n);
ATA_all=cell(1,n);
ATy_all=cell(1,n);
XX=zeros(n,n);
for i=1:n
    A=M0{i}(1:T,:);
    y=y0{i}(1:T,:);
    ATA = A'*A;
    ATy = A'*y;
    XX(:,i)=(ATA+(n/T)*eye(n)) \ (ATy);

    A_all{i}=A;
    y_all{i}=y;
    ATA_all{i}=ATA;
    ATy_all{i}=ATy;
end
sigmaG=std(XX(:));
sigma=((n/T)^0.5)*sigmaG;

if m==3 %符号网络
%     s=[min(XX(:));0;max(XX(:))];
    [~,s]=kmeans(XX(:),m,'Start',[min(XX(:));0;max(XX(:))]);
    s=sort(s);
end

if m==2 %一般网络
%     s=[min(XX(:));max(XX(:))];
    [~,s]=kmeans(XX(:),m,'Start',[0;max(XX(:))]);
    s=sort(s);
end

a=ones(1,m)/m;
Q=a.*exp((-(XX(:)-s').^2)/(2*sigmaG^2));
Q=Q./(sum(Q,2)+10^(-10));
%%%%%%%%%%重启参数%%%%%%%%%%%%%%%%%
s0=s;
Q0=Q;
a0=a;
sigma0=sigma;
sigmaG0=sigmaG;
XX0=XX;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%方法迭代部分%%%%%%%%%%%%%
for tttt=1:10
     kkk=1;
    for tt=1:t_max
        x_pre = XX;
        if  T<=1.5*n
            sd=max(sigma^2/sigmaG^2,n/T);
        else
            sd=max(sigma^2/sigmaG^2,0);
        end
        for k2=1:n
            ATA=ATA_all{k2};
            ATy=ATy_all{k2};
            ii=(k2-1)*n+1:(k2-1)*n+n;
            XX(:,k2)=(ATA+sd*eye(n))\(ATy+sd*Q(ii,:)*s);
        end

        sse=0;
        for k1=1:n
            A=A_all{k1};
            y=y_all{k1};
            x=XX(:,k1);
            sse=sse+sum((A*x-y).^2);
        end
        sigma=max(sqrt(sse/(n*n)),sig_1);

        s_p=s;
        s_pre=zeros(m,1);
        for ttt=1:t_max1
            Q=a.*exp((-(XX(:)-s').^2)/(2*sigmaG^2));
            Q=Q./(sum(Q,2)+10^(-10));
            a=sum(Q)./sum(sum(Q));
            s=sum(XX(:).*Q)./(sum(Q,1)+10^(-10));
            s=s';
            sigmaG = max(sqrt(sum(sum(Q.*((XX(:)-s').^2))))/sqrt(n*n),sig_0);
            if max(abs(s_pre(:)-s(:)))<MU
                break;
            end
            s_pre=s;
        end

        if max(abs(x_pre(:)-XX(:)))<MU && tt>=2
            break;
        end

        %%%%算法重启%%%%%%%%%%%%%%%
        if max(abs(sort(s_p)-sort(s)))>mean(abs(s0))
            s=s0;
            Q=Q0;
            a=a0;
            sigma=(0.8^kkk)*sigma0;
            kkk=kkk+1;
            sigmaG=sigmaG0;
            if kkk>=100
                break;
            end

        end
    end
    %%%%如果耦合强度重合，算法重启，加到初始耦合强度的距离%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ds=abs(s-s')+10*eye(m);
    if min(ds(:))<0.01
        sigma=sigma0;
        sigmaG=sigmaG0;
        a=a0;
        s=1.1.*s0;
        Q=a.*exp((-(XX0(:)-s').^2)/(2*sigmaG^2));
        Q=Q./(sum(Q,2)+10^(-10));
        Q0=Q;
        s0=s;
    else
        break;
    end
end

XX=XX-diag(diag(XX));
if min(ds(:))<0.01
    if m==3 %符号网络
        [~,s]=kmeans(XX(:),m,'Start',[min(XX(:));0;max(XX(:))]);
        s=sort(s);
    end
    if m==2 %一般网络
        [~,s]=kmeans(XX(:),m,'Start',[0;max(XX(:))]);
        s=sort(s);
    end
end



