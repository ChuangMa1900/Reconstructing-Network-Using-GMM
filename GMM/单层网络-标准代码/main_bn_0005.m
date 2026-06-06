function main_bn_0005
%噪音项
XI=0.0005;
D=5;%D次实验
i_r=[0.5,0.8]; %方程组的个数比例
res_all=cell(D,length(i_r));

%导入网络，设置权重
load bn_w
for d=1:D
    a=1;
    w=a*w;
    s=[0;a];
    m=length(s);%解的个数
    n=length(w);

    %生成数据
    [M0,y0]=Ge_LS_coupling_dynamic(w',2*n,10^-4,10,XI);

    %参数设置
    t_max=1000;
    t_max1=1000;
    MU1=0.001;

    for ii=1:length(i_r)      
        T=ceil(n*i_r(ii))
        res_A1=zeros(n,n);
        res_A2=zeros(n,n);
        for i=1:n
            C=M0{i}(1:T,:);
            y=y0{i}(1:T,:);

            %5次交叉验证的lasso
            [B,FitInfo]=lasso(C,y,'cv',5,'Intercept',false,'Standardize',false,'Lambda',linspace(0,1));
            idxLambda1SE = FitInfo.Index1SE;
            A1 = B(:,idxLambda1SE);
            res_A1(:,i)=A1;

            % %signal_lasso
            % [A2,~]=sol_cvx_signal_lasso_normal(C,y,n,a,0.6,0.4,0.1);
            % res_A2(:,i)=A2;
        end
        %lasso
        [~,s1]=kmeans(res_A1(:),m,'Start',[0;max(res_A1(:))]);
        F_err1=clc_F1_err(res_A1,w,s1);

        % %signal_lasso
        F_err2=[0,0];
        % F_err2=clc_F1_err(res_A2,w,s);

        %GMM
        [res_A3,s3]=RN_GMM_EM(M0,y0,m,T,t_max,t_max1,MU1);
        F_err3=clc_F1_err(res_A3,w,s3);
        [F_err1;F_err2;F_err3]
        res_all{d,ii}=[F_err1;F_err2;F_err3];

    end
    save res\bn_0005 res_all
end





