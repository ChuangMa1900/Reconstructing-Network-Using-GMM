function [C,y]=Ge_LS_coupling_dynamic(w,k,delta_t,S_L,XI)
%（耦合动力学）生成线性方程组数据
%输入：w 网络
%      k 方程组的个数
%delta_t  仿真间隔        delta_t=10^-4;
%delta_t*S_L 采样间隔%    S_L=10;
%     XI  数据采集误差%   本身带有误差，所以XI=0
%输出：元胞数组C与y，其中C{i}A_i=y{i}

    n=length(w);
    CC=eye(n);
    yy=zeros(n,1);
    for i=1:k
        [L,T]=Ge_data_coupling_dynamic(w,S_L,delta_t);
        %加噪音
        L=L+XI*randn(size(L));
    
        X_t=(L(:,:,1)+L(:,:,end))/2;
        D_X_t=(L(:,:,end)-L(:,:,1))/(T(end)-T(1));
        B=D_X_t(:,1)+X_t(:,2)+X_t(:,3);
        s=repmat(X_t(:,1),[1,n]);
        CC=[CC,s'-s];
        yy=[yy,B];
    end
    [n,m0]=size(yy);
    C=cellfun(@(x)reshape(x,n,m0)',mat2cell(CC,ones(1,n),n*m0)','UniformOutput',false);
    y=mat2cell(yy',m0,ones(1,n));
end