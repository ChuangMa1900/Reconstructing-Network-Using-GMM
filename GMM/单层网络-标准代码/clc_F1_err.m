function F_err=clc_F1_err(A,w,s)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A:预测值
%w：真实值
%s：真实类型值
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aa=mean(s);
    err=norm(A(:)-w(:));
    w1=A;
    w1(A>aa)=1;
    w1(A<aa)=0;
    w=w-diag(diag(w));
    w(w>0)=1;
    w1=w1-diag(diag(w1));
    F=clc_F1(w1,w);
    F_err=[F,err];
end