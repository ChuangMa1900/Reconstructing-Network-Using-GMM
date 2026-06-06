function w=Gen_Sig_Net(w,a)
n=length(w);
w=tril(w);
id=find(tril(w));
m=length(id);
w(id(randperm(m,ceil(m*a))))=-1;
w=w+w';