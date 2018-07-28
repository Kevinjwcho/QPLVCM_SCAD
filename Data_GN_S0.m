%%%QKLASSO settings(1)
function dat = Data_GN_S0(n, Co, p)

Mu = zeros(1,p);
%x = mvnrnd(Mu,Co,n);
x = [ones(n,1) mvnrnd(Mu,Co,n)];
%e = trnd(3,n,1);
e = normrnd(0,1,n,1);
u = unifrnd(0,1,n,1);
y = a1(u).*x(:,2) + a2(u).*x(:,3)+e;

dat{1} = y;
dat{2} = u;
dat{3} = x;

end