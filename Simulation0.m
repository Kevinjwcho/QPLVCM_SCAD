cvx_setup
cvx_quiet true
%-------------------------------------
%          Initial setting
%-------------------------------------

tau_list=(0.3:0.2:0.7);
M = length(tau_list);
N = 100;
p=3;
p0=p-1;
Co = 0.5.^abs(outerop(1:p0,1:p0,'-')); %covariance matrix
thre=10^-5*p;
iterbdd=20;
J=100; %discritization number
z=(1/J:1/J:1); %sequence from 0 to 1 by 0.01
hoptf=0.15; % bandwidth problem

%-----------------------------------
%          Data settings
%-----------------------------------

 MSE = zeros(M,1);
 ISE = zeros(M,p);
 hlam = zeros(M,1);
 Cn=1;
%Cn=sqrt(log(p));
%Cn=log(p);

tic %start stopwatch timer
for m=1:M
    
    tau = tau_list(m);




    
    txt = ['Processing tau=', num2str(tau), ', N=', num2str(N)];
    disp(txt)

    data = Data_GN_S0(N,Co,p0);
    U = data{2};
    X = data{3};
    Y=data{1};


lam_array= (0:0.005:0.5);
BIC_array=zeros(1,length(lam_array));
V_array=zeros(1,length(lam_array));
C_array=zeros(1,length(lam_array));

for lamk=1:length(lam_array) % select the best lambda with BIC
lam=lam_array(lamk);

txt = ['BIC lamda is ', num2str(lam)];
disp(txt)


[mhat, iter, error] = QPLVCM_SCAD(p, z, J,  N, X, Y, U, hoptf,lam, thre, iterbdd, tau); %jointly selection 


mY=(X*mhat(1:p,:))+(X*mhat((1:p)+p,:)).*(outerop(U,-z,'+')/hoptf);   
Yres=(Y*ones(1,J)-mY);
kw=norm_ker(outerop(U,z,'-'),hoptf);

V=nnz(mhat((1:p)+p,:).^2*ones(J,1));
C=nnz(mhat(1:p,:).^2*ones(J,1))-V;
V_array(lamk)=V;
C_array(lamk)=C;

% BIC_array(lamk)=log(sum(sum(qr_obj((Yres.*kw), tau)))/(N*J))+(V*log(N*hoptf)/(N*hoptf)+C*log(N)/N)*Cn;
 
if(tau == 0.3)
    mhatlm_name = ['result/mhatlm_tau3_', num2str(lam),'.csv'];
    csvwrite(mhatlm_name, mhat);
elseif(tau == 0.5)
    mhatlm_name = ['result/mhatlm_tau5_', num2str(lam),'.csv'];
    csvwrite(mhatlm_name, mhat);
elseif(tau == 0.7)
    mhatlm_name = ['result/mhatlm_tau7_', num2str(lam),'.csv'];
    csvwrite(mhatlm_name, mhat);
end
end

end

save Simul_0.mat

toc



