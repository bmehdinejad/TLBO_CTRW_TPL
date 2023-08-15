%% Descriptions:
%% This code is the forward solution of the CTRW-TPL. It uses the "u" variable as Laplace variable. 
%% As we know, the the solutions of the CTRW-TPL impose the Laplace transform on the time variable.
%% The solution used in this code is for Dirichlet boundary condition at the inlet and Neumann
%% boundary condition at the oulet.  
%%=========================================================================
tic
clc
clear all
close all
format short 
%%=========================================================================
%% The characterisctics of soil column and tracer tests:
L=100;             %% The length of the soil column (cm)
A=100;             %% The cross-sectional area of the soil column(cm^2)  
C0=1;              %% The concnentration of contaminant at the inlet
%%=========================================================================
%% The transport parameters of solute transport: 
v=0.57341;                      %% The average pore water velocity (cm/min)    
D=0.63176;                     %% The dispersion coefficent (cm^2/min)
%% ========================================================================
%% The transport parameters of CTRW model:
v_psi=v/L
D_psi=D/(L^2)
alpha_psi=D_psi/v_psi;
distance=L;
x=distance/L;
beta=2.88138;
t2=66643677601.28124;
%% ========================================================================
%% The experimental data at distance of L:
t_meas=[0.1
10
22.5
27
34.02
43.98
54
73.98
76.02
77.52
79.02
80.82
82.02
83.82
85.02
86.52
88.02
91.02
94.02
97.98
102
106.02
115.98
129
139.02
148.98
163.98
178.98
193.98
214.02
234
274.02
313.98];
Cr_meas=[0
0
0.0112
0.0108
0.0108
0.0108
0.0108
0.0633
0.2
0.275
0.3417
0.3767
0.4167
0.4333
0.4583
0.4833
0.5
0.5167
0.5667
0.6
0.65
0.6667
0.9667
1
1
1
1
1
1
1
1
1
1];
%%=========================================================================
%% The calculation of the TPL characteristic time, t1:
Tau = 1/t2; 
if (beta>1.4 && Tau<1e-5)
    t1=(D_psi*(beta-1))/(v_psi^2);
else
    t1=((beta*D_psi*((1/Tau).^(-1-beta)))*igamma(-beta,Tau))/((v_psi^2)*(-Tau +(Tau^beta)*(beta+Tau)*igamma(1-beta,Tau)));
    t2=t1/Tau;
end
% t1=1.01*10^(-1);
%% ========================================================================
%% The transport parameters of eipsilon algorithm 
alpha=0;
tol=1e-9;
T=2*max(t_meas);
M=10;
%% ========================================================================
%% The calculation of psi(u), memory, and structural functions and of C for Landa0:
Landa0=(alpha-log(tol)/(2*T));
psi_Landa0=((1+t2*Landa0)^beta)*(exp(t1*Landa0))*igamma(-beta,t1/t2+t1*Landa0)/igamma(-beta,t1/t2);
M_Landa0=t1*Landa0*psi_Landa0/(1-psi_Landa0);
w_Landa0=sqrt(1+4*Landa0*alpha_psi/(M_Landa0*v_psi));
z_Landa0=(1+w_Landa0)*x/(2*alpha_psi);
E_Landa0=(2*w_Landa0-x*w_Landa0+x)/(2*alpha_psi);
C_Landa0=(exp(E_Landa0)/Landa0)*(((w_Landa0+1)+(w_Landa0-1)*exp(w_Landa0*(x-1)/alpha_psi))/((w_Landa0+1)*exp(w_Landa0/alpha_psi)+(w_Landa0-1)));
%% ========================================================================
%% The calculation of psi(u), memory, and structural functions and of C for u:
for k=1:2*M
    Landa(k,1)=Landa0+i*k*pi/T;
    psi_u(k,1)=((1+t2*Landa(k,1))^beta)*(exp(t1*Landa(k,1)))*(igamma(-beta,t1/t2+t1*Landa(k,1)))/(igamma(-beta,t1/t2));
    M_u(k,1)=t1*Landa(k,1)*psi_u(k,1)/(1-psi_u(k,1));
    w(k,1)=sqrt(1+4*Landa(k,1)*alpha_psi/(M_u(k,1)*v_psi));
    z(k,1)=(1+w(k,1))*x/(2*alpha_psi);
    E(k,1)=(2*w(k,1)-x*w(k,1)+x)/(2*alpha_psi);
    C_u(k,1)=(exp(E(k,1))/Landa(k,1))*((w(k,1)+1)+(w(k,1)-1)*exp(w(k,1)*(x-1)/alpha_psi))/((w(k,1)+1)*exp(w(k,1)/alpha_psi)+(w(k,1)-1));
    a(k,1)=C_u(k,1);
end
Landa;
psi_u;
M_u;
w;
z;
E;
C_u;
a0=0.5*C_Landa0;
a=[a0;a];
for m=1:2*M
    q(m,1)=a(m+1,1)/a(m,1);
end
for m=1:2*M-1
    e(m,1)=q(m+1,1)-q(m,1);
    e(2*M,1)=0;
end
for r=2:M
    for m=1:2*M-2*r+1
        q(m,r)=q(m+1,r-1)*e(m+1,r-1)/e(m,r-1);
    end
       for m=1:2*M-2*r+1
        e(m,r)=q(m+1,r)-q(m,r)+e(m+1,r-1); 
    end
end
q;
e;
d0=a0;
for m=1:M
    d(1,2*m-1)=-q(1,m);
    d(1,2*m)=-e(1,m);
end
d;
for j=1:size(t_meas)
    z(j,1)=exp(i*pi*t_meas(j,1)/T);
    A(1,1)=d0;
    B(1,1)=1+d(1,1)*z(j,1);
    A(1,2)=A(1,1)+d(1,2)*z(j,1)*d0;
    B(1,2)=B(1,1)+d(1,2)*z(j,1);
for n=3:2*M
    A(1,n)=A(1,n-1)+d(1,n)*z(j,1)*A(1,n-2);
    B(1,n)=B(1,n-1)+d(1,n)*z(j,1)*B(1,n-2);
end
A;
B;
h=(1+(d(1,2*M-1)-d(1,2*M))*z(j,1))/2;
R=-h*(1-sqrt(1+d(1,2*M)*z(j,1)/((h)^2)));
A_2M=A(1,2*M-1)+R*A(1,2*M-2);
B_2M=B(1,2*M-1)+R*B(1,2*M-2);
Cr_t_calc(j,1)=(1/T)*exp(Landa0*t_meas(j,1))*real(A_2M/B_2M);
if Cr_t_calc(j,1)<0
    Cr_t_calc(j,1)=0;
end
end
Cr_t_calc
%% ========================================================================
%% Calculating the error value:
for ii=1:size(t_meas)
      df(ii,1)=(Cr_t_calc(ii,1)-Cr_meas(ii,1))^2; 
end
S=sum(df);
error=sqrt(S)
plot(t_meas,Cr_t_calc)
toc








