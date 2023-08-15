%% This code calculates the statistical criteria using the solute parameters estimated by DIM and CTRWMT at distances of L. 
clc
clear all
format long
%% ========================================================================
%% The experimental data at distance of L:
L=100;
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
N_1=size(t_meas);
N_2=N_1(1,1);
%%=========================================================================
%% The results of TLBO
Cr_calc_TLBO=[ 0
         0
    0.0000
    0.0000
         0
    0.0004
    0.0101
    0.2148
    0.2535
    0.2832
    0.3137
    0.3509
    0.3760
    0.4136
    0.4386
    0.4695
    0.4998
    0.5583
    0.6131
    0.6786
    0.7365
    0.7858
    0.8747
    0.9388
    0.9642
    0.9783
    0.9891
    0.9940
    0.9964
    0.9979
    0.9987
    0.9994
    0.9997];
%%=========================================================================
%% The results of CTRWMT
Cr_calc_CTRWMT=[ 0.0000
    0.0000
         0
    0.0001
    0.0015
    0.0214
    0.0932
    0.3695
    0.4002
    0.4224
    0.4443
    0.4700
    0.4869
    0.5115
    0.5275
    0.5470
    0.5659
    0.6020
    0.6357
    0.6766
    0.7138
    0.7471
    0.8144
    0.8760
    0.9086
    0.9321
    0.9560
    0.9709
    0.9804
    0.9880
    0.9925
    0.9967
    0.9984]; 
%%=========================================================================
%% Calculation of statistical indices
for i=1:size(t_meas)
    df_1(i,1)=Cr_meas(i,1)-mean(Cr_meas);
     df_2(i,1)=(Cr_meas(i,1)-mean(Cr_meas))^2;
    
    df_1_TLBO(i,1)=Cr_calc_TLBO(i,1)-Cr_meas(i,1);
    df_1_CTRWMT(i,1)=Cr_calc_CTRWMT(i,1)-Cr_meas(i,1);
    
    df_2_TLBO(i,1)=(Cr_calc_TLBO(i,1)-Cr_meas(i,1))^2;
    df_2_CTRWMT(i,1)=(Cr_calc_CTRWMT(i,1)-Cr_meas(i,1))^2;
    
    df_3_TLBO(i,1)=Cr_calc_TLBO(i,1)*Cr_meas(i,1);
    df_3_CTRWMT(i,1)=Cr_calc_CTRWMT(i,1)*Cr_meas(i,1);
    
    df_4_TLBO(i,1)=Cr_calc_TLBO(i,1)-mean(Cr_calc_TLBO);
    df_4_CTRWMT(i,1)=Cr_calc_CTRWMT(i,1)-mean(Cr_calc_CTRWMT);
    
    df_5_TLBO(i,1)=df_4_TLBO(i,1)*df_1(i,1);
    df_5_CTRWMT(i,1)=df_4_CTRWMT(i,1)*df_1(i,1);
    
    df_6_TLBO(i,1)=(Cr_calc_TLBO(i,1)-mean(Cr_calc_TLBO))^2;
    df_6_CTRWMT(i,1)=(Cr_calc_CTRWMT(i,1)-mean(Cr_calc_CTRWMT))^2;
    
    df_7_TLBO(i,1)=abs(Cr_calc_TLBO(i,1)-mean(Cr_meas))+abs(df_1(i,1));
    df_7_CTRWMT(i,1)=abs(Cr_calc_CTRWMT(i,1)-mean(Cr_meas))+abs(df_1(i,1));
    
end
HH_TLBO=sqrt(sum(df_2_TLBO)/sum(df_3_TLBO));
HH_CTRWMT=sqrt(sum(df_2_CTRWMT)/sum(df_3_CTRWMT));

MBE_TLBO=(1/N_2)*sum( df_1_TLBO);
MBE_CTRWMT=(1/N_2)*sum( df_1_CTRWMT);

RMSE_TLBO=sqrt((1/N_2)*sum(df_2_TLBO));
RMSE_CTRWMT=sqrt((1/N_2)*sum(df_2_CTRWMT));

SI_TLBO=RMSE_TLBO/mean(Cr_meas);
SI_CTRWMT=RMSE_CTRWMT/mean(Cr_meas);

TS_TLBO=MBE_TLBO*sqrt((N_2-1)/RMSE_TLBO^2-MBE_TLBO^2);
TS_CTRWMT=MBE_CTRWMT*sqrt((N_2-1)/RMSE_CTRWMT^2-MBE_CTRWMT^2);

SD_TLBO=(1/mean(Cr_meas))*sqrt(N_2*sum(df_2_TLBO)-(sum(df_1_TLBO))^2)/N_2;
SD_CTRWMT=(1/mean(Cr_meas))*sqrt(N_2*sum(df_2_CTRWMT)-(sum(df_1_CTRWMT))^2)/N_2;

U95_TLBO=1.96*sqrt(RMSE_TLBO^2+SD_TLBO^2);
U95_CTRWMT=1.96*sqrt(RMSE_CTRWMT^2+SD_CTRWMT^2);

R2_TLBO=(sum(df_5_TLBO))^2/(sum(df_2)*sum(df_6_TLBO));
R2_CTRWMT=(sum(df_5_CTRWMT))^2/(sum(df_2)*sum(df_6_CTRWMT));

WIA_TLBO=1-sum( df_2_TLBO)/sum(df_7_TLBO);
WIA_CTRWMT=1-sum( df_2_CTRWMT)/sum(df_7_CTRWMT);

GPI_TLBO=MBE_TLBO*RMSE_TLBO*U95_TLBO*TS_TLBO*(1-R2_TLBO);
GPI_CTRWMT=MBE_CTRWMT*RMSE_CTRWMT*U95_CTRWMT*TS_CTRWMT*(1-R2_CTRWMT);

Criteria_TLBO=[HH_TLBO;MBE_TLBO;RMSE_TLBO;SI_TLBO;TS_TLBO;SD_TLBO;U95_TLBO;R2_TLBO;WIA_TLBO;GPI_TLBO] 
Criteria_CTRWMT=[HH_CTRWMT;MBE_CTRWMT;RMSE_CTRWMT;SI_CTRWMT;TS_CTRWMT;SD_CTRWMT;U95_CTRWMT;R2_CTRWMT;WIA_CTRWMT;GPI_CTRWMT] 

