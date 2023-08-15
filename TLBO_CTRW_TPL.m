%% Descriptions:
%% This code otptimizes the CTRW-TPL parameters using TLBO algorithm.
%% Notice:This code is linked to the CTRW_TPL_fit.  
tic
clc
clear all
format long
%% Introducing the number and the range of variables
nVar=4;
x_min_1=0.00001;
x_max_1=0.01;
x_min_2=10^-5;
x_max_2=0.5;
x_min_3=0.8;
x_max_3=3;
x_min_4=50;
x_max_4=10^10;
%% Initialization  
nPop=15;                   %% Number of students 
nIt=300;                   %% Number of iterations
for i=1:nPop
    x(1,i)=x_min_1+rand(1,1)*(x_max_1-x_min_1); 
    x(2,i)=x_min_2+rand(1,1)*(x_max_2-x_min_2);   
    x(3,i)=x_min_3+rand(1,1)*(x_max_3-x_min_3); 
    x(4,i)=x_min_4+rand(1,1)*(x_max_4-x_min_4); 
    z(1,i)=CTRW_TPL_fit(x(1,i),x(2,i),x(3,i),x(4,i));
end
%% Producing the next generations 
Z=[];
I=[];
for It=1:nIt
%% Sorting the objective function
    for i=1:nPop       
    z(1,i)=CTRW_TPL_fit(x(1,i),x(2,i),x(3,i),x(4,i));
    k=2;
z_min_1=z(1,1);
x_min_1_1=x(1,1);
x_min_1_2=x(2,1);
x_min_1_3=x(3,1);
x_min_1_4=x(4,1);
    for j=k:nPop
        if z_min_1>z(1,j)
           a_z_1=z(1,j);
           z(1,j)=z_min_1;
           z(1,j-1)=a_z_1;
           a_x_1_1=x(1,j);
           x(1,j)=x_min_1_1;
           x(1,j-1)=a_x_1_1;
           a_x_1_2=x(2,j);
           x(2,j)=x_min_1_2;
           x(2,j-1)=a_x_1_2;
           a_x_1_3=x(3,j);
           x(3,j)=x_min_1_3;
           x(3,j-1)=a_x_1_3;
           a_x_1_4=x(4,j);
           x(4,j)=x_min_1_4;
           x(4,j-1)=a_x_1_4;
                   else
            z_min_1=z(1,j);
            x_min_1_1=x(1,j);
            x_min_1_2=x(2,j);
            x_min_1_3=x(3,j);
            x_min_1_4=x(4,j);
                    end
             k=k+1;
    end 
%% Determining the teacher and calculating the mean of each design variable
for j=1:nVar
    T(j,1)=x(j,1);
    M(j,1)=mean(x(j,:));
end
%% Calculationg the difference between the teacher and the mean
T_F=round(1+rand(1,1));
for j=1:nVar
Diff(j,1)=rand(1,1)*(T(j,1)-T_F*M(j,1));
end
x_new(1,i)=abs(x(1,i)+Diff(1,1));
x_new(2,i)=abs(x(2,i)+Diff(2,1));
x_new(3,i)=abs(x(3,i)+Diff(3,1));
x_new(4,i)=abs(x(4,i)+Diff(4,1));
z_new(1,i)=CTRW_TPL_fit(x_new(1,i),x_new(2,i),x_new(3,i),x_new(4,i));
if z_new(1,i)<z(1,i)
    x(1,i)=x_new(1,i);
    x(2,i)=x_new(2,i);
    x(3,i)=x_new(3,i); 
    x(4,i)=x_new(4,i); 
end 
    end
%% Learner phase
    for i=1:nPop
        zz(1,i)=CTRW_TPL_fit(x(1,i),x(2,i),x(3,i),x(4,i));        
    end
    for i=1:nPop
        AA=1000;
        BB=1000;
        while AA==BB
    AA=randi([1,nPop],1);
    BB=randi([1,nPop],1);
    if AA~=BB
       break 
    end
        end
        if zz(1,AA)<zz(1,BB)
       x_new_new(1,AA)=abs(x(1,AA)+rand(1,1)*(x(1,AA)-x(1,BB)));
       x_new_new(2,AA)=abs(x(2,AA)+rand(1,1)*(x(2,AA)-x(2,BB)));
       x_new_new(3,AA)=abs(x(3,AA)+rand(1,1)*(x(3,AA)-x(3,BB)));
       x_new_new(4,AA)=abs(x(4,AA)+rand(1,1)*(x(4,AA)-x(4,BB)));
       z_new_new(1,AA)=CTRW_TPL_fit(x_new_new(1,AA),x_new_new(2,AA),x_new_new(3,AA),x_new_new(4,AA));
       if z_new_new(1,AA)<zz(1,AA)
          x(1,AA)=x_new_new(1,AA);
          x(2,AA)=x_new_new(2,AA);
          x(3,AA)=x_new_new(3,AA);
          x(4,AA)=x_new_new(4,AA);
       end
    end
        if zz(1,BB)<zz(1,AA)
       x_new_new(1,AA)=abs(x(1,AA)+rand(1,1)*(x(1,BB)-x(1,AA)));
       x_new_new(2,AA)=abs(x(2,AA)+rand(1,1)*(x(2,BB)-x(2,AA)));
       x_new_new(3,AA)=abs(x(3,AA)+rand(1,1)*(x(3,BB)-x(3,AA)));
       x_new_new(4,AA)=abs(x(4,AA)+rand(1,1)*(x(4,BB)-x(4,AA)));
       z_new_new(1,AA)=CTRW_TPL_fit(x_new_new(1,AA),x_new_new(2,AA),x_new_new(3,AA),x_new_new(4,AA));
        end
        if z_new_new(1,AA)<zz(1,AA)
           x(1,AA)=x_new_new(1,AA);
           x(2,AA)=x_new_new(2,AA);
           x(3,AA)=x_new_new(3,AA);
           x(4,AA)=x_new_new(4,AA);
       end
    end
     I=[I;It];
    Z=[Z;z(1,1)];
fprintf('Iteration=%d\n',It);
fprintf('OF=%1.15f\n',z(1,1));
fprintf('v_psi=%1.15f\n',x(1,1));
fprintf('D_psi=%1.15f\n',x(2,1));
fprintf('beta=%1.15f\n',x(3,1));
fprintf('t2=%1.15f\n',x(4,1));
end
plot(I,Z)
toc




