clear all,close all,clc;
%%
% this is for SDPT3
cd('C:\Users\Artun\Desktop\BACK UP 2021-june 20\Documents 2020\MATLAB toolboxes for optimization\SDPT3-4.0');
run('Installmex.m')
run('startup.m')
% cd('C:\Users\Artun\Desktop\BACK UP 2021-june 20\Documents 2020\MATLAB toolboxes for optimization\SDPT3-4.0\Solver')
% run('sqlpdemo.m')

% this is for YALMIP
addpath(genpath('C:\Users\Artun\Desktop\BACK UP 2021-june 20\Documents 2020\MATLAB toolboxes for optimization\YALMIP-master'))
yalmip('clear');

cd('C:\Users\Artun\Desktop\ele 301 Lab 2021 summer\lec 04');
%% static full state feedback
A=[1,2;1,0];
B=[1;0];
eps1=1e-5;
nx=size(A,1);
nu=size(B,2);
P = sdpvar(nx,nx);
Z = sdpvar(nu,nx);
Constraints = [P>=eps1*eye(nx)];
Constraints =[Constraints ; A*P+[A*P]'+[B*Z]+[B*Z]'<=-eps1*eye(nx)];

ops = sdpsettings('solver','sdpt3');
sol = optimize(Constraints,[],ops);

if (sol.problem==1) % if there is a "problem" [i.e. infeasibility]
    disp('infeasible');
else % there is not a "problem" [i.e. infeasibility]
    disp('feasible');
end
P0=value(P);
Z0=value(Z);
if min(real(eig(P0)))>0
    disp('POS-DEF-MAT');
else
    disp('NOOOOOT POS-DEF-MAT----> BAD!!!!');
end

temp1=A*P0+[A*P0]'+[B*Z0]+[B*Z0]';
if max(real(eig(temp1)))<0
    disp('NEG-DEF-MAT');
else
    disp('NOOOOOT NEG-DEF-MAT----> BAD!!!!');
end

% K=Z0*inv(P0)
K=Z0/P0

eig(A+B*K)

%% example-1 [KEMIN ZHOU & J.DOYLE] "essentials of robust ctrl" [example 4.2]
k1=1;
k2=4;
b1=0.2;
b2=0.1;
m1=1;
m2=2;

A=[zeros(2,2),eye(2);
   -k1/m1,k1/m1,-b1/m1,b1/m1;
   k1/m2,-(k1+k2)/m2,b1/m2,-(b1+b2)/m2];

B=[zeros(2,2);
   diag([1/m1,1/m2])]

C=[eye(2),zeros(2,2)];

D=C*B*0; % zeros(size(C*B)) m,n,p --> dim of D: "p x m"
% hinfnorm(ss(A,B,C,D)) = 11.47
gamma=12


eps1=1e-5;
P = sdpvar(4,4);
Constraints = [P>=eps1*eye(4)];
Constraints =[Constraints ; [P*A+A'*P+C'*C,P*B+C'*D;B'*P+D'*C,D'*D-gamma^2*eye(2)]<=-eps1*eye(6)];

ops = sdpsettings('solver','sdpt3');
sol = optimize(Constraints,[],ops);
% optimize(Constraints, Objective,ops)

if (sol.problem==1) % if there is a "problem" [i.e. infeasibility]
    disp('infeasible');
else % there is not a "problem" [i.e. infeasibility]
    disp('feasible');
end
P0 = value(P);
%% example-2 bisection algo. using for loop

gamma_max=100;
gamma_min=0.1;

while(true)
    if [(gamma_max-gamma_min)/gamma_min]<1e-6
        gamma_test=(gamma_max+gamma_min)/2;
        gamma_test
        break;
    end
    
    gamma_test=(gamma_max+gamma_min)/2;
    
        yalmip('clear');
        eps1=1e-6;
        P = sdpvar(4,4);
        Constraints = [P>=eps1*eye(4)];
        Constraints =[Constraints ; [P*A+A'*P+C'*C,P*B+C'*D;B'*P+D'*C,D'*D-gamma_test^2*eye(2)]<=-eps1*eye(6)];

        ops = sdpsettings('solver','sdpt3');
        sol = optimize(Constraints,[],ops);

    if (sol.problem==1) % if there is a "problem" [i.e. infeasibility]
        disp('infeasible');
        % infeasible means that gamma_Real>gamma_test
        % gamma_test is a valid "lower-bound"
        gamma_min=gamma_test;
    else % there is not a "problem" [i.e. infeasibility]
        disp('feasible');
        % feasible means that gamma_Real<gamma_test
        % gamma_test is a valid "upper-bound"
        gamma_max=gamma_test;
    end
    
    
end







%% example-3 bisection algo. using "bisection command"

yalmip('clear');
eps1=1e-6;
P = sdpvar(4,4);
sdpvar gamma_test
Constraints = [P>=eps1*eye(4)];
Constraints =[Constraints ; ...
    [P*A+A'*P+C'*C,P*B+C'*D;B'*P+D'*C,D'*D-gamma_test^2*eye(2)]<=-eps1*eye(6)];

Objective = gamma_test;

ops=sdpsettings('solver','bisection','bisection.solver','sdpt3');
diagnostics = optimize(Constraints, Objective,ops)


value(gamma_test)

%%

























%