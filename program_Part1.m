clear
close all
clc
%% parameters 
md = 0.5;%mass of drone in kg
Iz = 0.05;% moment of inertia with respect to the z axis (kg*m^2)
%% state space reperesentation 
A = [0	1	0	0;
     0	0	0	0;
     0	0	0	1;
     0	0	0	0];
 
B = [0     0;
     1/md  0;
     0     0;
     0	   1/Iz];
     
C = [1 0 0 0;
     0 0 1 0];


%Dimensions of the system
[n, ~] = size(A);
[~, p] = size(B);
%% Check stability of linearized system

eig_Val_A = eig(A);

disp('Real part of the open loop eigenvalues of linearized system:');
disp(mat2str(real(eig_Val_A')));
%% Check controllability of the system

rank_R = rank(ctrb(A,B));
if( rank_R == n )
    disp('Controllability matrix is full rank and a system is controllable')
else
    disp('System is not controllable')
end

%% Stabilizer imposing convergence rate and minimizing actuator effort
yalmip('clear');

opts = sdpsettings('solver', 'mosek', 'verbose', 0);

W = sdpvar(n, n);
X = sdpvar(p, n);
fake_zero = 1e-3;
alpha_b =1;

k = sdpvar(1,1);

constr = [W >= eye(n);
         A*W + B*X + (A*W + B*X)' <= -2*alpha_b*W;
       [k*eye(n),     X';
              X,      k*eye(p)] >= fake_zero*eye(n+p)];
      
solvesdp(constr, k, opts);

W_c = value(W);
X_c = value(X);
P_c = inv(W_c);
k_c = value(k);
K_c = X_c/W_c;

[primal_res, dual_res] = check(constr);

if( all([primal_res; dual_res] > -fake_zero) )
    disp('-------------');
    disp('Controller imposing convergence rate and minimizing actuator effort');
    disp('Problem is feasible --> linearized system is asymptotically stable');
    disp('Gain matrix K');
    disp(mat2str(K_c));
    disp('Certified actuator bound k');
    disp(mat2str(k_c));
    disp('Real part of closed-loop eigenvalues ');
    disp(mat2str(real(eig(A+B*K_c)')));
else
    disp('Problem is infeasible');
end
                