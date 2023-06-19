clear
close all
clc

%% parameters 
md = 0.5;%mass of drone in kg
Iz = 0.05;% moment of inertia with respect to the z axis (kg*m^2)
%% state space reperesentation 
A = [0	1	0	0 0 0;
     0	0	1	0 0 0;
     0	0	0	0 0 0;
     0	0	0	0 1 0;
     0	0	0	0 0 1;
     0	0	0	0 0 0];
 

B = [0     0;
     0     0;
     1/md  0;
     0     0;
     0     0;
     0	   1/Iz];
     
C = [0 1 0 0 0 0;
    0 0 0 0 1 0];
D = [0 ,0;
    0 ,0];

% -- Define disturbance 
E = B; % The disturbance 
F =[0 ,0;
    0 ,0]; % Disturbance does not affect the measure


%% Plant Dimensions
%x[n x 1]
%u[p x 1]
%y[m x 1]
%w[d x 1]

[n,~] = size(A); %[n x n]
[~,p] = size(B); %[n x p]
[m,~] = size(C); %[m x n]
[~,d] = size(F); %[m x d]


%% Synthesis to ensure External Stability (|z|<gamma*|omega|)
yalmip('clear')

% -- Set the parameters of the LMI solver
opts = sdpsettings('solver', 'mosek', 'verbose', 0);

W = sdpvar(n,n);
X = sdpvar(p,n);
rho = sdpvar(1,1);
gamma = sdpvar(1,1);

fake_zero = 1e-3;
alpha_b = 1;

gamma_opt = [];
k_opt = [];

for k = 10:10:1400
%for k=1200
    constr = [W >= rho*eye(n);
              rho >= fake_zero;
              A*W + B*X + (A*W + B*X)' <= -2*alpha_b*W;
              [k*rho*eye(n),   X';
               X,          k*rho*eye(p)] >= fake_zero*eye(n+p);
              [A*W + B*X + (A*W + B*X)',        E,         (C*W + D*X)';
                          E',            -gamma*eye(d),         F';
                      C*W + D*X,                F, -gamma*eye(m)] <= -fake_zero*eye(n+d+m)];
                  
                  
     solvesdp(constr, gamma, opts);
     
     [primal_res, dual_res] = check(constr);
     
     if( all([primal_res; dual_res] > -fake_zero) )
         gamma_opt(end+1) = value(gamma);
         k_opt(end+1) = k;
     else
        disp('Problem is infeasible'); 
     end
end
% -- Extract results
W = value(W);
X_value = value(X);
K = X_value / W;
gamma = value(gamma);

%% OPtimality curve 


gamma_des = 0.1;

k_des = getMaxK_part2(gamma_des, k_opt, gamma_opt);
area(k_opt, gamma_opt, 'FaceColor',[0.3010 0.7450 0.9330], 'EdgeColor',[0 0.4470 0.7410]);
xlabel('Bound on |K|');
ylabel('Minimum value of \gamma');
axis([0 k_opt(end) 0 gamma_opt(1)+1/4]);
hold on;
figure (1);
plot([0 k_des], [gamma_des gamma_des], 'r', 'LineWidth', 1.5);
plot([k_des k_des], [0 gamma_des], 'r', 'LineWidth', 1.5);


%% Simulation of the state and output response 

tspan = [0 15];
init_con = [0 1 0.1 0 0.5 0]';
[t, x] = ode45(@part2, tspan, init_con);
z = C*x';
x_norm = vecnorm(x');

figure(2);
plot(t,x,'LineWidth',2)
title('state vectors')
figure(3)
plot(t,x(:,2),'LineWidth',2)
title('state vectors')

figure(4);
plot(t,z,'LineWidth',2)
title('output ez and eψ')
figure(5)
plot(t,x_norm,'LineWidth',2)
title('norm of state vector')



