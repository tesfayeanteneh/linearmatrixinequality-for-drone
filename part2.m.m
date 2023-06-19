function  dx = part2(t, x)
md = 0.5;%mass of drone in kg
Iz = 0.05;% moment of inertia with respect to the z axis (kg*m^2)
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
K1= [-138.7582 -134.2847   -4.0834   -0.0000   -0.0000   -0.0000
     -0.0000    0.0000    0.0000  -35.7527  -33.3610   -0.9383];

w = [2;0.2];

u = K1*x;
dx = A*x + B*u + E*w;
end