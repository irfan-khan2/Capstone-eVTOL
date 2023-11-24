close all;
clear all;


%% Define else where
syms phi(t) theta(t) si(t) t 
phi(t); % phi(t) roll %
theta(t); % pitch
si(t); % yaw

X1 = [];
X2 = [];
X3 = [];
X4 = [];
X5 = [];
X6 = [];


%%
%%Translational Conversion Matrix 
RbG=[ cos(si(t))*cos(theta(t)), sin(si(t))*cos(theta(t)), -sin(theta(t)); 
      cos(si(t))*sin(phi(t))*sin(theta(t))-cos(si(t))*sin(si(t)), sin(phi(t))*sin(si(t))*sin(theta(t))+cos(phi(t))*cos(si(t)), cos(theta(t))*sin(phi(t)); 
      cos(phi(t))*cos(si(t))*sin(theta(t))+sin(phi(t))*sin(si(t)), cos(phi(t))*sin(si(t))*sin(theta(t))-cos(si(t))*sin(phi(t)), cos(phi(t))*cos(theta(t))];
% body from global 

RGb= RbG' ;% global from body 

% Rotational Conversion Matrix 
%%
% Rotational Conversion Matrix 
S=[ 1, 0, -sin(theta(t)); 
    0, cos(phi(t)), sin(phi(t))*cos(theta(t));
    0, -sin(phi(t)), cos(phi(t))*cos(theta(t))];
    % global from body 

S_inv= S^-1 ; % body from global 
%% Translational EOM 
syms m g T1(t) T2(t) T3(t) T4(t) W1(t) W2(t) W3(t) W4(t) t
m; % mass in kg  

Fg=[0;0;m*g]; % airplane weight 
Fd=0; % drag =0 

% 4 speeds for each motor
W1(t) ;
W2(t)  ;
W3(t)  ;
W4(t)  ;

% Thrust as a function of time given motor speed
T1(t) ;
T2(t) ;
T3(t) ;
T4(t) ;

F_thrust= RGb*[0;0;T1(t)+T2(t)+T3(t)+T4(t)]; % Total thrust 

Acc_G=(1/m)*(Fg-F_thrust - Fd); % Translation EOM 


%% Rotational EOM 

syms  Ixx Ixy Ixz Iyx Iyy Iyz Izx Izy Izz 
syms P1x P1y P1z P2x P2y P2z P3x P3y P3z P4x P4y P4z

I=[1 8 0; 
    0 1 0; 
    6 0 1]; %Moment of inertia 

% I=[Ixx Ixy Ixz; 
%     Ixy Iyy Iyz; 
%     Ixz Iyz Izz]; %Moment of inertia 

% Position of each motor relative to body frame from CG
P1=[1;1;0];
P2=[-1;1;0];
P3= [-1;-1;0];
P4=[-1;1;0];


% P1=[P1x;P1y;P1z];
% P2=[P2x;P2y;P2z];
% P3= [P3x;P3y;P3z];
% P4=[P4x;P4y;P4z];

%Motors moment arm contribution in body axis frame
Tm1=cross(P1,[0;0;-T1(t)]);
Tm2=cross(P2,[0;0;-T2(t)]);
Tm3=cross(P3,[0;0;-T3(t)]);
Tm4=cross(P4,[0;0;-T4(t)]);

Tm_pitch_roll=Tm1+Tm2+Tm3+Tm4;

%kd is drag torque proportionality constant
kd=1;
% rotation direction of each motor. +ve=cw, -ve=ccw

syms M1 M2 M3 M4;

M1=1;
M2=-1;
M3=1;
M4=-1;

Tm_yaw=[0;0;kd*( M1*(W1(t))^2 + M2*(W2(t))^2 + M3*(W3(t))^2 + M4*(W4(t))^2)];

Tm_total= Tm_pitch_roll + Tm_yaw ;

%gyroscopic effect: 


% torque moment of inertia 
syms phi_derivative(t) theta_derivative(t) si_derivative(t) J_r

J_r=1;

T_g=[J_r*theta_derivative(t)*(M1*W1(t)+M2*W2(t))+M3*W3(t)+M4*W4(t); 
    -1*J_r*phi_derivative(t)*(M1*W1(t)+M2*W2(t))+M3*W3(t)+M4*W4(t);
    0];
T_g=0;
w=[phi_derivative(t);theta_derivative(t);si_derivative(t)];

% third term:
T_g_2=cross(w,I*w);

w_dot=I\((Tm_total)-(T_g)-(T_g_2));

w_dot_str=string(w_dot);
w_dot_str=strrep(w_dot_str,'theta_derivative(t)','X2');
w_dot_str=strrep(w_dot_str,'phi_derivative(t)','X4');
w_dot_str=strrep(w_dot_str,'si_derivative(t)','X6');
w_dot_str=strrep(w_dot_str,'si_derivative(t)','X6');
w_dot_str=strrep(w_dot_str,'si_derivative(t)','X6');
w_dot_str=strrep(w_dot_str,'si_derivative(t)','X6');
w_dot_str=strrep(w_dot_str,'W1(t)','t');
w_dot_str=strrep(w_dot_str,'W2(t)','t');
w_dot_str=strrep(w_dot_str,'W3(t)','t');
w_dot_str=strrep(w_dot_str,'W4(t)','t');

w_dot_str=strrep(w_dot_str,'T1(t)','t');
w_dot_str=strrep(w_dot_str,'T2(t)','t');
w_dot_str=strrep(w_dot_str,'T3(t)','t');
w_dot_str=strrep(w_dot_str,'T4(t)','t');


w_dot_str;
w_
%%
% fh1 = eval(['@(t,X2,X4,X6)',w_dot_str(1,1)]);
% fh3 = eval(w_dot_str(2,1));
% fh5 = eval(w_dot_str(3,1));
% 
% 
% % dvdt=@(t,X) [X(2);(2*t)/3 - t/2 - (2*t)/3 - t/2 - t^2/3 + t^2/3 - t^2/3 + t^2/3 + (X(4)*X(6))/12 + (9*X(4)*X(2))/4 - (9*X(6)*X(2))/4 + (13*X(4)^2)/6 - (11*X(6)^2)/6 - (5*X(2)^2)/12;
% %     X(4);t/2 - t/3 + t/3 + t/2 + (2*t^2)/3 - (2*t^2)/3 + (2*t^2)/3 - (2*t^2)/3 - (11*X(4)*X(6))/12 - (7*X(4)*X(2))/4 + (11*X(6)*X(2))/4 - (17*X(4)^2)/6 + (7*X(6)^2)/6 + (19*X(2)^2)/12;
% %     X(6);t/3 - t/2 - t/3 - t/2 - t^2/3 + t^2/3 - t^2/3 + t^2/3 + (19*X(4)*X(6))/12 - (X(4)*X(2))/4 - (3*X(6)*X(2))/4 + (7*X(4)^2)/6 + X(6)^2/6 - (23*X(2)^2)/12]
% 
% dwdt=@(t,X) [fh1;X2;fh3;X4;fh5;X6];
% 
% [X1,X2,X3,X4,X5,X6]=ode45(dwdt,[0 10],[0;0;0;0;0;0])



% Define symbolic variables and their derivatives
syms X2(t) X4(t) X6(t) X1(t) X3(t) X5(t)

% Convert w_dot_str to a function
w_dot_func = matlabFunction(w_dot_str, 'Vars', [t, X2, X4, X6, X1, X3, X5]);

% Define a time span
tspan = [0, 10]; % Adjust the time span as needed

% Define initial conditions for X2, X4, and X6

initial_conditions = [0,0,0,0,0,0];

% Solve the ODE using ode45
[t_sol, X_sol] = ode45(@(t, X) w_dot_func(t, X(1), X(2), X(3), X(4), X(5), X(6)), tspan, initial_conditions);

% Extract the solutions for X2, X4, and X6
X2_sol = X_sol(:, 1)

