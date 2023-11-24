
close all 
clear all 

t = 0;
%% Define else where
phi=@(t) t; % roll %
theta=@(t) t; % pitch
si=@(t) t; % yaw

%% Translational Conversion Matrix 
RbG=[ cos(si(t))*cos(theta(t)), sin(si(t))*cos(theta(t)), -sin(theta(t)); 
      cos(si(t))*sin(phi(t))*sin(theta(t))-cos(si(t))*sin(si(t)), sin(phi(t))*sin(si(t))*sin(theta(t))+cos(phi(t))*cos(si(t)), cos(theta(t))*sin(phi(t)); 
      cos(phi(t))*cos(si(t))*sin(theta(t))+sin(phi(t))*sin(si(t)), cos(phi(t))*sin(si(t))*sin(theta(t))-cos(si(t))*sin(phi(t)), cos(phi(t))*cos(theta(t))];
% body from global 

RGb= RbG' ;% global from body 

% Rotational Conversion Matrix 
%% Rotational Conversion Matrix 
S=[ 1, 0, -sin(theta(t)); 
    0, cos(phi(t)), sin(phi(t))*cos(theta(t));
    0, -sin(phi(t)), cos(phi(t))*cos(theta(t))];
    % global from body 

S_inv= S^-1 ; % body from global 
%% Translational EOM 
m= 1; % mass in kg 

Fg=[0;0;m*9.81]; % airplane weight 
Fd=0; % drag =0 

% 4 speeds for each motor
W1=@(t) t*1 ;
W2=@(t) t*1 ;
W3=@(t) t*1 ;
W4=@(t) t*1 ;

% Thrust as a function of time given motor speed
T1=@(t) W1(t)*1 ;
T2=@(t) W2(t)*1 ;
T3=@(t) W3(t)*1 ;
T4=@(t) W4(t)*1 ;

F_thrust= RGb*[0;0;T1(t)+T2(t)+T3(t)+T4(t)]; % Total thrust 

Acc_G=(1/m)*(Fg-F_thrust - Fd); % Translation EOM 


%% Rotational EOM 

I=[1 1 1; 
    1 1 1; 
    1 1 1]; %Moment of inertia 

% Position of each motor relative to body frame from CG
P1=[1;1;0];
P2=[1;-1;0];
P3=[-1;-1;0];
P4=[-1;1;0];

%Motors moment arm contribution in body axis frame
Tm1=cross(P1,[0;0;-T1(t)]);
Tm2=cross(P2,[0;0;-T2(t)]);
Tm3=cross(P3,[0;0;-T3(t)]);
Tm4=cross(P4,[0;0;-T4(t)]);

Tm_pitch_roll=Tm1+Tm2+Tm3+Tm4;

%kd is drag torque proportionality constant
kd=1;
% rotation direction of each motor. +ve=cw, -ve=ccw
M1=1;
M2=-1;
M3=1;
M4=-1;

Tm_yaw=[0;0;kd*( M1*(W1(t))^2 + M2*(W2(t))^2 + M3*(W3(t))^2 + M4*(W4(t))^2)];

Tm_total= Tm_pitch_roll + Tm_yaw ;

%% gyroscopic effect: 

phi_dot = @(t) t;
theta_dot = @(t) t;
si_dot = @(t) t;

%rotor moment of inertia

Jr = 0;
T_gyro =  [Jr*theta_dot*(W1-W2+W3-W4); -Jr*phi_dot*(W1-W2+W3-W4); 0];

