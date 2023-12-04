close all 
clear all 

%% Coordinate System 
%                            ^(+x)(+ve Roll)
%                            |
%                            |    
%                            | (+z)(into page)(+ve Yaw)
%(-ve Theta)(-y)<-----------(X)-------------->(+y)(+ve Theta) 
%                            |
%                            |
%                            |
%                            v (-x)(-ve Roll)

%% Declare all the symbols
syms phi(t) theta(t) si(t) t W1(t) W2(t) W3(t) W4(t) m g kt kd L


%% Translational Conversion Matrix 

RbG=[ cos(si(t))*cos(theta(t)), sin(si(t))*cos(theta(t)), -sin(theta(t)); 
          cos(si(t))*sin(phi(t))*sin(theta(t))-cos(phi(t))*sin(si(t)), sin(phi(t))*sin(si(t))*sin(theta(t))+cos(phi(t))*cos(si(t)), cos(theta(t))*sin(phi(t)); 
          cos(phi(t))*cos(si(t))*sin(theta(t))+sin(phi(t))*sin(si(t)), cos(phi(t))*sin(si(t))*sin(theta(t))-cos(si(t))*sin(phi(t)), cos(phi(t))*cos(theta(t))];
% body from global
RGb= RbG'; % global from body
%% Rotational Conversion Matrix 
S=[ 1, 0, -sin(theta(t)); 
    0, cos(phi(t)), sin(phi(t))*cos(theta(t));
    0, -sin(phi(t)), cos(phi(t))*cos(theta(t))];
% global from body 
S_inv= S^-1 ; % body from global
%% Translational EOM 
Fg=[0;0;m*g];% airplane weight 
Fd=0;% drag =0
% Thrust as a function of time given motor speed
T1=W1(t)^2*kt;
T2=W2(t)^2*kt;
T3=W3(t)^2*kt;
T4=W4(t)^2*kt;
F_thrust= RGb*[0;0;T1+T2+T3+T4]; % Total thrust
% Final Equation of motion
Acc_G=(1/m)*(Fg-F_thrust - Fd)
%% Rotational EOM 
syms Ixx Iyy Izz
I=[Ixx 0 0;0 Iyy 0;0 0 Izz];% moment of inertia tensor since drone is symmetrical, also Ixx = Iyy (square)
% Position of each motor relative to body frame from CG, Pi
P1=[L;0;0];    
P2=[0;L;0];
P3=[-L;0;0];
P4=[0;-L;0];
%Motors moment arm contribution in body axis frame
Tm1=cross(P1,[0;0;-T1]);
Tm2=cross(P2,[0;0;-T2]);
Tm3=cross(P3,[0;0;-T3]);
Tm4=cross(P4,[0;0;-T4]);
Tm_pitch_roll=Tm1+Tm2+Tm3+Tm4;
% rotation direction of each motor. +ve=cw, -ve=ccw
% Mi corresponds to Pi above
M1=1;
M2=-1;
M3=1;
M4=-1;
Tm_yaw=[0;0;kd*( M1*(W1(t))^2 + M2*(W2(t))^2 + M3*(W3(t))^2 + M4*(W4(t))^2)];
Tm_total= Tm_pitch_roll + Tm_yaw ;
% torque moment of inertia 
syms phi_derivative(t) theta_derivative(t) si_derivative(t) J_r p q r
w = [p; q; r]; 
%gyroscopic effect:
T_g=[J_r*q*(M1*W1(t)+M2*W2(t)+M3*W3(t)+M4*W4(t));
    -1*J_r*p*((M1*W1(t)+M2*W2(t))+M3*W3(t)+M4*W4(t));
    0];% rate of angualr acceleration in terms of p,q,r
% third term:
T_g_2=cross(w,I*w);
% Final Equation of motion
w_dot=I\((Tm_total)-(T_g)-(T_g_2)) % Rotational Equation of Motion (Phi_dot_dot; Theta_dot_dot; Si_dot_dot)