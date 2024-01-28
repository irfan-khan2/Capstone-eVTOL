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


%% Values of Arm Length
% Define a range of arm lengths L
% diameter of prop 0.13 m so L has to be greater than 0.13
L = 0.15;  % [m]distance of center of motor to the CoG of drone main body(user input - subject to change)

rpm = 1800*14.8*(2*pi/60); % max rpm output of each motor
g = 9.81;
m= 0.3; % mass of drone body without motors(user input - subject to change)
Fd=0; % drag =0 
kt=2.980e-6; % (user input - subject to change)
kd=1.140e-7; % kd is drag torque proportionality constant (user input - subject to change)
syms phi(t) theta(t) si(t) t 
% phi(t); % phi(t) roll %
% theta(t); % pitch
% si(t); % yaw

%% Translational Conversion Matrix 
RbG=[ cos(si(t))*cos(theta(t)), sin(si(t))*cos(theta(t)), -sin(theta(t));
    cos(si(t))*sin(phi(t))*sin(theta(t))-cos(phi(t))*sin(si(t)), sin(phi(t))*sin(si(t))*sin(theta(t))+cos(phi(t))*cos(si(t)), cos(theta(t))*sin(phi(t)); 
    cos(phi(t))*cos(si(t))*sin(theta(t))+sin(phi(t))*sin(si(t)), cos(phi(t))*sin(si(t))*sin(theta(t))-cos(si(t))*sin(phi(t)), cos(phi(t))*cos(theta(t))];
% body from global 
RGb= RbG' ;% global from body 

%% Rotational Conversion Matrix 
S=[ 1, 0, -sin(theta(t));
    0, cos(phi(t)), sin(phi(t))*cos(theta(t));
    0, -sin(phi(t)), cos(phi(t))*cos(theta(t))];
% global from body 
S_inv= S^-1 ; % body from global  

%% Translational EOM 
syms  T1(t) T2(t) T3(t) T4(t) W1(t) W2(t) W3(t) W4(t) t 
Fg=[0;0;m*g]; % airplane weight 

    
% Thrust as a function of time given motor speed
T1(t)=W1(t)^2*kt ;
T2(t)=W2(t)^2*kt ;
T3(t)=W3(t)^2*kt ;
T4(t)=W4(t)^2*kt;
F_thrust= RGb*[0;0;T1(t)+T2(t)+T3(t)+T4(t)]; % Total thrust 
% Final Equation of motion
Acc_G=(1/m)*(Fg-F_thrust - Fd); % Translation EOM (X_dot_dot; Y_dot_dot; Z_dot_dot)
    
%% Rotational EOM 
I=[1.92e-4 0 0; 
    0 1.92e-4 0;
    0 0 1.92e-4]; %Moment of inertia matrix(user input - subject to change)

% Position of each motor relative to body frame from CG
P1=[L;0;0];    
P2=[0;L;0];
P3=[-L;0;0];
P4=[0;-L;0];

I = update_inertia_tensor(I,P1,P2,P3,P4,0.0338);% 0.0338kg point mass of motors (user input - subject to change)
%Motors moment arm contribution in body axis frame
Tm1=cross(P1,[0;0;-T1(t)]);
Tm2=cross(P2,[0;0;-T2(t)]);
Tm3=cross(P3,[0;0;-T3(t)]);
Tm4=cross(P4,[0;0;-T4(t)]);
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
syms phi_derivative(t) theta_derivative(t) si_derivative(t)  
J_r=3.357e-5;
w=[phi_derivative(t) - sin(theta(t))*si_derivative(t);cos(phi(t))*theta_derivative(t) + cos(theta(t))*sin(phi(t))*si_derivative(t);cos(phi(t))*cos(theta(t))*si_derivative(t) - sin(phi(t))*theta_derivative(t)];
p=w(1,1);
q=w(2,1);
r=w(3,1);
%gyroscopic effect:
T_g=[J_r*q*(M1*W1(t)+M2*W2(t)+M3*W3(t)+M4*W4(t));
    -1*J_r*p*((M1*W1(t)+M2*W2(t))+M3*W3(t)+M4*W4(t));
    0];
% rate of angualr acceleration in terms of p,q,r

% third term:
T_g_2=cross(w,I*w);
% Final Equation of motion
w_dot=I\((Tm_total)-(T_g)-(T_g_2)); % Rotational Equation of Motion (Phi_dot_dot; Theta_dot_dot; Si_dot_dot)



%% Final form of equation conversion from symbols to numerics 

w_dot_str=string(w_dot);
w_dot_str=strrep(w_dot_str,'phi(t)','X(1)');
w_dot_str=strrep(w_dot_str,'theta(t)','X(3)');
w_dot_str=strrep(w_dot_str,'si(t)','X(5)');
w_dot_str=strrep(w_dot_str,'si_derivative(t)','X(6)');
w_dot_str=strrep(w_dot_str,'phi_derivative(t)','X(2)');
w_dot_str=strrep(w_dot_str,'theta_derivative(t)','X(4)');
w_dot_str=strrep(w_dot_str,'W1(t)',num2str(rpm*0.44));
w_dot_str=strrep(w_dot_str,'W2(t)',num2str(rpm*0.4403));
w_dot_str=strrep(w_dot_str,'W3(t)',num2str(rpm*0.44));
w_dot_str=strrep(w_dot_str,'W4(t)',num2str(rpm*0.4403));

w_dot_str;
% kt 

dvdt=@(t,X) [X(2);eval(w_dot_str(1,1));
    X(4);eval(w_dot_str(2,1));
    X(6);eval(w_dot_str(3,1))];

[t,X]=ode45(dvdt,0:0.5:5,[0.0;0;0;0;0;0]);



%% Evaluate numerical values of Acc_G using obtained X values

Acc_G = string(Acc_G);
Acc_G=strrep(Acc_G,'phi(t)','X(1)');
Acc_G=strrep(Acc_G,'theta(t)','X(3)');
Acc_G=strrep(Acc_G,'si(t)','X(5)');
Acc_G=strrep(Acc_G,'W1(t)',num2str(rpm*0.44));
Acc_G=strrep(Acc_G,'W2(t)',num2str(rpm*0.4403));
Acc_G=strrep(Acc_G,'W3(t)',num2str(rpm*0.44));
Acc_G=strrep(Acc_G,'W4(t)',num2str(rpm*0.4403));
acc_G_numeric = zeros(length(t), 3);
for i = 1:length(t)
    acc_G_numeric(i,:) = [eval(Acc_G(1,1)),eval(Acc_G(2,1)),eval(Acc_G(3,1))];
end
%% 
% Plot the results
figure;
% Adjust the plot layout as needed
hold on
plot(t, X(:,2));
plot(t, X(:,4));
plot(t, X(:,6));
title('Angular Velocity vs time');
xlabel('Time (sec)');
ylabel('Eulers Angualar Velocity (rad/s)');
legend('phi dot','theta dot','si dot');
hold off
grid on


figure;
% Adjust the plot layout as needed
hold on
plot(t, X(:,1));
plot(t, X(:,3));
plot(t, X(:,5));
title('Eulers angles vs time');
xlabel('Time (sec)');
ylabel('Eulers Angles (rad)');
legend('phi','theta','si');
hold off
grid on


%% Plot the results 

figure;
% Adjust the plot layout as needed
hold on
title('Transitional Acceleration vs time');
xlabel('Time (sec)');
ylabel('Acceleration (m/s^2)');
% Plot numerical values of Acc_G
plot(t, acc_G_numeric(:,1), '--');
plot(t, acc_G_numeric(:,2), '--');
plot(t, acc_G_numeric(:,3), '--');
legend('X Acc_G','Y Acc_G','Z Acc_G');
hold off
grid on
%% Function 
function inertia_tensor = update_inertia_tensor(inertia_tensor_, P1, P2, P3, P4, point_masses)
    % Check if the input parameters are valid
    if ~isequal(size(inertia_tensor_), [3, 3])
        error('Inertia tensor must be a 3x3 matrix');
    end
    m = point_masses;
    delta1 =  m * (-[0 -P1(3) P1(2); P1(3) 0 -P1(1); -P1(2) P1(1) 0] * [0 -P1(3) P1(2); P1(3) 0 -P1(1); -P1(2) P1(1) 0]);
    delta2 =  m * (-[0 -P2(3) P2(2); P2(3) 0 -P2(1); -P2(2) P2(1) 0] * [0 -P2(3) P2(2); P2(3) 0 -P2(1); -P2(2) P2(1) 0]);
    delta3 =  m * (-[0 -P3(3) P3(2); P3(3) 0 -P3(1); -P3(2) P3(1) 0] * [0 -P3(3) P3(2); P3(3) 0 -P3(1); -P3(2) P3(1) 0]);
    delta4 =  m * (-[0 -P4(3) P4(2); P4(3) 0 -P4(1); -P4(2) P4(1) 0] * [0 -P4(3) P4(2); P4(3) 0 -P4(1); -P4(2) P4(1) 0]);
    inertia_tensor = inertia_tensor_ + delta1 + delta2 + delta3 + delta4;
end



