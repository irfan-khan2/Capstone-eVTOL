
function [time,angles,size, accel] =Anglecalc(m,I,P,M,m_r)

% %% Define else where
% syms phi(t) theta(t) si(t) t 
% phi(t); % phi(t) roll %
% theta(t); % pitch
% si(t); % yaw
% %%
% %%Translational Conversion Matrix 
% RbG=[ cos(si(t))*cos(theta(t)), sin(si(t))*cos(theta(t)), -sin(theta(t)); 
%       cos(si(t))*sin(phi(t))*sin(theta(t))-cos(si(t))*sin(si(t)), sin(phi(t))*sin(si(t))*sin(theta(t))+cos(phi(t))*cos(si(t)), cos(theta(t))*sin(phi(t)); 
%       cos(phi(t))*cos(si(t))*sin(theta(t))+sin(phi(t))*sin(si(t)), cos(phi(t))*sin(si(t))*sin(theta(t))-cos(si(t))*sin(phi(t)), cos(phi(t))*cos(theta(t))];
% % body from global 
% 
% RGb= RbG' ;% global from body 
% 
% %%
% % Rotational Conversion Matrix 
% S=[ 1, 0, -sin(theta(t)); 
%     0, cos(phi(t)), sin(phi(t))*cos(theta(t));
%     0, -sin(phi(t)), cos(phi(t))*cos(theta(t))];
%     % global from body 
% 
% S_inv= S^-1 ; % body from global 
% %% Translational EOM 
% syms m g T1(t) T2(t) T3(t) T4(t) W1(t) W2(t) W3(t) W4(t) t
% m; % mass in kg  
% m=0.400;
% g=9.81;
% Fg=[0;0;m*g]; % airplane weight 
% Fd=0; % drag =0 
% % 4 speeds for each motor
% W1(t) ;
% W2(t)  ;
% W3(t)  ;
% W4(t)  ;
% % Thrust as a function of time given motor speed
% T1(t) ;
% T2(t) ;
% T3(t) ;
% T4(t) ;
% F_thrust= RGb*[0;0;T1(t)+T2(t)+T3(t)+T4(t)]; % Total thrust 
% Acc_G=(1/m)*(Fg-F_thrust - Fd); % Translation EOM 
% %% Rotational EOM 
% syms  Ixx Ixy Ixz Iyx Iyy Iyz Izx Izy Izz 
% syms P1x P1y P1z P2x P2y P2z P3x P3y P3z P4x P4y P4z
% I=[0.1 0.0 0.0; 
%    0.0 0.1 0.0; 
%    0.0 0.0 0.1]; %Moment of inertia 
% 
% % I=[Ixx Ixy Ixz; 
% %     Ixy Iyy Iyz; 
% %     Ixz Iyz Izz]; %Moment of inertia 
% % Position of each motor relative to body frame from CG
% P1=[0.1;0.1;0];
% P2=[0.1;-0.1;0];
% P3= [-0.1;-0.1;0];
% P4=[-0.1;0.1;0];
% % P1=[P1x;P1y;P1z];
% % P2=[P2x;P2y;P2z];
% % P3= [P3x;P3y;P3z];
% % P4=[P4x;P4y;P4z];
% %Motors moment arm contribution in body axis frame
% Tm1=cross(P1,[0;0;-T1(t)]);
% Tm2=cross(P2,[0;0;-T2(t)]);
% Tm3=cross(P3,[0;0;-T3(t)]);
% Tm4=cross(P4,[0;0;-T4(t)]);
% Tm_pitch_roll=Tm1+Tm2+Tm3+Tm4;
% 
% %kd is drag torque proportionality constant
% kd=1;
% % rotation direction of each motor. +ve=cw, -ve=ccw
% 
% syms M1 M2 M3 M4
% 
% M1=1;
% M2=-1;
% M3=1;
% M4=-1;
% 
% Tm_yaw=[0;0;kd*( M1*(W1(t))^2 + M2*(W2(t))^2 + M3*(W3(t))^2 + M4*(W4(t))^2)];
% 
% Tm_total= Tm_pitch_roll + Tm_yaw ;
% 
% %gyroscopic effect: 
% 
% 
%  % torque moment of inertia 
% syms phi_derivative(t) theta_derivative(t) si_derivative(t) J_r
% 
% J_r=0.1;
% 
% T_g=[J_r*theta_derivative(t)*(M1*W1(t)+M2*W2(t))+M3*W3(t)+M4*W4(t); 
%     -1*J_r*phi_derivative(t)*(M1*W1(t)+M2*W2(t))+M3*W3(t)+M4*W4(t);
%     0];
% 
% w=[phi_derivative(t);theta_derivative(t);si_derivative(t)];
% 
% % third term:
% T_g_2=cross(w,I*w);
% 
% w_dot=I\((Tm_total)-(T_g)-(T_g_2));
% 
% w_dot_str=string(w_dot);
% w_dot_str=strrep(w_dot_str,'theta_derivative(t)','X(2)');
% w_dot_str=strrep(w_dot_str,'phi_derivative(t)','X(4)');
% w_dot_str=strrep(w_dot_str,'si_derivative(t)','X(6)');
% w_dot_str=strrep(w_dot_str,'W1(t)','800*t');
% w_dot_str=strrep(w_dot_str,'W2(t)','700*t');
% w_dot_str=strrep(w_dot_str,'W3(t)','800*t');
% w_dot_str=strrep(w_dot_str,'W4(t)','700*t');
% 
% w_dot_str=strrep(w_dot_str,'T1(t)','0.5*t');
% w_dot_str=strrep(w_dot_str,'T2(t)','0.3*t');
% w_dot_str=strrep(w_dot_str,'T3(t)','0.5*t');
% w_dot_str=strrep(w_dot_str,'T4(t)','0.3*t');
% 
% 
% w_dot_str;
% 
% 
% 
% 
%  dvdt=@(t,X) [X(2);0.3*t - 0.5*t + 0.5*t - 0.3*t - 10*800*t + 10*700*t - (800*t*X(2))/100 + (700*t*X(2))/100;
%      X(4);0.5*t + 0.3*t - 0.5*t - 0.3*t - 10*800*t + 10*700*t + (800*t*X(4))/100 - (700*t*X(4))/100;
%      X(6);10*800^2 - 10*700^2 + 10*800^2 - 10*700^2];
% 
% % dvdt=@(t,X) [X(2);fh;X(4);fh2;X(6);fh3]
% 
% [t,X]=ode45(dvdt,[0:0.1:1],[0.0;0;0;0;0;0]);

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

%% Define else where
    syms phi(t) theta(t) si(t) t 
    phi(t); % phi(t) roll %
    theta(t); % pitch
    si(t); % yaw

    rpm = 1800*14.8*(2*pi/60); % max rpm output of each motor
    g = 9.81;
    Fd=0; % drag =0 
    kt=2.980e-6; % (user input - subject to change)
    kd=1.140e-7; % kd is drag torque proportionality constant (user input - subject to change)
    
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
%     I=[0.224 0 0; 
%         0 0.224 0; 
%         0 0 0.224]; %Moment of inertia matrix 
% % Position of each motor relative to body frame from CG
%     P1=[0.3;0;0];
%     P2=[0;0.3;0];
%     P3= [-0.3;0;0];
%     P4=[0;-0.3;0];
    P1 = P(:,1);
    P2 = P(:,2);
    P3 = P(:,3);
    P4 = P(:,4);

    I = update_inertia_tensor(I,P1,P2,P3,P4,m_r);
%Motors moment arm contribution in body axis frame
    Tm1=cross(P1,[0;0;-T1(t)]);
    Tm2=cross(P2,[0;0;-T2(t)]);
    Tm3=cross(P3,[0;0;-T3(t)]);
    Tm4=cross(P4,[0;0;-T4(t)]);
    Tm_pitch_roll=Tm1+Tm2+Tm3+Tm4;
 
% rotation direction of each motor. +ve=cw, -ve=ccw
% Mi corresponds to Pi above
    M1=M(1);
    M2=M(2);
    M3=M(3);
    M4=M(4);
    Tm_yaw=[0;0;kd*( M1*(W1(t))^2 + M2*(W2(t))^2 + M3*(W3(t))^2 + M4*(W4(t))^2)];
    Tm_total= Tm_pitch_roll + Tm_yaw ;
% torque moment of inertia 
syms phi_derivative(t) theta_derivative(t) si_derivative(t) 
    J_r=0.03;
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

%% User inputs
%     m=1; %mass in kg
%     m_r = m_r; % point mass of motor
%     I=I; % Drone moment of inertia 
%     kd=kd; % Drag torque proprotionality constant
%     kt=kt; % Motor thrust coefficient 
%     M1;M2;M3;M4; %Motros rotaions 
%     P1; P2; P3; P4; % Motors positions relative to CG
    
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

% dvdt=@(t,X) [X(2);fh;X(4);fh2;X(6);fh3]

[t,X]=ode45(dvdt,[0:0.5:5],[0.0;0;0;0;0;0])

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
time=t;
angles=X;
size=height(X);
accel = acc_G_numeric;

end
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