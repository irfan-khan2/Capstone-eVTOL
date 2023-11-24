close all 
clear all si


%% Define else where
syms phi theta si
phi; % phi roll %
theta; % pitch
si; % yaw



%%
%%Translational Conversion Matrix 
RbG=[ cos(si)*cos(theta), sin(si)*cos(theta), -sin(theta); 
      cos(si)*sin(phi)*sin(theta)-cos(si)*sin(si), sin(phi)*sin(si)*sin(theta)+cos(phi)*cos(si), cos(theta)*sin(phi); 
      cos(phi)*cos(si)*sin(theta)+sin(phi)*sin(si), cos(phi)*sin(si)*sin(theta)-cos(si)*sin(phi), cos(phi)*cos(theta)];
% body from global 

RGb= RbG' ;% global from body 

% Rotational Conversion Matrix 
%%
% Rotational Conversion Matrix 
S=[ 1, 0, -sin(theta); 
    0, cos(phi), sin(phi)*cos(theta);
    0, -sin(phi), cos(phi)*cos(theta)];
    % global from body 

S_inv= S^-1 ; % body from global 
%% Translational EOM 
syms m g T1 T2 T3 T4 W1 W2 W3 W4 
m; % mass in kg  

Fg=[0;0;m*g]; % airplane weight 
Fd=0; % drag =0 

% 4 speeds for each motor
W1= 10;
W2= 10;
W3= 10;
W4= 10;

% Thrust as a function of time given motor speed
T1= 5;
T2= 5;
T3= 5;
T4= 5;

F_thrust= RGb*[0;0;T1+T2+T3+T4]; % Total thrust 

Acc_G=(1/m)*(Fg-F_thrust - Fd); % Translation EOM 


%% Rotational EOM 

syms  Ixx Ixy Ixz Iyx Iyy Iyz Izx Izy Izz 
syms P1x P1y P1z P2x P2y P2z P3x P3y P3z P4x P4y P4z

% I=[Ixx Ixy Ixz; 
%     Ixy Iyy Iyz; 
%     Ixz Iyz Izz]; %Moment of inertia 
I=[1 8 0;
    0 1 0;
    0 0 1]; %Moment of inertia

% Position of each motor relative to body frame from CG
P1=[1;1;0];
P2=[-1;1;0];
P3= [-1;-1;0];
P4=[1;-1;0];

%Motors moment arm contribution in body axis frame
Tm1=cross(P1,[0;0;-T1]);
Tm2=cross(P2,[0;0;-T2]);
Tm3=cross(P3,[0;0;-T3]);
Tm4=cross(P4,[0;0;-T4]);

Tm_pitch_roll=Tm1+Tm2+Tm3+Tm4;

%kd is drag torque proportionality constant
kd=1;
% rotation direction of each motor. +ve=cw, -ve=ccw

syms M1 M2 M3 M4
M1 = 1;
M2 =-1;
M3 = 1;
M4 =-1;
Tm_yaw=[0;0;kd*( M1*(W1)^2 + M2*(W2)^2 + M3*(W3)^2 + M4*(W4)^2)];

Tm_total= Tm_pitch_roll + Tm_yaw ;

%gyroscopic effect: 


 % torque moment of inertia 
syms phi_derivative theta_derivative si_derivative J_r

T_g=[J_r*theta_derivative*(M1*W1+M2*W2)+M3*W3+M4*W4; 
    -1*J_r*phi_derivative*(M1*W1+M2*W2)+M3*W3+M4*W4;
    0];
T_g=0;
w=[phi_derivative;theta_derivative;si_derivative];

% third term:
T_g_2=cross(w,I*w);

J_r = 1;

w_dot = I\((Tm_total)-(T_g)-(T_g_2));


theta_=[];
theta_dot=[];
phi_=[];
phi_dot=[];
si_=[];
si_dot=[];

dwdt = zeros(6,1);
dwdt_theta_dot = matlabFunction(w_dot(1), 'Vars', [theta_derivative,theta,phi_derivative,phi,si_derivative,si]);
dwdt(4) = float_phi_derivative;
dwdt(5) = w_dot(3);
dwdt(6) = float_si_derivative;

[theta_dot,theta_,phi_dot,phi_,si_dot,si_] = ode45(@fun,[0 10],[0 0 0 0 0 0])

%% Function to convert symbols to floats 
function result_equation = convert_symbols_to_floats(input_equation)
    % Extract all the symbols from the input equation
    symbol_list = symvar(input_equation);
    
    % Create a cell array to store original symbols and undefined floats
    symbol_cells = cell(1, numel(symbol_list)*2);
    
    % Replace each symbol with an undefined float and store the mapping
    for i = 1:numel(symbol_list)
        symbol_cells{i} = symbol_list(i);
        symbol_cells{i+numel(symbol_list)} = sym(['_', char(symbol_list(i))]);
        input_equation = subs(input_equation, symbol_cells{i}, symbol_cells{i+numel(symbol_list)});
    end
    
    % Return the equation with undefined floats
    result_equation = input_equation;
end
%%
function dwdt = fun(theta_derivative,theta,phi_derivative,phi,si_derivative,si,t)
    dwdt = [dwdt_theta_dot,theta_derivative,dwdt_phi_dot,phi_derivative, dwdt_si_dot,si_derivative];
end