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
L_values = [0.1, 0.15, 0.2, 0.4];  % You can adjust these values as needed

rpm = 26280 * (2*pi/60);
% Initialize a cell array to store results for different L values
results = cell(length(L_values), 1);


for i = 1:length(L_values)
    L = L_values(i);


    syms phi(t) theta(t) si(t) t 
    phi(t); % phi(t) roll %
    theta(t); % pitch
    si(t); % yaw
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
        syms m g T1(t) T2(t) T3(t) T4(t) W1(t) W2(t) W3(t) W4(t) t kt
        Fg=[0;0;m*g]; % airplane weight 
        m= 0.300;
        Fd=0; % drag =0 
        kt=3.7102e-5;
    % 4 speeds for each motor
        W1(t) ;
        W2(t) ;
        W3(t) ;
        W4(t) ;
    % Thrust as a function of time given motor speed
        T1(t)=W1(t)^2*kt ;
        T2(t)=W2(t)^2*kt ;
        T3(t)=W3(t)^2*kt ;
        T4(t)=W4(t)^2*kt;
        F_thrust= RGb*[0;0;T1(t)+T2(t)+T3(t)+T4(t)]; % Total thrust 
    % Final Equation of motion
        Acc_G=(1/m)*(Fg-F_thrust - Fd); % Translation EOM (X_dot_dot; Y_dot_dot; Z_dot_dot)
        % Acc_G_str = string(Acc_G)
        % Acc_G_str=strrep(Acc_G_str,'phi(t)','X(1)');
        % Acc_G_str=strrep(Acc_G_str,'theta(t)','X(3)');
        % Acc_G_str=strrep(Acc_G_str,'si(t)','X(5)');
        % Acc_G_str=strrep(Acc_G_str,'W1(t)',num2str(rpm*0.502));
        % Acc_G_str=strrep(Acc_G_str,'W2(t)',num2str(rpm*0.502));
        % Acc_G_str=strrep(Acc_G_str,'W3(t)',num2str(rpm*0.502));
        % Acc_G_str=strrep(Acc_G_str,'W4(t)',num2str(rpm*0.502));
    %% Rotational EOM 
        I=[0.224 0 0; 
            0 0.224 0; 
            0 0 0.224]; %Moment of inertia matrix


    % Position of each motor relative to body frame from CG
        P1=[L;0;0];    
        P2=[0;L;0];
        P3=[-L;0;0];
        P4=[0;-L;0];

        I = update_inertia_tensor(I,[P1,P2,P3,P4],[0.05,0.05,0.05,0.05]);
    %Motors moment arm contribution in body axis frame
        Tm1=cross(P1,[0;0;-T1(t)]);
        Tm2=cross(P2,[0;0;-T2(t)]);
        Tm3=cross(P3,[0;0;-T3(t)]);
        Tm4=cross(P4,[0;0;-T4(t)]);
        Tm_pitch_roll=Tm1+Tm2+Tm3+Tm4;
    %kd is drag torque proportionality constant
        kd=7.6933e-07;

    % rotation direction of each motor. +ve=cw, -ve=ccw
    % Mi corresponds to Pi above
        syms M1 M2 M3 M4
        M1=1;
        M2=-1;
        M3=1;
        M4=-1;
        Tm_yaw=[0;0;kd*( M1*(W1(t))^2 + M2*(W2(t))^2 + M3*(W3(t))^2 + M4*(W4(t))^2)];
        Tm_total= Tm_pitch_roll + Tm_yaw ;
    % torque moment of inertia 
    syms phi_derivative(t) theta_derivative(t) si_derivative(t) J_r p q r 
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
%     J_r=J_r; % Motor moment of inertia assume point mass
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
        w_dot_str=strrep(w_dot_str,'W1(t)',num2str(rpm*0.502));
        w_dot_str=strrep(w_dot_str,'W2(t)',num2str(rpm*0.502));
        w_dot_str=strrep(w_dot_str,'W3(t)',num2str(rpm*0.5));
        w_dot_str=strrep(w_dot_str,'W4(t)',num2str(rpm*0.5));

    w_dot_str;
    % kt 

    dvdt=@(t,X) [X(2);eval(w_dot_str(1,1));
        X(4);eval(w_dot_str(2,1));
        X(6);eval(w_dot_str(3,1))];

% dvdt=@(t,X) [X(2);fh;X(4);fh2;X(6);fh3]

    [t,X]=ode45(dvdt,[0:0.1:1],[0.0;0;0;0;0;0]);
    results{i} = {L, t, X};
end


% Create arrays to store maximum values and corresponding L values
max_values = zeros(length(L_values), 1);

for i = 1:length(L_values)
    X = results{i}{3};  % Extract the state variable X
    max_values(i) = X(11, 4);  % Calculate the maximum value in the 3rd column
end

% Plot the relationship between L and the maximum value
figure;
plot(L_values, max_values, 'b-');
xlabel('Arm Length (L)');
ylabel('Maximum Value in 4th Column of X (theta dot)');
title('Relationship between L and Maximum Value in change in angular velocity:');
grid on;

% Plot the results for different L values
figure;
for i = 1:length(L_values)
    L = results{i}{1};
    t = results{i}{2};
    X = results{i}{3};

    subplot(2, 2, i);  % Adjust the subplot layout as needed
    hold on
    plot(t, X(:,2));
    plot(t, X(:,4));
    plot(t, X(:,6));
    title(['Arm Length L = ', num2str(L),'m']);
    xlabel('Time (sec)');
    ylabel('Eulers Angualar Velocity (rad/s)');
    legend('phi dot','theta dot','si dot');
    hold off
end


%% Function 
function inertia_tensor = update_inertia_tensor(inertia_tensor, point_mass_positions, point_masses)
    % Check if the input parameters are valid
    if ~isequal(size(inertia_tensor), [3, 3])
        error('Inertia tensor must be a 3x3 matrix');
    end
    if length(point_mass_positions) ~= length(point_masses)
        error('Number of point mass positions and masses must match');
    end

    % Loop through each point mass and update the inertia tensor
    for i = 1:length(point_mass_positions)
        r = point_mass_positions(i);  % Position vector of the point mass from the center of gravity
        m = point_masses(i);  % Mass of the point mass

        % Use the parallel axis theorem to update the inertia tensor
        r_squared = r * r';
        delta_I = m * [r_squared, 0, 0; 
                      0, r_squared, 0;
                      0, 0, r_squared];
        inertia_tensor = inertia_tensor + delta_I;
    end
end


