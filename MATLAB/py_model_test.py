import numpy as np
from scipy.integrate import odeint

# Define constants
m = 1.587
g = 9.81
kt = 3.7102e-5
kd = 7.6933e-07
J_r = 0.03

W1 = 400
W2 = 400
W3 = 400
W4 = 400

I = np.array([[0.224, 0, 0],
                  [0, 0.224, 0],
                  [0, 0, 0.224]])
P1 = np.array([0.3, 0, 0])
P2 = np.array([0, 0.3, 0])
P3 = np.array([-0.3, 0, 0])
P4 = np.array([0, -0.3, 0])

p_mass = [0.05,0.05,0.05,0.05]


Fd = 0

# Update moment of inertia tensor
def update_inertia_tensor(inertia_tensor, point_mass_positions, point_masses):
    # Check if the input parameters are valid
    if inertia_tensor.shape != (3, 3):
        raise ValueError("Inertia tensor must be a 3x3 matrix")
    if len(point_mass_positions) != len(point_masses):
        raise ValueError("Number of point mass positions and masses must match")

    # Loop through each point mass and update the inertia tensor
    for i in range(len(point_mass_positions)):
        r = point_mass_positions[i]  # Position vector of the point mass from the center of gravity
        m = point_masses[i]  # Mass of the point mass

        # Use the parallel axis theorem to update the inertia tensor
        r_squared = np.dot(r, r)
        delta_I = m * np.array([[r_squared, -np.dot(r, r), 0],
                               [-np.dot(r, r), r_squared, 0],
                               [0, 0, r_squared]])
        inertia_tensor += delta_I

    return inertia_tensor

# Define the system of ODEs
def drone_dynamics(X, t):
    phi = X[0] 
    phi_dot = X[1]
    theta =  X[2] 
    theta_dot =  X[3] 
    psi = X[4] 
    psi_dot = X[5] 
   
    
    RbG = np.array([[np.cos(psi) * np.cos(theta), np.sin(psi) * np.cos(theta), -np.sin(theta)],
                    [np.cos(psi) * np.sin(phi) * np.sin(theta) - np.cos(phi) * np.sin(psi),
                     np.sin(phi) * np.sin(psi) * np.sin(theta) + np.cos(phi) * np.cos(psi),
                     np.cos(theta) * np.sin(phi)],
                    [np.cos(phi) * np.cos(psi) * np.sin(theta) + np.sin(phi) * np.sin(psi),
                     np.cos(phi) * np.sin(psi) * np.sin(theta) - np.cos(psi) * np.sin(phi),
                     np.cos(phi) * np.cos(theta)]])
    
    Fg = np.array([0, 0, -m * g])
    Fd = 0
    
    W = np.array([W1*t, W2*t, W3*t, W4*t])
    T = np.array([W[0] ** 2 * kt, W[1] ** 2 * kt, W[2] ** 2 * kt, W[3] ** 2 * kt])
    
    F_thrust = np.dot(RbG, np.array([0, 0, np.sum(T)]))
    
    Acc_G = (1 / m) * (Fg - F_thrust - Fd)
    

    I_u = update_inertia_tensor(I,[P1,P2,P3,P4],p_mass)

    # print(I_u)

    Tm1 = np.cross(P1, np.array([0, 0, -T[0]]))
    Tm2 = np.cross(P2, np.array([0, 0, -T[1]]))
    Tm3 = np.cross(P3, np.array([0, 0, -T[2]]))
    Tm4 = np.cross(P4, np.array([0, 0, -T[3]]))
    Tm_pitch_roll = Tm1 + Tm2 + Tm3 + Tm4
    
    M = np.array([1, -1, 1, -1])
    Tm_yaw = np.array([0, 0, kd * np.sum(M * W**2)])
    Tm_total = Tm_pitch_roll + Tm_yaw


    w = np.array([phi_dot - np.sin(theta) * psi_dot,
                  np.cos(phi) * theta_dot + np.cos(theta) * np.sin(phi) * psi_dot,
                  np.cos(phi) * np.cos(theta) * psi_dot - np.sin(phi) * theta_dot])
    
    T_g = np.array([J_r * w[2] * (M[0] * W[0] + M[1] * W[1] + M[2] * W[2] + M[3] * W[3]),
                    -J_r * w[1] * ((M[0] * W[0] + M[1] * W[1]) + M[2] * W[2] + M[3] * W[3]),
                    0])
    
    T_g_2 = np.cross(w, np.dot(I_u, w))
    
    w_dot = np.dot(np.linalg.inv(I_u), Tm_total - T_g - T_g_2)
    
    
    return [phi, w_dot[0], theta, w_dot[1], psi, w_dot[2]]  # Only interested in phi, theta, psi accelerations

# Initial conditions
X0 = [0, 0, 0, 0, 0, 0]

# Time vector
t = np.linspace(0, 10, 10)  # 10 points from 0 to 1

# Solve the system of ODEs
X = odeint(drone_dynamics, y0 = X0, t = t)

# Access phi, theta, psi velocities:
phi_dot_ = X[:, 1]
theta_dot_ = X[:, 3]
psi_dot_ = X[:, 5]

print('theta_dot:',theta_dot_)
print('phi_dot:',phi_dot_)
print('psi_dot:',psi_dot_)
