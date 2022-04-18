clear; close all;
%Position = csvread('VICON_MASTER.csv');

%24-27 seconds deep squat
pos = csvread('VICON_MASTER.csv',4808,0,'A4809..AII5409');

Reaction = csvread('VICON_MASTER.csv',38863,0,'A38864..N41864');

Vicon_Moment = readmatrix('VICON_MOMENT_MAG_KG_ADJUSTED.xlsx','Range','A4:AG604');

Vicon_t = Vicon_Moment(:,1);

M_Vicon_Hip_L = Vicon_Moment(:,17);
M_Vicon_Hip_R = Vicon_Moment(:,22);

M_Vicon_Knee_L = Vicon_Moment(:,6);
M_Vicon_Knee_R = Vicon_Moment(:,11);

M_Vicon_Ankle_L = Vicon_Moment(:,28);
M_Vicon_Ankle_R = Vicon_Moment(:,33);

%10-13 seconds deep squat
%pos = csvread('VICON_MASTER.csv',2008,0,'A2009..AII2609');


% Remember that these parameters are zero-based, 
% so that column A maps to 0 and row 252 maps to 251.
t = pos(:, 1);

Ax = pos(:, 101); Ay = pos(:, 102); Az = pos(:, 103);   %Left Toe
Bx = pos(:, 380); By = pos(:, 381); Bz = pos(:, 382);   %Left Ankle
Cx = pos(:, 368); Cy = pos(:, 369); Cz = pos(:, 370);   %Left Knee
Dx = pos(:, 377); Dy = pos(:, 378); Dz = pos(:, 379);   %Left Hip
Ex = pos(:, 521); Ey = pos(:, 522); Ez = pos(:, 523);   %Left Shoulder

Fx = pos(:, 14); Fy = pos(:, 15); Fz = pos(:, 16);    %Neck
Gx = pos(:, 464); Gy = pos(:, 465); Gz = pos(:, 466);    %Head

Hx = pos(:, 125); Hy = pos(:, 126); Hz = pos(:, 127);   %Right Toe
Ix = pos(:, 428); Iy = pos(:, 429); Iz = pos(:, 430);   %Right Ankle
Jx = pos(:, 416); Jy = pos(:, 417); Jz = pos(:, 418);   %Right Knee
Kx = pos(:, 425); Ky = pos(:, 426); Kz = pos(:, 427);   %Right Hip
Lx = pos(:, 557); Ly = pos(:, 558); Lz = pos(:, 559);   %Right Shoulder
Mx = pos(:, 356); My = pos(:, 357); Mz = pos(:, 358);   %Pelvis
Nx = pos(:, 524); Ny = pos(:, 525); Nz = pos(:, 526);   %Left Wrist
Ox = pos(:, 560); Oy = pos(:, 561); Oz = pos(:, 562);   %Right Wrist
Px = pos(:, 98); Py = pos(:, 99); Pz = pos(:, 100);   %Left Heel
Qx = pos(:, 122); Qy = pos(:, 123); Qz = pos(:, 124);   %Right Heel

m_total = 80;
L_total = 1.86;

density = 1000;
g = [0, 0, -9.81];
g_z = -9.81;

head_scale = 1.1;

% Estimates for trunk are from Plagenhoef et al., 1983 
% https://exrx.net/Kinesiology/Segments

m_foot_L  = 0.5*0.0143*m_total;
m_leg_L   = 0.5*0.0475*m_total;
m_thigh_L = 0.5*0.2416*m_total;    %And pelvis
m_arm_L = 0.5*0.057*m_total;

m_torso = 0.551*m_total;     %not including arms

m_foot_R  = 0.5*0.0143*m_total;
m_leg_R   = 0.5*0.0475*m_total;
m_thigh_R = 0.5*0.2416*m_total;    %And pelvis
m_arm_R = 0.5*0.057*m_total;

m_head_neck  = 0.0826*m_total;

m_check = m_foot_L + m_leg_L + m_thigh_L + m_arm_L + m_foot_R + m_leg_R + m_thigh_R + m_torso + m_arm_R + m_head_neck;

% mass of body above hip
m_above_hip_L = m_arm_L + 0.5 * (m_torso +  m_head_neck);
m_above_hip_R = m_arm_R + 0.5 * (m_torso +  m_head_neck);

% mass of body above knee
m_above_knee_L = m_arm_L + m_thigh_L + 0.5 * (m_torso +  m_head_neck);

% mass of body above knee
m_above_knee_R = m_arm_R + m_thigh_R + 0.5 * (m_torso +  m_head_neck);

% mass of body above ankle
m_above_ankle_L = m_arm_L + m_leg_L + m_thigh_L + 0.5 * (m_torso +  m_head_neck);

% mass of body above ankle
m_above_ankle_R = m_arm_R + m_leg_R + m_thigh_R + 0.5 * (m_torso +  m_head_neck);

%vertical force through hip
F_matrix_hip_L_z = m_above_hip_L * g;
F_matrix_hip_R_z = m_above_hip_R * g;

%vertical force through knee
F_matrix_knee_L_z = m_above_knee_L * g;
F_matrix_knee_R_z = m_above_knee_R * g;

%vertical force through ankle
F_matrix_ankle_L_z = m_above_ankle_L * g;
F_matrix_ankle_R_z = m_above_ankle_R * g;
%preallocation for speed


num_rows = height(Ax);

for i = 1:num_rows

L_foot_L(i) = sqrt( (Ax(i) - Bx(i))^2 + (Ay(i) - By(i))^2 + (Az(i) - Bz(i))^2);
L_leg_L(i) = sqrt( (Cx(i) - Bx(i))^2 + (Cy(i) - By(i))^2 + (Cz(i) - Bz(i))^2);
L_thigh_L(i) = sqrt( (Dx(i) - Cx(i))^2 + (Dy(i) - Cy(i))^2 + (Dz(i) - Cz(i))^2);
L_torso_L(i) = sqrt( (Fx(i) - Mx(i))^2 + (Fy(i) - My(i))^2 + (Fz(i) - Mz(i))^2);
L_arm_L(i) = sqrt( (Nx(i) - Ex(i))^2 + (Ny(i) - Ey(i))^2 + (Nz(i) - Ez(i))^2);

L_foot_R(i) = sqrt( (Hx(i) - Ix(i))^2 + (Hy(i) - Iy(i))^2 + (Hz(i) - Iz(i))^2);
L_leg_R(i) = sqrt( (Jx(i) - Ix(i))^2 + (Jy(i) - Iy(i))^2 + (Jz(i) - Iz(i))^2);
L_thigh_R(i) = sqrt( (Kx(i) - Jx(i))^2 + (Ky(i) - Jy(i))^2 + (Kz(i) - Jz(i))^2);
L_torso_R(i) = sqrt( (Fx(i) - Mx(i))^2 + (Fy(i) - My(i))^2 + (Fz(i) - Mz(i))^2);
L_arm_R(i) = sqrt( (Ox(i) - Lx(i))^2 + (Oy(i) - Ly(i))^2 + (Oz(i) - Lz(i))^2);

L_head_neck(i)  = head_scale*sqrt( (Gx(i) - Fx(i))^2 + (Gy(i) - Fy(i))^2 + (Gz(i) - Fz(i))^2);

L_check(i) = L_leg_L(i) + L_thigh_L(i) + L_torso_L(i) + L_head_neck(i);

theta_foot_L(i)  = atand((By(i) - Ay(i))/sqrt(((Bx(i) - Ax(i))^2) + ((Bz(i) - Az(i))^2)));
theta_leg_L(i)   = atand((Cy(i) - By(i))/sqrt(((Cx(i) - Bx(i))^2) + ((Cz(i) - Bz(i))^2)));
theta_thigh_L(i) = atand((Dy(i) - Cy(i))/sqrt(((Dx(i) - Cx(i))^2) + ((Dz(i) - Cz(i))^2)));
theta_arm_L(i) = atand((Ey(i) - Ny(i))/sqrt(((Ex(i) - Nx(i))^2) + ((Ez(i) - Nz(i))^2)));

theta_torso(i) = atand((Fy(i) - My(i))/sqrt(((Fx(i) - Mx(i))^2) + ((Fz(i) - Mz(i))^2)));

theta_foot_R(i)  = atand((Iy(i) - Hy(i))/sqrt(((Ix(i) - Hx(i))^2) + ((Iz(i) - Hz(i))^2)));
theta_leg_R(i) = atand((Jy(i) - Iy(i))/sqrt(((Jx(i) - Ix(i))^2) + ((Jz(i) - Iz(i))^2)));
theta_thigh_R(i) = atand((Ky(i) - Jy(i))/sqrt(((Kx(i) - Jx(i))^2) + ((Kz(i) - Jz(i))^2)));
theta_arm_R(i) = atand((Ly(i) - Oy(i))/sqrt(((Lx(i) - Ox(i))^2) + ((Lz(i) - Oz(i))^2)));

theta_head(i)  = atand((Gy(i) - Fy(i))/sqrt(((Gx(i) - Fx(i))^2) + ((Gz(i) - Fz(i))^2)));

foot_L_cg_x(i) = (Ax(i) + Bx(i))/2;
foot_L_cg_y(i) = (Ay(i) + By(i))/2;
foot_L_cg_z(i) = (Az(i) + Bz(i))/2;

leg_L_cg_x(i) =  (Bx(i) + Cx(i))/2;
leg_L_cg_y(i) =  (By(i) + Cy(i))/2;
leg_L_cg_z(i) =  (Bz(i) + Cz(i))/2;

thigh_L_cg_x(i) =  (Cx(i) + Dx(i))/2;
thigh_L_cg_y(i) =  (Cy(i) + Dy(i))/2;
thigh_L_cg_z(i) =  (Cz(i) + Dz(i))/2;

arm_L_cg_x(i) =  (Nx(i) + Ex(i))/2;
arm_L_cg_y(i) =  (Ny(i) + Ey(i))/2;
arm_L_cg_z(i) =  (Nz(i) + Ez(i))/2;

torso_cg_x(i) =  (Fx(i) + Mx(i))/2;
torso_cg_y(i) =  (Fy(i) + My(i))/2;
torso_cg_z(i) =  (Fz(i) + Mz(i))/2;

head_neck_cg_x(i) =  (Fx(i) + Gx(i))/2;
head_neck_cg_y(i) =  (Fy(i) + Gy(i))/2;
head_neck_cg_z(i) =  (head_scale * (Fz(i) + Gz(i)))/2;

foot_R_cg_x(i) = (Ix(i) + Hx(i))/2;
foot_R_cg_y(i) = (Iy(i) + Hy(i))/2;
foot_R_cg_z(i) = (Iz(i) + Hz(i))/2;

leg_R_cg_x(i) =  (Jx(i) + Ix(i))/2;
leg_R_cg_y(i) =  (Jy(i) + Iy(i))/2;
leg_R_cg_z(i) =  (Jz(i) + Iz(i))/2;

thigh_R_cg_x(i) =  (Kx(i) + Jx(i))/2;
thigh_R_cg_y(i) =  (Ky(i) + Jy(i))/2;
thigh_R_cg_z(i) =  (Kz(i) + Jz(i))/2;

arm_R_cg_x(i) =  (Ox(i) + Lx(i))/2;
arm_R_cg_y(i) =  (Oy(i) + Ly(i))/2;
arm_R_cg_z(i) =  (Oz(i) + Lz(i))/2;


%Vector r direction (cg - joint centre)
%GRF acting though heel in this case
r_1_L(i,:) = [(foot_L_cg_x(i) - Px(i)), (foot_L_cg_y(i) - Py(i)), (foot_R_cg_z(i) - Pz(i))];
r_1_R(i,:) = [(foot_R_cg_x(i) - Qx(i)), (foot_R_cg_y(i) - Qy(i)), (foot_R_cg_z(i) - Qz(i))];

r_2_L(i,:) = [(foot_L_cg_x(i) - Bx(i)), (foot_L_cg_y(i) - By(i)), (foot_R_cg_z(i) - Bz(i))];
r_2_R(i,:) = [(foot_R_cg_x(i) - Ix(i)), (foot_R_cg_y(i) - Iy(i)), (foot_R_cg_z(i) - Iz(i))];

r_3_L(i,:) = [(leg_L_cg_x(i) - Bx(i)), (leg_L_cg_y(i) - By(i)), (leg_R_cg_z(i) - Bz(i))];
r_3_R(i,:) = [(leg_R_cg_x(i) - Ix(i)), (leg_R_cg_y(i) - Iy(i)), (leg_R_cg_z(i) - Iz(i))];

r_4_L(i,:) = [(leg_L_cg_x(i) - Cx(i)), (leg_L_cg_y(i) - Cy(i)), (leg_R_cg_z(i) - Cz(i))];
r_4_R(i,:) = [(leg_R_cg_x(i) - Jx(i)), (leg_R_cg_y(i) - Jy(i)), (leg_R_cg_z(i) - Jz(i))];

r_5_L(i,:) = [(thigh_L_cg_x(i) - Cx(i)), (thigh_L_cg_y(i) - Cy(i)), (thigh_R_cg_z(i) - Cz(i))];
r_5_R(i,:) = [(thigh_R_cg_x(i) - Jx(i)), (thigh_R_cg_y(i) - Jy(i)), (thigh_R_cg_z(i) - Jz(i))];

r_6_L(i,:) = [(thigh_L_cg_x(i) - Dx(i)), (thigh_L_cg_y(i) - Dy(i)), (thigh_R_cg_z(i) - Dz(i))];
r_6_R(i,:) = [(thigh_R_cg_x(i) - Kx(i)), (thigh_R_cg_y(i) - Ky(i)), (thigh_R_cg_z(i) - Kz(i))];

r_7_L(i,:) = [(Mx(i) - Dx(i)), (My(i) - Dy(i)), (Mz(i) - Dz(i))];
r_7_R(i,:) = [(Mx(i) - Kx(i)), (My(i) - Ky(i)), (Mz(i) - Kz(i))];

r_8(i,:) = [(torso_cg_x(i) - Mx(i)), (torso_cg_y(i) - My(i)), (torso_cg_z(i) - Mz(i))];
r_9(i,:) = [(torso_cg_x(i) - Fx(i)), (torso_cg_y(i) - Fy(i)), (torso_cg_z(i) - Fz(i))];

r_10(i,:) = [(head_neck_cg_x(i)- Fx(i)),(head_neck_cg_y(i)-Fy(i)), (head_neck_cg_z(i)-Fz(i))];

r_11_L(i,:) = [(foot_L_cg_x(i) - Ax(i)), (foot_L_cg_y(i) - Ay(i)), (foot_L_cg_z(i) - Az(i))];%Toe to cg
r_11_R(i,:) = [(foot_R_cg_x(i) - Hx(i)), (foot_R_cg_y(i) - Hy(i)), (foot_R_cg_z(i) - Hz(i))];

r_12_L(i,:) = [(Fx(i) - Ex(i)), (Fy(i) - Ey(i)), (Fz(i) - Ez(i))];%Shoulder to neck
r_12_R(i,:) = [(Fx(i) - Lx(i)), (Fy(i) - Ly(i)), (Fz(i) - Lz(i))];

r_13_L(i,:) = [(Ex(i) - Nx(i)), (Ey(i) - Ny(i)), (Ez(i) - Nz(i))];%Wrist to shoulder
r_13_R(i,:) = [(Lx(i) - Ox(i)), (Ly(i) - Oy(i)), (Lz(i) - Oz(i))];


% Whole body cg

whole_body_cg_x(i) =  (m_foot_L*foot_L_cg_x(i) + m_leg_L*leg_L_cg_x(i) ...
    + m_thigh_L*thigh_L_cg_x(i) + m_foot_R*foot_R_cg_x(i) ...
    + m_leg_R*leg_R_cg_x(i) + m_thigh_R*thigh_R_cg_x(i) ...
    + m_arm_L*arm_L_cg_x(i) + m_arm_R*arm_R_cg_x(i) ...
    + m_torso*torso_cg_x(i) +  m_head_neck*head_neck_cg_x(i))/m_total;

whole_body_cg_y(i) =  (m_foot_L*foot_L_cg_y(i) + m_leg_L*leg_L_cg_y(i) + m_thigh_L*thigh_L_cg_y(i) ...
    + m_foot_R*foot_R_cg_y(i) + m_leg_R*leg_R_cg_y(i) + m_thigh_R*thigh_R_cg_y(i) ...
    + m_arm_L*arm_L_cg_y(i) + m_arm_R*arm_R_cg_y(i) ...
    + m_torso*torso_cg_y(i) + m_head_neck*head_neck_cg_y(i))/m_total;

whole_body_cg_z(i) =  (m_foot_L*foot_L_cg_z(i) + m_leg_L*leg_L_cg_z(i) + m_thigh_L*thigh_L_cg_z(i) ...
    + m_foot_R*foot_R_cg_z(i) + m_leg_R*leg_R_cg_z(i) + m_thigh_R*thigh_R_cg_z(i) ...
    + m_arm_L*arm_L_cg_z(i) + m_arm_R*arm_R_cg_z(i) ...
    + m_torso*torso_cg_z(i) + m_head_neck*head_neck_cg_z(i))/m_total;

%Left Hip
% cg of body above left hip in pixels
cg_above_hip_L_x(i) =  (m_arm_L*arm_L_cg_x(i)  ...
        + 0.5 * m_torso*torso_cg_x(i) ...
    + 0.5 * m_head_neck*head_neck_cg_x(i))/m_above_hip_L;

cg_above_hip_L_y(i) =  (m_arm_L*arm_L_cg_y(i) ...
+ 0.5 * m_torso*torso_cg_y(i) ...
    +  0.5 * m_head_neck*head_neck_cg_y(i))/m_above_hip_L;

cg_above_hip_L_z(i) =  (m_arm_L*arm_L_cg_z(i) ...
+ 0.5 * m_torso*torso_cg_z(i) ...
    +  0.5 * m_head_neck*head_neck_cg_z(i))/m_above_hip_L;


%d is horizontal moment arm (meters) about the hip of the weight above hip 
d_hip_L_x(i) = (cg_above_hip_L_x(i)  - Dx(i));
d_hip_L_y(i) = (cg_above_hip_L_y(i) - Dy(i));
d_hip_L_z(i) = (cg_above_hip_L_z(i) - Dz(i));

d_matrix_hip_L(i,:) = [d_hip_L_x(i), d_hip_L_y(i), d_hip_L_z(i)];


%moment about hip (Nm)
M_hip_L(i,:) = cross(d_matrix_hip_L(i,:), F_matrix_hip_L_z(:));

M_hip_mag_L(i) = sqrt((M_hip_L(i, 1)^2) + (M_hip_L(i, 2)^2) + (M_hip_L(i, 3)^2));

%Right Hip
% cg of body above left hip in pixels
cg_above_hip_R_x(i) =  (m_arm_R*arm_R_cg_x(i)  ...
        + 0.5 * m_torso*torso_cg_x(i) ...
    + 0.5 * m_head_neck*head_neck_cg_x(i))/m_above_hip_R;

cg_above_hip_R_y(i) =  (m_arm_R*arm_R_cg_y(i) ...
+ 0.5 * m_torso*torso_cg_y(i) ...
    +  0.5 * m_head_neck*head_neck_cg_y(i))/m_above_hip_R;

cg_above_hip_R_z(i) =  (m_arm_R*arm_R_cg_z(i) ...
+ 0.5 * m_torso*torso_cg_z(i) ...
    +  0.5 * m_head_neck*head_neck_cg_z(i))/m_above_hip_R;

%d is horizontal moment arm (meters) about the hip of the weight above hip 
d_hip_R_x(i) = (cg_above_hip_R_x(i) - Kx(i));
d_hip_R_y(i) = (cg_above_hip_R_y(i) - Ky(i));
d_hip_R_z(i) = (cg_above_hip_R_z(i) - Kz(i));

d_matrix_hip_R(i,:) = [d_hip_R_x(i), d_hip_R_y(i), d_hip_R_z(i)];



%moment about hip (Nm)
M_hip_R(i,:) = cross(d_matrix_hip_R(i,:), F_matrix_hip_R_z(:));

M_hip_mag_R(i) = sqrt((M_hip_R(i, 1).^2) + (M_hip_R(i, 2).^2) + (M_hip_R(i, 3).^2));



%Left Knee

% cg of body above knee in pixels
cg_above_knee_L_x(i) =  (m_thigh_L*thigh_L_cg_x(i) + 0.5 * m_torso*torso_cg_x(i) ...
    + m_arm_L*arm_L_cg_x(i)  ...
    +  0.5 * m_head_neck*head_neck_cg_x(i))/m_above_knee_L;

cg_above_knee_L_y(i) =  (m_thigh_L*thigh_L_cg_y(i) + 0.5 * m_torso*torso_cg_y(i) ...
    + m_arm_L*arm_L_cg_y(i) ...
    +  0.5 * m_head_neck*head_neck_cg_y(i))/m_above_knee_L;

cg_above_knee_L_z(i) =  (m_thigh_L*thigh_L_cg_z(i) + 0.5 * m_torso*torso_cg_z(i) ...
    + m_arm_L*arm_L_cg_z(i) ...
    +  0.5 * m_head_neck*head_neck_cg_z(i))/m_above_knee_L;

%d is horizontal moment arm (meters) about the knee
d_knee_L_x(i) = (cg_above_knee_L_x(i) - Cx(i));
d_knee_L_y(i) = (cg_above_knee_L_y(i) - Cy(i));
d_knee_L_z(i) = (cg_above_knee_L_z(i) - Cz(i));

d_matrix_knee_L(i,:) = [d_knee_L_x(i), d_knee_L_y(i), d_knee_L_z(i)];

d_knee_L_horizontal(i) = sqrt((d_knee_L_x(i))^2 + (d_knee_L_y(i))^2);

%moment about knee (Nm)
M_knee_L(i,:) = cross(d_matrix_knee_L(i,:), F_matrix_knee_L_z(:));

M_knee_mag_L(i) = sqrt((M_knee_L(i, 1).^2) + (M_knee_L(i, 2).^2) + (M_knee_L(i, 3).^2));

%Right Knee

% cg of body above knee in pixels
cg_above_knee_R_x(i) =  (m_thigh_R*thigh_R_cg_x(i) + 0.5 * m_torso*torso_cg_x(i) ...
     + m_arm_R*arm_R_cg_x(i) ...
    +  0.5 * m_head_neck*head_neck_cg_x(i))/m_above_knee_R;

cg_above_knee_R_y(i) =  (m_thigh_R*thigh_R_cg_y(i) + 0.5 * m_torso*torso_cg_y(i) ...
    + m_arm_R*arm_R_cg_y(i) ...
    +  0.5 * m_head_neck*head_neck_cg_y(i))/m_above_knee_R;

cg_above_knee_R_z(i) =  (m_thigh_R*thigh_R_cg_z(i) + 0.5 * m_torso*torso_cg_z(i) ...
    + m_arm_R*arm_R_cg_z(i) ...
    +  0.5 * m_head_neck*head_neck_cg_z(i))/m_above_knee_R;

%d is horizontal moment arm (meters) about the knee
d_knee_R_x(i) = (cg_above_knee_R_x(i) - Jx(i));
d_knee_R_y(i) = (cg_above_knee_R_y(i) - Jy(i));
d_knee_R_z(i) = (cg_above_knee_R_z(i) - Jz(i));

d_matrix_knee_R(i,:) = [d_knee_R_x(i), d_knee_R_y(i), d_knee_R_z(i)];

d_knee_ankle_R_x(i) = (Jx(i) - Ix(i));
d_knee_ankle_R_y(i) = (Jy(i) - Iy(i));
d_knee_ankle_R_z(i) = (Jz(i) - Iz(i));

d_knee_ankle_matrix_R(i,:) = [d_knee_ankle_R_x(i), d_knee_ankle_R_y(i), d_knee_ankle_R_z(i)];

d_knee_R_horizontal(i) = sqrt((d_knee_R_x(i))^2 + (d_knee_R_y(i))^2);

%m = 601;
%n = 3;
%Average_reaction = zeros(m,n);
%Average_reaction(i,:) = [0, 0, 750];
%Av_react = [0, 0, -750];


%moment about knee (Nm)
M_knee_R(i,:) = cross(d_matrix_knee_R(i,:), F_matrix_knee_R_z(:)); %+ cross(d_knee_ankle_matrix_R(i,:), Av_react);

M_knee_mag_R(i) = sqrt((M_knee_R(i, 1).^2) + (M_knee_R(i, 2).^2) + (M_knee_R(i, 3).^2));


%Left Ankle

% cg of body above ankle in pixels
cg_above_ankle_L_x(i) =  (m_leg_L*leg_L_cg_x(i) + m_thigh_L*thigh_L_cg_x(i) + 0.5 * m_torso*torso_cg_x(i) ...
    + m_arm_L*arm_L_cg_x(i) ...
    +  0.5 * m_head_neck*head_neck_cg_x(i))/m_above_ankle_L;

cg_above_ankle_L_y(i) =  (m_leg_L*leg_L_cg_y(i) +m_thigh_L*thigh_L_cg_y(i) + 0.5 * m_torso*torso_cg_y(i) ...
    + m_arm_L*arm_L_cg_y(i)...
    +  0.5 * m_head_neck*head_neck_cg_y(i))/m_above_ankle_L;

cg_above_ankle_L_z(i) =  (m_leg_L*leg_L_cg_z(i) +m_thigh_L*thigh_L_cg_z(i) + 0.5 * m_torso*torso_cg_z(i) ...
    + m_arm_L*arm_L_cg_z(i) ...
    +  0.5 * m_head_neck*head_neck_cg_z(i))/m_above_ankle_L;

%d is horizontal moment arm (meters) about the ankle
d_ankle_L_x(i) = (cg_above_ankle_L_x(i) - Bx(i));
d_ankle_L_y(i) = (cg_above_ankle_L_y(i) - By(i));
d_ankle_L_z(i) = (cg_above_ankle_L_z(i) - Bz(i));

d_matrix_ankle_L(i,:) = [d_ankle_L_x(i), d_ankle_L_y(i), d_ankle_L_z(i)];



%moment about ankle (Nm)
M_ankle_L(i,:) = cross(d_matrix_ankle_L(i,:), F_matrix_ankle_L_z(:));

M_ankle_mag_L(i) = sqrt((M_ankle_L(i, 1).^2) + (M_ankle_L(i, 2).^2) + (M_ankle_L(i, 3).^2));

%Right Ankle

% cg of body above ankle in pixels
cg_above_ankle_R_x(i) =  (m_leg_R*leg_R_cg_x(i) + m_thigh_R*thigh_R_cg_x(i) + 0.5 * m_torso*torso_cg_x(i) ...
    + m_arm_R*arm_R_cg_x(i) ...
    +  0.5 * m_head_neck*head_neck_cg_x(i))/m_above_ankle_R;

cg_above_ankle_R_y(i) =  (m_leg_R*leg_R_cg_y(i) +m_thigh_R*thigh_R_cg_y(i) + 0.5 * m_torso*torso_cg_y(i) ...
    + m_arm_R*arm_R_cg_y(i) ...
    +  0.5 * m_head_neck*head_neck_cg_y(i))/m_above_ankle_R;

cg_above_ankle_R_z(i) =  (m_leg_R*leg_R_cg_z(i) +m_thigh_R*thigh_R_cg_z(i) + 0.5 * m_torso*torso_cg_z(i) ...
    + m_arm_R*arm_R_cg_z(i) ...
    +  0.5 * m_head_neck*head_neck_cg_z(i))/m_above_ankle_R;

%d is horizontal moment arm (meters) about the ankle
d_ankle_R_x(i) = (cg_above_ankle_R_x(i) - Ix(i));
d_ankle_R_y(i) = (cg_above_ankle_R_y(i) - Iy(i));
d_ankle_R_z(i) = (cg_above_ankle_R_z(i) - Iz(i));

d_matrix_ankle_R(i,:) = [d_ankle_R_x(i), d_ankle_R_y(i), d_ankle_R_z(i)];


%moment about ankle (Nm)
M_ankle_R(i,:) = cross(d_matrix_ankle_R(i,:), F_matrix_ankle_R_z(:));

M_ankle_mag_R(i) = sqrt((M_ankle_R(i, 1).^2) + (M_ankle_R(i, 2).^2) + (M_ankle_R(i, 3).^2));

%Vicon Ground Reaction
R_left = Reaction(:,4);
R_right = Reaction(:,10);
t_Reaction = Reaction(:,1);
R_total = -(R_left + R_right);

%Ground Reaction
theta(i) = pi/2 * sin((pi / 3) * t(i));
theta_dot(i) = ((pi^2)/6) * cos((pi / 3) * t(i));
theta_ddot(i) = -((pi^3)/18) * sin((pi / 3) * t(i));  %rad/s^2
r(i) = L_thigh_L(i)/2;

a(i) = r(i) * theta_ddot(i);     %m/s^2
R(i) = m_total*g_z - m_total*a(i);
end

%transpose matrices
M_hip_mag_L_tranpose = transpose(M_hip_mag_L);
M_hip_mag_R_tranpose = transpose(M_hip_mag_R);

M_knee_mag_L_transpose = transpose(M_knee_mag_L);
M_knee_mag_R_transpose = transpose(M_knee_mag_R);

M_ankle_mag_L_transpose = transpose(M_ankle_mag_L);
M_ankle_mag_R_transpose = transpose(M_ankle_mag_R);


Sx = [Ax, Bx, Cx, Dx, Ex, Nx, Ox, Fx, Lx, Kx, Mx, Gx, Hx, Ix, Jx];
Sy = [Ay, By, Cy, Dy, Ey, Ny, My, Fy, Ly, Ky, My, Gy, Hy, Iy, Jy];
Sz = [Az, Bz, Cz, Dz, Ez, Nz, Oz, Fz, Lz, Kz, Mz, Gz, Hz, Iz, Jz];

CGx = [foot_L_cg_x;leg_L_cg_x;thigh_L_cg_x;torso_cg_x;head_neck_cg_x;thigh_R_cg_x;leg_R_cg_x;foot_R_cg_x];
CGy = [foot_L_cg_y;leg_L_cg_y;thigh_L_cg_y;torso_cg_y;head_neck_cg_y;thigh_R_cg_y;leg_R_cg_y;foot_R_cg_y];
CGz = [foot_L_cg_z;leg_L_cg_z;thigh_L_cg_z;torso_cg_z;head_neck_cg_z;thigh_R_cg_z;leg_R_cg_z;foot_R_cg_z];

Above_CG_x = [cg_above_hip_L_x; cg_above_hip_R_x; cg_above_knee_R_x; cg_above_knee_L_x; cg_above_ankle_L_x; cg_above_ankle_R_x];
Above_CG_y = [cg_above_hip_L_y; cg_above_hip_R_y; cg_above_knee_R_y; cg_above_knee_L_y; cg_above_ankle_L_y; cg_above_ankle_R_y];
Above_CG_z = [cg_above_hip_L_z; cg_above_hip_R_z; cg_above_knee_R_z; cg_above_knee_L_z; cg_above_ankle_L_z; cg_above_ankle_R_z];





figure (1)
plot(t,M_knee_mag_L)
title('Torque-Time plot of Left Knee')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Knee_L)

figure (2)
plot(t,M_knee_mag_R)
title('Torque-Time plot of Right Knee')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Knee_R)

figure (3)
plot(t,M_hip_mag_L)
title('Torque-Time plot of Left Hip')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Hip_L)

figure (4)
plot(t,M_hip_mag_R)
title('Torque-Time plot of Right Hip')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Hip_R)

figure (5)
plot(t,M_ankle_mag_L)
title('Torque-Time plot of Left Ankle')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Ankle_L)


figure (6)
plot(t,M_ankle_mag_R)
title('Torque-Time plot of Right Ankle')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Ankle_R)

figure (7)
plot(t,R)
title('Ground Reaction Force')
xlabel('Time (s)')
ylabel('Force (N)')
hold on;
plot(t_Reaction, R_total)

figure(8)
plot3(Sx(1,:), Sy(1,:), Sz(1,:), '-o', 'Color', 'b', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on
plot3(CGx(:,1), CGy(:,1), CGz(:,1), '-o', 'Color', 'g', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

figure(9)
plot3(Sx(150,:), Sy(150,:), Sz(150,:), '-o', 'Color', 'b', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on
plot3(CGx(:,150), CGy(:,150), CGz(:,150), '-o', 'Color', 'g', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

figure(10)
plot3(Sx(300,:), Sy(300,:), Sz(300,:), '-o', 'Color', 'b', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on
plot3(CGx(:,300), CGy(:,300), CGz(:,300), '-o', 'Color', 'g', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

figure(11)
plot3(Sx(450,:), Sy(450,:), Sz(450,:), '-o', 'Color', 'b', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on
plot3(CGx(:,450), CGy(:,450), CGz(:,450), '-o', 'Color', 'g', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

figure(12)
plot3(Sx(600,:), Sy(600,:), Sz(600,:), '-o', 'Color', 'b', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on
plot3(CGx(:,600), CGy(:,600), CGz(:,600), '-o', 'Color', 'g', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

%title('Vicon Total Ground Reaction Force')
%xlabel('Time (s)')
%ylabel('Force (N)')


figure(13)
plot3(Sx(1,:), Sy(1,:), Sz(1,:), '-o', 'Color', 'b', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on
plot3(Above_CG_x(:,1), Above_CG_y(:,1), Above_CG_z(:,1), 'o', 'Color', 'r', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on
plot3(CGx(:,1), CGy(:,1), CGz(:,1), 'o', 'Color', 'g', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

figure(14)
plot3(Sx(300,:), Sy(300,:), Sz(300,:), '-o', 'Color', 'b', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on
plot3(Above_CG_x(:,300), Above_CG_y(:,300), Above_CG_z(:,300), 'o', 'Color', 'r', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on
plot3(CGx(:,300), CGy(:,300), CGz(:,300), 'o', 'Color', 'g', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

figure(15)
plot3(Sx(1,:), Sy(1,:), Sz(1,:), 'o', 'Color', 'b', 'MarkerSize',4,'MarkerFaceColor','#D9FFFF')
hold on
plot3(CGx(:,1), CGy(:,1), CGz(:,1), 'o', 'Color', 'g', 'MarkerSize',4,'MarkerFaceColor','#D9FFFF')
xlabel('x')
ylabel('y')
zlabel('z')
hold on
quiver3(Px(1), Py(1), Pz(1), r_1_L(1,1), r_1_L(1,2), r_1_L(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Bx(1), By(1), Bz(1), r_2_L(1,1), r_2_L(1,2), r_2_L(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Bx(1), By(1), Bz(1), r_3_L(1,1), r_3_L(1,2), r_3_L(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Cx(1), Cy(1), Cz(1), r_4_L(1,1), r_4_L(1,2), r_4_L(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Cx(1), Cy(1), Cz(1), r_5_L(1,1), r_5_L(1,2), r_5_L(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Dx(1), Dy(1), Dz(1), r_6_L(1,1), r_6_L(1,2), r_6_L(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Dx(1), Dy(1), Dz(1), r_7_L(1,1), r_7_L(1,2), r_7_L(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Qx(1), Qy(1), Qz(1), r_1_R(1,1), r_1_R(1,2), r_1_R(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ix(1), Iy(1), Iz(1), r_2_R(1,1), r_2_R(1,2), r_2_R(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ix(1), Iy(1), Iz(1), r_3_R(1,1), r_3_R(1,2), r_3_R(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Jx(1), Jy(1), Jz(1), r_4_R(1,1), r_4_R(1,2), r_4_R(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Jx(1), Jy(1), Jz(1), r_5_R(1,1), r_5_R(1,2), r_5_R(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Kx(1), Ky(1), Kz(1), r_6_R(1,1), r_6_R(1,2), r_6_R(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Kx(1), Ky(1), Kz(1), r_7_R(1,1), r_7_R(1,2), r_7_R(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Mx(1), My(1), Mz(1), r_8(1,1), r_8(1,2), r_8(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Fx(1), Fy(1), Fz(1), r_9(1,1), r_9(1,2), r_9(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Fx(1), Fy(1), Fz(1), r_10(1,1), r_10(1,2), r_10(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ax(1), Ay(1), Az(1), r_11_L(1,1), r_11_L(1,2), r_11_L(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Hx(1), Hy(1), Hz(1), r_11_R(1,1), r_11_R(1,2), r_11_R(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ex(1), Ey(1), Ez(1), r_12_L(1,1), r_12_L(1,2), r_12_L(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Lx(1), Ly(1), Lz(1), r_12_R(1,1), r_12_R(1,2), r_12_R(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Nx(1), Ny(1), Nz(1), r_13_L(1,1), r_13_L(1,2), r_13_L(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ox(1), Oy(1), Oz(1), r_13_R(1,1), r_13_R(1,2), r_13_R(1,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
axis equal


figure(16)
plot3(Sx(350,:), Sy(350,:), Sz(350,:), 'o', 'Color', 'b', 'MarkerSize',4,'MarkerFaceColor','#D9FFFF')
hold on
plot3(CGx(:,350), CGy(:,350), CGz(:,350), 'o', 'Color', 'g', 'MarkerSize',4,'MarkerFaceColor','#D9FFFF')
xlabel('x')
ylabel('y')
zlabel('z')
hold on
quiver3(Px(350), Py(350), Pz(350), r_1_L(350,1), r_1_L(350,2), r_1_L(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Bx(350), By(350), Bz(350), r_2_L(350,1), r_2_L(350,2), r_2_L(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Bx(350), By(350), Bz(350), r_3_L(350,1), r_3_L(350,2), r_3_L(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Cx(350), Cy(350), Cz(350), r_4_L(350,1), r_4_L(350,2), r_4_L(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Cx(350), Cy(350), Cz(350), r_5_L(350,1), r_5_L(350,2), r_5_L(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Dx(350), Dy(350), Dz(350), r_6_L(350,1), r_6_L(350,2), r_6_L(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Dx(350), Dy(350), Dz(350), r_7_L(350,1), r_7_L(350,2), r_7_L(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Qx(350), Qy(350), Qz(350), r_1_R(350,1), r_1_R(350,2), r_1_R(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ix(350), Iy(350), Iz(350), r_2_R(350,1), r_2_R(350,2), r_2_R(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ix(350), Iy(350), Iz(350), r_3_R(350,1), r_3_R(350,2), r_3_R(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Jx(350), Jy(350), Jz(350), r_4_R(350,1), r_4_R(350,2), r_4_R(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Jx(350), Jy(350), Jz(350), r_5_R(350,1), r_5_R(350,2), r_5_R(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Kx(350), Ky(350), Kz(350), r_6_R(350,1), r_6_R(350,2), r_6_R(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Kx(350), Ky(350), Kz(350), r_7_R(350,1), r_7_R(350,2), r_7_R(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Mx(350), My(350), Mz(350), r_8(350,1), r_8(350,2), r_8(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Fx(350), Fy(350), Fz(350), r_9(350,1), r_9(350,2), r_9(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Fx(350), Fy(350), Fz(350), r_10(350,1), r_10(350,2), r_10(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ax(350), Ay(350), Az(350), r_11_L(350,1), r_11_L(350,2), r_11_L(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Hx(350), Hy(350), Hz(350), r_11_R(350,1), r_11_R(350,2), r_11_R(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ex(350), Ey(350), Ez(350), r_12_L(350,1), r_12_L(350,2), r_12_L(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Lx(350), Ly(350), Lz(350), r_12_R(350,1), r_12_R(350,2), r_12_R(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Nx(350), Ny(350), Nz(350), r_13_L(350,1), r_13_L(350,2), r_13_L(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ox(350), Oy(350), Oz(350), r_13_R(350,1), r_13_R(350,2), r_13_R(350,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
axis equal

figure(17)
plot3(Sx(250,:), Sy(250,:), Sz(250,:), 'o', 'Color', 'b', 'MarkerSize',4,'MarkerFaceColor','#D9FFFF')
hold on
plot3(CGx(:,250), CGy(:,250), CGz(:,250), 'o', 'Color', 'g', 'MarkerSize',4,'MarkerFaceColor','#D9FFFF')
xlabel('x')
ylabel('y')
zlabel('z')
hold on
quiver3(Px(250), Py(250), Pz(250), r_1_L(250,1), r_1_L(250,2), r_1_L(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Bx(250), By(250), Bz(250), r_2_L(250,1), r_2_L(250,2), r_2_L(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Bx(250), By(250), Bz(250), r_3_L(250,1), r_3_L(250,2), r_3_L(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Cx(250), Cy(250), Cz(250), r_4_L(250,1), r_4_L(250,2), r_4_L(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Cx(250), Cy(250), Cz(250), r_5_L(250,1), r_5_L(250,2), r_5_L(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Dx(250), Dy(250), Dz(250), r_6_L(250,1), r_6_L(250,2), r_6_L(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Dx(250), Dy(250), Dz(250), r_7_L(250,1), r_7_L(250,2), r_7_L(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Qx(250), Qy(250), Qz(250), r_1_R(250,1), r_1_R(250,2), r_1_R(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ix(250), Iy(250), Iz(250), r_2_R(250,1), r_2_R(250,2), r_2_R(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ix(250), Iy(250), Iz(250), r_3_R(250,1), r_3_R(250,2), r_3_R(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Jx(250), Jy(250), Jz(250), r_4_R(250,1), r_4_R(250,2), r_4_R(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Jx(250), Jy(250), Jz(250), r_5_R(250,1), r_5_R(250,2), r_5_R(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Kx(250), Ky(250), Kz(250), r_6_R(250,1), r_6_R(250,2), r_6_R(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Kx(250), Ky(250), Kz(250), r_7_R(250,1), r_7_R(250,2), r_7_R(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Mx(250), My(250), Mz(250), r_8(250,1), r_8(250,2), r_8(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Fx(250), Fy(250), Fz(250), r_9(250,1), r_9(250,2), r_9(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Fx(250), Fy(250), Fz(250), r_10(250,1), r_10(250,2), r_10(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ax(250), Ay(250), Az(250), r_11_L(250,1), r_11_L(250,2), r_11_L(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Hx(250), Hy(250), Hz(250), r_11_R(250,1), r_11_R(250,2), r_11_R(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ex(250), Ey(250), Ez(250), r_12_L(250,1), r_12_L(250,2), r_12_L(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Lx(250), Ly(250), Lz(250), r_12_R(250,1), r_12_R(250,2), r_12_R(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Nx(250), Ny(250), Nz(250), r_13_L(250,1), r_13_L(250,2), r_13_L(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
quiver3(Ox(250), Oy(250), Oz(250), r_13_R(250,1), r_13_R(250,2), r_13_R(250,3), 'LineWidth',3, 'Color','b', 'Marker', '*', 'MarkerSize',6)
axis equal


figure (18)
subplot(3, 2, 1)
plot(t,M_hip_mag_L)
title('Torque-Time plot of Left Hip')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Hip_L)
ylim([0 8e4])

subplot(3,2,2)
plot(t,M_hip_mag_R)
title('Torque-Time plot of Right Hip')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Hip_R)
ylim([0 8e4])

subplot(3,2,3)
plot(t,M_knee_mag_L)
title('Torque-Time plot of Left Knee')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Knee_L)
ylim([0 12e4])

subplot(3,2,4)
plot(t,M_knee_mag_R)
title('Torque-Time plot of Right Knee')
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Knee_R)
ylim([0 12e4])

subplot(3,2,5)
plot(t,M_ankle_mag_L)
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Ankle_L)
title('Torque-Time plot of Left Ankle')
ylim([0 8e4])

subplot(3,2,6)
plot(t, M_ankle_mag_R)
xlabel('Time (s)')
ylabel('Torque (Nmm)')
hold on
plot(t, M_Vicon_Ankle_R)
title('Torque-Time plot of Right Ankle')
ylim([0 8e4])

figure(19)
plot(t,d_knee_L_horizontal)
title('Horizontal distance/moment arm for the left and right knee')
xlabel('Time (s)')
ylabel('Distance (mm)')
hold on
plot(t, d_knee_R_horizontal)
legend('Left','Right')
ylim([0 320])

