# https://ahrs.readthedocs.io/en/latest/filters/ekf.html
import numpy as np
from ahrs.filters import EKF
from ahrs.common.orientation import acc2q, q2euler, RAD2DEG, q2R
from ahrs.common.frames import enu2ecef
from ahrs.common.quaternion import Quaternion
import serial
import open3d as o3d
import keyboard

# Rotation from visualizer to NED
R_vis_ned = np.array([
    [ 0,  1,  0],
    [ 0,  0, -1],
    [-1,  0,  0]
])

# Setup visualization
vis = o3d.visualization.Visualizer()
vis.create_window(width=800, height=600)

# XYZ = RGB
axis_ned = o3d.geometry.TriangleMesh.create_coordinate_frame(size=1.0, origin=[0,0,0])
axis_ned.rotate(R_vis_ned, center=[0,0,0])
vis.add_geometry(axis_ned)

arrow = o3d.io.read_triangle_mesh("simple_arrow.obj")
arrow.compute_vertex_normals()
arrow.paint_uniform_color([1,0,0])
arrow.translate(-arrow.get_center())
arrow.scale(0.07, center=[0,0,0])
vis.add_geometry(arrow)

# Initialize EKF
var_acc = 1e-8
var_gyr = 1e-1
alpha = 0.1
ekf = EKF(frame='NED', var_acc=var_acc, var_gyr=var_gyr)

# Setup serial connection to IMU
S = serial.Serial('COM3', 115200)

# Initial values
Q = Quaternion().to_array()
acc_last = None

def get_imu_data():
    # Returns acceleration and gyroscope data from IMU in NED frame
    # units of returned data: acc[m/s^2], gyr[rad/s]
    line = S.readline().decode('utf-8').strip()
    columns = line.split(' ')
    if len(columns) != 6:
        return None
    [ax, ay, az, wx, wy, wz] = [float(x) for x in columns]
    acc = np.array([ax, ay, az])
    gyr = np.array([wx, wy, wz]) * np.pi / 180  # Convert to rads
    return acc, gyr

def estimate_gyro_bias(alpha=0.01, num_samples=50):
    # Estimate gyro bias by averaging gyro data
    B = None
    for i in range(num_samples):
        data = get_imu_data()
        if data is None:
            continue
        _ , gyr = data
        # simple low-pass filter(iir)
        B = alpha * gyr + (1 - alpha) * (B if B is not None else gyr)
    return B

B = estimate_gyro_bias()
print(f"Gyro bias: {B}")

while not (keyboard.is_pressed('q') or keyboard.is_pressed('esc')):
    data = get_imu_data()
    if data is None:
        continue

    # Data filter and bias correction
    acc, gyr = data
    acc = alpha * acc + (1 - alpha) * (acc_last if acc_last is not None else acc)
    acc_last = acc
    gyr = gyr - B

    # Confirm that acceleration data is correct(with respect to NED frame)
    # print(f"ax: {acc[0]:4.1f}, ay: {acc[1]:4.1f}, az: {acc[2]:4.1f}") # [0,0,-9.8] is expected when flat

    # Update EKF, data is in NED frame(both sensors and EKF)
    Q = ekf.update(q=Q, gyr=gyr, acc=acc)

    # Confirm that the Quaternion is correct(with respect to NED frame)
    # euler = q2euler(Q) * RAD2DEG
    # print(f"q2euler-> Roll: {euler[0]:4.0f}, Pitch: {euler[1]:4.0f}, Yaw: {euler[2]:4.0f}")

    # Update visualization
    R_ned = q2R(Q) # Rotation in NED frame
    R_vis = R_vis_ned @ q2R(Q) @ R_vis_ned.T # Convert from NED to visualizer frame, rotate, then convert back to NED frame
    arrow.rotate(R_vis, center=[0,0,0]) # Rotate arrow
    vis.update_geometry(arrow)
    vis.poll_events()
    vis.update_renderer()
    arrow.rotate(np.linalg.inv(R_vis), center=[0,0,0]) # Reset rotation

vis.destroy_window()
S.close()