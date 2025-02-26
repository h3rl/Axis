<!-- Logo Section  -->
<div align="center">
    <h1>Axis</h1>
    <img src="./assets/logo.svg">
</div>

<div align="center">
    <p>STM32-based IMU orientation estimation with Kalman filtering</p>
</div>

## About
Axis is an experimental project exploring IMU sensor filtering and orientation estimation. It uses an MPU6050 IMU sensor connected to an STM32F103 microcontroller to provide real-time attitude determination through Kalman filtering. This project was created as a learning exercise to explore IMU sensor fusion and filtering techniques.

<!-- 
## Demo
Comming soon?

<div align="center">
    <img alt="demo" src="demo.gif">
</div>
-->

## Features
- STM32F103RBTx microcontroller
- MPU6050 IMU integration
- Extended Kalman filter(EKF) for orientation estimation
- Built with CMake and STM32CubeMX

## Getting started

1. Clone the repository:

```bash
    git clone https://github.com/h3rl/Axis
    cd Axis
    code .
```
2. Install recommended VS Code extensions when prompted

## Project Roadmap
- [x] Basic IMU readings implementation
- [x] Extended Kalman filter python POC
- [ ] Extended Kalman filter c implementation
- [ ] DMP integration and comparison

## Credits

- [RcdMathLib](https://git.imp.fu-berlin.de/zkasmi/RcdMathLib/-/tree/master) - Matrix math and vector operations library
- [libdriver/mpu6050](https://github.com/libdriver/mpu6050) - MPU6050 sensor driver implementation
