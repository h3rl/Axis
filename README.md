<!-- Logo Section  -->
<div align="center">
    <svg width="290" height="290" viewBox="-50 -50 100 100" xmlns="http://www.w3.org/2000/svg">
        <defs>
            <marker id="arrowhead" markerWidth="10" markerHeight="7" refX="6" refY="3.5" orient="auto">
                <polygon points="0 0, 10 3.5, 0 7" fill="black"/>
            </marker>
        </defs>
        <line x1="0" y1="0" x2="40" y2="0" stroke="black" stroke-width="3" marker-end="url(#arrowhead)" transform="rotate(10)"/>
        <line x1="0" y1="0" x2="40" y2="0" stroke="black" stroke-width="3" marker-end="url(#arrowhead)" transform="rotate(130)"/>
        <line x1="0" y1="0" x2="40" y2="0" stroke="black" stroke-width="3" marker-end="url(#arrowhead)" transform="rotate(250)"/>
    </svg>
</div>

<div align="center">
    <p style="font-size:3rem; font-weight: bold;">Axis</p>
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
- Kalman filter for orientation estimation
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
<input checked="" disabled="" type="checkbox"> Basic IMU readings implementation<br>
<input disabled="" type="checkbox"> Kalman filter orientation estimation<br>
<input disabled="" type="checkbox"> DMP integration and comparison

## Credits

- [RcdMathLib](https://git.imp.fu-berlin.de/zkasmi/RcdMathLib/-/tree/master) - Matrix math and vector operations library
- [libdriver/mpu6050](https://github.com/libdriver/mpu6050) - MPU6050 sensor driver implementation
