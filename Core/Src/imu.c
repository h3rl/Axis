/*
 * imu.c
 *
 *  Created on: Feb 18, 2025
 *      Author: halvard
 * 
 *  Description: imu functions
 */

#include "imu.h"
#include "driver_mpu6050_basic.h"
#include "util.h"

int imu_init(imu_t *imu)
{
    zeromem(imu, sizeof(imu_t));

    if(mpu6050_basic_init(MPU6050_ADDRESS_AD0_LOW) != 0)
    {
        print("MPU6050 init failed!\r\n");
        return 1;
    }
    print("MPU6050 ok\r\n");
    return 0;
}

int imu_process(imu_t *imu)
{
    float* g = imu->g;
    float* dps = imu->dps;
    if(mpu6050_basic_read(g, dps) != 0)
    {
        print("MPU6050 read failed!\r\n");
        return 1;
    }

    return 0;
}

void imu_deinit(void)
{
    mpu6050_basic_deinit();
}