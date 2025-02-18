/*
 * imu.c
 *
 *  Created on: Feb 18, 2025
 *      Author: halvard
 * 
 *  Description: imu functions
 */

#ifndef _IMU_H__
#define _IMU_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct imu_t
{
    float g[3];
    float dps[3];
} imu_t;

int imu_init(imu_t *imu);
int imu_process(imu_t *imu);
//void imu_deinit(void); /* Can be ignored as library only sets mpu6050 to sleep on deinit */

#ifdef __cplusplus
}
#endif

#endif /* _IMU_H__ */
