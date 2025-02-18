/*
 * retarget.h
 *
 *  Created on: Jan 28, 2025
 *      Author: halvard
 * 
 *  Description: redirect printf to UART
 */

#ifndef _RETARGET_H__
#define _RETARGET_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stm32f1xx_hal.h>

void RetargetInit(UART_HandleTypeDef *huart);

#ifdef __cplusplus
}
#endif

#endif /* _RETARGET_H__ */
