/*
 * Copyright (C) 2014 Beat Küng <beat-kueng@gmx.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 */

#ifndef _TIMER_H_
#define _TIMER_H_

#include <stdint.h>

/*!
  The function returns the number of ticks since the certain event (e.g. when the machine was turned on).
  It can be used to initialize rand or to measure a function execution time by reading the tick count
  before and after the function call. The granularity of ticks depends on the hardware and OS used. Use
  getTickFrequency() to convert ticks to seconds.
*/
int64_t getTickCount(void);

/*!
  Returns the number of ticks per seconds.

  The function returns the number of ticks (as returned by getTickCount()) per second.
  The following code computes the execution time in milliseconds:

  \code
  double exec_time = (double)getTickCount();
  // do something ...
  exec_time = ((double)getTickCount() - exec_time)*1000./getTickFrequency();
  \endcode
*/
double getTickFrequency(void);

/*!
 * convert tick difference to seconds
 */
double getTickSeconds(int64_t dt);

#endif /* _TIMER_H_ */
