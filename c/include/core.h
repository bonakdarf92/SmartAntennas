//
// Created by faridb on 26.12.20.
//

#ifndef SMART_ANTENNAS_CORE_H
#define SMART_ANTENNAS_CORE_H
#include "osqp/osqp.h"



int init(OSQPWorkspace *work, OSQPSettings  *settings, OSQPData * data);

void run(OSQPWorkspace * w);


#endif //SMART_ANTENNAS_CORE_H
