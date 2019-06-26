//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_CONFIG_H
#define TRACER_FINAL_CONFIG_H
const int MAX_DEPTH = 10;
const int h = 1024, w = 1536;
const int ITER = 20;

//const int SAMPLE = 1024;
int SAMPLE = 40;
int SPP_NUM = SAMPLE * 4;

const int PHOTON_SAMPLE = 100000;
const int PHOTON_GLOBAL  = 50000000;

const int    PHOTON_GLOBAL_N = 64;
const double PHOTON_GLOBAL_R = 32.0;
const int    PHOTON_CAUSTIC_N = 64;
const double PHOTON_CAUSTIC_R = 0.5;
const double PHOTON_GLOBAL_MUL = 1.0;

#endif //TRACER_FINAL_CONFIG_H
