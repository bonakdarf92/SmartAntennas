//
// Created by Farid Bonakdar on 26.12.20.
//
#include "core.h"

c_float P_x[3] = {4.0, 1.0, 2.0, };
const c_int P_nnz = 3;
c_int P_i[3] = {0, 0, 1, };
c_int P_p[3] = {0, 1, 3, };
const c_float q[2] = {1.0, 1.0, };
c_float A_x[4] = {1.0, 1.0, 1.0, 1.0, };
const c_int A_nnz = 4;
c_int A_i[4] = {0, 1, 0, 2, };
c_int A_p[3] = {0, 2, 4, };
const c_float l[3] = {1.0, 0.0, 0.0, };
const c_float u[3] = {1.0, 0.7, 0.7, };
const c_int n = 2;
const c_int m = 3;

// Exitflag
c_int exitflag;// = 0;


int init(OSQPWorkspace *work, OSQPSettings  *settings, OSQPData * data){
    // Workspace structures

    //OSQPData      *data     = (OSQPData *)c_malloc(sizeof(OSQPData));

    // Populate data
    if (data) {
        data->n = n;
        data->m = m;
        data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
        data->q = q;
        data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
        data->l = l;
        data->u = u;
    }
    // Define solver settings as default
    if (settings) {
        osqp_set_default_settings(settings);
        settings->alpha = 1.0; // Change alpha parameter
    }

    // Setup workspace
    exitflag = osqp_setup(&work, data, settings);
    return exitflag;
}


void run(OSQPWorkspace * w) {
    osqp_solve(w);
    //return 0;
}