//// Written by Johan De Bock <johan.debock@gmail.com>
////
//// memory requirements : td_nb * (Rw * Rh * 4 + bb_bail * 24) + (Rlrmax[0] + Rlrmax[1] + Rlrmax[2] + 3) * 4 + (Ww * Wh + Tw * Th) * 12
////   R2000x2000 bb1000 th3 :    48 MB
////   R4000x4000 bb1000 th3 :   192 MB
////   R8000x8000 bb1000 th3 :   768 MB
//// R16000x16000 bb1000 th3 :  3072 MB
//// R32000x32000 bb1000 th3 : 12288 MB
//// R40000x40000 bb1000 th3 : 19200 MB
//// R50000x50000 bb1000 th3 : 30000 MB tested
#define __USE_MINGW_ANSI_STDIO 1
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <omp.h>

#include "dSFMT.h"
#include "png.h"

#ifdef WINDOWS
#include <windows.h>
#define wait_ms(x) Sleep(x)
#endif
#ifdef LINUX
#include <unistd.h>
#define wait_ms(x) usleep((x) * 1000)
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//// daemon
int daemon_mode = 0;

//// strings
char filename[512];
char dirname[512];
char commandname[512];

//// keep total number of paths plotted for a render
long long unsigned int Ppsum = 0;

//// auto write png functionality
long long unsigned int Ppsum_autoPNG_last = 0; // Ppsum of last written auto png
long long unsigned int Ppsum_autoPNG_delta = 1e9; // Ppsum difference between each auto png write, default 1e9

//// td threads
#define TD_MAX 18 // maximum number of calc threads
int td_nb = 3; // number of calc threads
int td_stop = 0; // signal td_stop to all threads
int td_pause = 0; // signal td_pause to all calc threads
int td_paused[TD_MAX]; // which calc threads are paused

//// lr layer
//// suppose 3 calc threads
//// layer mode 0 : sum 3 threads -> 1 index in 1 RGB/Grey color table -> pixels RGB/Grey
//// layer mode 1 : sum 3 threads -> 1 index in 3 Grey color tables -> pixels RGB
//// layer mode 2 : 3 threads -> 3 indexes in 3 Grey color tables -> pixels RGB
int lr_mode; // layer mode
#define LR_MODE_NB 3 // number of layer modes
#define LR_NB 3 // number of layers

//// C seeding complex rectangle
#define C_BATCH_NB 100 // number of c to try in one batch, without interruption
double Cr_lo = -2.0, Cr_up = 2.0; // real lower and upper bound of rectangle C in complex plane
double Ci_lo = 0.0, Ci_up = 2.0; // imaginary lower and upper bound of rectangle C in complex plane : only positive : symmetric
double Cr_ra = Cr_up - Cr_lo; // real range of rectangle C in complex plane
double Ci_ra = Ci_up - Ci_lo; // imaginary range of rectangle C in complex plane

//// P path (in complex plane), Mandelbrot
//// P[n] = z_{n+1} = z_n2+c
//// P = [c, c2+c, (c2+c)2+c, ...]
typedef struct {
    double r; // real
    double i; // imaginary
} complex; // complex number
complex* P[TD_MAX]; // 1 path , per thread
unsigned int* RPo[TD_MAX]; // offsets of a path in render , per thread

long long unsigned int Pp[TD_MAX]; // number of paths plotted , per thread

//// bb BuddhaBrot types and parameters
//// BuddhaBrot type 0 : the Buddhabrot
//// BuddhaBrot type 1 : the Anti-Buddhabrot
//// BuddhaBrot type 2 : the Anti-Buddhabrot with some lobes cut
//// Mandelbrot set = Mset = set of c for which : z_n bounded : lim_{n->inf} z_n != inf
//// notMset = set of c for which : z_n unbounded : lim_{n->inf} z_n == inf
////
//// Buddhabrot : notMset : unbounded
//// plot the unbounded paths = P[n_inf]2 > 4 for n_inf < bailout
//// -> n_inf = core_mandelbrot(c)
//// skip the certainly bounded paths
////
//// Anti-Buddhabrot : Mset : bounded
//// plot the bounded paths = P[n_inf]2 <= 4 for n_inf == bailout
//// -> n_inf = core_mandelbrot(c)
int bb_type[TD_MAX]; // BuddhaBrot type , per thread
#define BB_TYPE_NB 3 // number of BuddhaBrot types
int bb_bail[TD_MAX]; // bailout , per thread
int bb_pps[TD_MAX]; // path plot start : start plotting for n = bb_pps , per thread
int bb_ppe[TD_MAX]; // path plot end : end plotting for n = n_inf - bb_ppe , per thread
int bb_minn[TD_MAX]; // path minimum n_inf : plot path only if n_inf >= bb_minn , per thread

//// R render (count matrix (of pixels))
//// registers the number of times a path passed in a pixel
//// mapped in a complex rectangle
unsigned int* R[TD_MAX]; // 1 render , per thread
int Rw = 0; // render width
int Rh = 0; // render height

unsigned int Rlrmax[LR_NB]; // maximum count sum , per layer

double Rr_lo, Rr_up; // real lower and upper bound of rectangle R in complex plane
double Rr_ra; // real range of rectangle R in complex plane
double Ri_lo, Ri_up; // imaginary lower and upper bound of rectangle R in complex plane
double Ri_ra; // imaginary range of rectangle R in complex plane

double Rhdivr; // Rhdivr = Rh / Rr_ra
double Rwdivi; // Rwdivi = Rw / Ri_ra

double Rr_f = 0.5; // factor to change Rr_ra
double Ri_f = 0.5; // factor to change Ri_ra TODO
double Rw_f = 0.5; // factor to change render width
double Rh_f = 0.5; // factor to change render height

//// H histogram
//// a histogram of the values in a render, bins 0, 1, 2, ..., Rlrmax
unsigned int* H[LR_NB]; // 1 histogram , per layer
unsigned int Hl[LR_NB] = {0, 0, 0}; // histogram length > Rlrmax , per layer

//// W window (count matrix (of pixels))
//// for the pixels shown in the BuddhaBrot window
//// it is located somewhere in the render
unsigned int* W[LR_NB]; // 1 window , per layer
int Ww = 1000; // window width
int Wh = 1000; // window height
int RWow = 0; // x offset of window in render
int RWoh = 0; // y offset of window in render

double Wr_lo, Wr_up; // real lower and upper bound of rectangle W in complex plane
double Wi_lo, Wi_up; // imaginary lower and upper bound of rectangle W in complex plane

int recalc_WH_if_paused = 0; // indicate that H and W need to be recalculated when changing layer mode/coloring method/panning window when paused

//// T tile (count matrix (of pixels))
//// for the pixels outputted to a png image
//// it is located somewhere in the render
unsigned int* T[LR_NB]; // 1 tile , per layer
int Tw = 4000; // tile width
int Th = 4000; // tile height

//// cm coloring methods
//// coloring method 0 : rank-order mapping : each unique sum get a unique color : independent of sum frequency
//// coloring method 1 : histogram equalization : frequent sums get spread more in color table : ideally resulting in flat histogram or each color evenly represented
//// coloring method 2 : linear
//// coloring method 0 : problematic : when lots of low freq unique sums : 25% of unique -> spread over 25% of color table
//// coloring method 1 : problematic : when a few high freq sums : 25% of pixels = sum -> 25% jump in color table
//// coloring method 0 : effect on Buddhabrot 0->255 : more detail in mid sums : low sums lost in black
//// coloring method 1 : effect on Buddhabrot 0->255 : more detail in low sums : high sums lost in white
//// coloring method 0 : if cm0n > 1 : i = e * Rank[sum] / (cm0n - 1)
//// coloring method 0 : else          i = 0
//// coloring method 1 : if cm1n > 0 : i = e * CumulativeH[sum] / cm1n :
//// coloring method 1 : else          i = 0                               TODO check e * () in code optim
int cm[LR_NB]; // coloring method , per layer
#define CM_NB 3 // number of coloring methods
unsigned int cm0n[LR_NB]; // normalization value for coloring method 0 = number of filled bins in histogram = number of unique values in render , per layer
unsigned int cm1n[LR_NB]; // normalization value for coloring method 1 = (number of pixels in render - number of 0 pixels in render) , per layer

//// csf coloring sum function
//// sum function 0 : none
//// sum function 1 : log
//// sum function 0 : sum -> sum
//// sum function 1 : if sum >= (unsigned int)csfp1 : sum -> csfp1 + log(sum - csfp1 + 1)
//// sum function 1 : else              sum -> sum
int csf[LR_NB]; // coloring sum function , per layer
#define CSF_NB 2 // number of coloring sum functions
int csfp1[LR_NB]; // parameter 1 , per layer

//// CT color tables
//// index -> Grey value
unsigned char* CT[LR_NB]; // 1 color table , per layer
int CTe[LR_NB]; // color table end (= length-1) , per layer
int ct_type[LR_NB]; // 1 color table type, per layer
#define CT_TYPE_NB 6 // number of color table types

double ct_f[LR_NB] = {1.0, 1.0, 1.0}; // scale index in color table with this factor (with clipping) i = MIN(e, ct_f * i) , per layer
int ct_o[LR_NB]; // start offset in color table (0 index gets mapped on this color) , per layer
int ct_v[LR_NB] = {0, 0, 0}; // cycle speed : diff per frame in index , per layer

//// plc permutate layer color functionality
int plc = 0;
#define PLC_NB 6




int ct_cycle(int i, int e)
{
    return ((i > e) ? (i - e - 1) : ((i < 0) ? (e + i + 1) : i));
}

void ct_load(int load_selectedlayer, int load_ct_type1, int load_ct_type2, int load_ct_type3)
{
    ct_type[0] = load_ct_type1;
    ct_type[1] = load_ct_type2;
    ct_type[2] = load_ct_type3;
    int layer_start = 0;
    int layer_end = LR_NB - 1;

    if (load_selectedlayer != -1) {
        layer_start = load_selectedlayer;
        layer_end = load_selectedlayer;
    }

    for (int layer_iter = layer_start; layer_iter <= layer_end; layer_iter++) {
        if (ct_type[layer_iter] == 0) {
            CTe[layer_iter] = 255;
            free(CT[layer_iter]);
            CT[layer_iter] = (unsigned char*)calloc(CTe[layer_iter] + 1, sizeof(unsigned char));

            for (int i = 0; i <= 255; i++) {
                CT[layer_iter][i] = i;
            }
        }

        if (ct_type[layer_iter] == 1) {
            CTe[layer_iter] = 255;
            free(CT[layer_iter]);
            CT[layer_iter] = (unsigned char*)calloc(CTe[layer_iter] + 1, sizeof(unsigned char));

            for (int i = 0, j = 255; j >= 0; i++, j--) {
                CT[layer_iter][i] = j;
            }
        }

        if (ct_type[layer_iter] == 2) {
            CTe[layer_iter] = 509;
            free(CT[layer_iter]);
            CT[layer_iter] = (unsigned char*)calloc(CTe[layer_iter] + 1, sizeof(unsigned char));

            for (int i = 0; i < 255; i++) {
                CT[layer_iter][i] = i;
            }

            for (int i = 255, j = 255; j > 0; i++, j--) {
                CT[layer_iter][i] = j;
            }
        }

        if (ct_type[layer_iter] == 3) {
            CTe[layer_iter] = 509;
            free(CT[layer_iter]);
            CT[layer_iter] = (unsigned char*)calloc(CTe[layer_iter] + 1, sizeof(unsigned char));

            for (int i = 0, j = 255; j > 0; i++, j--) {
                CT[layer_iter][i] = j;
            }

            for (int i = 255, j = 0; j < 255; i++, j++) {
                CT[layer_iter][i] = j;
            }
        }

        if (ct_type[layer_iter] == 4) {
            CTe[layer_iter] = 255;
            free(CT[layer_iter]);
            CT[layer_iter] = (unsigned char*)calloc(CTe[layer_iter] + 1, sizeof(unsigned char));

            for (int i = 0; i <= 255; i++) {
                CT[layer_iter][i] = 0;
            }
        }

        if (ct_type[layer_iter] == 5) {
            CTe[layer_iter] = 255;
            free(CT[layer_iter]);
            CT[layer_iter] = (unsigned char*)calloc(CTe[layer_iter] + 1, sizeof(unsigned char));

            for (int i = 0; i <= 255; i++) {
                CT[layer_iter][i] = 255;
            }
        }
    }
}




void pause_calcthreads_and_wait()
{
    td_pause = 1;
    #pragma omp flush(td_pause)

    while (1) {
        int td_paused_nb = 0;

        for (int td_i = 0; td_i < td_nb; td_i += 1) {
            td_paused_nb += td_paused[td_i];
        }

        if (td_paused_nb == td_nb) {
            break;
        }

        wait_ms(50);
    }
}

void reset_R(int td_i)
{
    Pp[td_i] = 0;
    memset(R[td_i], 0, (unsigned int)Rw * Rh * sizeof(unsigned int));
    Hl[0] = 0;
    Hl[1] = 0;
    Hl[2] = 0;
    Ppsum_autoPNG_last = 0;
}

void realloc_R(int td_i)
{
    Pp[td_i] = 0;

    if (R[td_i] != NULL) {
        free(R[td_i]);
        R[td_i] = NULL;
    }

    R[td_i] = (unsigned int*)calloc((unsigned int)Rw * Rh, sizeof(unsigned int));
    Hl[0] = 0;
    Hl[1] = 0;
    Hl[2] = 0;
    Ppsum_autoPNG_last = 0;
}

void realloc_P(int td_i)
{
    if (P[td_i] != NULL) {
        free(P[td_i]);
        P[td_i] = NULL;
    }

    P[td_i] = (complex*)calloc(bb_bail[td_i], sizeof(complex));

    if (RPo[td_i] != NULL) {
        free(RPo[td_i]);
        RPo[td_i] = NULL;
    }

    RPo[td_i] = (unsigned int*)calloc(2 * bb_bail[td_i], sizeof(unsigned int));
}

void free_calcthread(int td_i)
{
    Pp[td_i] = 0;

    if (R[td_i] != NULL) {
        free(R[td_i]);
        R[td_i] = NULL;
    }

    if (P[td_i] != NULL) {
        free(P[td_i]);
        P[td_i] = NULL;
    }

    if (RPo[td_i] != NULL) {
        free(RPo[td_i]);
        RPo[td_i] = NULL;
    }
}

void realloc_calcthread(int td_i)
{
    realloc_R(td_i);
    realloc_P(td_i);
}

void decrease_num_calcthreads()
{
    pause_calcthreads_and_wait();
    int new_num_calcthreads = MAX(td_nb - 3, 3);

    for (int td_i = new_num_calcthreads; td_i < td_nb; td_i++) {
        for (unsigned int Ri = 0; Ri < (unsigned int)Rw * Rh; Ri++) {
            R[td_i % 3][Ri] += R[td_i][Ri];
        }

        Pp[td_i % 3] += Pp[td_i];
        free_calcthread(td_i);
    }

    td_nb = new_num_calcthreads;
    td_pause = 0;
    #pragma omp flush(td_pause, td_nb)
}

void increase_num_calcthreads()
{
    int new_num_calcthreads = MIN(td_nb + 3, TD_MAX);

    for (int td_i = td_nb; td_i < new_num_calcthreads; td_i++) {
        bb_type[td_i] = bb_type[td_i - 3];
        bb_bail[td_i] = bb_bail[td_i - 3];
        bb_pps[td_i] = bb_pps[td_i - 3];
        bb_ppe[td_i] = bb_ppe[td_i - 3];
        bb_minn[td_i] = bb_minn[td_i - 3];
        realloc_calcthread(td_i);
    }

    td_nb = new_num_calcthreads;
    #pragma omp flush(td_nb)
}

void load_location(int pause_calcthreads, double load_zoom, double load_centerx, double load_centery, int load_renderw, int load_renderh)
{
    if (pause_calcthreads) {
        pause_calcthreads_and_wait();
    }

    double rangediv2x = 2.0 / load_zoom;
    double rangediv2y = 2.0 / load_zoom;
    Rr_lo = load_centerx - rangediv2x;
    Rr_up = load_centerx + rangediv2x;
    Ri_lo = load_centery - rangediv2y;
    Ri_up = load_centery + rangediv2y;
    Rr_ra = Rr_up - Rr_lo;
    Ri_ra = Ri_up - Ri_lo;

    if (Rw == load_renderw && Rh == load_renderh) {
        for (int td_i = 0; td_i < td_nb; td_i += 1) {
            reset_R(td_i);
        }
    } else {
        if (Rw == 0) {
            RWow = 0.5 * (load_renderw - Ww);
        } else if (Rw < load_renderw) {
            RWow = RWow / Rw_f + 0.5 * Ww * (1.0 / Rw_f - 1.0);
        } else {
            RWow = RWow * Rw_f + 0.5 * Ww * (Rw_f - 1.0);
        }

        if (Rh == 0) {
            RWoh = 0.5 * (load_renderh - Wh);
        } else if (Rh < load_renderh) {
            RWoh = RWoh / Rh_f + 0.5 * Wh * (1.0 / Rh_f - 1.0);
        } else {
            RWoh = RWoh * Rh_f + 0.5 * Wh * (Rh_f - 1.0);
        }

        Rw = load_renderw;
        Rh = load_renderh;

        if (Rw < Ww) {
            Rw = Ww;
            RWow = 0;
        } else if (RWow < 0) {
            RWow = 0;
        } else if (RWow + Ww > Rw) {
            RWow = Rw - Ww;
        }

        if (Rh < Wh) {
            Rh = Wh;
            RWoh = 0;
        } else if (RWoh < 0) {
            RWoh = 0;
        } else if (RWoh + Wh > Rh) {
            RWoh = Rh - Wh;
        }

        for (int td_i = 0; td_i < td_nb; td_i += 1) {
            realloc_R(td_i);
        }
    }

    Rhdivr = Rh / Rr_ra;
    Rwdivi = Rw / Ri_ra;
    Wr_lo = Rr_lo + RWoh / Rhdivr;
    Wr_up = Wr_lo + Wh / Rhdivr;
    Wi_lo = Ri_lo + RWow / Rwdivi;
    Wi_up = Wi_lo + Ww / Rwdivi;

    if (pause_calcthreads) {
        td_pause = 0;
        #pragma omp flush(td_pause)
    }
}

void load_bb_param(int load_selectedlayer, int load_bb_type, int load_bb_bail, int load_bb_pps, int load_bb_ppe, int load_bb_minn)
{
    pause_calcthreads_and_wait();
    load_bb_bail = MAX(load_bb_bail, 0);
    load_bb_pps = MIN(MAX(load_bb_pps, 0), load_bb_bail);
    load_bb_ppe = MIN(MAX(load_bb_ppe, 0), load_bb_bail);
    load_bb_minn = MIN(MAX(load_bb_minn, 0), load_bb_bail);
    int layer_start = 0;
    int layer_step = 1;

    if (load_selectedlayer != -1) {
        layer_start = load_selectedlayer;
        layer_step = LR_NB;
    }

    for (int td_i = layer_start; td_i < td_nb; td_i += layer_step) {
        bb_type[td_i] = load_bb_type;

        if (bb_bail[td_i] != load_bb_bail) {
            bb_bail[td_i] = load_bb_bail;
            realloc_P(td_i);
        }

        bb_pps[td_i] = load_bb_pps;
        bb_ppe[td_i] = load_bb_ppe;
        bb_minn[td_i] = load_bb_minn;
        reset_R(td_i);
    }

    td_pause = 0;
    #pragma omp flush(td_pause)
}

void load_location_bb_color_param(int pause_calcthreads, double load_zoom, double load_centerx, double load_centery, int load_renderw, int load_renderh, int load_lr_mode, int load_bb_type1, int load_bb_bail1, int load_bb_pps1, int load_bb_ppe1, int load_bb_minn1, int load_bb_type2, int load_bb_bail2, int load_bb_pps2, int load_bb_ppe2, int load_bb_minn2, int load_bb_type3, int load_bb_bail3, int load_bb_pps3, int load_bb_ppe3, int load_bb_minn3, int load_ct_type1, int load_ct_type2, int load_ct_type3, int load_cm1, int load_csf1, int load_csfp11, double load_ct_f1, int load_ct_o1, int load_cm2, int load_csf2, int load_csfp12, double load_ct_f2, int load_ct_o2, int load_cm3, int load_csf3, int load_csfp13, double load_ct_f3, int load_ct_o3)
{
    if (pause_calcthreads) {
        pause_calcthreads_and_wait();
    }

    load_bb_bail1 = MAX(load_bb_bail1, 0);
    load_bb_bail2 = MAX(load_bb_bail2, 0);
    load_bb_bail3 = MAX(load_bb_bail3, 0);
    load_bb_pps1 = MIN(MAX(load_bb_pps1, 0), load_bb_bail1);
    load_bb_pps2 = MIN(MAX(load_bb_pps2, 0), load_bb_bail2);
    load_bb_pps3 = MIN(MAX(load_bb_pps3, 0), load_bb_bail3);
    load_bb_ppe1 = MIN(MAX(load_bb_ppe1, 0), load_bb_bail1);
    load_bb_ppe2 = MIN(MAX(load_bb_ppe2, 0), load_bb_bail2);
    load_bb_ppe3 = MIN(MAX(load_bb_ppe3, 0), load_bb_bail3);
    load_bb_minn1 = MIN(MAX(load_bb_minn1, 0), load_bb_bail1);
    load_bb_minn2 = MIN(MAX(load_bb_minn2, 0), load_bb_bail2);
    load_bb_minn3 = MIN(MAX(load_bb_minn3, 0), load_bb_bail3);
    lr_mode = load_lr_mode;
    ct_load(-1, load_ct_type1, load_ct_type2, load_ct_type3);
    cm[0] = load_cm1;
    cm[1] = load_cm2;
    cm[2] = load_cm3;
    csf[0] = load_csf1;
    csf[1] = load_csf2;
    csf[2] = load_csf3;
    csfp1[0] = load_csfp11;
    csfp1[1] = load_csfp12;
    csfp1[2] = load_csfp13;
    ct_f[0] = load_ct_f1;
    ct_f[1] = load_ct_f2;
    ct_f[2] = load_ct_f3;
    ct_o[0] = load_ct_o1;
    ct_o[1] = load_ct_o2;
    ct_o[2] = load_ct_o3;

    for (int td_i = 0; td_i < td_nb; td_i += LR_NB) {
        bb_type[td_i] = load_bb_type1;

        if (bb_bail[td_i] != load_bb_bail1) {
            bb_bail[td_i] = load_bb_bail1;
            realloc_P(td_i);
        }

        bb_pps[td_i] = load_bb_pps1;
        bb_ppe[td_i] = load_bb_ppe1;
        bb_minn[td_i] = load_bb_minn1;
    }

    for (int td_i = 1; td_i < td_nb; td_i += LR_NB) {
        bb_type[td_i] = load_bb_type2;

        if (bb_bail[td_i] != load_bb_bail2) {
            bb_bail[td_i] = load_bb_bail2;
            realloc_P(td_i);
        }

        bb_pps[td_i] = load_bb_pps2;
        bb_ppe[td_i] = load_bb_ppe2;
        bb_minn[td_i] = load_bb_minn2;
    }

    for (int td_i = 2; td_i < td_nb; td_i += LR_NB) {
        bb_type[td_i] = load_bb_type3;

        if (bb_bail[td_i] != load_bb_bail3) {
            bb_bail[td_i] = load_bb_bail3;
            realloc_P(td_i);
        }

        bb_pps[td_i] = load_bb_pps3;
        bb_ppe[td_i] = load_bb_ppe3;
        bb_minn[td_i] = load_bb_minn3;
    }

    load_location(0, load_zoom, load_centerx, load_centery, load_renderw, load_renderh);

    if (pause_calcthreads) {
        td_pause = 0;
        #pragma omp flush(td_pause)
    }
}

int core_mandelbrot(int td_i, complex c)
{
    int n;
    complex z = {0.0, 0.0};

    for (n = 0; n < bb_bail[td_i]; n++) {
        P[td_i][n].r = z.r * z.r - z.i * z.i + c.r;
        P[td_i][n].i = z.r * z.i + z.r * z.i + c.i;
        z = P[td_i][n];

        if (z.r * z.r + z.i * z.i > 4.0) {
            break;
        }
    }

    return (n);
}

void calculation_thread(int td_i)
{
    while (1) {
        #pragma omp flush(td_pause)

        if (td_pause == 1) {
            td_paused[td_i] = 1;

            while (1) {
                #pragma omp flush(td_pause, td_stop)

                if (td_pause == 0 || td_stop == 1) {
                    break;
                }

                wait_ms(1000);
            }

            td_paused[td_i] = 0;
        }

        while (1) {
            #pragma omp flush(td_nb, td_stop)

            if (td_i < td_nb || td_stop == 1) {
                break;
            }

            wait_ms(1000);
        }

        #pragma omp flush(td_stop)

        if (td_stop == 1) {
            break;
        }

        if (bb_type[td_i] == 0) {
            for (int c_iter = 0; c_iter < C_BATCH_NB; c_iter++) {
                complex c;

                while (1) {
                    c.r = Cr_lo + Cr_ra * dsfmt_gv_genrand_close_open();
                    c.i = Ci_lo + Ci_ra * dsfmt_gv_genrand_close_open();
                    double c2 = c.r * c.r + c.i * c.i;

                    //// skip the certainly bounded paths
                    //// skip period 1 cardioid period 2 circle period 3 circle period 4 circle
                    //// period 1 cardioid =   exact = 3/32 = 0.09375
                    //// http://www.mrob.com/pub/muency/r2a.html
                    //// period 2 circle   =   exact = -15/16 = -0.9375
                    //// http://www.mrob.com/pub/muency/r212a.html
                    //// period 3 circle   = inexact = estimated by slightly smaller circle c (-0.125,0.744) r 0.092
                    //// http://www.mrob.com/pub/muency/r213a.html
                    //// http://www.mrob.com/pub/muency/muatomsizes.html
                    //// < -0.560697 --> corrected by testing (bail = 1000000, see fix! below) ever larger circles to < -0.5603 --> r 0.0941328847
                    //// period 4 circle   = inexact = estimated by slightly smaller circle c (-1.308,0) r 0.058
                    if (!(c2 * (8.0 * c2 - 3.0) + c.r < 0.09375 || c2 + c.r + c.r < -0.9375 || c2 + 0.25 * c.r - 1.488 * c.i < -0.5603 || c2 + 2.616 * c.r < -1.7075)) {
                        break;
                    } else if (0) {
                        //// certainly bounded
                        int n_inf = core_mandelbrot(td_i, c);

                        if (n_inf < bb_bail[td_i]) {
                            printf("\r                                                                                ");
                            printf("\rfix! notMset %f %f\n", c.r, c.i);
                            fflush(stdout);
                        }
                    }
                }

                //// plot the unbounded paths = P[n_inf]2 > 4 for n_inf < bailout
                int n_inf = core_mandelbrot(td_i, c);

                if (n_inf < bb_bail[td_i] && n_inf >= bb_minn[td_i]) {
                    int minimum_one_point_inside_window = 0;
                    int offset_count = 0;

                    //// |P[td_i][n_inf]| > 2.0 but is handy for getting rid of radius 2.0 circle in plot
                    for (int xy_p = bb_pps[td_i]; xy_p <= n_inf - bb_ppe[td_i]; xy_p++) {
                        if ((P[td_i][xy_p].r >= Rr_up) || (P[td_i][xy_p].r < Rr_lo)) {
                            continue;
                        }

                        if ((P[td_i][xy_p].i >= Ri_lo) && (P[td_i][xy_p].i < Ri_up)) {
                            minimum_one_point_inside_window = 1;
                            int ix = Rhdivr * (P[td_i][xy_p].r - Rr_lo);
                            int iy = Rwdivi * (P[td_i][xy_p].i - Ri_lo);
                            RPo[td_i][offset_count++] = (unsigned int)ix * Rw + iy;

                            if ((-P[td_i][xy_p].i >= Ri_lo) && (-P[td_i][xy_p].i < Ri_up)) {
                                int iyc = Rwdivi * (-P[td_i][xy_p].i - Ri_lo);
                                RPo[td_i][offset_count++] = (unsigned int)ix * Rw + iyc;
                            }
                        } else {
                            if ((-P[td_i][xy_p].i >= Ri_lo) && (-P[td_i][xy_p].i < Ri_up)) {
                                minimum_one_point_inside_window = 1;
                                int ix = Rhdivr * (P[td_i][xy_p].r - Rr_lo);
                                int iyc = Rwdivi * (-P[td_i][xy_p].i - Ri_lo);
                                RPo[td_i][offset_count++] = (unsigned int)ix * Rw + iyc;
                            }
                        }
                    }

                    if (minimum_one_point_inside_window) {
                        Pp[td_i] += 2;
                    }

                    for (int offset_iter = 0; offset_iter < offset_count; offset_iter++) {
                        R[td_i][RPo[td_i][offset_iter]]++;
                    }
                }
            }
        }

        if (bb_type[td_i] == 1) {
            for (int c_iter = 0; c_iter < C_BATCH_NB; c_iter++) {
                complex c;
                c.r = Cr_lo + Cr_ra * dsfmt_gv_genrand_close_open();
                c.i = Ci_lo + Ci_ra * dsfmt_gv_genrand_close_open();
                int n_inf = core_mandelbrot(td_i, c);

                //// plot the bounded paths = P[n_inf]2 <= 4 for n_inf == bailout
                if (n_inf == bb_bail[td_i]) {
                    int minimum_one_point_inside_window = 0;
                    int offset_count = 0;

                    for (int xy_p = bb_pps[td_i]; xy_p < n_inf - bb_ppe[td_i]; xy_p++) {
                        if ((P[td_i][xy_p].r >= Rr_up) || (P[td_i][xy_p].r < Rr_lo)) {
                            continue;
                        }

                        if ((P[td_i][xy_p].i >= Ri_lo) && (P[td_i][xy_p].i < Ri_up)) {
                            minimum_one_point_inside_window = 1;
                            int ix = Rhdivr * (P[td_i][xy_p].r - Rr_lo);
                            int iy = Rwdivi * (P[td_i][xy_p].i - Ri_lo);
                            RPo[td_i][offset_count++] = (unsigned int)ix * Rw + iy;

                            if ((-P[td_i][xy_p].i >= Ri_lo) && (-P[td_i][xy_p].i < Ri_up)) {
                                int iyc = Rwdivi * (-P[td_i][xy_p].i - Ri_lo);
                                RPo[td_i][offset_count++] = (unsigned int)ix * Rw + iyc;
                            }
                        } else {
                            if ((-P[td_i][xy_p].i >= Ri_lo) && (-P[td_i][xy_p].i < Ri_up)) {
                                minimum_one_point_inside_window = 1;
                                int ix = Rhdivr * (P[td_i][xy_p].r - Rr_lo);
                                int iyc = Rwdivi * (-P[td_i][xy_p].i - Ri_lo);
                                RPo[td_i][offset_count++] = (unsigned int)ix * Rw + iyc;
                            }
                        }
                    }

                    if (minimum_one_point_inside_window) {
                        Pp[td_i] += 2;
                    }

                    for (int offset_iter = 0; offset_iter < offset_count; offset_iter++) {
                        R[td_i][RPo[td_i][offset_iter]]++;
                    }
                }
            }
        }

        if (bb_type[td_i] == 2) {
            for (int c_iter = 0; c_iter < C_BATCH_NB; c_iter++) {
                complex c;

                while (1) {
                    c.r = Cr_lo + Cr_ra * dsfmt_gv_genrand_close_open();
                    c.i = Ci_lo + Ci_ra * dsfmt_gv_genrand_close_open();
                    double c2 = c.r * c.r + c.i * c.i;

                    //// skip the bounded paths with large footprint = small periods
                    //// period 3 circle   = inexact = estimated by slightly larger circle c (-0.125,0.744) r 0.092
                    //// < -0.560697 --> corrected by testing (bail = 1000000, look at output png) ever smaller circles to < -0.5601 --> r 0.0951892851
                    if (!(c2 * (8.0 * c2 - 3.0) + c.r < 0.09375 || c2 + c.r + c.r < -0.9375 || c2 + 0.25 * c.r - 1.488 * c.i < -0.5601 || c2 + 2.616 * c.r < -1.7072)) {
                        break;
                    }
                }

                int n_inf = core_mandelbrot(td_i, c);

                if (n_inf == bb_bail[td_i]) {
                    int minimum_one_point_inside_window = 0;
                    int offset_count = 0;

                    for (int xy_p = bb_pps[td_i]; xy_p < n_inf - bb_ppe[td_i]; xy_p++) {
                        if ((P[td_i][xy_p].r >= Rr_up) || (P[td_i][xy_p].r < Rr_lo)) {
                            continue;
                        }

                        if ((P[td_i][xy_p].i >= Ri_lo) && (P[td_i][xy_p].i < Ri_up)) {
                            minimum_one_point_inside_window = 1;
                            int ix = Rhdivr * (P[td_i][xy_p].r - Rr_lo);
                            int iy = Rwdivi * (P[td_i][xy_p].i - Ri_lo);
                            RPo[td_i][offset_count++] = (unsigned int)ix * Rw + iy;

                            if ((-P[td_i][xy_p].i >= Ri_lo) && (-P[td_i][xy_p].i < Ri_up)) {
                                int iyc = Rwdivi * (-P[td_i][xy_p].i - Ri_lo);
                                RPo[td_i][offset_count++] = (unsigned int)ix * Rw + iyc;
                            }
                        } else {
                            if ((-P[td_i][xy_p].i >= Ri_lo) && (-P[td_i][xy_p].i < Ri_up)) {
                                minimum_one_point_inside_window = 1;
                                int ix = Rhdivr * (P[td_i][xy_p].r - Rr_lo);
                                int iyc = Rwdivi * (-P[td_i][xy_p].i - Ri_lo);
                                RPo[td_i][offset_count++] = (unsigned int)ix * Rw + iyc;
                            }
                        }
                    }

                    if (minimum_one_point_inside_window) {
                        Pp[td_i] += 2;
                    }

                    for (int offset_iter = 0; offset_iter < offset_count; offset_iter++) {
                        R[td_i][RPo[td_i][offset_iter]]++;
                    }
                }
            }
        }
    }
}

int save_param_file()
{
    FILE* parameters_file;

    if ((parameters_file = fopen("BuddhaBrot-MT-daemon-parameters.txt", "wb")) == NULL) {
        return (0);
    }

    fprintf(parameters_file, "centerx %lf\ncentery %lf\nzoom %lf\nwindowcenterx %lf\nwindowcentery %lf\nwindowzoom %lf\nrenderw %d\nrenderh %d\nwindowoffsetw %d\nwindowoffseth %d\n\nlayermode %d\n\nbuddhatype1 %d\nbailout1 %d\npathplotstart1 %d\npathplotend1 %d\npathminninf1 %d\n\nbuddhatype2 %d\nbailout2 %d\npathplotstart2 %d\npathplotend2 %d\npathminninf2 %d\n\nbuddhatype3 %d\nbailout3 %d\npathplotstart3 %d\npathplotend3 %d\npathminninf3 %d\n\ncolortabletype1 %d\ncolortabletype2 %d\ncolortabletype3 %d\n\ncoloringmethod1 %d\ncoloringsumfunc1 %d\ncoloringsumfuncparam11 %d\ncolortablescale1 %d\ncolortableoffset1 %d\n\ncoloringmethod2 %d\ncoloringsumfunc2 %d\ncoloringsumfuncparam12 %d\ncolortablescale2 %d\ncolortableoffset2 %d\n\ncoloringmethod3 %d\ncoloringsumfunc3 %d\ncoloringsumfuncparam13 %d\ncolortablescale3 %d\ncolortableoffset3 %d\n\nnumberofcalculationthreads %d\n", 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), 0.5 * (Wr_lo + Wr_up), 0.5 * (Wi_lo + Wi_up), 4.0 / (Wr_up - Wr_lo), Rw, Rh, RWow, RWoh, lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], ct_type[0], ct_type[1], ct_type[2], cm[0], csf[0], csfp1[0], (int)(10.0 * (ct_f[0] - 1.0)), ct_o[0], cm[1], csf[1], csfp1[1], (int)(10.0 * (ct_f[1] - 1.0)), ct_o[1], cm[2], csf[2], csfp1[2], (int)(10.0 * (ct_f[2] - 1.0)), ct_o[2], td_nb);
    fprintf(parameters_file, "\n");

    for (int td_i = 0; td_i < TD_MAX; td_i++) {
        fprintf(parameters_file, "pathsplottedthread%02d %llu\n", td_i, Pp[td_i]);
    }

    fprintf(parameters_file, "\n");
    fprintf(parameters_file, "pathsplottedsum %llu\n", Ppsum);
    fprintf(parameters_file, "\n");
    fprintf(parameters_file, "rendercountmatrixmax1 %d\n", Rlrmax[0]);
    fprintf(parameters_file, "rendercountmatrixmax2 %d\n", Rlrmax[1]);
    fprintf(parameters_file, "rendercountmatrixmax3 %d\n", Rlrmax[2]);
    fclose(parameters_file);
    return (1);
}

int load_param_file(int pause_calcthreads, int load_status_files, int minmem)
{
    FILE* parameters_file;

    if ((parameters_file = fopen("BuddhaBrot-MT-daemon-parameters.txt", "rb")) == NULL) {
        return (0);
    }

    double centerx = 0.0, centery = 0.0, zoom = 0.0;
    double windowcenterx = 0.0, windowcentery = 0.0, windowzoom = 0.0;
    int renderw = 0, renderh = 0;
    int windowoffsetw = 0, windowoffseth = 0;
    int layermode = 0;
    int buddhatype1 = 0, bailout1 = 0, pathplotstart1 = 0, pathplotend1 = 0, pathminninf1 = 0;
    int buddhatype2 = 0, bailout2 = 0, pathplotstart2 = 0, pathplotend2 = 0, pathminninf2 = 0;
    int buddhatype3 = 0, bailout3 = 0, pathplotstart3 = 0, pathplotend3 = 0, pathminninf3 = 0;
    int colortabletype1 = 0;
    int colortabletype2 = 0;
    int colortabletype3 = 0;
    int coloringmethod1 = 0, coloringsumfunc1 = 0, coloringsumfuncparam11 = 0, colortablescale1 = 0, colortableoffset1 = 0;
    int coloringmethod2 = 0, coloringsumfunc2 = 0, coloringsumfuncparam12 = 0, colortablescale2 = 0, colortableoffset2 = 0;
    int coloringmethod3 = 0, coloringsumfunc3 = 0, coloringsumfuncparam13 = 0, colortablescale3 = 0, colortableoffset3 = 0;
    int numberofcalculationthreads = 0;

    if (fscanf(parameters_file, "centerx %lf\ncentery %lf\nzoom %lf\nwindowcenterx %lf\nwindowcentery %lf\nwindowzoom %lf\nrenderw %d\nrenderh %d\nwindowoffsetw %d\nwindowoffseth %d\n\nlayermode %d\n\nbuddhatype1 %d\nbailout1 %d\npathplotstart1 %d\npathplotend1 %d\npathminninf1 %d\n\nbuddhatype2 %d\nbailout2 %d\npathplotstart2 %d\npathplotend2 %d\npathminninf2 %d\n\nbuddhatype3 %d\nbailout3 %d\npathplotstart3 %d\npathplotend3 %d\npathminninf3 %d\n\ncolortabletype1 %d\ncolortabletype2 %d\ncolortabletype3 %d\n\ncoloringmethod1 %d\ncoloringsumfunc1 %d\ncoloringsumfuncparam11 %d\ncolortablescale1 %d\ncolortableoffset1 %d\n\ncoloringmethod2 %d\ncoloringsumfunc2 %d\ncoloringsumfuncparam12 %d\ncolortablescale2 %d\ncolortableoffset2 %d\n\ncoloringmethod3 %d\ncoloringsumfunc3 %d\ncoloringsumfuncparam13 %d\ncolortablescale3 %d\ncolortableoffset3 %d\n\nnumberofcalculationthreads %d\n", &centerx, &centery, &zoom, &windowcenterx, &windowcentery, &windowzoom, &renderw, &renderh, &windowoffsetw, &windowoffseth, &layermode, &buddhatype1, &bailout1, &pathplotstart1, &pathplotend1, &pathminninf1, &buddhatype2, &bailout2, &pathplotstart2, &pathplotend2, &pathminninf2, &buddhatype3, &bailout3, &pathplotstart3, &pathplotend3, &pathminninf3, &colortabletype1, &colortabletype2, &colortabletype3, &coloringmethod1, &coloringsumfunc1, &coloringsumfuncparam11, &colortablescale1, &colortableoffset1, &coloringmethod2, &coloringsumfunc2, &coloringsumfuncparam12, &colortablescale2, &colortableoffset2, &coloringmethod3, &coloringsumfunc3, &coloringsumfuncparam13, &colortablescale3, &colortableoffset3, &numberofcalculationthreads)) {}

    if (minmem == 0) {
        while (numberofcalculationthreads > td_nb) {
            increase_num_calcthreads();
        }

        while (numberofcalculationthreads < td_nb) {
            decrease_num_calcthreads();
        }
    } else {
        while (3 < td_nb) {
            decrease_num_calcthreads();
        }
    }

    load_location_bb_color_param(pause_calcthreads, zoom, centerx, centery, renderw, renderh, layermode, buddhatype1, bailout1, pathplotstart1, pathplotend1, pathminninf1, buddhatype2, bailout2, pathplotstart2, pathplotend2, pathminninf2, buddhatype3, bailout3, pathplotstart3, pathplotend3, pathminninf3, colortabletype1, colortabletype2, colortabletype3, coloringmethod1, coloringsumfunc1, coloringsumfuncparam11, (double)colortablescale1 / 10.0 + 1.0, colortableoffset1, coloringmethod2, coloringsumfunc2, coloringsumfuncparam12, (double)colortablescale2 / 10.0 + 1.0, colortableoffset2, coloringmethod3, coloringsumfunc3, coloringsumfuncparam13, (double)colortablescale3 / 10.0 + 1.0, colortableoffset3);
    RWow = windowoffsetw;
    Wi_lo = Ri_lo + RWow / Rwdivi;
    Wi_up = Wi_lo + Ww / Rwdivi;
    RWoh = windowoffseth;
    Wr_lo = Rr_lo + RWoh / Rhdivr;
    Wr_up = Wr_lo + Wh / Rhdivr;

    if (load_status_files == 1) {
        if (fscanf(parameters_file, "\n")) {}

        for (int td_i = 0; td_i < TD_MAX; td_i++) {
            int temp;

            if (fscanf(parameters_file, "pathsplottedthread%d %llu\n", &temp, &Pp[td_i])) {}
        }
    }

    Ppsum = 0;

    for (int td_i = 0; td_i < td_nb; td_i += 1) {
        Ppsum += Pp[td_i];
    }

    Ppsum_autoPNG_last = (Ppsum / Ppsum_autoPNG_delta) * Ppsum_autoPNG_delta;
    fclose(parameters_file);
    return (numberofcalculationthreads);
}

int save_status_files()
{
    printf("\r                                                                                ");
    printf("\rsaving status...");
    fflush(stdout);
    pause_calcthreads_and_wait();
    Ppsum = 0;

    for (int td_i = 0; td_i < td_nb; td_i += 1) {
        Ppsum += Pp[td_i];
    }

    if (lr_mode == 0) {
        Rlrmax[0] = 0;
        Rlrmax[1] = 0;
        Rlrmax[2] = 0;

        for (int Wy = 0, Ry = RWoh; Wy < Wh; Wy++, Ry++) {
            for (int Wx = 0, Rx = RWow; Wx < Ww; Wx++, Rx++) {
                unsigned int Ri = (unsigned int)Ry * Rw + Rx;
                unsigned int sum = 0;

                for (int td_i = 0; td_i < td_nb; td_i += 1) {
                    sum += R[td_i][Ri];
                }

                if (csf[0] == 1) {
                    if (sum >= (unsigned int)csfp1[0]) {
                        sum = csfp1[0] + (unsigned int) log((double)(sum - csfp1[0] + 1));
                    }
                }

                if (sum > Rlrmax[0]) {
                    Rlrmax[0] = sum;
                }
            }
        }
    }

    if (lr_mode == 1) {
        for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
            Rlrmax[layer_iter] = 0;

            for (int Wy = 0, Ry = RWoh; Wy < Wh; Wy++, Ry++) {
                for (int Wx = 0, Rx = RWow; Wx < Ww; Wx++, Rx++) {
                    unsigned int Ri = (unsigned int)Ry * Rw + Rx;
                    unsigned int sum = 0;

                    for (int td_i = 0; td_i < td_nb; td_i += 1) {
                        sum += R[td_i][Ri];
                    }

                    if (csf[layer_iter] == 1) {
                        if (sum >= (unsigned int)csfp1[layer_iter]) {
                            sum = csfp1[layer_iter] + (unsigned int) log((double)(sum - csfp1[layer_iter] + 1));
                        }
                    }

                    if (sum > Rlrmax[layer_iter]) {
                        Rlrmax[layer_iter] = sum;
                    }
                }
            }
        }
    }

    if (lr_mode == 2) {
        for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
            Rlrmax[layer_iter] = 0;

            for (int Wy = 0, Ry = RWoh; Wy < Wh; Wy++, Ry++) {
                for (int Wx = 0, Rx = RWow; Wx < Ww; Wx++, Rx++) {
                    unsigned int Ri = (unsigned int)Ry * Rw + Rx;
                    unsigned int sum = 0;

                    for (int td_i = layer_iter; td_i < td_nb; td_i += LR_NB) {
                        sum += R[td_i][Ri];
                    }

                    if (csf[layer_iter] == 1) {
                        if (sum >= (unsigned int)csfp1[layer_iter]) {
                            sum = csfp1[layer_iter] + (unsigned int) log((double)(sum - csfp1[layer_iter] + 1));
                        }
                    }

                    if (sum > Rlrmax[layer_iter]) {
                        Rlrmax[layer_iter] = sum;
                    }
                }
            }
        }
    }

    save_param_file();

    for (int td_i = 0; td_i < td_nb; td_i += 1) {
        sprintf(filename, "BuddhaBrot-MT-status-thread%02i.bin", td_i);
        FILE* status_file;

        if ((status_file = fopen(filename, "wb")) == NULL) {
            return (0);
        }

        fwrite(R[td_i], sizeof(unsigned int), (unsigned int)Rw * Rh, status_file);
        fclose(status_file);
    }

    td_pause = 0;
    #pragma omp flush(td_pause)
    printf("done");
    fflush(stdout);
    return (1);
}

int load_status_files()
{
    printf("\r                                                                                ");
    printf("\rloading status...");
    fflush(stdout);
    pause_calcthreads_and_wait();
    load_param_file(0, 1, 0);

    for (int td_i = 0; td_i < td_nb; td_i += 1) {
        sprintf(filename, "BuddhaBrot-MT-status-thread%02i.bin", td_i);
        FILE* status_file;

        if ((status_file = fopen(filename, "rb")) == NULL) {
            return (0);
        }

        if (fread(R[td_i], sizeof(unsigned int), (unsigned int)Rw * Rh, status_file)) {}

        fclose(status_file);
    }

    td_pause = 0;
    #pragma omp flush(td_pause)
    printf("done");
    fflush(stdout);
    return (1);
}

int load_status_files_minmem()
{
    printf("\r                                                                                ");
    printf("\rloading status...");
    fflush(stdout);
    pause_calcthreads_and_wait();
    int original_num_calcthreads = load_param_file(0, 1, 1);
    unsigned int* tmp_R = (unsigned int*)calloc((unsigned int)Rw * Rh, sizeof(unsigned int));

    for (int td_i = 0; td_i < 3; td_i++) {
        sprintf(filename, "BuddhaBrot-MT-status-thread%02i.bin", td_i);
        FILE* status_file;

        if ((status_file = fopen(filename, "rb")) == NULL) {
            return (0);
        }

        if (fread(R[td_i], sizeof(unsigned int), (unsigned int)Rw * Rh, status_file)) {}

        fclose(status_file);
    }

    for (int td_i = 3; td_i < original_num_calcthreads; td_i++) {
        sprintf(filename, "BuddhaBrot-MT-status-thread%02i.bin", td_i);
        FILE* status_file;

        if ((status_file = fopen(filename, "rb")) == NULL) {
            return (0);
        }

        if (fread(tmp_R, sizeof(unsigned int), (unsigned int)Rw * Rh, status_file)) {}

        for (unsigned int Ri = 0; Ri < (unsigned int)Rw * Rh; Ri++) {
            R[td_i % 3][Ri] += tmp_R[Ri];
        }

        Pp[td_i % 3] += Pp[td_i];
        Pp[td_i] = 0;
        fclose(status_file);
    }

    free(tmp_R);
    td_pause = 0;
    #pragma omp flush(td_pause)
    printf("done");
    fflush(stdout);
    return (1);
}

void writeRtoPNG(const char* filename)
{
    printf("\r                                                                                ");
    printf("\rwriting to png...");
    fflush(stdout);
    pause_calcthreads_and_wait();

    if (lr_mode == 0) {
        Rlrmax[0] = 0;

        for (unsigned int Ri = 0; Ri < (unsigned int)Rw * Rh; Ri++) {
            unsigned int sum = 0;

            for (int td_i = 0; td_i < td_nb; td_i += 1) {
                sum += R[td_i][Ri];
            }

            if (csf[0] == 1) {
                if (sum >= (unsigned int)csfp1[0]) {
                    sum = csfp1[0] + (unsigned int) log((double)(sum - csfp1[0] + 1));
                }
            }

            if (sum > Rlrmax[0]) {
                Rlrmax[0] = sum;
            }
        }

        if (Rlrmax[0] >= Hl[0]) {
            if (H[0] != NULL) {
                free(H[0]);
                H[0] = NULL;
            }

            Hl[0] = Rlrmax[0] + 1;
            H[0] = (unsigned int*)calloc(Hl[0], sizeof(unsigned int));
        } else {
            memset(H[0], 0, Hl[0] * sizeof(unsigned int));
        }

        Rlrmax[1] = 0;
        Rlrmax[2] = 0;
        Hl[1] = 0;
        Hl[2] = 0;

        if (H[1] != NULL) {
            free(H[1]);
            H[1] = NULL;
        }

        if (H[2] != NULL) {
            free(H[2]);
            H[2] = NULL;
        }

        for (unsigned int Ri = 0; Ri < (unsigned int)Rw * Rh; Ri++) {
            unsigned int sum = 0;

            for (int td_i = 0; td_i < td_nb; td_i += 1) {
                sum += R[td_i][Ri];
            }

            if (csf[0] == 1) {
                if (sum >= (unsigned int)csfp1[0]) {
                    sum = csfp1[0] + (unsigned int) log((double)(sum - csfp1[0] + 1));
                }
            }

            H[0][sum]++;
        }

        if (cm[0] == 0) {
            cm0n[0] = 0;

            for (unsigned int i = 0; i <= Rlrmax[0]; i++) {
                if (H[0][i] > 0) {
                    H[0][i] = cm0n[0]++;
                }
            }
        }

        if (cm[0] == 1) {
            cm1n[0] = (unsigned int)Rw * Rh - H[0][0];
            H[0][0] = 0;

            for (unsigned int i = 2; i <= Rlrmax[0]; i++) {
                H[0][i] += H[0][i - 1];
            }
        }
    }

    if (lr_mode == 1) {
        for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
            Rlrmax[layer_iter] = 0;

            for (unsigned int Ri = 0; Ri < (unsigned int)Rw * Rh; Ri++) {
                unsigned int sum = 0;

                for (int td_i = 0; td_i < td_nb; td_i += 1) {
                    sum += R[td_i][Ri];
                }

                if (csf[layer_iter] == 1) {
                    if (sum >= (unsigned int)csfp1[layer_iter]) {
                        sum = csfp1[layer_iter] + (unsigned int) log((double)(sum - csfp1[layer_iter] + 1));
                    }
                }

                if (sum > Rlrmax[layer_iter]) {
                    Rlrmax[layer_iter] = sum;
                }
            }

            if (Rlrmax[layer_iter] >= Hl[layer_iter]) {
                if (H[layer_iter] != NULL) {
                    free(H[layer_iter]);
                    H[layer_iter] = NULL;
                }

                Hl[layer_iter] = Rlrmax[layer_iter] + 1;
                H[layer_iter] = (unsigned int*)calloc(Hl[layer_iter], sizeof(unsigned int));
            } else {
                memset(H[layer_iter], 0, Hl[layer_iter] * sizeof(unsigned int));
            }

            for (unsigned int Ri = 0; Ri < (unsigned int)Rw * Rh; Ri++) {
                unsigned int sum = 0;

                for (int td_i = 0; td_i < td_nb; td_i += 1) {
                    sum += R[td_i][Ri];
                }

                if (csf[layer_iter] == 1) {
                    if (sum >= (unsigned int)csfp1[layer_iter]) {
                        sum = csfp1[layer_iter] + (unsigned int) log((double)(sum - csfp1[layer_iter] + 1));
                    }
                }

                H[layer_iter][sum]++;
            }

            if (cm[layer_iter] == 0) {
                cm0n[layer_iter] = 0;

                for (unsigned int i = 0; i <= Rlrmax[layer_iter]; i++) {
                    if (H[layer_iter][i] > 0) {
                        H[layer_iter][i] = cm0n[layer_iter]++;
                    }
                }
            }

            if (cm[layer_iter] == 1) {
                cm1n[layer_iter] = (unsigned int)Rw * Rh - H[layer_iter][0];
                H[layer_iter][0] = 0;

                for (unsigned int i = 2; i <= Rlrmax[layer_iter]; i++) {
                    H[layer_iter][i] += H[layer_iter][i - 1];
                }
            }
        }
    }

    if (lr_mode == 2) {
        for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
            Rlrmax[layer_iter] = 0;

            for (unsigned int Ri = 0; Ri < (unsigned int)Rw * Rh; Ri++) {
                unsigned int sum = 0;

                for (int td_i = layer_iter; td_i < td_nb; td_i += LR_NB) {
                    sum += R[td_i][Ri];
                }

                if (csf[layer_iter] == 1) {
                    if (sum >= (unsigned int)csfp1[layer_iter]) {
                        sum = csfp1[layer_iter] + (unsigned int) log((double)(sum - csfp1[layer_iter] + 1));
                    }
                }

                if (sum > Rlrmax[layer_iter]) {
                    Rlrmax[layer_iter] = sum;
                }
            }

            if (Rlrmax[layer_iter] >= Hl[layer_iter]) {
                if (H[layer_iter] != NULL) {
                    free(H[layer_iter]);
                    H[layer_iter] = NULL;
                }

                Hl[layer_iter] = Rlrmax[layer_iter] + 1;
                H[layer_iter] = (unsigned int*)calloc(Hl[layer_iter], sizeof(unsigned int));
            } else {
                memset(H[layer_iter], 0, Hl[layer_iter] * sizeof(unsigned int));
            }

            for (unsigned int Ri = 0; Ri < (unsigned int)Rw * Rh; Ri++) {
                unsigned int sum = 0;

                for (int td_i = layer_iter; td_i < td_nb; td_i += LR_NB) {
                    sum += R[td_i][Ri];
                }

                if (csf[layer_iter] == 1) {
                    if (sum >= (unsigned int)csfp1[layer_iter]) {
                        sum = csfp1[layer_iter] + (unsigned int) log((double)(sum - csfp1[layer_iter] + 1));
                    }
                }

                H[layer_iter][sum]++;
            }

            if (cm[layer_iter] == 0) {
                cm0n[layer_iter] = 0;

                for (unsigned int i = 0; i <= Rlrmax[layer_iter]; i++) {
                    if (H[layer_iter][i] > 0) {
                        H[layer_iter][i] = cm0n[layer_iter]++;
                    }
                }
            }

            if (cm[layer_iter] == 1) {
                cm1n[layer_iter] = (unsigned int)Rw * Rh - H[layer_iter][0];
                H[layer_iter][0] = 0;

                for (unsigned int i = 2; i <= Rlrmax[layer_iter]; i++) {
                    H[layer_iter][i] += H[layer_iter][i - 1];
                }
            }
        }
    }

    FILE* outfile = fopen(filename, "wb");
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    png_init_io(png_ptr, outfile);
    png_set_compression_level(png_ptr, 1);
    png_set_IHDR(png_ptr, info_ptr, Rw, Rh, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png_ptr, info_ptr);
    png_bytep row_buffer = (png_bytep)calloc(Rw * 3, sizeof(png_byte));

    if (lr_mode == 0) {
        if (cm[0] == 0) {
            for (int png_y = 0; png_y < Rh; png_y++) {
                for (int png_x = 0; png_x < Rw; png_x++) {
                    unsigned int Ri = (unsigned int)png_y * Rw + png_x;
                    unsigned int sum = 0;

                    for (int td_i = 0; td_i < td_nb; td_i += 1) {
                        sum += R[td_i][Ri];
                    }

                    if (csf[0] == 1) {
                        if (sum >= (unsigned int)csfp1[0]) {
                            sum = csfp1[0] + (unsigned int) log((double)(sum - csfp1[0] + 1));
                        }
                    }

                    int ct_i = 0;

                    if (cm0n[0] > 1) {
                        ct_i = MIN(CTe[0], ct_f[0] * CTe[0] * ((double)H[0][sum] / (cm0n[0] - 1)));
                    }

                    ct_i = ct_cycle(ct_i + ct_o[0], CTe[0]);
                    row_buffer[png_x * 3 + 0] = CT[0][ct_i];
                    row_buffer[png_x * 3 + 1] = CT[0][ct_i];
                    row_buffer[png_x * 3 + 2] = CT[0][ct_i];
                }

                png_write_row(png_ptr, row_buffer);
            }
        }

        if (cm[0] == 1) {
            for (int png_y = 0; png_y < Rh; png_y++) {
                for (int png_x = 0; png_x < Rw; png_x++) {
                    unsigned int Ri = (unsigned int)png_y * Rw + png_x;
                    unsigned int sum = 0;

                    for (int td_i = 0; td_i < td_nb; td_i += 1) {
                        sum += R[td_i][Ri];
                    }

                    if (csf[0] == 1) {
                        if (sum >= (unsigned int)csfp1[0]) {
                            sum = csfp1[0] + (unsigned int) log((double)(sum - csfp1[0] + 1));
                        }
                    }

                    int ct_i = 0;

                    if (cm1n[0] > 0) {
                        ct_i = MIN(CTe[0], ct_f[0] * CTe[0] * ((double)H[0][sum] / cm1n[0]));
                    }

                    ct_i = ct_cycle(ct_i + ct_o[0], CTe[0]);
                    row_buffer[png_x * 3 + 0] = CT[0][ct_i];
                    row_buffer[png_x * 3 + 1] = CT[0][ct_i];
                    row_buffer[png_x * 3 + 2] = CT[0][ct_i];
                }

                png_write_row(png_ptr, row_buffer);
            }
        }

        if (cm[0] == 2) {
            for (int png_y = 0; png_y < Rh; png_y++) {
                for (int png_x = 0; png_x < Rw; png_x++) {
                    unsigned int Ri = (unsigned int)png_y * Rw + png_x;
                    unsigned int sum = 0;

                    for (int td_i = 0; td_i < td_nb; td_i += 1) {
                        sum += R[td_i][Ri];
                    }

                    if (csf[0] == 1) {
                        if (sum >= (unsigned int)csfp1[0]) {
                            sum = csfp1[0] + (unsigned int) log((double)(sum - csfp1[0] + 1));
                        }
                    }

                    int ct_i = 0;

                    if (Rlrmax[0] > 0) {
                        ct_i = MIN(CTe[0], ct_f[0] * CTe[0] * ((double)sum / Rlrmax[0]));
                    }

                    ct_i = ct_cycle(ct_i + ct_o[0], CTe[0]);
                    row_buffer[png_x * 3 + 0] = CT[0][ct_i];
                    row_buffer[png_x * 3 + 1] = CT[0][ct_i];
                    row_buffer[png_x * 3 + 2] = CT[0][ct_i];
                }

                png_write_row(png_ptr, row_buffer);
            }
        }
    }

    if (lr_mode == 1) {
        for (int png_y = 0; png_y < Rh; png_y++) {
            for (int png_x = 0; png_x < Rw; png_x++) {
                unsigned int Ri = (unsigned int)png_y * Rw + png_x;
                int ct_i[LR_NB] = {0};

                for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
                    unsigned int sum = 0;

                    for (int td_i = 0; td_i < td_nb; td_i += 1) {
                        sum += R[td_i][Ri];
                    }

                    if (csf[layer_iter] == 1) {
                        if (sum >= (unsigned int)csfp1[layer_iter]) {
                            sum = csfp1[layer_iter] + (unsigned int) log((double)(sum - csfp1[layer_iter] + 1));
                        }
                    }

                    if (cm[layer_iter] == 0) {
                        if (cm0n[layer_iter] > 1) {
                            ct_i[layer_iter] = MIN(CTe[layer_iter], ct_f[layer_iter] * CTe[layer_iter] * ((double)H[layer_iter][sum] / (cm0n[layer_iter] - 1)));
                        }
                    }

                    if (cm[layer_iter] == 1) {
                        if (cm1n[layer_iter] > 0) {
                            ct_i[layer_iter] = MIN(CTe[layer_iter], ct_f[layer_iter] * CTe[layer_iter] * ((double)H[layer_iter][sum] / cm1n[layer_iter]));
                        }
                    }

                    if (cm[layer_iter] == 2) {
                        if (Rlrmax[layer_iter] > 0) {
                            ct_i[layer_iter] = MIN(CTe[layer_iter], ct_f[layer_iter] * CTe[layer_iter] * ((double)sum / Rlrmax[layer_iter]));
                        }
                    }

                    ct_i[layer_iter] = ct_cycle(ct_i[layer_iter] + ct_o[layer_iter], CTe[layer_iter]);
                }

                row_buffer[png_x * 3 + 0] = CT[0][ct_i[0]];
                row_buffer[png_x * 3 + 1] = CT[1][ct_i[1]];
                row_buffer[png_x * 3 + 2] = CT[2][ct_i[2]];
            }

            png_write_row(png_ptr, row_buffer);
        }
    }

    if (lr_mode == 2) {
        for (int png_y = 0; png_y < Rh; png_y++) {
            for (int png_x = 0; png_x < Rw; png_x++) {
                unsigned int Ri = (unsigned int)png_y * Rw + png_x;
                int ct_i[LR_NB] = {0};

                for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
                    unsigned int sum = 0;

                    for (int td_i = layer_iter; td_i < td_nb; td_i += LR_NB) {
                        sum += R[td_i][Ri];
                    }

                    if (csf[layer_iter] == 1) {
                        if (sum >= (unsigned int)csfp1[layer_iter]) {
                            sum = csfp1[layer_iter] + (unsigned int) log((double)(sum - csfp1[layer_iter] + 1));
                        }
                    }

                    if (cm[layer_iter] == 0) {
                        if (cm0n[layer_iter] > 1) {
                            ct_i[layer_iter] = MIN(CTe[layer_iter], ct_f[layer_iter] * CTe[layer_iter] * ((double)H[layer_iter][sum] / (cm0n[layer_iter] - 1)));
                        }
                    }

                    if (cm[layer_iter] == 1) {
                        if (cm1n[layer_iter] > 0) {
                            ct_i[layer_iter] = MIN(CTe[layer_iter], ct_f[layer_iter] * CTe[layer_iter] * ((double)H[layer_iter][sum] / cm1n[layer_iter]));
                        }
                    }

                    if (cm[layer_iter] == 2) {
                        if (Rlrmax[layer_iter] > 0) {
                            ct_i[layer_iter] = MIN(CTe[layer_iter], ct_f[layer_iter] * CTe[layer_iter] * ((double)sum / Rlrmax[layer_iter]));
                        }
                    }

                    ct_i[layer_iter] = ct_cycle(ct_i[layer_iter] + ct_o[layer_iter], CTe[layer_iter]);
                }

                row_buffer[png_x * 3 + 0] = CT[0][ct_i[0]];
                row_buffer[png_x * 3 + 1] = CT[1][ct_i[1]];
                row_buffer[png_x * 3 + 2] = CT[2][ct_i[2]];
            }

            png_write_row(png_ptr, row_buffer);
        }
    }

    png_write_end(png_ptr, NULL);

    if (row_buffer != NULL) {
        free(row_buffer);
        row_buffer = NULL;
    }

    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(outfile);
    td_pause = 0;
    #pragma omp flush(td_pause)
    printf("done");
    fflush(stdout);
}

void writeTtoPNG(const char* filename, int local_png_offset_x, int local_png_offset_y, int local_png_width, int local_png_height)
{
    printf("\r                                                                                ");
    printf("\rwriting to png...");
    fflush(stdout);
    pause_calcthreads_and_wait();

    for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
        T[layer_iter] = (unsigned int*)calloc((unsigned int)local_png_width * local_png_height, sizeof(unsigned int));
    }

    if (lr_mode == 0) {
        Rlrmax[0] = 0;

        for (int png_y = 0, Ry = local_png_offset_y; png_y < local_png_height; png_y++, Ry++) {
            for (int png_x = 0, Rx = local_png_offset_x; png_x < local_png_width; png_x++, Rx++) {
                unsigned int Ri = (unsigned int)Ry * Rw + Rx;
                unsigned int sum = 0;

                for (int td_i = 0; td_i < td_nb; td_i += 1) {
                    sum += R[td_i][Ri];
                }

                if (csf[0] == 1) {
                    if (sum >= (unsigned int)csfp1[0]) {
                        sum = csfp1[0] + (unsigned int) log((double)(sum - csfp1[0] + 1));
                    }
                }

                if (sum > Rlrmax[0]) {
                    Rlrmax[0] = sum;
                }

                T[0][png_y * local_png_width + png_x] = sum;
            }
        }

        if (Rlrmax[0] >= Hl[0]) {
            if (H[0] != NULL) {
                free(H[0]);
                H[0] = NULL;
            }

            Hl[0] = Rlrmax[0] + 1;
            H[0] = (unsigned int*)calloc(Hl[0], sizeof(unsigned int));
        } else {
            memset(H[0], 0, Hl[0] * sizeof(unsigned int));
        }

        Rlrmax[1] = 0;
        Rlrmax[2] = 0;
        Hl[1] = 0;
        Hl[2] = 0;

        if (H[1] != NULL) {
            free(H[1]);
            H[1] = NULL;
        }

        if (H[2] != NULL) {
            free(H[2]);
            H[2] = NULL;
        }

        for (int png_y = 0; png_y < local_png_height; png_y++) {
            for (int png_x = 0; png_x < local_png_width; png_x++) {
                H[0][T[0][png_y * local_png_width + png_x]]++;
            }
        }

        if (cm[0] == 0) {
            cm0n[0] = 0;

            for (unsigned int i = 0; i <= Rlrmax[0]; i++) {
                if (H[0][i] > 0) {
                    H[0][i] = cm0n[0]++;
                }
            }
        }

        if (cm[0] == 1) {
            cm1n[0] = (unsigned int)local_png_width * local_png_height - H[0][0];
            H[0][0] = 0;

            for (unsigned int i = 2; i <= Rlrmax[0]; i++) {
                H[0][i] += H[0][i - 1];
            }
        }
    }

    if (lr_mode == 1) {
        for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
            Rlrmax[layer_iter] = 0;

            for (int png_y = 0, Ry = local_png_offset_y; png_y < local_png_height; png_y++, Ry++) {
                for (int png_x = 0, Rx = local_png_offset_x; png_x < local_png_width; png_x++, Rx++) {
                    unsigned int Ri = (unsigned int)Ry * Rw + Rx;
                    unsigned int sum = 0;

                    for (int td_i = 0; td_i < td_nb; td_i += 1) {
                        sum += R[td_i][Ri];
                    }

                    if (csf[layer_iter] == 1) {
                        if (sum >= (unsigned int)csfp1[layer_iter]) {
                            sum = csfp1[layer_iter] + (unsigned int) log((double)(sum - csfp1[layer_iter] + 1));
                        }
                    }

                    if (sum > Rlrmax[layer_iter]) {
                        Rlrmax[layer_iter] = sum;
                    }

                    T[layer_iter][png_y * local_png_width + png_x] = sum;
                }
            }

            if (Rlrmax[layer_iter] >= Hl[layer_iter]) {
                if (H[layer_iter] != NULL) {
                    free(H[layer_iter]);
                    H[layer_iter] = NULL;
                }

                Hl[layer_iter] = Rlrmax[layer_iter] + 1;
                H[layer_iter] = (unsigned int*)calloc(Hl[layer_iter], sizeof(unsigned int));
            } else {
                memset(H[layer_iter], 0, Hl[layer_iter] * sizeof(unsigned int));
            }

            for (int png_y = 0; png_y < local_png_height; png_y++) {
                for (int png_x = 0; png_x < local_png_width; png_x++) {
                    H[layer_iter][T[layer_iter][png_y * local_png_width + png_x]]++;
                }
            }

            if (cm[layer_iter] == 0) {
                cm0n[layer_iter] = 0;

                for (unsigned int i = 0; i <= Rlrmax[layer_iter]; i++) {
                    if (H[layer_iter][i] > 0) {
                        H[layer_iter][i] = cm0n[layer_iter]++;
                    }
                }
            }

            if (cm[layer_iter] == 1) {
                cm1n[layer_iter] = (unsigned int)local_png_width * local_png_height - H[layer_iter][0];
                H[layer_iter][0] = 0;

                for (unsigned int i = 2; i <= Rlrmax[layer_iter]; i++) {
                    H[layer_iter][i] += H[layer_iter][i - 1];
                }
            }
        }
    }

    if (lr_mode == 2) {
        for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
            Rlrmax[layer_iter] = 0;

            for (int png_y = 0, Ry = local_png_offset_y; png_y < local_png_height; png_y++, Ry++) {
                for (int png_x = 0, Rx = local_png_offset_x; png_x < local_png_width; png_x++, Rx++) {
                    unsigned int Ri = (unsigned int)Ry * Rw + Rx;
                    unsigned int sum = 0;

                    for (int td_i = layer_iter; td_i < td_nb; td_i += LR_NB) {
                        sum += R[td_i][Ri];
                    }

                    if (csf[layer_iter] == 1) {
                        if (sum >= (unsigned int)csfp1[layer_iter]) {
                            sum = csfp1[layer_iter] + (unsigned int) log((double)(sum - csfp1[layer_iter] + 1));
                        }
                    }

                    if (sum > Rlrmax[layer_iter]) {
                        Rlrmax[layer_iter] = sum;
                    }

                    T[layer_iter][png_y * local_png_width + png_x] = sum;
                }
            }

            if (Rlrmax[layer_iter] >= Hl[layer_iter]) {
                if (H[layer_iter] != NULL) {
                    free(H[layer_iter]);
                    H[layer_iter] = NULL;
                }

                Hl[layer_iter] = Rlrmax[layer_iter] + 1;
                H[layer_iter] = (unsigned int*)calloc(Hl[layer_iter], sizeof(unsigned int));
            } else {
                memset(H[layer_iter], 0, Hl[layer_iter] * sizeof(unsigned int));
            }

            for (int png_y = 0; png_y < local_png_height; png_y++) {
                for (int png_x = 0; png_x < local_png_width; png_x++) {
                    H[layer_iter][T[layer_iter][png_y * local_png_width + png_x]]++;
                }
            }

            if (cm[layer_iter] == 0) {
                cm0n[layer_iter] = 0;

                for (unsigned int i = 0; i <= Rlrmax[layer_iter]; i++) {
                    if (H[layer_iter][i] > 0) {
                        H[layer_iter][i] = cm0n[layer_iter]++;
                    }
                }
            }

            if (cm[layer_iter] == 1) {
                cm1n[layer_iter] = (unsigned int)local_png_width * local_png_height - H[layer_iter][0];
                H[layer_iter][0] = 0;

                for (unsigned int i = 2; i <= Rlrmax[layer_iter]; i++) {
                    H[layer_iter][i] += H[layer_iter][i - 1];
                }
            }
        }
    }

    FILE* outfile = fopen(filename, "wb");
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    png_init_io(png_ptr, outfile);
    png_set_compression_level(png_ptr, 1);
    png_set_IHDR(png_ptr, info_ptr, local_png_width, local_png_height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png_ptr, info_ptr);
    png_bytep row_buffer = (png_bytep)calloc(local_png_width * 3, sizeof(png_byte));

    if (lr_mode == 0) {
        if (cm[0] == 0) {
            for (int png_y = 0; png_y < local_png_height; png_y++) {
                for (int png_x = 0; png_x < local_png_width; png_x++) {
                    int ct_i = 0;

                    if (cm0n[0] > 1) {
                        ct_i = MIN(CTe[0], ct_f[0] * CTe[0] * ((double)H[0][T[0][png_y * local_png_width + png_x]] / (cm0n[0] - 1)));
                    }

                    ct_i = ct_cycle(ct_i + ct_o[0], CTe[0]);
                    row_buffer[png_x * 3 + 0] = CT[0][ct_i];
                    row_buffer[png_x * 3 + 1] = CT[0][ct_i];
                    row_buffer[png_x * 3 + 2] = CT[0][ct_i];
                }

                png_write_row(png_ptr, row_buffer);
            }
        }

        if (cm[0] == 1) {
            for (int png_y = 0; png_y < local_png_height; png_y++) {
                for (int png_x = 0; png_x < local_png_width; png_x++) {
                    int ct_i = 0;

                    if (cm1n[0] > 0) {
                        ct_i = MIN(CTe[0], ct_f[0] * CTe[0] * ((double)H[0][T[0][png_y * local_png_width + png_x]] / cm1n[0]));
                    }

                    ct_i = ct_cycle(ct_i + ct_o[0], CTe[0]);
                    row_buffer[png_x * 3 + 0] = CT[0][ct_i];
                    row_buffer[png_x * 3 + 1] = CT[0][ct_i];
                    row_buffer[png_x * 3 + 2] = CT[0][ct_i];
                }

                png_write_row(png_ptr, row_buffer);
            }
        }

        if (cm[0] == 2) {
            for (int png_y = 0; png_y < local_png_height; png_y++) {
                for (int png_x = 0; png_x < local_png_width; png_x++) {
                    int ct_i = 0;

                    if (Rlrmax[0] > 0) {
                        ct_i = MIN(CTe[0], ct_f[0] * CTe[0] * ((double)T[0][png_y * local_png_width + png_x] / Rlrmax[0]));
                    }

                    ct_i = ct_cycle(ct_i + ct_o[0], CTe[0]);
                    row_buffer[png_x * 3 + 0] = CT[0][ct_i];
                    row_buffer[png_x * 3 + 1] = CT[0][ct_i];
                    row_buffer[png_x * 3 + 2] = CT[0][ct_i];
                }

                png_write_row(png_ptr, row_buffer);
            }
        }
    }

    if (lr_mode == 1) {
        for (int png_y = 0; png_y < local_png_height; png_y++) {
            for (int png_x = 0; png_x < local_png_width; png_x++) {
                int ct_i[LR_NB] = {0};

                for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
                    if (cm[layer_iter] == 0) {
                        if (cm0n[layer_iter] > 1) {
                            ct_i[layer_iter] = MIN(CTe[layer_iter], ct_f[layer_iter] * CTe[layer_iter] * ((double)H[layer_iter][T[layer_iter][png_y * local_png_width + png_x]] / (cm0n[layer_iter] - 1)));
                        }
                    }

                    if (cm[layer_iter] == 1) {
                        if (cm1n[layer_iter] > 0) {
                            ct_i[layer_iter] = MIN(CTe[layer_iter], ct_f[layer_iter] * CTe[layer_iter] * ((double)H[layer_iter][T[layer_iter][png_y * local_png_width + png_x]] / cm1n[layer_iter]));
                        }
                    }

                    if (cm[layer_iter] == 2) {
                        if (Rlrmax[layer_iter] > 0) {
                            ct_i[layer_iter] = MIN(CTe[layer_iter], ct_f[layer_iter] * CTe[layer_iter] * ((double)T[layer_iter][png_y * local_png_width + png_x] / Rlrmax[layer_iter]));
                        }
                    }

                    ct_i[layer_iter] = ct_cycle(ct_i[layer_iter] + ct_o[layer_iter], CTe[layer_iter]);
                }

                row_buffer[png_x * 3 + 0] = CT[0][ct_i[0]];
                row_buffer[png_x * 3 + 1] = CT[1][ct_i[1]];
                row_buffer[png_x * 3 + 2] = CT[2][ct_i[2]];
            }

            png_write_row(png_ptr, row_buffer);
        }
    }

    if (lr_mode == 2) {
        for (int png_y = 0; png_y < local_png_height; png_y++) {
            for (int png_x = 0; png_x < local_png_width; png_x++) {
                int ct_i[LR_NB] = {0};

                for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
                    if (cm[layer_iter] == 0) {
                        if (cm0n[layer_iter] > 1) {
                            ct_i[layer_iter] = MIN(CTe[layer_iter], ct_f[layer_iter] * CTe[layer_iter] * ((double)H[layer_iter][T[layer_iter][png_y * local_png_width + png_x]] / (cm0n[layer_iter] - 1)));
                        }
                    }

                    if (cm[layer_iter] == 1) {
                        if (cm1n[layer_iter] > 0) {
                            ct_i[layer_iter] = MIN(CTe[layer_iter], ct_f[layer_iter] * CTe[layer_iter] * ((double)H[layer_iter][T[layer_iter][png_y * local_png_width + png_x]] / cm1n[layer_iter]));
                        }
                    }

                    if (cm[layer_iter] == 2) {
                        if (Rlrmax[layer_iter] > 0) {
                            ct_i[layer_iter] = MIN(CTe[layer_iter], ct_f[layer_iter] * CTe[layer_iter] * ((double)T[layer_iter][png_y * local_png_width + png_x] / Rlrmax[layer_iter]));
                        }
                    }

                    ct_i[layer_iter] = ct_cycle(ct_i[layer_iter] + ct_o[layer_iter], CTe[layer_iter]);
                }

                row_buffer[png_x * 3 + 0] = CT[0][ct_i[0]];
                row_buffer[png_x * 3 + 1] = CT[1][ct_i[1]];
                row_buffer[png_x * 3 + 2] = CT[2][ct_i[2]];
            }

            png_write_row(png_ptr, row_buffer);
        }
    }

    png_write_end(png_ptr, NULL);

    if (row_buffer != NULL) {
        free(row_buffer);
        row_buffer = NULL;
    }

    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(outfile);

    for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
        if (T[layer_iter] != NULL) {
            free(T[layer_iter]);
            T[layer_iter] = NULL;
        }
    }

    td_pause = 0;
    #pragma omp flush(td_pause)
    printf("done");
    fflush(stdout);
}

void writeRtoPNG_and_generate_filename()
{
    if (lr_mode == 0) {
        sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i%i%i cm%i%i.%i.%i.%i R%ix%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type[0], ct_type[1], ct_type[2], cm[0], csf[0], csfp1[0], (int)(10.0 * (ct_f[0] - 1.0)), ct_o[0], Rw, Rh, (double)Ppsum);
    }

    if (lr_mode == 1) {
        sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i%i%i cm%i%i.%i.%i.%i cm%i%i.%i.%i.%i cm%i%i.%i.%i.%i R%ix%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type[0], ct_type[1], ct_type[2], cm[0], csf[0], csfp1[0], (int)(10.0 * (ct_f[0] - 1.0)), ct_o[0], cm[1], csf[1], csfp1[1], (int)(10.0 * (ct_f[1] - 1.0)), ct_o[1], cm[2], csf[2], csfp1[2], (int)(10.0 * (ct_f[2] - 1.0)), ct_o[2], Rw, Rh, (double)Ppsum);
    }

    if (lr_mode == 2) {
        sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i%i%i cm%i%i.%i.%i.%i cm%i%i.%i.%i.%i cm%i%i.%i.%i.%i R%ix%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type[0], ct_type[1], ct_type[2], cm[0], csf[0], csfp1[0], (int)(10.0 * (ct_f[0] - 1.0)), ct_o[0], cm[1], csf[1], csfp1[1], (int)(10.0 * (ct_f[1] - 1.0)), ct_o[1], cm[2], csf[2], csfp1[2], (int)(10.0 * (ct_f[2] - 1.0)), ct_o[2], Rw, Rh, (double)Ppsum);
    }

    writeRtoPNG(filename);
}

void writeRTtoPNG_and_generate_filenames()
{
    if (lr_mode == 0) {
        sprintf(dirname, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i%i%i cm%i%i.%i.%i.%i %g", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type[0], ct_type[1], ct_type[2], cm[0], csf[0], csfp1[0], (int)(10.0 * (ct_f[0] - 1.0)), ct_o[0], (double)Ppsum);
    }

    if (lr_mode == 1) {
        sprintf(dirname, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i%i%i cm%i%i.%i.%i.%i cm%i%i.%i.%i.%i cm%i%i.%i.%i.%i %g", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type[0], ct_type[1], ct_type[2], cm[0], csf[0], csfp1[0], (int)(10.0 * (ct_f[0] - 1.0)), ct_o[0], cm[1], csf[1], csfp1[1], (int)(10.0 * (ct_f[1] - 1.0)), ct_o[1], cm[2], csf[2], csfp1[2], (int)(10.0 * (ct_f[2] - 1.0)), ct_o[2], (double)Ppsum);
    }

    if (lr_mode == 2) {
        sprintf(dirname, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i%i%i cm%i%i.%i.%i.%i cm%i%i.%i.%i.%i cm%i%i.%i.%i.%i %g", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type[0], ct_type[1], ct_type[2], cm[0], csf[0], csfp1[0], (int)(10.0 * (ct_f[0] - 1.0)), ct_o[0], cm[1], csf[1], csfp1[1], (int)(10.0 * (ct_f[1] - 1.0)), ct_o[1], cm[2], csf[2], csfp1[2], (int)(10.0 * (ct_f[2] - 1.0)), ct_o[2], (double)Ppsum);
    }

    sprintf(commandname, "mkdir \"%s\"", dirname);

    if (system(commandname)) {
    }

    for (int png_offset_x = 0; png_offset_x + Tw <= Rw; png_offset_x += Tw) {
        for (int png_offset_y = 0; png_offset_y + Th <= Rh; png_offset_y += Th) {
            double Tr_lo = Rr_lo + png_offset_y / Rhdivr;
            double Tr_up = Tr_lo + Th / Rhdivr;
            double Ti_lo = Ri_lo + png_offset_x / Rwdivi;
            double Ti_up = Ti_lo + Tw / Rwdivi;
            sprintf(filename, "%s/T%ix%i %05i %05i W(%.6f %.6f %.1f).png", dirname, Tw, Th, png_offset_y, png_offset_x, 0.5 * (Tr_lo + Tr_up), 0.5 * (Ti_lo + Ti_up), 4.0 / (Tr_up - Tr_lo));
            writeTtoPNG(filename, png_offset_x, png_offset_y, Tw, Th);
        }
    }
}

void load_status_files_thread()
{
    load_status_files();
    double memory_usage = 0;
    memory_usage += (double)Rw * Rh * 4 * td_nb;
    memory_usage += (double)bb_bail[0] * 24 * td_nb / 3;
    memory_usage += (double)bb_bail[1] * 24 * td_nb / 3;
    memory_usage += (double)bb_bail[2] * 24 * td_nb / 3;
    memory_usage += ((double)Rlrmax[0] + Rlrmax[1] + Rlrmax[2] + 3) * 4;
    memory_usage += ((double)Ww * Wh + Tw * Th) * 12;
    memory_usage /= 1024.0 * 1024.0 * 1024.0;
    printf("\r                                                                                ");
    printf("\rstarted: calcthreads %i   mem %.3f\n", td_nb, memory_usage);
    fflush(stdout);

    while (1) {
        while (Ppsum < Ppsum_autoPNG_last + Ppsum_autoPNG_delta) {
            wait_ms(100);
            Ppsum = 0;

            for (int td_i = 0; td_i < td_nb; td_i += 1) {
                Ppsum += Pp[td_i];
            }

            printf("\r                                                                                ");
            printf("\r#paths plotted = %g", (double)Ppsum);
            fflush(stdout);
        }

        cm[0] = 0;
        cm[1] = 0;
        cm[2] = 0;
        writeRtoPNG_and_generate_filename();
        cm[0] = 1;
        cm[1] = 1;
        cm[2] = 1;
        writeRtoPNG_and_generate_filename();
        save_status_files();
        Ppsum_autoPNG_last = (Ppsum / Ppsum_autoPNG_delta) * Ppsum_autoPNG_delta;
    }
}

void load_param_file_thread()
{
    load_param_file(1, 0, 0);
    double memory_usage = 0;
    memory_usage += (double)Rw * Rh * 4 * td_nb;
    memory_usage += (double)bb_bail[0] * 24 * td_nb / 3;
    memory_usage += (double)bb_bail[1] * 24 * td_nb / 3;
    memory_usage += (double)bb_bail[2] * 24 * td_nb / 3;
    memory_usage += ((double)Rlrmax[0] + Rlrmax[1] + Rlrmax[2] + 3) * 4;
    memory_usage += ((double)Ww * Wh + Tw * Th) * 12;
    memory_usage /= 1024.0 * 1024.0 * 1024.0;
    printf("\r                                                                                ");
    printf("\rstarted: calcthreads %i   mem %.3f\n", td_nb, memory_usage);
    fflush(stdout);

    while (1) {
        while (Ppsum < Ppsum_autoPNG_last + Ppsum_autoPNG_delta) {
            wait_ms(100);
            Ppsum = 0;

            for (int td_i = 0; td_i < td_nb; td_i += 1) {
                Ppsum += Pp[td_i];
            }

            printf("\r                                                                                ");
            printf("\r#paths plotted = %g", (double)Ppsum);
            fflush(stdout);
        }

        cm[0] = 0;
        cm[1] = 0;
        cm[2] = 0;
        writeRtoPNG_and_generate_filename();
        cm[0] = 1;
        cm[1] = 1;
        cm[2] = 1;
        writeRtoPNG_and_generate_filename();
        save_status_files();
        Ppsum_autoPNG_last = (Ppsum / Ppsum_autoPNG_delta) * Ppsum_autoPNG_delta;
    }
}

void batch_render()
{
    for (int bail = 100; bail <= 100000000; bail *= 10) {
        int minn = 0;
        printf("\r                                                                                ");
        printf("\rbail %i   minn %i\n", bail, minn);
        fflush(stdout);
        load_location_bb_color_param(1, 1.6, -0.4, 0.0, 1000, 1000, 0, 0, bail, 0, 0, minn, 0, bail, 0, 0, minn, 0, bail, 0, 0, minn, 0, 0, 0, 0, 0, 5, 1.0, 0, 0, 0, 5, 1.0, 0, 0, 0, 5, 1.0, 0);
        Ppsum = 0;

        while (Ppsum < Ppsum_autoPNG_delta) {
            wait_ms(100);
            Ppsum = 0;

            for (int td_i = 0; td_i < td_nb; td_i += 1) {
                Ppsum += Pp[td_i];
            }

            printf("\r                                                                                ");
            printf("\r#paths plotted = %g", (double)Ppsum);
            fflush(stdout);
        }

        cm[0] = 0;
        cm[1] = 0;
        cm[2] = 0;
        writeRtoPNG_and_generate_filename();
        cm[0] = 1;
        cm[1] = 1;
        cm[2] = 1;
        writeRtoPNG_and_generate_filename();
    }

    td_stop = 1;
    #pragma omp flush(td_stop)
}

void daemon_thread()
{
    if (daemon_mode == 0) {
        load_status_files_thread();
    }

    if (daemon_mode == 1) {
        load_param_file_thread();
    }

    if (daemon_mode == 2) {
        batch_render();
    }
}

int main(int argc, char* argv[])
{
    if (argc >= 4) {
        daemon_mode = atoi(argv[1]);
        td_nb = MIN(atoi(argv[2]), TD_MAX);
        Ppsum_autoPNG_delta = atof(argv[3]);
    }

    //// init
    dsfmt_gv_init_gen_rand(0);
    load_location_bb_color_param(0, 1.6, -0.4, 0.0, 1000, 1000, 0, 0, 1000, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 0, 0, 0, 0, 5, 1.0, 0, 0, 0, 5, 1.0, 0, 0, 0, 5, 1.0, 0);

    for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
        W[layer_iter] = (unsigned int*)calloc(Ww * Wh, sizeof(unsigned int));
        H[layer_iter] = NULL;
    }

    #pragma omp parallel sections num_threads(TD_MAX + 1)
    {
        #pragma omp section
        daemon_thread();
        #pragma omp section
        calculation_thread(0);
        #pragma omp section
        calculation_thread(1);
        #pragma omp section
        calculation_thread(2);
        #pragma omp section
        calculation_thread(3);
        #pragma omp section
        calculation_thread(4);
        #pragma omp section
        calculation_thread(5);
        #pragma omp section
        calculation_thread(6);
        #pragma omp section
        calculation_thread(7);
        #pragma omp section
        calculation_thread(8);
        #pragma omp section
        calculation_thread(9);
        #pragma omp section
        calculation_thread(10);
        #pragma omp section
        calculation_thread(11);
        #pragma omp section
        calculation_thread(12);
        #pragma omp section
        calculation_thread(13);
        #pragma omp section
        calculation_thread(14);
        #pragma omp section
        calculation_thread(15);
        #pragma omp section
        calculation_thread(16);
        #pragma omp section
        calculation_thread(17);
    }

    for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
        if (W[layer_iter] != NULL) {
            free(W[layer_iter]);
            W[layer_iter] = NULL;
        }

        if (H[layer_iter] != NULL) {
            free(H[layer_iter]);
            H[layer_iter] = NULL;
        }
    }

    for (int td_i = 0; td_i < td_nb; td_i += 1) {
        if (R[td_i] != NULL) {
            free(R[td_i]);
            R[td_i] = NULL;
        }

        if (P[td_i] != NULL) {
            free(P[td_i]);
            P[td_i] = NULL;
        }

        if (RPo[td_i] != NULL) {
            free(RPo[td_i]);
            RPo[td_i] = NULL;
        }
    }

    return (0);
}
