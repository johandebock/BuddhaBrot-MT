//// Written by Johan De Bock <johan.debock@gmail.com>
////
//// memory requirements : td_nb * (Rw * Rh * 4 + bb_bail * 24) + (Rlrmax[0] + Rlrmax[1] + Rlrmax[2] + 3) * 4 + (Ww * Wh + Tw * Th) * 12
////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <omp.h>

#include "SDL.h"
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

//// strings
char titlebar[512];
int titlebarswitch = 0;
char filename[512];
char dirname[512];
char commandname[512];

//// SDL
SDL_Window* sdl_window;
SDL_Renderer* sdl_renderer;
SDL_Surface* sdl_surface;
SDL_Texture* sdl_texture;
SDL_Event sdl_event;


//// keep total number of paths plotted for a render
long long unsigned int Ppsum;

//// fps frames per second regulation
double fps = 1.0;
Uint32 t_last_frame = 0; // time last frame was rendered

//// auto write PNG regulation
//// auto write mode 0 : no auto write
//// auto write mode 1 : auto write based on elapsed time
//// auto write mode 2 : auto write based on number of paths plotted
int autoWtoPNG = 0; // auto write window to one png
int autoRTtoPNG = 0; // auto write tiled render to multiple pngs
int autoRtoPNG = 0; // auto write render to one png
#define AUTOWRITE_MODE_NB 3 // number of auto write PNG modes

//// auto write mode 1
Uint32 t_autoPNG_last = 0; // time last auto PNG was written
Uint32 t_autoPNG_delta = 6e5; // time between each auto PNG write, default 10min

//// auto write mode 2
long long unsigned int Ppsum_autoPNG_last = 0; // Ppsum for last written auto PNG
long long unsigned int Ppsum_autoPNG_delta = 1e9; // Ppsum difference between each written auto PNG, default 1e9

//// td threads
#define TD_MAX 18 // maximum number of calc threads
int td_nb; // number of calc threads
int td_vc_nb = 2; // number of virtual cores
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
//// P[n] = z_{n+1} = z_n�+c
//// P = [c, c�+c, (c�+c)�+c, ...]
typedef struct {
    double r; // real
    double i; // imaginary
} complex; // complex number
complex* P[TD_MAX]; // 1 path , per thread
int* RPo[TD_MAX]; // offsets of a path in render , per thread

long long unsigned int Pp[TD_MAX]; // number of paths plotted , per thread

//// bb BuddhaBrot types and parameters
//// BuddhaBrot type 0 : the Buddhabrot
//// BuddhaBrot type 1 : the Anti-Buddhabrot
//// BuddhaBrot type 2 : the Anti-Buddhabrot with some lobes cut
//// Mandelbrot set = Mset = set of c for which : z_n bounded : lim_{n->inf} z_n != inf
//// notMset = set of c for which : z_n unbounded : lim_{n->inf} z_n == inf
////
//// Buddhabrot : notMset : unbounded
//// plot the unbounded paths = P[n_inf]� > 4 for n_inf < bailout
//// -> n_inf = core_mandelbrot(c)
//// skip the certainly bounded paths
////
//// Anti-Buddhabrot : Mset : bounded
//// plot the bounded paths = P[n_inf]� <= 4 for n_inf == bailout
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
int Rw; // render width
int Rh; // render height

unsigned int Rlrmax[LR_NB]; // maximum count , per layer

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
int RWox; // x offset of window in render
int RWoy; // y offset of window in render

double Wr_lo, Wr_up; // real lower and upper bound of rectangle W in complex plane
double Wi_lo, Wi_up; // imaginary lower and upper bound of rectangle W in complex plane

int recalc_WH_if_paused = 0; // indicate that H and W need to be recalculated when changing layer mode/coloring method/panning window when paused

//// T tile (count matrix (of pixels))
//// for the pixels outputted to a png image
//// it is located somewhere in the render
unsigned int* T[LR_NB]; // 1 tile , per layer
int Tw = 1000; // tile width
int Th = 1000; // tile height

//// cm coloring methods
//// coloring method 0 : rank-order mapping
//// coloring method 1 : histogram mapping
//// coloring method 0 : if cm0n > 1 ? ci = ct_e * H[count] / (cm0n - 1) : ci = ct_s TODO
//// coloring method 1 : if cm1n > 0 ? ci = ct_e * H[count] / cm1n : ci = ct_s TODO check ct_e * () in code
//// coloring method 2 : 0 , log
//// coloring method 3 : 1 , log
//// log : if count >= cm_log ? count = cm_log + log(count - cm_log + 1) : count = count
int cm[LR_NB]; // coloring method , per layer
#define CM_NORMAL_NB 2 // number of normal coloring methods
#define CM_LOG_NB 2 // number of normal coloring methods
unsigned int cm_log[LR_NB]; // logarithmic when >= cm_log
int cm0n[LR_NB]; // normalization value for coloring method 0 2 = number of filled bins in histogram = number of unique values in render
int cm1n[LR_NB]; // normalization value for coloring method 1 3 = (number of pixels in render - number of 0 pixels in render)

//// CT color tables
//// index -> RGB/Grey value
//// color table type 0 : Grey : 0->255->0 0->255->0 0->255->0
//// color table type 1 : Grey : 0->255 0->255 0->255
//// color table type 2 : Grey : 255->0 255->0 255->0
//// color table type 3 : RGB : 0->R->RG->RGB->RG->R->0
typedef unsigned char rgb[3]; // an RGB value
rgb* CT; // 1 color table
int ct_type; // color table type
#define CT_NONCYCLE_NB 8 // number of noncycle color table types : 1X
#define CT_CYCLE_NB 3 // number of cycle color table types : 2X
int ct_e; // last index in color table = len - 1
int ct_o[LR_NB]; // start offset in ct (0 index gets mapped on this color) , per layer
int ct_v[LR_NB] = {0, 0, 0}; // cycle speed : diff per frame in index, per layer




int ct_cycle(int i)
{
    return ((i > ct_e) ? (i - ct_e - 1) : ((i < 0) ? (ct_e + i + 1) : i));
}

void ct_load(int ct_type)
{
    //// 1X : noncycle
    if (ct_type == 10) {
        ct_e = 255;
        free(CT);
        CT = (rgb*)calloc(ct_e + 1, sizeof(rgb));

        for (int i = 0, j = 0; j <= 255; i++, j++) {
            CT[i][0] = j;
            CT[i][1] = j;
            CT[i][2] = j;
        }
    }

    if (ct_type == 11) {
        ct_e = 255;
        free(CT);
        CT = (rgb*)calloc(ct_e + 1, sizeof(rgb));

        for (int i = 0, j = 255; j >= 0; i++, j--) {
            CT[i][0] = j;
            CT[i][1] = j;
            CT[i][2] = j;
        }
    }

    if (ct_type == 12) {
        ct_e = 255;
        free(CT);
        CT = (rgb*)calloc(ct_e + 1, sizeof(rgb));

        for (int i = 0, j = 0; j <= 255; i++, j++) {
            CT[i][0] = 255 - j;
            CT[i][1] = j;
            CT[i][2] = j;
        }
    }

    if (ct_type == 13) {
        ct_e = 255;
        free(CT);
        CT = (rgb*)calloc(ct_e + 1, sizeof(rgb));

        for (int i = 0, j = 0; j <= 255; i++, j++) {
            CT[i][0] = j;
            CT[i][1] = 255 - j;
            CT[i][2] = j;
        }
    }

    if (ct_type == 14) {
        ct_e = 255;
        free(CT);
        CT = (rgb*)calloc(ct_e + 1, sizeof(rgb));

        for (int i = 0, j = 0; j <= 255; i++, j++) {
            CT[i][0] = j;
            CT[i][1] = j;
            CT[i][2] = 255 - j;
        }
    }

    if (ct_type == 15) {
        ct_e = 255;
        free(CT);
        CT = (rgb*)calloc(ct_e + 1, sizeof(rgb));

        for (int i = 0, j = 0; j <= 255; i++, j++) {
            CT[i][0] = j;
            CT[i][1] = 255 - j;
            CT[i][2] = 255 - j;
        }
    }

    if (ct_type == 16) {
        ct_e = 255;
        free(CT);
        CT = (rgb*)calloc(ct_e + 1, sizeof(rgb));

        for (int i = 0, j = 0; j <= 255; i++, j++) {
            CT[i][0] = 255 - j;
            CT[i][1] = j;
            CT[i][2] = 255 - j;
        }
    }

    if (ct_type == 17) {
        ct_e = 255;
        free(CT);
        CT = (rgb*)calloc(ct_e + 1, sizeof(rgb));

        for (int i = 0, j = 0; j <= 255; i++, j++) {
            CT[i][0] = 255 - j;
            CT[i][1] = 255 - j;
            CT[i][2] = j;
        }
    }

    //// 2X : cycle
    if (ct_type == 20) {
        ct_e = 509;
        free(CT);
        CT = (rgb*)calloc(ct_e + 1, sizeof(rgb));

        for (int i = 0, j = 0; j < 255; i++, j++) {
            CT[i][0] = j;
            CT[i][1] = j;
            CT[i][2] = j;
        }

        for (int i = 255, j = 255; j > 0; i++, j--) {
            CT[i][0] = j;
            CT[i][1] = j;
            CT[i][2] = j;
        }
    }

    if (ct_type == 21) {
        ct_e = 509;
        free(CT);
        CT = (rgb*)calloc(ct_e + 1, sizeof(rgb));

        for (int i = 0, j = 255; j > 0; i++, j--) {
            CT[i][0] = j;
            CT[i][1] = j;
            CT[i][2] = j;
        }

        for (int i = 255, j = 0; j < 255; i++, j++) {
            CT[i][0] = j;
            CT[i][1] = j;
            CT[i][2] = j;
        }
    }

    if (ct_type == 22) {
        ct_e = 1529;
        free(CT);
        CT = (rgb*)calloc(ct_e + 1, sizeof(rgb));

        for (int i = 0, j = 0; j < 255; i++, j++) {
            CT[i][0] = j;
            CT[i][1] = 0;
            CT[i][2] = 0;
        }

        for (int i = 255, j = 0; j < 255; i++, j++) {
            CT[i][0] = 255;
            CT[i][1] = j;
            CT[i][2] = 0;
        }

        for (int i = 510, j = 0; j < 255; i++, j++) {
            CT[i][0] = 255;
            CT[i][1] = 255;
            CT[i][2] = j;
        }

        for (int i = 765, j = 255; j > 0; i++, j--) {
            CT[i][0] = 255;
            CT[i][1] = 255;
            CT[i][2] = j;
        }

        for (int i = 1020, j = 255; j > 0; i++, j--) {
            CT[i][0] = 255;
            CT[i][1] = j;
            CT[i][2] = 0;
        }

        for (int i = 1275, j = 255; j > 0; i++, j--) {
            CT[i][0] = j;
            CT[i][1] = 0;
            CT[i][2] = 0;
        }
    }
}




void reponsive_caption_update(const char* string)
{
    SDL_SetWindowTitle(sdl_window, string);
    SDL_PollEvent(&sdl_event);
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
    memset(R[td_i], 0, Rw * Rh * sizeof(unsigned int));
}

void realloc_R(int td_i)
{
    Pp[td_i] = 0;

    if (R[td_i] != NULL) {
        free(R[td_i]);
        R[td_i] = NULL;
    }

    R[td_i] = (unsigned int*)calloc(Rw * Rh, sizeof(unsigned int));
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

    RPo[td_i] = (int*)calloc(2 * bb_bail[td_i], sizeof(int));
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
        for (int Ri = 0; Ri < Rw * Rh; Ri++) {
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

void load_location(int pause_calcthreads, double zoom, double x_l, double x_u, double y_l, double y_u, int cmw, int cmh)
{
    if (pause_calcthreads) {
        pause_calcthreads_and_wait();
    }

    if (zoom > 0.0) {
        double rangediv2 = 2.0 / zoom;
        x_u = x_l + rangediv2;
        x_l = x_l - rangediv2;
        y_u = y_l + rangediv2;
        y_l = y_l - rangediv2;
    }

    Rr_lo = x_l;
    Rr_up = x_u;
    Ri_lo = y_l;
    Ri_up = y_u;
    Rr_ra = Rr_up - Rr_lo;
    Ri_ra = Ri_up - Ri_lo;

    if (Rw == cmw && Rh == cmh) {
        RWox = 0;
        RWoy = 0;

        for (int td_i = 0; td_i < td_nb; td_i += 1) {
            reset_R(td_i);
        }
    } else {
        if (Rw == 0) {
            RWox = 0;
        } else {
            if (Rw < cmw) {
                RWox = RWox / Rw_f + 0.5 * Ww * (1.0 / Rw_f - 1.0);
            } else {
                RWox = RWox * Rw_f + 0.5 * Ww * (Rw_f - 1.0);
            }
        }

        if (Rh == 0) {
            RWoy = 0;
        } else {
            if (Rh < cmh) {
                RWoy = RWoy / Rh_f + 0.5 * Wh * (1.0 / Rh_f - 1.0);
            } else {
                RWoy = RWoy * Rh_f + 0.5 * Wh * (Rh_f - 1.0);
            }
        }

        Rw = cmw;
        Rh = cmh;

        if (Rw < Ww) {
            Rw = Ww;
            RWox = 0;
        } else if (RWox < 0) {
            RWox = 0;
        } else if (RWox + Ww > Rw) {
            RWox = Rw - Ww;
        }

        if (Rh < Wh) {
            Rh = Wh;
            RWoy = 0;
        } else if (RWoy < 0) {
            RWoy = 0;
        } else if (RWoy + Wh > Rh) {
            RWoy = Rh - Wh;
        }

        for (int td_i = 0; td_i < td_nb; td_i += 1) {
            realloc_R(td_i);
        }
    }

    Rhdivr = Rh / Rr_ra;
    Rwdivi = Rw / Ri_ra;
    Wr_lo = Rr_lo + RWoy / Rhdivr;
    Wr_up = Wr_lo + Wh / Rhdivr;
    Wi_lo = Ri_lo + RWox / Rwdivi;
    Wi_up = Wi_lo + Ww / Rwdivi;
    Hl[0] = 0;
    Hl[1] = 0;
    Hl[2] = 0;

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

    if (load_selectedlayer != -1) {
        for (int td_i = load_selectedlayer; td_i < td_nb; td_i += LR_NB) {
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
    }

    if (load_selectedlayer == -1) {
        for (int td_i = 0; td_i < td_nb; td_i += 1) {
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
    }

    td_pause = 0;
    #pragma omp flush(td_pause)
}

void load_location_bb_color_param(int pause_calcthreads, double zoom, double x_l, double x_u, double y_l, double y_u, int cmw, int cmh, int load_lr_mode, int load_bb_type1, int load_bb_bail1, int load_bb_pps1, int load_bb_ppe1, int load_bb_minn1, int load_bb_type2, int load_bb_bail2, int load_bb_pps2, int load_bb_ppe2, int load_bb_minn2, int load_bb_type3, int load_bb_bail3, int load_bb_pps3, int load_bb_ppe3, int load_bb_minn3, int load_ct_type, int load_cm1, int load_cm_log1, int load_ct_o1, int load_cm2, int load_cm_log2, int load_ct_o2, int load_cm3, int load_cm_log3, int load_ct_o3)
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
    ct_type = load_ct_type;
    ct_load(ct_type);
    cm[0] = load_cm1;
    cm[1] = load_cm2;
    cm[2] = load_cm3;
    cm_log[0] = load_cm_log1;
    cm_log[1] = load_cm_log2;
    cm_log[2] = load_cm_log3;
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

    load_location(0, zoom, x_l, x_u, y_l, y_u, cmw, cmh);

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
    printf("thread %i start\n", td_i);

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
                            printf("fix! notMset %f %f\n", c.r, c.i);
                        }
                    }
                }

                //// plot the unbounded paths = P[n_inf]� > 4 for n_inf < bailout
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
                            RPo[td_i][offset_count++] = ix * Rw + iy;

                            if ((-P[td_i][xy_p].i >= Ri_lo) && (-P[td_i][xy_p].i < Ri_up)) {
                                int iyc = Rwdivi * (-P[td_i][xy_p].i - Ri_lo);
                                RPo[td_i][offset_count++] = ix * Rw + iyc;
                            }
                        } else {
                            if ((-P[td_i][xy_p].i >= Ri_lo) && (-P[td_i][xy_p].i < Ri_up)) {
                                minimum_one_point_inside_window = 1;
                                int ix = Rhdivr * (P[td_i][xy_p].r - Rr_lo);
                                int iyc = Rwdivi * (-P[td_i][xy_p].i - Ri_lo);
                                RPo[td_i][offset_count++] = ix * Rw + iyc;
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

                //// plot the bounded paths = P[n_inf]� <= 4 for n_inf == bailout
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
                            RPo[td_i][offset_count++] = ix * Rw + iy;

                            if ((-P[td_i][xy_p].i >= Ri_lo) && (-P[td_i][xy_p].i < Ri_up)) {
                                int iyc = Rwdivi * (-P[td_i][xy_p].i - Ri_lo);
                                RPo[td_i][offset_count++] = ix * Rw + iyc;
                            }
                        } else {
                            if ((-P[td_i][xy_p].i >= Ri_lo) && (-P[td_i][xy_p].i < Ri_up)) {
                                minimum_one_point_inside_window = 1;
                                int ix = Rhdivr * (P[td_i][xy_p].r - Rr_lo);
                                int iyc = Rwdivi * (-P[td_i][xy_p].i - Ri_lo);
                                RPo[td_i][offset_count++] = ix * Rw + iyc;
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
                            RPo[td_i][offset_count++] = ix * Rw + iy;

                            if ((-P[td_i][xy_p].i >= Ri_lo) && (-P[td_i][xy_p].i < Ri_up)) {
                                int iyc = Rwdivi * (-P[td_i][xy_p].i - Ri_lo);
                                RPo[td_i][offset_count++] = ix * Rw + iyc;
                            }
                        } else {
                            if ((-P[td_i][xy_p].i >= Ri_lo) && (-P[td_i][xy_p].i < Ri_up)) {
                                minimum_one_point_inside_window = 1;
                                int ix = Rhdivr * (P[td_i][xy_p].r - Rr_lo);
                                int iyc = Rwdivi * (-P[td_i][xy_p].i - Ri_lo);
                                RPo[td_i][offset_count++] = ix * Rw + iyc;
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

    printf("thread %i end\n", td_i);
}

int save_param_file()
{
    FILE* parameters_file;

    if ((parameters_file = fopen("BuddhaBrot-MT-parameters.txt", "wb")) == NULL) {
        return (0);
    }

    fprintf(parameters_file, "centerx %lf\ncentery %lf\nzoom %lf\nwindowcenterx %lf\nwindowcentery %lf\nwindowzoom %lf\nrendersizex %d\nrendersizey %d\nwindowoffsetx %d\nwindowoffsety %d\n\nlayermode %d\n\nbuddhatype1 %d\nbailout1 %d\npathplotstart1 %d\npathplotend1 %d\npathminninf1 %d\n\nbuddhatype2 %d\nbailout2 %d\npathplotstart2 %d\npathplotend2 %d\npathminninf2 %d\n\nbuddhatype3 %d\nbailout3 %d\npathplotstart3 %d\npathplotend3 %d\npathminninf3 %d\n\ncolortabletype %d\n\ncoloringmethod1 %d\ncoloringmethodlog1 %d\ncolortableoffset1 %d\n\ncoloringmethod2 %d\ncoloringmethodlog2 %d\ncolortableoffset2 %d\n\ncoloringmethod3 %d\ncoloringmethodlog3 %d\ncolortableoffset3 %d\n\nnumberofcalculationthreads %d\n", 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), 0.5 * (Wr_lo + Wr_up), 0.5 * (Wi_lo + Wi_up), 4.0 / (Wr_up - Wr_lo), Rw, Rh, RWox, RWoy, lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], td_nb);
    fprintf(parameters_file, "\n");

    for (int td_i = 0; td_i < TD_MAX; td_i++) {
        fprintf(parameters_file, "pathsplottedthread%02d %I64u\n", td_i, Pp[td_i]);
    }

    fprintf(parameters_file, "\n");
    fprintf(parameters_file, "pathsplottedsum %I64u\n", Ppsum);
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

    if ((parameters_file = fopen("BuddhaBrot-MT-parameters.txt", "rb")) == NULL) {
        return (0);
    }

    double centerx = 0.0, centery = 0.0, zoom = 0.0;
    double windowcenterx = 0.0, windowcentery = 0.0, windowzoom = 0.0;
    int rendersizex = 0, rendersizey = 0;
    int windowoffsetx = 0, windowoffsety = 0;
    int layermode = 0;
    int buddhatype1 = 0, bailout1 = 0, pathplotstart1 = 0, pathplotend1 = 0, pathminninf1 = 0;
    int buddhatype2 = 0, bailout2 = 0, pathplotstart2 = 0, pathplotend2 = 0, pathminninf2 = 0;
    int buddhatype3 = 0, bailout3 = 0, pathplotstart3 = 0, pathplotend3 = 0, pathminninf3 = 0;
    int colortabletype = 0;
    int coloringmethod1 = 0, coloringmethodlog1 = 0, colortableoffset1 = 0;
    int coloringmethod2 = 0, coloringmethodlog2 = 0, colortableoffset2 = 0;
    int coloringmethod3 = 0, coloringmethodlog3 = 0, colortableoffset3 = 0;
    int numberofcalculationthreads = 0;

    if (fscanf(parameters_file, "centerx %lf\ncentery %lf\nzoom %lf\nwindowcenterx %lf\nwindowcentery %lf\nwindowzoom %lf\nrendersizex %d\nrendersizey %d\nwindowoffsetx %d\nwindowoffsety %d\n\nlayermode %d\n\nbuddhatype1 %d\nbailout1 %d\npathplotstart1 %d\npathplotend1 %d\npathminninf1 %d\n\nbuddhatype2 %d\nbailout2 %d\npathplotstart2 %d\npathplotend2 %d\npathminninf2 %d\n\nbuddhatype3 %d\nbailout3 %d\npathplotstart3 %d\npathplotend3 %d\npathminninf3 %d\n\ncolortabletype %d\n\ncoloringmethod1 %d\ncoloringmethodlog1 %d\ncolortableoffset1 %d\n\ncoloringmethod2 %d\ncoloringmethodlog2 %d\ncolortableoffset2 %d\n\ncoloringmethod3 %d\ncoloringmethodlog3 %d\ncolortableoffset3 %d\n\nnumberofcalculationthreads %d\n", &centerx, &centery, &zoom, &windowcenterx, &windowcentery, &windowzoom, &rendersizex, &rendersizey, &windowoffsetx, &windowoffsety, &layermode, &buddhatype1, &bailout1, &pathplotstart1, &pathplotend1, &pathminninf1, &buddhatype2, &bailout2, &pathplotstart2, &pathplotend2, &pathminninf2, &buddhatype3, &bailout3, &pathplotstart3, &pathplotend3, &pathminninf3, &colortabletype, &coloringmethod1, &coloringmethodlog1, &colortableoffset1, &coloringmethod2, &coloringmethodlog2, &colortableoffset2, &coloringmethod3, &coloringmethodlog3, &colortableoffset3, &numberofcalculationthreads)) {}

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

    load_location_bb_color_param(pause_calcthreads, zoom, centerx, 0.0, centery, 0.0, rendersizex, rendersizey, layermode, buddhatype1, bailout1, pathplotstart1, pathplotend1, pathminninf1, buddhatype2, bailout2, pathplotstart2, pathplotend2, pathminninf2, buddhatype3, bailout3, pathplotstart3, pathplotend3, pathminninf3, colortabletype, coloringmethod1, coloringmethodlog1, colortableoffset1, coloringmethod2, coloringmethodlog2, colortableoffset2, coloringmethod3, coloringmethodlog3, colortableoffset3);
    RWox = windowoffsetx;
    Wi_lo = Ri_lo + RWox / Rwdivi;
    Wi_up = Wi_lo + Ww / Rwdivi;
    RWoy = windowoffsety;
    Wr_lo = Rr_lo + RWoy / Rhdivr;
    Wr_up = Wr_lo + Wh / Rhdivr;

    if (load_status_files == 1) {
        if (fscanf(parameters_file, "\n")) {}

        for (int td_i = 0; td_i < TD_MAX; td_i++) {
            int temp;

            if (fscanf(parameters_file, "pathsplottedthread%d %I64u\n", &temp, &Pp[td_i])) {}
        }
    }

    Ppsum = 0;

    for (int td_i = 0; td_i < td_nb; td_i += 1) {
        Ppsum += Pp[td_i];
    }

    fclose(parameters_file);
    return (numberofcalculationthreads);
}

int save_status_files()
{
    reponsive_caption_update("BuddhaBrot-MT: saving status: pausing calculation threads...");
    pause_calcthreads_and_wait();
    Ppsum = 0;

    for (int td_i = 0; td_i < td_nb; td_i += 1) {
        Ppsum += Pp[td_i];
    }

    if (lr_mode == 0) {
        Rlrmax[0] = 0;
        Rlrmax[1] = 0;
        Rlrmax[2] = 0;

        for (int Wy = 0, Ry = RWoy; Wy < Wh; Wy++, Ry++) {
            for (int Wx = 0, Rx = RWox; Wx < Ww; Wx++, Rx++) {
                int Ri = Ry * Rw + Rx;
                unsigned int sum = 0;

                for (int td_i = 0; td_i < td_nb; td_i += 1) {
                    sum += R[td_i][Ri];
                }

                if (cm[0] == 2 || cm[0] == 3) {
                    if (sum >= cm_log[0]) {
                        sum = cm_log[0] + (unsigned int) log((double)(sum - cm_log[0] + 1));
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

            for (int Wy = 0, Ry = RWoy; Wy < Wh; Wy++, Ry++) {
                for (int Wx = 0, Rx = RWox; Wx < Ww; Wx++, Rx++) {
                    int Ri = Ry * Rw + Rx;
                    unsigned int sum = 0;

                    for (int td_i = 0; td_i < td_nb; td_i += 1) {
                        sum += R[td_i][Ri];
                    }

                    if (cm[layer_iter] == 2 || cm[layer_iter] == 3) {
                        if (sum >= cm_log[layer_iter]) {
                            sum = cm_log[layer_iter] + (unsigned int) log((double)(sum - cm_log[layer_iter] + 1));
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

            for (int Wy = 0, Ry = RWoy; Wy < Wh; Wy++, Ry++) {
                for (int Wx = 0, Rx = RWox; Wx < Ww; Wx++, Rx++) {
                    int Ri = Ry * Rw + Rx;
                    unsigned int sum = 0;

                    for (int td_i = layer_iter; td_i < td_nb; td_i += LR_NB) {
                        sum += R[td_i][Ri];
                    }

                    if (cm[layer_iter] == 2 || cm[layer_iter] == 3) {
                        if (sum >= cm_log[layer_iter]) {
                            sum = cm_log[layer_iter] + (unsigned int) log((double)(sum - cm_log[layer_iter] + 1));
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
        sprintf(titlebar, "BuddhaBrot-MT: saving status: saving render count matrix of thread %i...", td_i);
        reponsive_caption_update(titlebar);
        sprintf(filename, "BuddhaBrot-MT-status-thread%02i.bin", td_i);
        FILE* status_file;

        if ((status_file = fopen(filename, "wb")) == NULL) {
            return (0);
        }

        fwrite(R[td_i], sizeof(unsigned int), Rw * Rh, status_file);
        fclose(status_file);
    }

    td_pause = 0;
    #pragma omp flush(td_pause)
    return (1);
}

int load_status_files()
{
    reponsive_caption_update("BuddhaBrot-MT: loading status: pausing calculation threads...");
    pause_calcthreads_and_wait();
    load_param_file(0, 1, 0);

    for (int td_i = 0; td_i < td_nb; td_i += 1) {
        sprintf(titlebar, "BuddhaBrot-MT: loading status: loading render count matrix of thread %i...", td_i);
        reponsive_caption_update(titlebar);
        sprintf(filename, "BuddhaBrot-MT-status-thread%02i.bin", td_i);
        FILE* status_file;

        if ((status_file = fopen(filename, "rb")) == NULL) {
            return (0);
        }

        if (fread(R[td_i], sizeof(unsigned int), Rw * Rh, status_file)) {}

        fclose(status_file);
    }

    td_pause = 0;
    #pragma omp flush(td_pause)
    return (1);
}

int load_status_files_minmem()
{
    reponsive_caption_update("BuddhaBrot-MT: loading status: pausing calculation threads...");
    pause_calcthreads_and_wait();
    int original_num_calcthreads = load_param_file(0, 1, 1);
    unsigned int* tmp_R = (unsigned int*)calloc(Rw * Rh, sizeof(unsigned int));

    for (int td_i = 0; td_i < 3; td_i++) {
        sprintf(titlebar, "BuddhaBrot-MT: loading status: loading render count matrix of thread %i...", td_i);
        reponsive_caption_update(titlebar);
        sprintf(filename, "BuddhaBrot-MT-status-thread%02i.bin", td_i);
        FILE* status_file;

        if ((status_file = fopen(filename, "rb")) == NULL) {
            return (0);
        }

        if (fread(R[td_i], sizeof(unsigned int), Rw * Rh, status_file)) {}

        fclose(status_file);
    }

    for (int td_i = 3; td_i < original_num_calcthreads; td_i++) {
        sprintf(titlebar, "BuddhaBrot-MT: loading status: loading render count matrix of thread %i...", td_i);
        reponsive_caption_update(titlebar);
        sprintf(filename, "BuddhaBrot-MT-status-thread%02i.bin", td_i);
        FILE* status_file;

        if ((status_file = fopen(filename, "rb")) == NULL) {
            return (0);
        }

        if (fread(tmp_R, sizeof(unsigned int), Rw * Rh, status_file)) {}

        for (int Ri = 0; Ri < Rw * Rh; Ri++) {
            R[td_i % 3][Ri] += tmp_R[Ri];
        }

        Pp[td_i % 3] += Pp[td_i];
        Pp[td_i] = 0;
        fclose(status_file);
    }

    free(tmp_R);
    td_pause = 0;
    #pragma omp flush(td_pause)
    return (1);
}

void writeRtoPNG(const char* filename)
{
    reponsive_caption_update("BuddhaBrot-MT: writing to 8-bit png: pausing calculation threads...");
    pause_calcthreads_and_wait();

    if (lr_mode == 0) {
        Rlrmax[0] = 0;

        for (int Ri = 0; Ri < Rw * Rh; Ri++) {
            unsigned int sum = 0;

            for (int td_i = 0; td_i < td_nb; td_i += 1) {
                sum += R[td_i][Ri];
            }

            if (cm[0] == 2 || cm[0] == 3) {
                if (sum >= cm_log[0]) {
                    sum = cm_log[0] + (unsigned int) log((double)(sum - cm_log[0] + 1));
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

        for (int Ri = 0; Ri < Rw * Rh; Ri++) {
            unsigned int sum = 0;

            for (int td_i = 0; td_i < td_nb; td_i += 1) {
                sum += R[td_i][Ri];
            }

            if (cm[0] == 2 || cm[0] == 3) {
                if (sum >= cm_log[0]) {
                    sum = cm_log[0] + (unsigned int) log((double)(sum - cm_log[0] + 1));
                }
            }

            H[0][sum]++;
        }

        if (cm[0] == 0 || cm[0] == 2) {
            cm0n[0] = 0;

            for (unsigned int i = 0; i <= Rlrmax[0]; i++) {
                if (H[0][i] > 0) {
                    H[0][i] = cm0n[0]++;
                }
            }
        }

        if (cm[0] == 1 || cm[0] == 3) {
            cm1n[0] = Rw * Rh - H[0][0];
            H[0][0] = 0;

            for (unsigned int i = 2; i <= Rlrmax[0]; i++) {
                H[0][i] += H[0][i - 1];
            }
        }
    }

    if (lr_mode == 1) {
        for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
            Rlrmax[layer_iter] = 0;

            for (int Ri = 0; Ri < Rw * Rh; Ri++) {
                unsigned int sum = 0;

                for (int td_i = 0; td_i < td_nb; td_i += 1) {
                    sum += R[td_i][Ri];
                }

                if (cm[layer_iter] == 2 || cm[layer_iter] == 3) {
                    if (sum >= cm_log[layer_iter]) {
                        sum = cm_log[layer_iter] + (unsigned int) log((double)(sum - cm_log[layer_iter] + 1));
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

            for (int Ri = 0; Ri < Rw * Rh; Ri++) {
                unsigned int sum = 0;

                for (int td_i = 0; td_i < td_nb; td_i += 1) {
                    sum += R[td_i][Ri];
                }

                if (cm[layer_iter] == 2 || cm[layer_iter] == 3) {
                    if (sum >= cm_log[layer_iter]) {
                        sum = cm_log[layer_iter] + (unsigned int) log((double)(sum - cm_log[layer_iter] + 1));
                    }
                }

                H[layer_iter][sum]++;
            }

            if (cm[layer_iter] == 0 || cm[layer_iter] == 2) {
                cm0n[layer_iter] = 0;

                for (unsigned int i = 0; i <= Rlrmax[layer_iter]; i++) {
                    if (H[layer_iter][i] > 0) {
                        H[layer_iter][i] = cm0n[layer_iter]++;
                    }
                }
            }

            if (cm[layer_iter] == 1 || cm[layer_iter] == 3) {
                cm1n[layer_iter] = Rw * Rh - H[layer_iter][0];
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

            for (int Ri = 0; Ri < Rw * Rh; Ri++) {
                unsigned int sum = 0;

                for (int td_i = layer_iter; td_i < td_nb; td_i += LR_NB) {
                    sum += R[td_i][Ri];
                }

                if (cm[layer_iter] == 2 || cm[layer_iter] == 3) {
                    if (sum >= cm_log[layer_iter]) {
                        sum = cm_log[layer_iter] + (unsigned int) log((double)(sum - cm_log[layer_iter] + 1));
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

            for (int Ri = 0; Ri < Rw * Rh; Ri++) {
                unsigned int sum = 0;

                for (int td_i = layer_iter; td_i < td_nb; td_i += LR_NB) {
                    sum += R[td_i][Ri];
                }

                if (cm[layer_iter] == 2 || cm[layer_iter] == 3) {
                    if (sum >= cm_log[layer_iter]) {
                        sum = cm_log[layer_iter] + (unsigned int) log((double)(sum - cm_log[layer_iter] + 1));
                    }
                }

                H[layer_iter][sum]++;
            }

            if (cm[layer_iter] == 0 || cm[layer_iter] == 2) {
                cm0n[layer_iter] = 0;

                for (unsigned int i = 0; i <= Rlrmax[layer_iter]; i++) {
                    if (H[layer_iter][i] > 0) {
                        H[layer_iter][i] = cm0n[layer_iter]++;
                    }
                }
            }

            if (cm[layer_iter] == 1 || cm[layer_iter] == 3) {
                cm1n[layer_iter] = Rw * Rh - H[layer_iter][0];
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
        if (cm[0] == 0 || cm[0] == 2) {
            for (int y = 0; y < Rh; y++) {
                sprintf(titlebar, "BuddhaBrot-MT: writing to 8-bit png: %.0f%% done...", (double)y / Rh * 100);
                reponsive_caption_update(titlebar);

                for (int x = 0; x < Rw; x++) {
                    int Ri = y * Rw + x;
                    unsigned int sum = 0;

                    for (int td_i = 0; td_i < td_nb; td_i += 1) {
                        sum += R[td_i][Ri];
                    }

                    if (cm[0] == 2) {
                        if (sum >= cm_log[0]) {
                            sum = cm_log[0] + (unsigned int) log((double)(sum - cm_log[0] + 1));
                        }
                    }

                    int ct_i = 0;

                    if (cm0n[0] > 1) {
                        ct_i = ct_e * ((double)H[0][sum] / (cm0n[0] - 1));
                    }

                    ct_i = ct_cycle(ct_i + ct_o[0]);
                    row_buffer[x * 3 + 0] = CT[ct_i][0];
                    row_buffer[x * 3 + 1] = CT[ct_i][1];
                    row_buffer[x * 3 + 2] = CT[ct_i][2];
                }

                png_write_row(png_ptr, row_buffer);
            }
        }

        if (cm[0] == 1 || cm[0] == 3) {
            for (int y = 0; y < Rh; y++) {
                sprintf(titlebar, "BuddhaBrot-MT: writing to 8-bit png: %.0f%% done...", (double)y / Rh * 100);
                reponsive_caption_update(titlebar);

                for (int x = 0; x < Rw; x++) {
                    int Ri = y * Rw + x;
                    unsigned int sum = 0;

                    for (int td_i = 0; td_i < td_nb; td_i += 1) {
                        sum += R[td_i][Ri];
                    }

                    if (cm[0] == 3) {
                        if (sum >= cm_log[0]) {
                            sum = cm_log[0] + (unsigned int) log((double)(sum - cm_log[0] + 1));
                        }
                    }

                    int ct_i = 0;

                    if (cm1n[0] > 0) {
                        ct_i = ct_e * ((double)H[0][sum] / cm1n[0]);
                    }

                    ct_i = ct_cycle(ct_i + ct_o[0]);
                    row_buffer[x * 3 + 0] = CT[ct_i][0];
                    row_buffer[x * 3 + 1] = CT[ct_i][1];
                    row_buffer[x * 3 + 2] = CT[ct_i][2];
                }

                png_write_row(png_ptr, row_buffer);
            }
        }
    }

    if (lr_mode == 1) {
        for (int y = 0; y < Rh; y++) {
            sprintf(titlebar, "BuddhaBrot-MT: writing to 8-bit png: %.0f%% done...", (double)y / Rh * 100);
            reponsive_caption_update(titlebar);

            for (int x = 0; x < Rw; x++) {
                int Ri = y * Rw + x;
                int ct_i[LR_NB] = {0};

                for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
                    unsigned int sum = 0;

                    for (int td_i = 0; td_i < td_nb; td_i += 1) {
                        sum += R[td_i][Ri];
                    }

                    if (cm[layer_iter] == 2 || cm[layer_iter] == 3) {
                        if (sum >= cm_log[layer_iter]) {
                            sum = cm_log[layer_iter] + (unsigned int) log((double)(sum - cm_log[layer_iter] + 1));
                        }
                    }

                    if (cm[layer_iter] == 0 || cm[layer_iter] == 2) {
                        if (cm0n[layer_iter] > 1) {
                            ct_i[layer_iter] = ct_e * ((double)H[layer_iter][sum] / (cm0n[layer_iter] - 1));
                        }
                    }

                    if (cm[layer_iter] == 1 || cm[layer_iter] == 3) {
                        if (cm1n[layer_iter] > 0) {
                            ct_i[layer_iter] = ct_e * ((double)H[layer_iter][sum] / cm1n[layer_iter]);
                        }
                    }

                    ct_i[layer_iter] = ct_cycle(ct_i[layer_iter] + ct_o[layer_iter]);
                }

                row_buffer[x * 3 + 0] = CT[ct_i[0]][0];
                row_buffer[x * 3 + 1] = CT[ct_i[1]][1];
                row_buffer[x * 3 + 2] = CT[ct_i[2]][2];
            }

            png_write_row(png_ptr, row_buffer);
        }
    }

    if (lr_mode == 2) {
        for (int y = 0; y < Rh; y++) {
            sprintf(titlebar, "BuddhaBrot-MT: writing to 8-bit png: %.0f%% done...", (double)y / Rh * 100);
            reponsive_caption_update(titlebar);

            for (int x = 0; x < Rw; x++) {
                int Ri = y * Rw + x;
                int ct_i[LR_NB] = {0};

                for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
                    unsigned int sum = 0;

                    for (int td_i = layer_iter; td_i < td_nb; td_i += LR_NB) {
                        sum += R[td_i][Ri];
                    }

                    if (cm[layer_iter] == 2 || cm[layer_iter] == 3) {
                        if (sum >= cm_log[layer_iter]) {
                            sum = cm_log[layer_iter] + (unsigned int) log((double)(sum - cm_log[layer_iter] + 1));
                        }
                    }

                    if (cm[layer_iter] == 0 || cm[layer_iter] == 2) {
                        if (cm0n[layer_iter] > 1) {
                            ct_i[layer_iter] = ct_e * ((double)H[layer_iter][sum] / (cm0n[layer_iter] - 1));
                        }
                    }

                    if (cm[layer_iter] == 1 || cm[layer_iter] == 3) {
                        if (cm1n[layer_iter] > 0) {
                            ct_i[layer_iter] = ct_e * ((double)H[layer_iter][sum] / cm1n[layer_iter]);
                        }
                    }

                    ct_i[layer_iter] = ct_cycle(ct_i[layer_iter] + ct_o[layer_iter]);
                }

                row_buffer[x * 3 + 0] = CT[ct_i[0]][0];
                row_buffer[x * 3 + 1] = CT[ct_i[1]][1];
                row_buffer[x * 3 + 2] = CT[ct_i[2]][2];
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
}

void writeTtoPNG(const char* filename, int local_png_offset_x, int local_png_offset_y, int local_png_width, int local_png_height)
{
    reponsive_caption_update("BuddhaBrot-MT: writing to 8-bit png: pausing calculation threads...");
    pause_calcthreads_and_wait();

    for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
        T[layer_iter] = (unsigned int*)calloc(local_png_width * local_png_height, sizeof(unsigned int));
    }

    if (lr_mode == 0) {
        Rlrmax[0] = 0;

        for (int png_y = 0, Ry = local_png_offset_y; png_y < local_png_height; png_y++, Ry++) {
            for (int png_x = 0, Rx = local_png_offset_x; png_x < local_png_width; png_x++, Rx++) {
                int Ri = Ry * Rw + Rx;
                unsigned int sum = 0;

                for (int td_i = 0; td_i < td_nb; td_i += 1) {
                    sum += R[td_i][Ri];
                }

                if (cm[0] == 2 || cm[0] == 3) {
                    if (sum >= cm_log[0]) {
                        sum = cm_log[0] + (unsigned int) log((double)(sum - cm_log[0] + 1));
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

        if (cm[0] == 0 || cm[0] == 2) {
            cm0n[0] = 0;

            for (unsigned int i = 0; i <= Rlrmax[0]; i++) {
                if (H[0][i] > 0) {
                    H[0][i] = cm0n[0]++;
                }
            }
        }

        if (cm[0] == 1 || cm[0] == 3) {
            cm1n[0] = local_png_width * local_png_height - H[0][0];
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
                    int Ri = Ry * Rw + Rx;
                    unsigned int sum = 0;

                    for (int td_i = 0; td_i < td_nb; td_i += 1) {
                        sum += R[td_i][Ri];
                    }

                    if (cm[layer_iter] == 2 || cm[layer_iter] == 3) {
                        if (sum >= cm_log[layer_iter]) {
                            sum = cm_log[layer_iter] + (unsigned int) log((double)(sum - cm_log[layer_iter] + 1));
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

            if (cm[layer_iter] == 0 || cm[layer_iter] == 2) {
                cm0n[layer_iter] = 0;

                for (unsigned int i = 0; i <= Rlrmax[layer_iter]; i++) {
                    if (H[layer_iter][i] > 0) {
                        H[layer_iter][i] = cm0n[layer_iter]++;
                    }
                }
            }

            if (cm[layer_iter] == 1 || cm[layer_iter] == 3) {
                cm1n[layer_iter] = local_png_width * local_png_height - H[layer_iter][0];
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
                    int Ri = Ry * Rw + Rx;
                    unsigned int sum = 0;

                    for (int td_i = layer_iter; td_i < td_nb; td_i += LR_NB) {
                        sum += R[td_i][Ri];
                    }

                    if (cm[layer_iter] == 2 || cm[layer_iter] == 3) {
                        if (sum >= cm_log[layer_iter]) {
                            sum = cm_log[layer_iter] + (unsigned int) log((double)(sum - cm_log[layer_iter] + 1));
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

            if (cm[layer_iter] == 0 || cm[layer_iter] == 2) {
                cm0n[layer_iter] = 0;

                for (unsigned int i = 0; i <= Rlrmax[layer_iter]; i++) {
                    if (H[layer_iter][i] > 0) {
                        H[layer_iter][i] = cm0n[layer_iter]++;
                    }
                }
            }

            if (cm[layer_iter] == 1 || cm[layer_iter] == 3) {
                cm1n[layer_iter] = local_png_width * local_png_height - H[layer_iter][0];
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
        if (cm[0] == 0 || cm[0] == 2) {
            for (int png_y = 0; png_y < local_png_height; png_y++) {
                sprintf(titlebar, "BuddhaBrot-MT: writing to 8-bit png: %.0f%% done...", (double)png_y / local_png_height * 100);
                reponsive_caption_update(titlebar);

                for (int png_x = 0; png_x < local_png_width; png_x++) {
                    int ct_i = 0;

                    if (cm0n[0] > 1) {
                        ct_i = ct_e * ((double)H[0][T[0][png_y * local_png_width + png_x]] / (cm0n[0] - 1));
                    }

                    ct_i = ct_cycle(ct_i + ct_o[0]);
                    row_buffer[png_x * 3 + 0] = CT[ct_i][0];
                    row_buffer[png_x * 3 + 1] = CT[ct_i][1];
                    row_buffer[png_x * 3 + 2] = CT[ct_i][2];
                }

                png_write_row(png_ptr, row_buffer);
            }
        }

        if (cm[0] == 1 || cm[0] == 3) {
            for (int png_y = 0; png_y < local_png_height; png_y++) {
                sprintf(titlebar, "BuddhaBrot-MT: writing to 8-bit png: %.0f%% done...", (double)png_y / Rh * 100);
                reponsive_caption_update(titlebar);

                for (int png_x = 0; png_x < local_png_width; png_x++) {
                    int ct_i = 0;

                    if (cm1n[0] > 0) {
                        ct_i = ct_e * ((double)H[0][T[0][png_y * local_png_width + png_x]] / cm1n[0]);
                    }

                    ct_i = ct_cycle(ct_i + ct_o[0]);
                    row_buffer[png_x * 3 + 0] = CT[ct_i][0];
                    row_buffer[png_x * 3 + 1] = CT[ct_i][1];
                    row_buffer[png_x * 3 + 2] = CT[ct_i][2];
                }

                png_write_row(png_ptr, row_buffer);
            }
        }
    }

    if (lr_mode == 1) {
        for (int png_y = 0; png_y < local_png_height; png_y++) {
            sprintf(titlebar, "BuddhaBrot-MT: writing to 8-bit png: %.0f%% done...", (double)png_y / Rh * 100);
            reponsive_caption_update(titlebar);

            for (int png_x = 0; png_x < local_png_width; png_x++) {
                int ct_i[LR_NB] = {0};

                for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
                    if (cm[layer_iter] == 0 || cm[layer_iter] == 2) {
                        if (cm0n[layer_iter] > 1) {
                            ct_i[layer_iter] = ct_e * ((double)H[layer_iter][T[layer_iter][png_y * local_png_width + png_x]] / (cm0n[layer_iter] - 1));
                        }
                    }

                    if (cm[layer_iter] == 1 || cm[layer_iter] == 3) {
                        if (cm1n[layer_iter] > 0) {
                            ct_i[layer_iter] = ct_e * ((double)H[layer_iter][T[layer_iter][png_y * local_png_width + png_x]] / cm1n[layer_iter]);
                        }
                    }

                    ct_i[layer_iter] = ct_cycle(ct_i[layer_iter] + ct_o[layer_iter]);
                }

                row_buffer[png_x * 3 + 0] = CT[ct_i[0]][0];
                row_buffer[png_x * 3 + 1] = CT[ct_i[1]][1];
                row_buffer[png_x * 3 + 2] = CT[ct_i[2]][2];
            }

            png_write_row(png_ptr, row_buffer);
        }
    }

    if (lr_mode == 2) {
        for (int png_y = 0; png_y < local_png_height; png_y++) {
            sprintf(titlebar, "BuddhaBrot-MT: writing to 8-bit png: %.0f%% done...", (double)png_y / Rh * 100);
            reponsive_caption_update(titlebar);

            for (int png_x = 0; png_x < local_png_width; png_x++) {
                int ct_i[LR_NB] = {0};

                for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
                    if (cm[layer_iter] == 0 || cm[layer_iter] == 2) {
                        if (cm0n[layer_iter] > 1) {
                            ct_i[layer_iter] = ct_e * ((double)H[layer_iter][T[layer_iter][png_y * local_png_width + png_x]] / (cm0n[layer_iter] - 1));
                        }
                    }

                    if (cm[layer_iter] == 1 || cm[layer_iter] == 3) {
                        if (cm1n[layer_iter] > 0) {
                            ct_i[layer_iter] = ct_e * ((double)H[layer_iter][T[layer_iter][png_y * local_png_width + png_x]] / cm1n[layer_iter]);
                        }
                    }

                    ct_i[layer_iter] = ct_cycle(ct_i[layer_iter] + ct_o[layer_iter]);
                }

                row_buffer[png_x * 3 + 0] = CT[ct_i[0]][0];
                row_buffer[png_x * 3 + 1] = CT[ct_i[1]][1];
                row_buffer[png_x * 3 + 2] = CT[ct_i[2]][2];
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
}

void sdl_message_check()
{
    recalc_WH_if_paused = 0;

    if (sdl_event.type == SDL_QUIT) {
        td_stop = 1;
        #pragma omp flush(td_stop)
    }

    if (sdl_event.type == SDL_KEYDOWN) {
        //// switch title bar display
        if (sdl_event.key.keysym.sym == SDLK_ESCAPE && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            titlebarswitch = 1 - titlebarswitch;
        }

        //// saving status
        if (sdl_event.key.keysym.sym == SDLK_F9 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            save_status_files();
        }

        //// loading status
        if (sdl_event.key.keysym.sym == SDLK_F10 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_status_files();
        }

        if (sdl_event.key.keysym.sym == SDLK_F10 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_param_file(1, 0, 0);
        }

        if (sdl_event.key.keysym.sym == SDLK_F10 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_status_files_minmem();
        }

        //// td_pause
        if (sdl_event.key.keysym.sym == SDLK_F11 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            td_pause = 1 - td_pause;
            #pragma omp flush(td_pause)
        }

        //// increase threads
        if (sdl_event.key.keysym.sym == SDLK_F11 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            increase_num_calcthreads();
        }

        //// decrease threads
        if (sdl_event.key.keysym.sym == SDLK_F11 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            decrease_num_calcthreads();
        }

        //// set fps at 1 frame per second
        if (sdl_event.key.keysym.sym == SDLK_F12 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            fps = 1.0;
        }

        //// set fps at 10 frames per second
        if (sdl_event.key.keysym.sym == SDLK_F12 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            fps = 10.0;
        }

        //// set fps at 30 frames per second
        if (sdl_event.key.keysym.sym == SDLK_F12 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            fps = 30.0;
        }

        //// change lr_mode
        if (sdl_event.key.keysym.sym == SDLK_F1 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            lr_mode = (lr_mode + 1) % LR_MODE_NB;
            recalc_WH_if_paused = 1;
        }

        //// change ct_type noncycle
        if (sdl_event.key.keysym.sym == SDLK_F2 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_type = 10 + MAX(ct_type - 9, 0) % CT_NONCYCLE_NB;
            ct_load(ct_type);
        }

        //// change ct_type cycle
        if (sdl_event.key.keysym.sym == SDLK_F3 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_type = 20 +  MAX(ct_type - 19, 0) % CT_CYCLE_NB;
            ct_load(ct_type);
        }

        //// change bb_type layer 123
        if (sdl_event.key.keysym.sym == SDLK_BACKQUOTE && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, (bb_type[0] + 1) % BB_TYPE_NB, bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0]);
        }

        //// change bb_bail layer 123
        if (sdl_event.key.keysym.sym == SDLK_1 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0] + 1, bb_pps[0], bb_ppe[0], bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_1 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0] * 10, bb_pps[0], bb_ppe[0], bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_1 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0] - 1, bb_pps[0], bb_ppe[0], bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_1 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0] / 10, bb_pps[0], bb_ppe[0], bb_minn[0]);
        }

        //// change bb_bail layer 1
        if (sdl_event.key.keysym.sym == SDLK_q && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0] + 1, bb_pps[0], bb_ppe[0], bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_q && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0] * 10, bb_pps[0], bb_ppe[0], bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_q && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0] - 1, bb_pps[0], bb_ppe[0], bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_q && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0] / 10, bb_pps[0], bb_ppe[0], bb_minn[0]);
        }

        //// change bb_bail layer 2
        if (sdl_event.key.keysym.sym == SDLK_a && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1] + 1, bb_pps[1], bb_ppe[1], bb_minn[1]);
        }

        if (sdl_event.key.keysym.sym == SDLK_a && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1] * 10, bb_pps[1], bb_ppe[1], bb_minn[1]);
        }

        if (sdl_event.key.keysym.sym == SDLK_a && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1] - 1, bb_pps[1], bb_ppe[1], bb_minn[1]);
        }

        if (sdl_event.key.keysym.sym == SDLK_a && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1] / 10, bb_pps[1], bb_ppe[1], bb_minn[1]);
        }

        //// change bb_bail layer 3
        if (sdl_event.key.keysym.sym == SDLK_z && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2] + 1, bb_pps[2], bb_ppe[2], bb_minn[2]);
        }

        if (sdl_event.key.keysym.sym == SDLK_z && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2] * 10, bb_pps[2], bb_ppe[2], bb_minn[2]);
        }

        if (sdl_event.key.keysym.sym == SDLK_z && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2] - 1, bb_pps[2], bb_ppe[2], bb_minn[2]);
        }

        if (sdl_event.key.keysym.sym == SDLK_z && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2] / 10, bb_pps[2], bb_ppe[2], bb_minn[2]);
        }

        //// change bb_pps layer 123
        if (sdl_event.key.keysym.sym == SDLK_2 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0], bb_pps[0] + 1, bb_ppe[0], bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_2 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0], bb_pps[0] * 10, bb_ppe[0], bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_2 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0], bb_pps[0] - 1, bb_ppe[0], bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_2 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0], bb_pps[0] / 10, bb_ppe[0], bb_minn[0]);
        }

        //// change bb_pps layer 1
        if (sdl_event.key.keysym.sym == SDLK_w && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0], bb_pps[0] + 1, bb_ppe[0], bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_w && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0], bb_pps[0] * 10, bb_ppe[0], bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_w && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0], bb_pps[0] - 1, bb_ppe[0], bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_w && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0], bb_pps[0] / 10, bb_ppe[0], bb_minn[0]);
        }

        //// change bb_pps layer 2
        if (sdl_event.key.keysym.sym == SDLK_s && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1], bb_pps[1] + 1, bb_ppe[1], bb_minn[1]);
        }

        if (sdl_event.key.keysym.sym == SDLK_s && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1], bb_pps[1] * 10, bb_ppe[1], bb_minn[1]);
        }

        if (sdl_event.key.keysym.sym == SDLK_s && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1], bb_pps[1] - 1, bb_ppe[1], bb_minn[1]);
        }

        if (sdl_event.key.keysym.sym == SDLK_s && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1], bb_pps[1] / 10, bb_ppe[1], bb_minn[1]);
        }

        //// change bb_pps layer 3
        if (sdl_event.key.keysym.sym == SDLK_x && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2], bb_pps[2] + 1, bb_ppe[2], bb_minn[2]);
        }

        if (sdl_event.key.keysym.sym == SDLK_x && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2], bb_pps[2] * 10, bb_ppe[2], bb_minn[2]);
        }

        if (sdl_event.key.keysym.sym == SDLK_x && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2], bb_pps[2] - 1, bb_ppe[2], bb_minn[2]);
        }

        if (sdl_event.key.keysym.sym == SDLK_x && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2], bb_pps[2] / 10, bb_ppe[2], bb_minn[2]);
        }

        //// change bb_ppe layer 123
        if (sdl_event.key.keysym.sym == SDLK_3 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0] + 1, bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_3 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0] * 10, bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_3 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0] - 1, bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_3 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0] / 10, bb_minn[0]);
        }

        //// change bb_ppe layer 1
        if (sdl_event.key.keysym.sym == SDLK_e && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0] + 1, bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_e && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0] * 10, bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_e && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0] - 1, bb_minn[0]);
        }

        if (sdl_event.key.keysym.sym == SDLK_e && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0] / 10, bb_minn[0]);
        }

        //// change bb_ppe layer 2
        if (sdl_event.key.keysym.sym == SDLK_d && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1] + 1, bb_minn[1]);
        }

        if (sdl_event.key.keysym.sym == SDLK_d && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1] * 10, bb_minn[1]);
        }

        if (sdl_event.key.keysym.sym == SDLK_d && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1] - 1, bb_minn[1]);
        }

        if (sdl_event.key.keysym.sym == SDLK_d && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1] / 10, bb_minn[1]);
        }

        //// change bb_ppe layer 3
        if (sdl_event.key.keysym.sym == SDLK_c && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2] + 1, bb_minn[2]);
        }

        if (sdl_event.key.keysym.sym == SDLK_c && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2] * 10, bb_minn[2]);
        }

        if (sdl_event.key.keysym.sym == SDLK_c && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2] - 1, bb_minn[2]);
        }

        if (sdl_event.key.keysym.sym == SDLK_c && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2] / 10, bb_minn[2]);
        }

        //// change bb_minn layer 123
        if (sdl_event.key.keysym.sym == SDLK_4 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0] + 1);
        }

        if (sdl_event.key.keysym.sym == SDLK_4 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0] * 10);
        }

        if (sdl_event.key.keysym.sym == SDLK_4 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0] - 1);
        }

        if (sdl_event.key.keysym.sym == SDLK_4 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(-1, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0] / 10);
        }

        //// change bb_minn layer 1
        if (sdl_event.key.keysym.sym == SDLK_r && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0] + 1);
        }

        if (sdl_event.key.keysym.sym == SDLK_r && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0] * 10);
        }

        if (sdl_event.key.keysym.sym == SDLK_r && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0] - 1);
        }

        if (sdl_event.key.keysym.sym == SDLK_r && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(0, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0] / 10);
        }

        //// change bb_minn layer 2
        if (sdl_event.key.keysym.sym == SDLK_f && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1] + 1);
        }

        if (sdl_event.key.keysym.sym == SDLK_f && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1] * 10);
        }

        if (sdl_event.key.keysym.sym == SDLK_f && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1] - 1);
        }

        if (sdl_event.key.keysym.sym == SDLK_f && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(1, bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1] , bb_minn[1] / 10);
        }

        //// change bb_minn layer 3
        if (sdl_event.key.keysym.sym == SDLK_v && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2] + 1);
        }

        if (sdl_event.key.keysym.sym == SDLK_v && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2] * 10);
        }

        if (sdl_event.key.keysym.sym == SDLK_v && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2] - 1);
        }

        if (sdl_event.key.keysym.sym == SDLK_v && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_bb_param(2, bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2] / 10);
        }

        //// change normal cm layer 123
        if (sdl_event.key.keysym.sym == SDLK_5 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            cm[0] = (cm[0] + 1) % CM_NORMAL_NB;
            cm[1] = cm[0];
            cm[2] = cm[0];
            recalc_WH_if_paused = 1;
        }

        //// change normal cm layer 1
        if (sdl_event.key.keysym.sym == SDLK_t&& !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            cm[0] = (cm[0] + 1) % CM_NORMAL_NB;
            recalc_WH_if_paused = 1;
        }

        //// change normal cm layer 2
        if (sdl_event.key.keysym.sym == SDLK_g && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            cm[1] = (cm[1] + 1) % CM_NORMAL_NB;
            recalc_WH_if_paused = 1;
        }

        //// change normal cm layer 3
        if (sdl_event.key.keysym.sym == SDLK_b && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            cm[2] = (cm[2] + 1) % CM_NORMAL_NB;
            recalc_WH_if_paused = 1;
        }

        //// change log cm layer 123
        if (sdl_event.key.keysym.sym == SDLK_5 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            cm[0] = 2 + MAX(cm[0] - 1, 0) % CM_LOG_NB;
            cm[1] = cm[0];
            cm[2] = cm[0];
            recalc_WH_if_paused = 1;
        }

        //// change log cm layer 1
        if (sdl_event.key.keysym.sym == SDLK_t&& (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            cm[0] = 2 + MAX(cm[0] - 1, 0) % CM_LOG_NB;
            recalc_WH_if_paused = 1;
        }

        //// change log cm layer 2
        if (sdl_event.key.keysym.sym == SDLK_g && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            cm[1] = 2 + MAX(cm[1] - 1, 0) % CM_LOG_NB;
            recalc_WH_if_paused = 1;
        }

        //// change log cm layer 3
        if (sdl_event.key.keysym.sym == SDLK_b && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            cm[2] = 2 + MAX(cm[2] - 1, 0) % CM_LOG_NB;
            recalc_WH_if_paused = 1;
        }

        //// change cm_log layer 123
        if (sdl_event.key.keysym.sym == SDLK_5 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            cm_log[0] += 1;
            cm_log[1] = cm_log[0];
            cm_log[2] = cm_log[0];
            recalc_WH_if_paused = 1;
        }

        if (sdl_event.key.keysym.sym == SDLK_5 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            if (cm_log[0] != 0) {
                cm_log[0] -= 1;
            }

            cm_log[1] = cm_log[0];
            cm_log[2] = cm_log[0];
            recalc_WH_if_paused = 1;
        }

        //// change cm_log layer 1
        if (sdl_event.key.keysym.sym == SDLK_t&& !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            cm_log[0] += 1;
            recalc_WH_if_paused = 1;
        }

        if (sdl_event.key.keysym.sym == SDLK_t&& (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            if (cm_log[0] != 0) {
                cm_log[0] -= 1;
            }

            recalc_WH_if_paused = 1;
        }

        //// change cm_log layer 2
        if (sdl_event.key.keysym.sym == SDLK_g && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            cm_log[1] += 1;
            recalc_WH_if_paused = 1;
        }

        if (sdl_event.key.keysym.sym == SDLK_g && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            if (cm_log[1] != 0) {
                cm_log[1] -= 1;
            }

            recalc_WH_if_paused = 1;
        }

        //// change cm_log layer 3
        if (sdl_event.key.keysym.sym == SDLK_b && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            cm_log[2] += 1;
            recalc_WH_if_paused = 1;
        }

        if (sdl_event.key.keysym.sym == SDLK_b && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            if (cm_log[2] != 0) {
                cm_log[2] -= 1;
            }

            recalc_WH_if_paused = 1;
        }

        //// change ct_o layer 123
        if (sdl_event.key.keysym.sym == SDLK_6 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_o[0] = ct_cycle(ct_o[0] + 1);
            ct_o[1] = ct_o[0];
            ct_o[2] = ct_o[0];
        }

        if (sdl_event.key.keysym.sym == SDLK_6 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_o[0] = ct_cycle(ct_o[0] + 10);
            ct_o[1] = ct_o[0];
            ct_o[2] = ct_o[0];
        }

        if (sdl_event.key.keysym.sym == SDLK_6 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_o[0] = 0;
            ct_o[1] = 0;
            ct_o[2] = 0;
        }

        //// change ct_o layer 1
        if (sdl_event.key.keysym.sym == SDLK_y && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_o[0] = ct_cycle(ct_o[0] + 1);
        }

        if (sdl_event.key.keysym.sym == SDLK_y && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_o[0] = ct_cycle(ct_o[0] + 10);
        }

        if (sdl_event.key.keysym.sym == SDLK_y && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_o[0] = 0;
        }

        //// change ct_o layer 2
        if (sdl_event.key.keysym.sym == SDLK_h && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_o[1] = ct_cycle(ct_o[1] + 1);
        }

        if (sdl_event.key.keysym.sym == SDLK_h && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_o[1] = ct_cycle(ct_o[1] + 10);
        }

        if (sdl_event.key.keysym.sym == SDLK_h && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_o[1] = 0;
        }

        //// change ct_o layer 3
        if (sdl_event.key.keysym.sym == SDLK_n && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_o[2] = ct_cycle(ct_o[2] + 1);
        }

        if (sdl_event.key.keysym.sym == SDLK_n && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_o[2] = ct_cycle(ct_o[2] + 10);
        }

        if (sdl_event.key.keysym.sym == SDLK_n && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_o[2] = 0;
        }

        //// change ct_v layer 123
        if (sdl_event.key.keysym.sym == SDLK_7 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_v[0]++;
            ct_v[1] = ct_v[0];
            ct_v[2] = ct_v[0];
        }

        if (sdl_event.key.keysym.sym == SDLK_7 && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_v[0]--;
            ct_v[1] = ct_v[0];
            ct_v[2] = ct_v[0];
        }

        if (sdl_event.key.keysym.sym == SDLK_7 && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_v[0] = 0;
            ct_v[1] = 0;
            ct_v[2] = 0;
        }

        //// change ct_v layer 1
        if (sdl_event.key.keysym.sym == SDLK_u && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_v[0]++;
        }

        if (sdl_event.key.keysym.sym == SDLK_u && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_v[0]--;
        }

        if (sdl_event.key.keysym.sym == SDLK_u && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_v[0] = 0;
        }

        //// change ct_v layer 2
        if (sdl_event.key.keysym.sym == SDLK_j && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_v[1]++;
        }

        if (sdl_event.key.keysym.sym == SDLK_j && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_v[1]--;
        }

        if (sdl_event.key.keysym.sym == SDLK_j && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_v[1] = 0;
        }

        //// change ct_v layer 3
        if (sdl_event.key.keysym.sym == SDLK_m && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_v[2]++;
        }

        if (sdl_event.key.keysym.sym == SDLK_m && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_v[2]--;
        }

        if (sdl_event.key.keysym.sym == SDLK_m && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            ct_v[2] = 0;
        }

        //// increase render size
        if (sdl_event.key.keysym.sym == SDLK_PAGEUP && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_location(1, 4.0 / (Rr_up - Rr_lo), 0.5 * (Rr_lo + Rr_up), 0, 0.5 * (Ri_lo + Ri_up), 0, Rw / Rw_f, Rh / Rh_f);
        }

        //// decrease render size
        if (sdl_event.key.keysym.sym == SDLK_PAGEDOWN && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            if ((RWox > 0) || (RWoy > 0) || (Rw > Ww) || (Rh > Wh)) {
                load_location(1, 4.0 / (Rr_up - Rr_lo), 0.5 * (Rr_lo + Rr_up), 0, 0.5 * (Ri_lo + Ri_up), 0, Rw * Rw_f, Rh * Rh_f);
            }
        }

        //// pan window left in render, 10% of window
        if (sdl_event.key.keysym.sym == SDLK_LEFT && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            RWox -= Ww * 0.1;

            if (RWox < 0) {
                RWox = 0;
            }

            Wi_lo = Ri_lo + RWox / Rwdivi;
            Wi_up = Wi_lo + Ww / Rwdivi;
            recalc_WH_if_paused = 1;
        }

        //// pan window left in render, 1% of window
        if (sdl_event.key.keysym.sym == SDLK_LEFT && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            RWox -= Ww * 0.01;

            if (RWox < 0) {
                RWox = 0;
            }

            Wi_lo = Ri_lo + RWox / Rwdivi;
            Wi_up = Wi_lo + Ww / Rwdivi;
            recalc_WH_if_paused = 1;
        }

        //// pan window right in render, 10% of window
        if (sdl_event.key.keysym.sym == SDLK_RIGHT && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            RWox += Ww * 0.1;

            if (RWox + Ww > Rw) {
                RWox = Rw - Ww;
            }

            Wi_lo = Ri_lo + RWox / Rwdivi;
            Wi_up = Wi_lo + Ww / Rwdivi;
            recalc_WH_if_paused = 1;
        }

        //// pan window right in render, 1% of window
        if (sdl_event.key.keysym.sym == SDLK_RIGHT && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            RWox += Ww * 0.01;

            if (RWox + Ww > Rw) {
                RWox = Rw - Ww;
            }

            Wi_lo = Ri_lo + RWox / Rwdivi;
            Wi_up = Wi_lo + Ww / Rwdivi;
            recalc_WH_if_paused = 1;
        }

        //// pan window up in render, 10% of window
        if (sdl_event.key.keysym.sym == SDLK_UP && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            RWoy -= Wh * 0.1;

            if (RWoy < 0) {
                RWoy = 0;
            }

            Wr_lo = Rr_lo + RWoy / Rhdivr;
            Wr_up = Wr_lo + Wh / Rhdivr;
            recalc_WH_if_paused = 1;
        }

        //// pan window up in render, 1% of window
        if (sdl_event.key.keysym.sym == SDLK_UP && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            RWoy -= Wh * 0.01;

            if (RWoy < 0) {
                RWoy = 0;
            }

            Wr_lo = Rr_lo + RWoy / Rhdivr;
            Wr_up = Wr_lo + Wh / Rhdivr;
            recalc_WH_if_paused = 1;
        }

        //// pan window down in render, 10% of window
        if (sdl_event.key.keysym.sym == SDLK_DOWN && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            RWoy += Wh * 0.1;

            if (RWoy + Wh > Rh) {
                RWoy = Rh - Wh;
            }

            Wr_lo = Rr_lo + RWoy / Rhdivr;
            Wr_up = Wr_lo + Wh / Rhdivr;
            recalc_WH_if_paused = 1;
        }

        //// pan window down in render, 1% of window
        if (sdl_event.key.keysym.sym == SDLK_DOWN && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            RWoy += Wh * 0.01;

            if (RWoy + Wh > Rh) {
                RWoy = Rh - Wh;
            }

            Wr_lo = Rr_lo + RWoy / Rhdivr;
            Wr_up = Wr_lo + Wh / Rhdivr;
            recalc_WH_if_paused = 1;
        }

        //// zoom in fractal
        if (sdl_event.key.keysym.sym == SDLK_PAGEUP && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_location(1, 4.0 / (Rr_up - Rr_lo) / Rr_f, 0.5 * (Rr_lo + Rr_up), 0, 0.5 * (Ri_lo + Ri_up), 0, Rw, Rh);
        }

        //// zoom out fractal
        if (sdl_event.key.keysym.sym == SDLK_PAGEDOWN && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_location(1, 4.0 / (Rr_up - Rr_lo) * Rr_f, 0.5 * (Rr_lo + Rr_up), 0, 0.5 * (Ri_lo + Ri_up), 0, Rw, Rh);
        }

        //// pan fractal location left, 10% of render size
        if (sdl_event.key.keysym.sym == SDLK_LEFT && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_location(1, 4.0 / (Rr_up - Rr_lo), 0.5 * (Rr_lo + Rr_up), 0, 0.5 * (Ri_lo + Ri_up) - Ri_ra * 0.1, 0, Rw, Rh);
        }

        //// pan fractal location left, 1% of render size
        if (sdl_event.key.keysym.sym == SDLK_LEFT && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_location(1, 4.0 / (Rr_up - Rr_lo), 0.5 * (Rr_lo + Rr_up), 0, 0.5 * (Ri_lo + Ri_up) - Ri_ra * 0.01, 0, Rw, Rh);
        }

        //// pan fractal location right, 10% of render size
        if (sdl_event.key.keysym.sym == SDLK_RIGHT && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_location(1, 4.0 / (Rr_up - Rr_lo), 0.5 * (Rr_lo + Rr_up), 0, 0.5 * (Ri_lo + Ri_up) + Ri_ra * 0.1, 0, Rw, Rh);
        }

        //// pan fractal location right, 1% of render size
        if (sdl_event.key.keysym.sym == SDLK_RIGHT && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_location(1, 4.0 / (Rr_up - Rr_lo), 0.5 * (Rr_lo + Rr_up), 0, 0.5 * (Ri_lo + Ri_up) + Ri_ra * 0.01, 0, Rw, Rh);
        }

        //// pan fractal location up, 10% of render size
        if (sdl_event.key.keysym.sym == SDLK_UP && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_location(1, 4.0 / (Rr_up - Rr_lo), 0.5 * (Rr_lo + Rr_up) - Rr_ra * 0.1, 0, 0.5 * (Ri_lo + Ri_up), 0, Rw, Rh);
        }

        //// pan fractal location up, 1% of render size
        if (sdl_event.key.keysym.sym == SDLK_UP && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_location(1, 4.0 / (Rr_up - Rr_lo), 0.5 * (Rr_lo + Rr_up) - Rr_ra * 0.01, 0, 0.5 * (Ri_lo + Ri_up), 0, Rw, Rh);
        }

        //// pan fractal location down, 10% of render size
        if (sdl_event.key.keysym.sym == SDLK_DOWN && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_location(1, 4.0 / (Rr_up - Rr_lo), 0.5 * (Rr_lo + Rr_up) + Rr_ra * 0.1, 0, 0.5 * (Ri_lo + Ri_up), 0, Rw, Rh);
        }

        //// pan fractal location down, 1% of render size
        if (sdl_event.key.keysym.sym == SDLK_DOWN && (sdl_event.key.keysym.mod & KMOD_SHIFT) && (sdl_event.key.keysym.mod & KMOD_CTRL)) {
            load_location(1, 4.0 / (Rr_up - Rr_lo), 0.5 * (Rr_lo + Rr_up) + Rr_ra * 0.01, 0, 0.5 * (Ri_lo + Ri_up), 0, Rw, Rh);
        }

        //// writing window to one png
        if (sdl_event.key.keysym.sym == SDLK_BACKSPACE && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            if (lr_mode == 0) {
                sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i W(%.6f %.6f %.1f) ct%i cm%i.%i.%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Wr_lo + Wr_up), 0.5 * (Wi_lo + Wi_up), 4.0 / (Wr_up - Wr_lo), ct_type, cm[0], cm_log[0], ct_o[0], (double)Ppsum);
            }

            if (lr_mode == 1) {
                sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i W(%.6f %.6f %.1f) ct%i cm%i.%i.%i cm.%i.%i.%i cm%i.%i.%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Wr_lo + Wr_up), 0.5 * (Wi_lo + Wi_up), 4.0 / (Wr_up - Wr_lo), ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], (double)Ppsum);
            }

            if (lr_mode == 2) {
                sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i W(%.6f %.6f %.1f) ct%i cm%i.%i.%i cm.%i.%i.%i cm%i.%i.%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Wr_lo + Wr_up), 0.5 * (Wi_lo + Wi_up), 4.0 / (Wr_up - Wr_lo), ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], (double)Ppsum);
            }

            writeTtoPNG(filename, RWox, RWoy, Ww, Wh);
        }

        if (sdl_event.key.keysym.sym == SDLK_BACKSPACE && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            autoWtoPNG = (autoWtoPNG + 1) % AUTOWRITE_MODE_NB;
        }

        //// writing tiled render to multiple pngs
        if (sdl_event.key.keysym.sym == SDLK_BACKSLASH && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            if (lr_mode == 0) {
                sprintf(dirname, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i cm%i.%i.%i R%ix%i T%ix%i %g", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type, cm[0], cm_log[0], ct_o[0], Rw, Rh, Tw, Th, (double)Ppsum);
            }

            if (lr_mode == 1) {
                sprintf(dirname, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i cm%i.%i.%i cm.%i.%i.%i cm%i.%i.%i R%ix%i T%ix%i %g", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], Rw, Rh, Tw, Th, (double)Ppsum);
            }

            if (lr_mode == 2) {
                sprintf(dirname, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i cm%i.%i.%i cm.%i.%i.%i cm%i.%i.%i R%ix%i T%ix%i %g", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], Rw, Rh, Tw, Th, (double)Ppsum);
            }

            sprintf(commandname, "mkdir \"%s\"", dirname);
            system(commandname);

            for (int png_offset_x = 0; png_offset_x < Rw; png_offset_x += Tw) {
                for (int png_offset_y = 0; png_offset_y < Rh; png_offset_y += Th) {
                    double Tr_lo = Rr_lo + png_offset_y / Rhdivr;
                    double Tr_up = Tr_lo + Th / Rhdivr;
                    double Ti_lo = Ri_lo + png_offset_x / Rwdivi;
                    double Ti_up = Ti_lo + Tw / Rwdivi;
                    sprintf(filename, "%s/T%ix%i %05i %05i W(%.6f %.6f %.1f).png", dirname, Tw, Th, png_offset_y, png_offset_x, 0.5 * (Tr_lo + Tr_up), 0.5 * (Ti_lo + Ti_up), 4.0 / (Tr_up - Tr_lo));
                    writeTtoPNG(filename, png_offset_x, png_offset_y, Tw, Th);
                }
            }
        }

        if (sdl_event.key.keysym.sym == SDLK_BACKSLASH && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            autoRTtoPNG = (autoRTtoPNG + 1) % AUTOWRITE_MODE_NB;
        }

        //// writing render to one png
        if (sdl_event.key.keysym.sym == SDLK_RETURN && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            if (lr_mode == 0) {
                sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i cm%i.%i.%i R%ix%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type, cm[0], cm_log[0], ct_o[0], Rw, Rh, (double)Ppsum);
            }

            if (lr_mode == 1) {
                sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i cm%i.%i.%i cm.%i.%i.%i cm%i.%i.%i R%ix%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], Rw, Rh, (double)Ppsum);
            }

            if (lr_mode == 2) {
                sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i cm%i.%i.%i cm.%i.%i.%i cm%i.%i.%i R%ix%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], Rw, Rh, (double)Ppsum);
            }

            writeRtoPNG(filename);
        }

        if (sdl_event.key.keysym.sym == SDLK_RETURN && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            autoRtoPNG = (autoRtoPNG + 1) % AUTOWRITE_MODE_NB;
        }

        if (sdl_event.key.keysym.sym == SDLK_LEFTBRACKET && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            Tw = MIN(Tw / Rw_f, Rw);
        }

        if (sdl_event.key.keysym.sym == SDLK_LEFTBRACKET && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            Tw = MAX(Tw * Rw_f, 1000);
        }

        if (sdl_event.key.keysym.sym == SDLK_RIGHTBRACKET && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            Th = MIN(Th / Rh_f, Rh);
        }

        if (sdl_event.key.keysym.sym == SDLK_RIGHTBRACKET && (sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            Th = MAX(Th * Rh_f, 1000);
        }

        if (sdl_event.key.keysym.sym == SDLK_MINUS && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            t_autoPNG_delta = MAX(t_autoPNG_delta / 10, 6e4);
        }

        if (sdl_event.key.keysym.sym == SDLK_EQUALS && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            t_autoPNG_delta = MIN(t_autoPNG_delta * 10, 6e7);
        }

        if (sdl_event.key.keysym.sym == SDLK_SEMICOLON && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            Ppsum_autoPNG_delta = MAX(Ppsum_autoPNG_delta / 10, 1e4);
        }

        if (sdl_event.key.keysym.sym == SDLK_QUOTE && !(sdl_event.key.keysym.mod & KMOD_SHIFT) && !(sdl_event.key.keysym.mod & KMOD_CTRL)) {
            Ppsum_autoPNG_delta = Ppsum_autoPNG_delta * 10;
        }
    }
}

void visualisation_thread()
{
    printf("thread v start\n");

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        td_stop = 1;
        #pragma omp flush(td_stop)
    } else {
        sdl_window = SDL_CreateWindow("BuddhaBrot-MT", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, Ww, Wh, 0);
        sdl_renderer = SDL_CreateRenderer(sdl_window, -1, SDL_RENDERER_ACCELERATED);
        sdl_surface = SDL_CreateRGBSurface(0, Ww, Wh, 32, 0xff000000, 0x00ff0000, 0x0000ff00, 0x000000ff);
        sdl_texture = SDL_CreateTexture(sdl_renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, Ww, Wh);
        dsfmt_gv_init_gen_rand(0);

        while (1) {
            #pragma omp flush(td_stop)

            if (td_stop == 1) {
                sprintf(titlebar, "BuddhaBrot-MT: exiting: closing down %i calculation threads...", td_nb);
                reponsive_caption_update(titlebar);
                break;
            }

            Ppsum = 0;

            for (int td_i = 0; td_i < td_nb; td_i += 1) {
                Ppsum += Pp[td_i];
            }

            wait_ms(MAX(1000.0 / fps - (SDL_GetTicks() - t_last_frame), 0));
            t_last_frame = SDL_GetTicks();

            while (SDL_PollEvent(&sdl_event)) {
                sdl_message_check();
            }

            if ((autoWtoPNG == 1 && (SDL_GetTicks() > t_autoPNG_last + t_autoPNG_delta)) || (autoWtoPNG == 2 && (Ppsum > Ppsum_autoPNG_last + Ppsum_autoPNG_delta))) {
                if (lr_mode == 0) {
                    sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i W(%.6f %.6f %.1f) ct%i cm%i.%i.%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Wr_lo + Wr_up), 0.5 * (Wi_lo + Wi_up), 4.0 / (Wr_up - Wr_lo), ct_type, cm[0], cm_log[0], ct_o[0], (double)Ppsum);
                }

                if (lr_mode == 1) {
                    sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i W(%.6f %.6f %.1f) ct%i cm%i.%i.%i cm.%i.%i.%i cm%i.%i.%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Wr_lo + Wr_up), 0.5 * (Wi_lo + Wi_up), 4.0 / (Wr_up - Wr_lo), ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], (double)Ppsum);
                }

                if (lr_mode == 2) {
                    sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i W(%.6f %.6f %.1f) ct%i cm%i.%i.%i cm.%i.%i.%i cm%i.%i.%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Wr_lo + Wr_up), 0.5 * (Wi_lo + Wi_up), 4.0 / (Wr_up - Wr_lo), ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], (double)Ppsum);
                }

                writeTtoPNG(filename, RWox, RWoy, Ww, Wh);
            }

            if ((autoRTtoPNG == 1 && (SDL_GetTicks() > t_autoPNG_last + t_autoPNG_delta)) || (autoRTtoPNG == 2 && (Ppsum > Ppsum_autoPNG_last + Ppsum_autoPNG_delta))) {
                if (lr_mode == 0) {
                    sprintf(dirname, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i cm%i.%i.%i R%ix%i T%ix%i %g", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type, cm[0], cm_log[0], ct_o[0], Rw, Rh, Tw, Th, (double)Ppsum);
                }

                if (lr_mode == 1) {
                    sprintf(dirname, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i cm%i.%i.%i cm.%i.%i.%i cm%i.%i.%i R%ix%i T%ix%i %g", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], Rw, Rh, Tw, Th, (double)Ppsum);
                }

                if (lr_mode == 2) {
                    sprintf(dirname, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i cm%i.%i.%i cm.%i.%i.%i cm%i.%i.%i R%ix%i T%ix%i %g", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], Rw, Rh, Tw, Th, (double)Ppsum);
                }

                sprintf(commandname, "mkdir \"%s\"", dirname);
                system(commandname);

                for (int png_offset_x = 0; png_offset_x < Rw; png_offset_x += Tw) {
                    for (int png_offset_y = 0; png_offset_y < Rh; png_offset_y += Th) {
                        double Tr_lo = Rr_lo + png_offset_y / Rhdivr;
                        double Tr_up = Tr_lo + Th / Rhdivr;
                        double Ti_lo = Ri_lo + png_offset_x / Rwdivi;
                        double Ti_up = Ti_lo + Tw / Rwdivi;
                        sprintf(filename, "%s/T%ix%i %05i %05i W(%.6f %.6f %.1f).png", dirname, Tw, Th, png_offset_y, png_offset_x, 0.5 * (Tr_lo + Tr_up), 0.5 * (Ti_lo + Ti_up), 4.0 / (Tr_up - Tr_lo));
                        writeTtoPNG(filename, png_offset_x, png_offset_y, Tw, Th);
                    }
                }
            }

            if ((autoRtoPNG && (SDL_GetTicks() > t_autoPNG_last + t_autoPNG_delta)) || (autoRtoPNG == 2 && (Ppsum > Ppsum_autoPNG_last + Ppsum_autoPNG_delta))) {
                if (lr_mode == 0) {
                    sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i cm%i.%i.%i R%ix%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type, cm[0], cm_log[0], ct_o[0], Rw, Rh, (double)Ppsum);
                }

                if (lr_mode == 1) {
                    sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i cm%i.%i.%i cm.%i.%i.%i cm%i.%i.%i R%ix%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], Rw, Rh, (double)Ppsum);
                }

                if (lr_mode == 2) {
                    sprintf(filename, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) ct%i cm%i.%i.%i cm.%i.%i.%i cm%i.%i.%i R%ix%i %g.png", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], Rw, Rh, (double)Ppsum);
                }

                writeRtoPNG(filename);
            }

            if ((autoWtoPNG == 1 || autoRTtoPNG == 1 || autoRtoPNG == 1) && (SDL_GetTicks() > t_autoPNG_last + t_autoPNG_delta)) {
                t_autoPNG_last = SDL_GetTicks();
            }

            if ((autoWtoPNG == 2 || autoRTtoPNG == 2 || autoRtoPNG == 2) && (Ppsum > Ppsum_autoPNG_last + Ppsum_autoPNG_delta)) {
                Ppsum_autoPNG_last = Ppsum;
            }

            #pragma omp flush(td_pause)

            if (td_pause == 0 || recalc_WH_if_paused == 1) {
                if (lr_mode == 0) {
                    Rlrmax[0] = 0;

                    for (int Wy = 0, Ry = RWoy; Wy < Wh; Wy++, Ry++) {
                        for (int Wx = 0, Rx = RWox; Wx < Ww; Wx++, Rx++) {
                            int Ri = Ry * Rw + Rx;
                            unsigned int sum = 0;

                            for (int td_i = 0; td_i < td_nb; td_i += 1) {
                                sum += R[td_i][Ri];
                            }

                            if (cm[0] == 2 || cm[0] == 3) {
                                if (sum >= cm_log[0]) {
                                    sum = cm_log[0] + (unsigned int) log((double)(sum - cm_log[0] + 1));
                                }
                            }

                            if (sum > Rlrmax[0]) {
                                Rlrmax[0] = sum;
                            }

                            W[0][Wy * Ww + Wx] = sum;
                        }
                    }

                    if (Rlrmax[0] >= Hl[0]) {
                        if (H[0] != NULL) {
                            free(H[0]);
                            H[0] = NULL;
                        }

                        Hl[0] = MAX(Rlrmax[0] + 1, 2 * Hl[0]);
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

                    for (int Wy = 0; Wy < Wh; Wy++) {
                        for (int Wx = 0; Wx < Ww; Wx++) {
                            H[0][W[0][Wy * Ww + Wx]]++;
                        }
                    }

                    if (cm[0] == 0 || cm[0] == 2) {
                        cm0n[0] = 0;

                        for (unsigned int i = 0; i <= Rlrmax[0]; i++) {
                            if (H[0][i] > 0) {
                                H[0][i] = cm0n[0]++;
                            }
                        }
                    }

                    if (cm[0] == 1 || cm[0] == 3) {
                        cm1n[0] = Ww * Wh - H[0][0];
                        H[0][0] = 0;

                        for (unsigned int i = 2; i <= Rlrmax[0]; i++) {
                            H[0][i] += H[0][i - 1];
                        }
                    }
                }

                if (lr_mode == 1) {
                    for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
                        Rlrmax[layer_iter] = 0;

                        for (int Wy = 0, Ry = RWoy; Wy < Wh; Wy++, Ry++) {
                            for (int Wx = 0, Rx = RWox; Wx < Ww; Wx++, Rx++) {
                                int Ri = Ry * Rw + Rx;
                                unsigned int sum = 0;

                                for (int td_i = 0; td_i < td_nb; td_i += 1) {
                                    sum += R[td_i][Ri];
                                }

                                if (cm[layer_iter] == 2 || cm[layer_iter] == 3) {
                                    if (sum >= cm_log[layer_iter]) {
                                        sum = cm_log[layer_iter] + (unsigned int) log((double)(sum - cm_log[layer_iter] + 1));
                                    }
                                }

                                if (sum > Rlrmax[layer_iter]) {
                                    Rlrmax[layer_iter] = sum;
                                }

                                W[layer_iter][Wy * Ww + Wx] = sum;
                            }
                        }

                        if (Rlrmax[layer_iter] >= Hl[layer_iter]) {
                            if (H[layer_iter] != NULL) {
                                free(H[layer_iter]);
                                H[layer_iter] = NULL;
                            }

                            Hl[layer_iter] = MAX(Rlrmax[layer_iter] + 1, 2 * Hl[layer_iter]);
                            H[layer_iter] = (unsigned int*)calloc(Hl[layer_iter], sizeof(unsigned int));
                        } else {
                            memset(H[layer_iter], 0, Hl[layer_iter] * sizeof(unsigned int));
                        }

                        for (int Wy = 0; Wy < Wh; Wy++) {
                            for (int Wx = 0; Wx < Ww; Wx++) {
                                H[layer_iter][W[layer_iter][Wy * Ww + Wx]]++;
                            }
                        }

                        if (cm[layer_iter] == 0 || cm[layer_iter] == 2) {
                            cm0n[layer_iter] = 0;

                            for (unsigned int i = 0; i <= Rlrmax[layer_iter]; i++) {
                                if (H[layer_iter][i] > 0) {
                                    H[layer_iter][i] = cm0n[layer_iter]++;
                                }
                            }
                        }

                        if (cm[layer_iter] == 1 || cm[layer_iter] == 3) {
                            cm1n[layer_iter] = Ww * Wh - H[layer_iter][0];
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

                        for (int Wy = 0, Ry = RWoy; Wy < Wh; Wy++, Ry++) {
                            for (int Wx = 0, Rx = RWox; Wx < Ww; Wx++, Rx++) {
                                int Ri = Ry * Rw + Rx;
                                unsigned int sum = 0;

                                for (int td_i = layer_iter; td_i < td_nb; td_i += LR_NB) {
                                    sum += R[td_i][Ri];
                                }

                                if (cm[layer_iter] == 2 || cm[layer_iter] == 3) {
                                    if (sum >= cm_log[layer_iter]) {
                                        sum = cm_log[layer_iter] + (unsigned int) log((double)(sum - cm_log[layer_iter] + 1));
                                    }
                                }

                                if (sum > Rlrmax[layer_iter]) {
                                    Rlrmax[layer_iter] = sum;
                                }

                                W[layer_iter][Wy * Ww + Wx] = sum;
                            }
                        }

                        if (Rlrmax[layer_iter] >= Hl[layer_iter]) {
                            if (H[layer_iter] != NULL) {
                                free(H[layer_iter]);
                                H[layer_iter] = NULL;
                            }

                            Hl[layer_iter] = MAX(Rlrmax[layer_iter] + 1, 2 * Hl[layer_iter]);
                            H[layer_iter] = (unsigned int*)calloc(Hl[layer_iter], sizeof(unsigned int));
                        } else {
                            memset(H[layer_iter], 0, Hl[layer_iter] * sizeof(unsigned int));
                        }

                        for (int Wy = 0; Wy < Wh; Wy++) {
                            for (int Wx = 0; Wx < Ww; Wx++) {
                                H[layer_iter][W[layer_iter][Wy * Ww + Wx]]++;
                            }
                        }

                        if (cm[layer_iter] == 0 || cm[layer_iter] == 2) {
                            cm0n[layer_iter] = 0;

                            for (unsigned int i = 0; i <= Rlrmax[layer_iter]; i++) {
                                if (H[layer_iter][i] > 0) {
                                    H[layer_iter][i] = cm0n[layer_iter]++;
                                }
                            }
                        }

                        if (cm[layer_iter] == 1 || cm[layer_iter] == 3) {
                            cm1n[layer_iter] = Ww * Wh - H[layer_iter][0];
                            H[layer_iter][0] = 0;

                            for (unsigned int i = 2; i <= Rlrmax[layer_iter]; i++) {
                                H[layer_iter][i] += H[layer_iter][i - 1];
                            }
                        }
                    }
                }
            }

            Uint32* pixel_p;

            if (lr_mode == 0) {
                ct_o[0] = ct_cycle(ct_o[0] + ct_v[0]);

                if (cm[0] == 0 || cm[0] == 2) {
                    for (int Wy = 0; Wy < Wh; Wy++) {
                        for (int Wx = 0; Wx < Ww; Wx++) {
                            int ct_i = 0;

                            if (cm0n[0] > 1) {
                                ct_i = ct_e * ((double)H[0][W[0][Wy * Ww + Wx]] / (cm0n[0] - 1));
                            }

                            ct_i = ct_cycle(ct_i + ct_o[0]);
                            pixel_p = (Uint32*)sdl_surface->pixels + Wy * Ww + Wx;
                            *pixel_p = SDL_MapRGBA(sdl_surface->format, CT[ct_i][0], CT[ct_i][1], CT[ct_i][2], SDL_ALPHA_OPAQUE);
                        }
                    }
                }

                if (cm[0] == 1 || cm[0] == 3) {
                    for (int Wy = 0; Wy < Wh; Wy++) {
                        for (int Wx = 0; Wx < Ww; Wx++) {
                            int ct_i = 0;

                            if (cm1n[0] > 0) {
                                ct_i = ct_e * ((double)H[0][W[0][Wy * Ww + Wx]] / cm1n[0]);
                            }

                            ct_i = ct_cycle(ct_i + ct_o[0]);
                            pixel_p = (Uint32*)sdl_surface->pixels + Wy * Ww + Wx;
                            *pixel_p = SDL_MapRGBA(sdl_surface->format, CT[ct_i][0], CT[ct_i][1], CT[ct_i][2], SDL_ALPHA_OPAQUE);
                        }
                    }
                }
            }

            if (lr_mode == 1) {
                for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
                    ct_o[layer_iter] = ct_cycle(ct_o[layer_iter] + ct_v[layer_iter]);
                }

                for (int Wy = 0; Wy < Wh; Wy++) {
                    for (int Wx = 0; Wx < Ww; Wx++) {
                        int ct_i[LR_NB] = {0};

                        for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
                            if (cm[layer_iter] == 0 || cm[layer_iter] == 2) {
                                if (cm0n[layer_iter] > 1) {
                                    ct_i[layer_iter] = ct_e * ((double)H[layer_iter][W[layer_iter][Wy * Ww + Wx]] / (cm0n[layer_iter] - 1));
                                }
                            }

                            if (cm[layer_iter] == 1 || cm[layer_iter] == 3) {
                                if (cm1n[layer_iter] > 0) {
                                    ct_i[layer_iter] = ct_e * ((double)H[layer_iter][W[layer_iter][Wy * Ww + Wx]] / cm1n[layer_iter]);
                                }
                            }

                            ct_i[layer_iter] = ct_cycle(ct_i[layer_iter] + ct_o[layer_iter]);
                        }

                        pixel_p = (Uint32*)sdl_surface->pixels + Wy * Ww + Wx;
                        *pixel_p = SDL_MapRGBA(sdl_surface->format, CT[ct_i[0]][0], CT[ct_i[1]][1], CT[ct_i[2]][2], SDL_ALPHA_OPAQUE);
                    }
                }
            }

            if (lr_mode == 2) {
                for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
                    ct_o[layer_iter] = ct_cycle(ct_o[layer_iter] + ct_v[layer_iter]);
                }

                for (int Wy = 0; Wy < Wh; Wy++) {
                    for (int Wx = 0; Wx < Ww; Wx++) {
                        int ct_i[LR_NB] = {0};

                        for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
                            if (cm[layer_iter] == 0 || cm[layer_iter] == 2) {
                                if (cm0n[layer_iter] > 1) {
                                    ct_i[layer_iter] = ct_e * ((double)H[layer_iter][W[layer_iter][Wy * Ww + Wx]] / (cm0n[layer_iter] - 1));
                                }
                            }

                            if (cm[layer_iter] == 1 || cm[layer_iter] == 3) {
                                if (cm1n[layer_iter] > 0) {
                                    ct_i[layer_iter] = ct_e * ((double)H[layer_iter][W[layer_iter][Wy * Ww + Wx]] / cm1n[layer_iter]);
                                }
                            }

                            ct_i[layer_iter] = ct_cycle(ct_i[layer_iter] + ct_o[layer_iter]);
                        }

                        pixel_p = (Uint32*)sdl_surface->pixels + Wy * Ww + Wx;
                        *pixel_p = SDL_MapRGBA(sdl_surface->format, CT[ct_i[0]][0], CT[ct_i[1]][1], CT[ct_i[2]][2], SDL_ALPHA_OPAQUE);
                    }
                }
            }

            if (titlebarswitch == 0) {
                sprintf(titlebar, "lm%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i bb%i.%i.%i.%i.%i R(%.6f %.6f %.1f) W(%.6f %.6f %.1f) ct%i cm%i.%i.%i cm.%i.%i.%i cm%i.%i.%i th%i %g", lr_mode, bb_type[0], bb_bail[0], bb_pps[0], bb_ppe[0], bb_minn[0], bb_type[1], bb_bail[1], bb_pps[1], bb_ppe[1], bb_minn[1], bb_type[2], bb_bail[2], bb_pps[2], bb_ppe[2], bb_minn[2], 0.5 * (Rr_lo + Rr_up), 0.5 * (Ri_lo + Ri_up), 4.0 / (Rr_up - Rr_lo), 0.5 * (Wr_lo + Wr_up), 0.5 * (Wi_lo + Wi_up), 4.0 / (Wr_up - Wr_lo), ct_type, cm[0], cm_log[0], ct_o[0], cm[1], cm_log[1], ct_o[1], cm[2], cm_log[2], ct_o[2], td_nb, (double)Ppsum);
            } else {
                sprintf(titlebar, "lm%i R%ix%i T%ix%i autopng%i%i%i delta %.0e %.0e th%i Rmax %g Pp %g", lr_mode, Rw, Rh, Tw, Th, autoWtoPNG, autoRTtoPNG, autoRtoPNG, (double)t_autoPNG_delta, (double)Ppsum_autoPNG_delta, td_nb, (double)(Rlrmax[0] + Rlrmax[1] + Rlrmax[2]), (double)Ppsum);
            }

            SDL_SetWindowTitle(sdl_window, titlebar);
            SDL_UpdateTexture(sdl_texture, NULL, sdl_surface->pixels, Ww * sizeof(Uint32));
            SDL_RenderClear(sdl_renderer);
            SDL_RenderCopy(sdl_renderer, sdl_texture, NULL, NULL);
            SDL_RenderPresent(sdl_renderer);
        }

        SDL_DestroyTexture(sdl_texture);
        SDL_FreeSurface(sdl_surface);
        SDL_DestroyRenderer(sdl_renderer);
        SDL_DestroyWindow(sdl_window);
        SDL_Quit();
    }

    printf("thread v end\n");
}

int main(int argc, char* argv[])
{
    printf("main start\n");
    #pragma omp parallel
    {
        td_vc_nb = omp_get_num_threads();
    }

    if (argc >= 2) {
        td_nb = MIN(atoi(argv[1]), TD_MAX);
    } else {
        td_nb = MIN(ceil((double)td_vc_nb / 3) * 3, TD_MAX);
    }

    load_location_bb_color_param(0, 1.0, 0.0, 0.0, 0.0, 0.0, 1000, 1000, 0, 0, 1000, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 1000, 0, 0, 0, 10, 0, 5, 0, 0, 5, 0, 0, 5, 0);

    for (int layer_iter = 0; layer_iter < LR_NB; layer_iter++) {
        W[layer_iter] = (unsigned int*)calloc(Ww * Wh, sizeof(unsigned int));
        H[layer_iter] = NULL;
    }

    #pragma omp parallel sections num_threads(TD_MAX + 1)
    {
        #pragma omp section
        visualisation_thread();
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

    printf("main end\n");
    return (0);
}
