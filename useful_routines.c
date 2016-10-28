/////////////////////
// HEADERS //
////////////////////
#include<stdio.h>
#include<stdlib.h>  
#include<string.h>  
#include<malloc.h>
#include<math.h>
#include<time.h>

//CABECERAS DE GSL                                                                  
#include<gsl/gsl_errno.h>
#include<gsl/gsl_sort.h>
#include <gsl/gsl_permutation.h>
#include<gsl/gsl_sort_uint.h>
#include<gsl/gsl_fit.h>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>

////////////////////
// MACROS //
///////////////////

#define X 0
#define Y 1
#define Z 2

////////////////////////////////////
// DATA STRUCTURES //
///////////////////////////////////

//******************
// For Gadget data
typedef struct 
{
  int type;
#ifdef LONGIDS
  unsigned long long id;
#else
  unsigned int id;
#endif
  float pos[3];
  float vel[3];
  float mass;
  //float Z;
  float pot;
  float acce;
  float timestep;
}particulas;

typedef struct
{
  float U;
  float rho;
  float Ne;
  float Nh;
  float h;
  float sfr;
  //float stellar_age;
  float ecr;
} gas_properties;

typedef struct
{
  int Npart[6];
  double mass[6];
  double time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned int npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHighWord[6];
  int  flag_entropy_instead_u;
  char     fill[256 - 6*4 - 6*8 - 2*8 - 2*4 - 6*4 - 2*4 - 4*8 - 2*4 - 6*4 - 1*4];  /* fills to 256 Bytes */
} io_header;

//*********************

typedef struct
{
  int id;
  double cm[3];
  double vcm[3];
  double lcm[3];
  double M;
} CM;

typedef struct
{
  double xini;
  double xfin;
  double x;
}Bin;

////////////////////////////////////
// GLOBAL VARIABLES //
///////////////////////////////////
particulas *particles;
gas_properties *gaspro;
io_header Header;
CM cmdummy;
Bin *bins;

int N_part_total,N_part[6],N_min,N_max;
int n0,n1,n2,n3,n4,n5;
int returnRead;

double Mass_total,Mass_tot[6];

double G = 43007.1;

float *U;

//////////////////////
// ROUTINES //
/////////////////////

#ifdef READ_GADGET1
#include<../useful_routines/input_output_data.c>
#endif

#ifdef TRANSLATIONS_ROTATIOS
#include<../useful_routines/translatios_rotations.c>
#endif

#ifdef BINS
#include<../useful_routines/makes_bins.c>
#endif

#ifdef MEDLEY
#include<../useful_routines/medley.c>
#endif

