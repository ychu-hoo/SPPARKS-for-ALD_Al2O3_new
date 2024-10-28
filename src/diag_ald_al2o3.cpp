/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "diag_ald_al2o3.h"
#include "app.h"
#include "app_ald_al2o3.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

enum{VACANCY,O,OH,OH2,
AlX3O,AlX3OH,AlX3OH2,AlX2O,AlX2OH,AlX2OH2,
AlXO,AlXOH,AlXOH2,AlO,AlOH,AlOH2,
Al,AlX,AlX2,
OH2Al,OH2AlX,OH2AlX2,OHAl,OHAlX,OHAlX2,
OAl,OAlX,OAlX2,
QCM,LIGANDS,EVENTS,ONE,TWO,THREE};


// enum{VACANCY,O,OH,//2
// OH2, ZnX2O, ZnX2OH, ZnX2OH2, // 6 DEZ adsorption 
// ZnXO, ZnXOH, ZnO, ZnOH, Zn, ZnX, // 13 DEZ surface species
// OH2Zn, OH2ZnX, OHZn, OHZnX, OZn, // 18 water pulse
// Zn_i, ZnX_i, // 20 inert sites
// QCM, OXYGEN, ZINC, ADS_DEZ, HYDROGEN, DEZ, MEZ, LIGANDS, EVENTS,ONE,TWO,THREE
// };       // same as DiagAld


/* ---------------------------------------------------------------------- */

DiagAldAl2o3::DiagAldAl2o3(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"ald/al2o3") != 0)
    error->all(FLERR,"Diag style incompatible with app style");

  nlist = 0;

  int iarg = iarg_child;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"list") == 0) {
      nlist = narg - iarg - 1;
      list = new char*[nlist];
      int j = 0;
      for (int i = iarg+1; i < narg; i++) {
	int n = strlen(arg[i]) + 1;
	list[j] = new char[n];
	strcpy(list[j],arg[i]);
	j++;
      }
      iarg = narg;
    } else error->all(FLERR,"Diag_style aldal2o3 requires app_style aldalo3");                                      
  }

  if (nlist == 0) error->all(FLERR,"Diag_style aldalo3 requires app_style aldal2o3");
  which = new int[nlist];
  index = new int[nlist];
  ivector = new int[nlist];
}

/* ---------------------------------------------------------------------- */

DiagAldAl2o3::~DiagAldAl2o3()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] ivector;
}

/* ---------------------------------------------------------------------- */

void DiagAldAl2o3::init()
{
  appaldal2o3 = (AppAldAl2o3 *) app;
  
  int none = appaldal2o3->none;
  int ntwo = appaldal2o3->ntwo;
  int nthree = appaldal2o3->nthree;
  for (int i = 0; i < nlist; i++) {
      if (strcmp(list[i],"O") == 0) which[i] = O;
      else if (strcmp(list[i],"OH") == 0) which[i] = OH;
      else if (strcmp(list[i],"OH2") == 0) which[i] = OH2;
      else if (strcmp(list[i],"VAC") == 0) which[i] = VACANCY;
      else if (strcmp(list[i],"QCM") == 0) which[i] = QCM;
      else if (strcmp(list[i],"LIGANDS") == 0) which[i] = LIGANDS;
      else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;
      else if (strcmp(list[i],"AlX3O") == 0) which[i] = AlX3O;
      else if (strcmp(list[i],"AlX3OH") == 0) which[i] = AlX3OH;
      else if (strcmp(list[i],"AlX3OH2") == 0) which[i] = AlX3OH2;
      else if (strcmp(list[i],"AlX2O") == 0) which[i] = AlX2O;
      else if (strcmp(list[i],"AlX2OH") == 0) which[i] = AlX2OH;
      else if (strcmp(list[i],"AlX2OH2") == 0) which[i] = AlX2OH2;
      else if (strcmp(list[i],"AlXO") == 0) which[i] = AlXO;
      else if (strcmp(list[i],"AlXOH") == 0) which[i] = AlXOH;
      else if (strcmp(list[i],"AlXOH2") == 0) which[i] = AlXOH2;
      else if (strcmp(list[i],"AlO") == 0) which[i] = AlO;
      else if (strcmp(list[i],"AlOH") == 0) which[i] = AlOH;
      else if (strcmp(list[i],"AlOH2") == 0) which[i] = AlOH2;
      else if (strcmp(list[i],"Al") == 0) which[i] = Al;
      else if (strcmp(list[i],"AlX") == 0) which[i] = AlX;
      else if (strcmp(list[i],"AlX2") == 0) which[i] = AlX2;
      else if (strcmp(list[i],"OH2Al") == 0) which[i] = OH2Al;
      else if (strcmp(list[i],"OH2AlX") == 0) which[i] = OH2AlX;
      else if (strcmp(list[i],"OH2AlX2") == 0) which[i] = OH2AlX2;
      else if (strcmp(list[i],"OHAl") == 0) which[i] = OHAl;
      else if (strcmp(list[i],"OHAlX") == 0) which[i] = OHAlX;
      else if (strcmp(list[i],"OHAlX2") == 0) which[i] = OHAlX2;
      else if (strcmp(list[i],"OAl") == 0) which[i] = OAl;
      else if (strcmp(list[i],"OAlX") == 0) which[i] = OAlX;
      else if (strcmp(list[i],"OAlX2") == 0) which[i] = OAlX2;

    else if (list[i][0] == 's') {
      which[i] = ONE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > none) 
	error->all(FLERR,"Diag_style aldal2o3 requires app_style aldal2o3");
      index[i] = n - 1;
    } else if (list[i][0] == 'd') {
      which[i] = TWO;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ntwo) 
	error->all(FLERR,"Diag_style aldal2o3 requires app_style aldal2o3");
      index[i] = n - 1;
    } else if (list[i][0] == 'v') {
      which[i] = THREE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nthree) 
	error->all(FLERR,"Diag_style aldal2o3 requires app_style aldal2o3");
      index[i] = n - 1;
    } else error->all(FLERR,"Diag_style aldal2o3 requires app_style aldal2o3");
  }

  siteflag = 1; 

  for (int i = 0; i < nlist; i++) ivector[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DiagAldAl2o3::compute()
{
  int sites[800],ivalue;
// here as well we have to consider some modification, generally it does not seem so difficult
  if (siteflag) {
    sites[O] = 0; sites[OH] = 0; sites[VACANCY] = 0; sites[OH2] = 0;
    sites[AlX3O] = 0; sites[AlX3OH] = 0;sites[AlX3OH2] = 0; sites[AlX2O] = 0; sites[AlX2OH] = 0;   
    sites[AlX2OH2] = 0; sites[AlXO] = 0;sites[AlXOH] = 0; sites[AlXOH2] = 0; sites[AlO] = 0; sites[AlOH] = 0;   
    sites[AlOH2] = 0; sites[Al] = 0;sites[AlX] = 0; sites[AlX2] = 0;
    sites[OH2Al] = 0; sites[OH2AlX] = 0;sites[OH2AlX2] = 0; sites[OHAl] = 0; sites[OHAlX] = 0;   
    sites[OHAlX2] = 0; sites[OAl] = 0;sites[OAlX] = 0; sites[OAlX2] = 0;  


    int *element = appaldal2o3->element;
    int nlocal = appaldal2o3->nlocal;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == OH) ivalue = sites[OH];
    else if (which[i] == O) ivalue = sites[O];
    else if (which[i] == VACANCY) ivalue = sites[VACANCY];
    else if (which[i] == OH2) ivalue = sites[OH2];
    else if (which[i] == AlX3O) ivalue = sites[AlX3O];
    else if (which[i] == AlX3OH) ivalue = sites[AlX3OH];
    else if (which[i] == AlX3OH2) ivalue = sites[AlX3OH2];
    else if (which[i] == AlX2O) ivalue = sites[AlX2O];
    else if (which[i] == AlX2OH) ivalue = sites[AlX2OH];
    else if (which[i] == AlX2OH2) ivalue = sites[AlX2OH2];
    else if (which[i] == AlXO) ivalue = sites[AlXO];
    else if (which[i] == AlXOH) ivalue = sites[AlXOH];
    else if (which[i] == AlXOH2) ivalue = sites[AlXOH2];
    else if (which[i] == AlO) ivalue = sites[AlO];
    else if (which[i] == AlOH) ivalue = sites[AlOH];
    else if (which[i] == AlOH2) ivalue = sites[AlOH2];
    else if (which[i] == Al) ivalue = sites[Al];
    else if (which[i] == AlX) ivalue = sites[AlX];
    else if (which[i] == AlX2) ivalue = sites[AlX2];
    else if (which[i] == OH2Al) ivalue = sites[OH2Al];
    else if (which[i] == OH2AlX) ivalue = sites[OH2AlX];
    else if (which[i] == OH2AlX2) ivalue = sites[OH2AlX2];
    else if (which[i] == OHAl) ivalue = sites[OHAl];
    else if (which[i] == OHAlX) ivalue = sites[OHAlX];
    else if (which[i] == OHAlX2) ivalue = sites[OHAlX2];
    else if (which[i] == OAl) ivalue = sites[OAl];
    else if (which[i] == OAlX) ivalue = sites[OAlX];
    else if (which[i] == OAlX2) ivalue = sites[OAlX2];
    else if (which[i] == QCM) ivalue = 18.02*sites[OH2] + 17.01*sites[OH] + 15.99*sites[O] + 88.09*sites[AlX3O] + 89.09*sites[AlX3OH] + 90.10*sites[AlX3OH2] + 73.05*sites[AlX2O] + 74.06*sites[AlX2OH] + 75.07*sites[AlX2OH2] + 58.02*sites[AlXO] + 59.02*sites[AlXOH] + 60.03*sites[AlXOH2] + 42.98*sites[AlO] + 43.99*sites[AlOH] + 44.99*sites[AlOH2] + 26.98*sites[Al] + 42.02*sites[AlX] + 57.05*sites[AlX2] + 44.99*sites[OH2Al] + 60.03*sites[OH2AlX] + 75.07*sites[OH2AlX2] + 43.99*sites[OHAl] + 59.02*sites[OHAlX] + 74.06*sites[OHAlX2] + 42.98*sites[OAl] + 58.02*sites[OAlX] + 73.05*sites[OAlX2];
    else if (which[i] == LIGANDS) ivalue = 3*(sites[AlX3O]+sites[AlX3OH]+sites[AlX3OH2])+2*(sites[AlX2O]+sites[AlX2OH]+sites[AlX2OH2])+(sites[AlXO]+sites[AlXOH]+sites[AlXOH2])+(sites[AlX]+sites[OH2AlX]+sites[OHAlX]+sites[OAlX])+2*(sites[AlX2]+sites[OH2AlX2]+sites[OHAlX2]+sites[OAlX2]);
    else if (which[i] == EVENTS) ivalue = appaldal2o3->nevents;
    else if (which[i] == ONE) ivalue = appaldal2o3->scount[index[i]];
    else if (which[i] == TWO) ivalue = appaldal2o3->dcount[index[i]];
    else if (which[i] == THREE) ivalue = appaldal2o3->vcount[index[i]];
    
    MPI_Allreduce(&ivalue,&ivector[i],1,MPI_INT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagAldAl2o3::stats(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str,"%6d ",ivector[i]);
    str += strlen(str);
  }
}

/* ---------------------------------------------------------------------- */

void DiagAldAl2o3::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str,"%6s ",list[i]);
    str += strlen(str);
  }
}
