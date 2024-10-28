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

#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "app_ald_al2o3.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{VACANCY,O,OH,OH2,
AlX3O,AlX3OH,AlX3OH2,AlX2O,AlX2OH,AlX2OH2,
AlXO,AlXOH,AlXOH2,AlO,AlOH,AlOH2,
Al,AlX,AlX2,
OH2Al,OH2AlX,OH2AlX2,OHAl,OHAlX,OHAlX2,
OAl,OAlX,OAlX2};

// enum{VACANCY,O,OH,// 3
// OH2, ZnX2O, ZnX2OH, ZnX2OH2, // 7 DEZ adsorption 
// ZnXO, ZnXOH, ZnO, ZnOH, Zn, ZnX, // 13 DEZ surface species
// OH2Zn, OH2ZnX, OHZn, OHZnX, OZn}; // 18 water pulse



#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppAldAl2o3::AppAldAl2o3(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = 2;
  ndouble = 0;
  delpropensity = 1;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;
  allow_masking = 0;
  

  create_arrays();

  if (narg != 1) error->all(FLERR,"Illegal app_style command");

  cycle = 0;
  pressureOn = 1;
  hello = 1;
  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;

  // reaction lists

  none = ntwo = nthree = 0;
  srate = drate = vrate = NULL;
  spropensity = dpropensity = vpropensity = NULL;
  sinput = soutput = NULL;
  dinput = doutput = NULL;
  vinput = voutput = NULL;
  comneigh = NULL;
  scount = dcount = vcount = NULL;
  sA = dA = vA = NULL;
  scoord = dcoord = vcoord = NULL;
  sexpon = dexpon = vexpon = NULL;
  spresson = dpresson = vpresson = NULL;
}

/* ---------------------------------------------------------------------- */

AppAldAl2o3::~AppAldAl2o3()
{
  delete [] esites;
  delete [] echeck;
  memory->sfree(events);
  memory->sfree(firstevent);
  memory->sfree(srate);
  memory->sfree(drate);
  memory->sfree(vrate);
  memory->sfree(spropensity);
  memory->sfree(dpropensity);
  memory->sfree(vpropensity);
  memory->sfree(sinput);
  memory->sfree(soutput);
  memory->sfree(dinput);
  memory->sfree(doutput);
  memory->sfree(vinput);
  memory->sfree(voutput);
  memory->sfree(comneigh);
  memory->sfree(scount);
  memory->sfree(dcount);
  memory->sfree(vcount);
  memory->sfree(sA);
  memory->sfree(dA);
  memory->sfree(vA);
  memory->sfree(scoord);
  memory->sfree(dcoord);
  memory->sfree(vcoord);
  memory->sfree(sexpon);
  memory->sfree(dexpon);
  memory->sfree(vexpon);
  memory->sfree(spresson);
  memory->sfree(dpresson);
  memory->sfree(vpresson);
}

/* ---------------------------------------------------------------------- */

void AppAldAl2o3::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"event") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal event command E1");
    int rstyle = atoi(arg[0]);
    grow_reactions(rstyle);

    if (rstyle == 1) {
      if (narg != 9) error->all(FLERR,"Illegal event arg command");
//type I
      if (strcmp(arg[1],"O") == 0) sinput[none] = O;
      else if (strcmp(arg[1],"OH") == 0) sinput[none] = OH; 
      else if (strcmp(arg[1],"OH2") == 0) sinput[none] = OH2;
      else if (strcmp(arg[1],"AlX3O") == 0) sinput[none] = AlX3O;
      else if (strcmp(arg[1],"AlX3OH") == 0) sinput[none] = AlX3OH;
      else if (strcmp(arg[1],"AlX3OH2") == 0) sinput[none] = AlX3OH2;
      else if (strcmp(arg[1],"AlX2O") == 0) sinput[none] = AlX2O;
      else if (strcmp(arg[1],"AlX2OH") == 0) sinput[none] = AlX2OH;
      else if (strcmp(arg[1],"AlX2OH2") == 0) sinput[none] = AlX2OH2;
      else if (strcmp(arg[1],"AlXO") == 0) sinput[none] = AlXO;
      else if (strcmp(arg[1],"AlXOH") == 0) sinput[none] = AlXOH;
      else if (strcmp(arg[1],"AlXOH2") == 0) sinput[none] = AlXOH2;
      else if (strcmp(arg[1],"AlO") == 0) sinput[none] = AlO;
      else if (strcmp(arg[1],"AlOH") == 0) sinput[none] = AlOH;
      else if (strcmp(arg[1],"AlOH2") == 0) sinput[none] = AlOH2;
      else if (strcmp(arg[1],"Al") == 0) sinput[none] = Al;
      else if (strcmp(arg[1],"AlX") == 0) sinput[none] = AlX;
      else if (strcmp(arg[1],"AlX2") == 0) sinput[none] = AlX2;
      else if (strcmp(arg[1],"OH2Al") == 0) sinput[none] = OH2Al;
      else if (strcmp(arg[1],"OH2AlX") == 0) sinput[none] = OH2AlX;
      else if (strcmp(arg[1],"OH2AlX2") == 0) sinput[none] = OH2AlX2;
      else if (strcmp(arg[1],"OHAl") == 0) sinput[none] = OHAl;
      else if (strcmp(arg[1],"OHAlX") == 0) sinput[none] = OHAlX;
      else if (strcmp(arg[1],"OHAlX2") == 0) sinput[none] = OHAlX2;
      else if (strcmp(arg[1],"OAl") == 0) sinput[none] = OAl;
      else if (strcmp(arg[1],"OAlX") == 0) sinput[none] = OAlX;
      else if (strcmp(arg[1],"OAlX2") == 0) sinput[none] = OAlX2;
  
      else error->all(FLERR,"Illegal event arg1 command");

      if (strcmp(arg[2],"O") == 0) soutput[none] = O;
      else if (strcmp(arg[2],"VAC") == 0) soutput[none] = VACANCY; 
      else if (strcmp(arg[2],"OH") == 0) soutput[none] = OH; 
      else if (strcmp(arg[2],"OH2") == 0) soutput[none] = OH2;
      else if (strcmp(arg[2],"AlX3O") == 0) soutput[none] = AlX3O;
      else if (strcmp(arg[2],"AlX3OH") == 0) soutput[none] = AlX3OH;
      else if (strcmp(arg[2],"AlX3OH2") == 0) soutput[none] = AlX3OH2;
      else if (strcmp(arg[2],"AlX2O") == 0) soutput[none] = AlX2O;
      else if (strcmp(arg[2],"AlX2OH") == 0) soutput[none] = AlX2OH;
      else if (strcmp(arg[2],"AlX2OH2") == 0) soutput[none] = AlX2OH2;
      else if (strcmp(arg[2],"AlXO") == 0) soutput[none] = AlXO;
      else if (strcmp(arg[2],"AlXOH") == 0) soutput[none] = AlXOH;
      else if (strcmp(arg[2],"AlXOH2") == 0) soutput[none] = AlXOH2;
      else if (strcmp(arg[2],"AlO") == 0) soutput[none] = AlO;
      else if (strcmp(arg[2],"AlOH") == 0) soutput[none] = AlOH;
      else if (strcmp(arg[2],"AlOH2") == 0) soutput[none] = AlOH2;
      else if (strcmp(arg[2],"Al") == 0) soutput[none] = Al;
      else if (strcmp(arg[2],"AlX") == 0) soutput[none] = AlX;
      else if (strcmp(arg[2],"AlX2") == 0) soutput[none] = AlX2;
      else if (strcmp(arg[2],"OH2Al") == 0) soutput[none] = OH2Al;
      else if (strcmp(arg[2],"OH2AlX") == 0) soutput[none] = OH2AlX;
      else if (strcmp(arg[2],"OH2AlX2") == 0) soutput[none] = OH2AlX2;
      else if (strcmp(arg[2],"OHAl") == 0) soutput[none] = OHAl;
      else if (strcmp(arg[2],"OHAlX") == 0) soutput[none] = OHAlX;
      else if (strcmp(arg[2],"OHAlX2") == 0) soutput[none] = OHAlX2;
      else if (strcmp(arg[2],"OAl") == 0) soutput[none] = OAl;
      else if (strcmp(arg[2],"OAlX") == 0) soutput[none] = OAlX;
      else if (strcmp(arg[2],"OAlX2") == 0) soutput[none] = OAlX2;

      else error->all(FLERR,"Illegal event command E2");
      
      sA[none] = atof(arg[3]);
      if (sA[none] == 0.0) error->warning(FLERR,"Illegal coef during reading command");
      sexpon[none] = atoi(arg[4]);
      srate[none] = atof(arg[5]);
      scoord[none] = atoi(arg[6]);
      spresson[none] = atoi(arg[7]);

      none++;

//type II 
    } else if (rstyle == 2) {
      if (narg != 11) error->all(FLERR,"Illegal event command E3");

      if (strcmp(arg[1],"O") == 0) dinput[ntwo][0] = O;
      else if (strcmp(arg[1],"VAC") == 0) dinput[ntwo][0] = VACANCY; 
      else if (strcmp(arg[1],"OH") == 0) dinput[ntwo][0] = OH; 
      else if (strcmp(arg[1],"OH2") == 0) dinput[ntwo][0] = OH2;
      else if (strcmp(arg[1],"AlX3O") == 0) dinput[ntwo][0] = AlX3O;
      else if (strcmp(arg[1],"AlX3OH") == 0) dinput[ntwo][0] = AlX3OH;
      else if (strcmp(arg[1],"AlX3OH2") == 0) dinput[ntwo][0] = AlX3OH2;
      else if (strcmp(arg[1],"AlX2O") == 0) dinput[ntwo][0] = AlX2O;
      else if (strcmp(arg[1],"AlX2OH") == 0) dinput[ntwo][0] = AlX2OH;
      else if (strcmp(arg[1],"AlX2OH2") == 0) dinput[ntwo][0] = AlX2OH2;
      else if (strcmp(arg[1],"AlXO") == 0) dinput[ntwo][0] = AlXO;
      else if (strcmp(arg[1],"AlXOH") == 0) dinput[ntwo][0] = AlXOH;
      else if (strcmp(arg[1],"AlXOH2") == 0) dinput[ntwo][0] = AlXOH2;
      else if (strcmp(arg[1],"AlO") == 0) dinput[ntwo][0] = AlO;
      else if (strcmp(arg[1],"AlOH") == 0) dinput[ntwo][0] = AlOH;
      else if (strcmp(arg[1],"AlOH2") == 0) dinput[ntwo][0] = AlOH2;
      else if (strcmp(arg[1],"Al") == 0) dinput[ntwo][0] = Al;
      else if (strcmp(arg[1],"AlX") == 0) dinput[ntwo][0] = AlX;
      else if (strcmp(arg[1],"AlX2") == 0) dinput[ntwo][0] = AlX2;
      else if (strcmp(arg[1],"OH2Al") == 0) dinput[ntwo][0] = OH2Al;
      else if (strcmp(arg[1],"OH2AlX") == 0) dinput[ntwo][0] = OH2AlX;
      else if (strcmp(arg[1],"OH2AlX2") == 0) dinput[ntwo][0] = OH2AlX2;
      else if (strcmp(arg[1],"OHAl") == 0) dinput[ntwo][0] = OHAl;
      else if (strcmp(arg[1],"OHAlX") == 0) dinput[ntwo][0] = OHAlX;
      else if (strcmp(arg[1],"OHAlX2") == 0) dinput[ntwo][0] = OHAlX2;
      else if (strcmp(arg[1],"OAl") == 0) dinput[ntwo][0] = OAl;
      else if (strcmp(arg[1],"OAlX") == 0) dinput[ntwo][0] = OAlX;
      else if (strcmp(arg[1],"OAlX2") == 0) dinput[ntwo][0] = OAlX2;

      else error->all(FLERR,"Illegal event command E4");

      if (strcmp(arg[2],"O") == 0) doutput[ntwo][0] = O;
      else if (strcmp(arg[2],"VAC") == 0) doutput[ntwo][0] = VACANCY; 
      else if (strcmp(arg[2],"OH") == 0) doutput[ntwo][0] = OH; 
      else if (strcmp(arg[2],"OH2") == 0) doutput[ntwo][0] = OH2;
      else if (strcmp(arg[2],"AlX3O") == 0) doutput[ntwo][0] = AlX3O;
      else if (strcmp(arg[2],"AlX3OH") == 0) doutput[ntwo][0] = AlX3OH;
      else if (strcmp(arg[2],"AlX3OH2") == 0) doutput[ntwo][0] = AlX3OH2;
      else if (strcmp(arg[2],"AlX2O") == 0) doutput[ntwo][0] = AlX2O;
      else if (strcmp(arg[2],"AlX2OH") == 0) doutput[ntwo][0] = AlX2OH;
      else if (strcmp(arg[2],"AlX2OH2") == 0) doutput[ntwo][0] = AlX2OH2;
      else if (strcmp(arg[2],"AlXO") == 0) doutput[ntwo][0] = AlXO;
      else if (strcmp(arg[2],"AlXOH") == 0) doutput[ntwo][0] = AlXOH;
      else if (strcmp(arg[2],"AlXOH2") == 0) doutput[ntwo][0] = AlXOH2;
      else if (strcmp(arg[2],"AlO") == 0) doutput[ntwo][0] = AlO;
      else if (strcmp(arg[2],"AlOH") == 0) doutput[ntwo][0] = AlOH;
      else if (strcmp(arg[2],"AlOH2") == 0) doutput[ntwo][0] = AlOH2;
      else if (strcmp(arg[2],"Al") == 0) doutput[ntwo][0] = Al;
      else if (strcmp(arg[2],"AlX") == 0) doutput[ntwo][0] = AlX;
      else if (strcmp(arg[2],"AlX2") == 0) doutput[ntwo][0] = AlX2;
      else if (strcmp(arg[2],"OH2Al") == 0) doutput[ntwo][0] = OH2Al;
      else if (strcmp(arg[2],"OH2AlX") == 0) doutput[ntwo][0] = OH2AlX;
      else if (strcmp(arg[2],"OH2AlX2") == 0) doutput[ntwo][0] = OH2AlX2;
      else if (strcmp(arg[2],"OHAl") == 0) doutput[ntwo][0] = OHAl;
      else if (strcmp(arg[2],"OHAlX") == 0) doutput[ntwo][0] = OHAlX;
      else if (strcmp(arg[2],"OHAlX2") == 0) doutput[ntwo][0] = OHAlX2;
      else if (strcmp(arg[2],"OAl") == 0) doutput[ntwo][0] = OAl;
      else if (strcmp(arg[2],"OAlX") == 0) doutput[ntwo][0] = OAlX;
      else if (strcmp(arg[2],"OAlX2") == 0) doutput[ntwo][0] = OAlX2;
  
      else error->all(FLERR,"Illegal event command2");

      if (strcmp(arg[3],"O") == 0) dinput[ntwo][1] = O;
      else if (strcmp(arg[3],"VAC") == 0) dinput[ntwo][1] = VACANCY; 
      else if (strcmp(arg[3],"OH") == 0) dinput[ntwo][1] = OH; 
      else if (strcmp(arg[3],"OH2") == 0) dinput[ntwo][1] = OH2;
      else if (strcmp(arg[3],"AlX3O") == 0) dinput[ntwo][1] = AlX3O;
      else if (strcmp(arg[3],"AlX3OH") == 0) dinput[ntwo][1] = AlX3OH;
      else if (strcmp(arg[3],"AlX3OH2") == 0) dinput[ntwo][1] = AlX3OH2;
      else if (strcmp(arg[3],"AlX2O") == 0) dinput[ntwo][1] = AlX2O;
      else if (strcmp(arg[3],"AlX2OH") == 0) dinput[ntwo][1] = AlX2OH;
      else if (strcmp(arg[3],"AlX2OH2") == 0) dinput[ntwo][1] = AlX2OH2;
      else if (strcmp(arg[3],"AlXO") == 0) dinput[ntwo][1] = AlXO;
      else if (strcmp(arg[3],"AlXOH") == 0) dinput[ntwo][1] = AlXOH;
      else if (strcmp(arg[3],"AlXOH2") == 0) dinput[ntwo][1] = AlXOH2;
      else if (strcmp(arg[3],"AlO") == 0) dinput[ntwo][1] = AlO;
      else if (strcmp(arg[3],"AlOH") == 0) dinput[ntwo][1] = AlOH;
      else if (strcmp(arg[3],"AlOH2") == 0) dinput[ntwo][1] = AlOH2;
      else if (strcmp(arg[3],"Al") == 0) dinput[ntwo][1] = Al;
      else if (strcmp(arg[3],"AlX") == 0) dinput[ntwo][1] = AlX;
      else if (strcmp(arg[3],"AlX2") == 0) dinput[ntwo][1] = AlX2;
      else if (strcmp(arg[3],"OH2Al") == 0) dinput[ntwo][1] = OH2Al;
      else if (strcmp(arg[3],"OH2AlX") == 0) dinput[ntwo][1] = OH2AlX;
      else if (strcmp(arg[3],"OH2AlX2") == 0) dinput[ntwo][1] = OH2AlX2;
      else if (strcmp(arg[3],"OHAl") == 0) dinput[ntwo][1] = OHAl;
      else if (strcmp(arg[3],"OHAlX") == 0) dinput[ntwo][1] = OHAlX;
      else if (strcmp(arg[3],"OHAlX2") == 0) dinput[ntwo][1] = OHAlX2;
      else if (strcmp(arg[3],"OAl") == 0) dinput[ntwo][1] = OAl;
      else if (strcmp(arg[3],"OAlX") == 0) dinput[ntwo][1] = OAlX;
      else if (strcmp(arg[3],"OAlX2") == 0) dinput[ntwo][1] = OAlX2;
  
      else error->all(FLERR,"Illegal event command2");

      if (strcmp(arg[4],"O") == 0) doutput[ntwo][1] = O;
      else if (strcmp(arg[4],"VAC") == 0) doutput[ntwo][1] = VACANCY; 
      else if (strcmp(arg[4],"OH") == 0) doutput[ntwo][1] = OH; 
      else if (strcmp(arg[4],"OH2") == 0) doutput[ntwo][1] = OH2;
      else if (strcmp(arg[4],"AlX3O") == 0) doutput[ntwo][1] = AlX3O;
      else if (strcmp(arg[4],"AlX3OH") == 0) doutput[ntwo][1] = AlX3OH;
      else if (strcmp(arg[4],"AlX3OH2") == 0) doutput[ntwo][1] = AlX3OH2;
      else if (strcmp(arg[4],"AlX2O") == 0) doutput[ntwo][1] = AlX2O;
      else if (strcmp(arg[4],"AlX2OH") == 0) doutput[ntwo][1] = AlX2OH;
      else if (strcmp(arg[4],"AlX2OH2") == 0) doutput[ntwo][1] = AlX2OH2;
      else if (strcmp(arg[4],"AlXO") == 0) doutput[ntwo][1] = AlXO;
      else if (strcmp(arg[4],"AlXOH") == 0) doutput[ntwo][1] = AlXOH;
      else if (strcmp(arg[4],"AlXOH2") == 0) doutput[ntwo][1] = AlXOH2;
      else if (strcmp(arg[4],"AlO") == 0) doutput[ntwo][1] = AlO;
      else if (strcmp(arg[4],"AlOH") == 0) doutput[ntwo][1] = AlOH;
      else if (strcmp(arg[4],"AlOH2") == 0) doutput[ntwo][1] = AlOH2;
      else if (strcmp(arg[4],"Al") == 0) doutput[ntwo][1] = Al;
      else if (strcmp(arg[4],"AlX") == 0) doutput[ntwo][1] = AlX;
      else if (strcmp(arg[4],"AlX2") == 0) doutput[ntwo][1] = AlX2;
      else if (strcmp(arg[4],"OH2Al") == 0) doutput[ntwo][1] = OH2Al;
      else if (strcmp(arg[4],"OH2AlX") == 0) doutput[ntwo][1] = OH2AlX;
      else if (strcmp(arg[4],"OH2AlX2") == 0) doutput[ntwo][1] = OH2AlX2;
      else if (strcmp(arg[4],"OHAl") == 0) doutput[ntwo][1] = OHAl;
      else if (strcmp(arg[4],"OHAlX") == 0) doutput[ntwo][1] = OHAlX;
      else if (strcmp(arg[4],"OHAlX2") == 0) doutput[ntwo][1] = OHAlX2;
      else if (strcmp(arg[4],"OAl") == 0) doutput[ntwo][1] = OAl;
      else if (strcmp(arg[4],"OAlX") == 0) doutput[ntwo][1] = OAlX;
      else if (strcmp(arg[4],"OAlX2") == 0) doutput[ntwo][1] = OAlX2;

      else error->all(FLERR,"Illegal event command2");

      dA[ntwo] = atof(arg[5]);
      dexpon[ntwo] = atoi(arg[6]);
      if (dexpon[ntwo] != 0.0) error->warning(FLERR,"Illegal expon command2");
      drate[ntwo] = atof(arg[7]);
      dcoord[ntwo] = atoi(arg[8]);
      dpresson[ntwo] = atoi(arg[9]);
      ntwo++;
// type III
    }else if (rstyle == 3) {
      if (narg != 11) error->all(FLERR,"Illegal event command31");

      if (strcmp(arg[1],"O") == 0) vinput[nthree][0] = O;
      else if (strcmp(arg[1],"VAC") == 0) vinput[nthree][0] = VACANCY; 
      else if (strcmp(arg[1],"OH") == 0) vinput[nthree][0] = OH; 
      else if (strcmp(arg[1],"OH2") == 0) vinput[nthree][0] = OH2;
      else if (strcmp(arg[1],"AlX3O") == 0) vinput[nthree][0] = AlX3O;
      else if (strcmp(arg[1],"AlX3OH") == 0) vinput[nthree][0] = AlX3OH;
      else if (strcmp(arg[1],"AlX3OH2") == 0) vinput[nthree][0] = AlX3OH2;
      else if (strcmp(arg[1],"AlX2O") == 0) vinput[nthree][0] = AlX2O;
      else if (strcmp(arg[1],"AlX2OH") == 0) vinput[nthree][0] = AlX2OH;
      else if (strcmp(arg[1],"AlX2OH2") == 0) vinput[nthree][0] = AlX2OH2;
      else if (strcmp(arg[1],"AlXO") == 0) vinput[nthree][0] = AlXO;
      else if (strcmp(arg[1],"AlXOH") == 0) vinput[nthree][0] = AlXOH;
      else if (strcmp(arg[1],"AlXOH2") == 0) vinput[nthree][0] = AlXOH2;
      else if (strcmp(arg[1],"AlO") == 0) vinput[nthree][0] = AlO;
      else if (strcmp(arg[1],"AlOH") == 0) vinput[nthree][0] = AlOH;
      else if (strcmp(arg[1],"AlOH2") == 0) vinput[nthree][0] = AlOH2;
      else if (strcmp(arg[1],"Al") == 0) vinput[nthree][0] = Al;
      else if (strcmp(arg[1],"AlX") == 0) vinput[nthree][0] = AlX;
      else if (strcmp(arg[1],"AlX2") == 0) vinput[nthree][0] = AlX2;
      else if (strcmp(arg[1],"OH2Al") == 0) vinput[nthree][0] = OH2Al;
      else if (strcmp(arg[1],"OH2AlX") == 0) vinput[nthree][0] = OH2AlX;
      else if (strcmp(arg[1],"OH2AlX2") == 0) vinput[nthree][0] = OH2AlX2;
      else if (strcmp(arg[1],"OHAl") == 0) vinput[nthree][0] = OHAl;
      else if (strcmp(arg[1],"OHAlX") == 0) vinput[nthree][0] = OHAlX;
      else if (strcmp(arg[1],"OHAlX2") == 0) vinput[nthree][0] = OHAlX2;
      else if (strcmp(arg[1],"OAl") == 0) vinput[nthree][0] = OAl;
      else if (strcmp(arg[1],"OAlX") == 0) vinput[nthree][0] = OAlX;
      else if (strcmp(arg[1],"OAlX2") == 0) vinput[nthree][0] = OAlX2;

      else error->all(FLERR,"Illegal event command32");

      if (strcmp(arg[2],"O") == 0) voutput[nthree][0] = O;
      else if (strcmp(arg[2],"VAC") == 0) voutput[nthree][0] = VACANCY; 
      else if (strcmp(arg[2],"OH") == 0) voutput[nthree][0] = OH; 
      else if (strcmp(arg[2],"OH2") == 0) voutput[nthree][0] = OH2;
      else if (strcmp(arg[2],"AlX3O") == 0) voutput[nthree][0] = AlX3O;
      else if (strcmp(arg[2],"AlX3OH") == 0) voutput[nthree][0] = AlX3OH;
      else if (strcmp(arg[2],"AlX3OH2") == 0) voutput[nthree][0] = AlX3OH2;
      else if (strcmp(arg[2],"AlX2O") == 0) voutput[nthree][0] = AlX2O;
      else if (strcmp(arg[2],"AlX2OH") == 0) voutput[nthree][0] = AlX2OH;
      else if (strcmp(arg[2],"AlX2OH2") == 0) voutput[nthree][0] = AlX2OH2;
      else if (strcmp(arg[2],"AlXO") == 0) voutput[nthree][0] = AlXO;
      else if (strcmp(arg[2],"AlXOH") == 0) voutput[nthree][0] = AlXOH;
      else if (strcmp(arg[2],"AlXOH2") == 0) voutput[nthree][0] = AlXOH2;
      else if (strcmp(arg[2],"AlO") == 0) voutput[nthree][0] = AlO;
      else if (strcmp(arg[2],"AlOH") == 0) voutput[nthree][0] = AlOH;
      else if (strcmp(arg[2],"AlOH2") == 0) voutput[nthree][0] = AlOH2;
      else if (strcmp(arg[2],"Al") == 0) voutput[nthree][0] = Al;
      else if (strcmp(arg[2],"AlX") == 0) voutput[nthree][0] = AlX;
      else if (strcmp(arg[2],"AlX2") == 0) voutput[nthree][0] = AlX2;
      else if (strcmp(arg[2],"OH2Al") == 0) voutput[nthree][0] = OH2Al;
      else if (strcmp(arg[2],"OH2AlX") == 0) voutput[nthree][0] = OH2AlX;
      else if (strcmp(arg[2],"OH2AlX2") == 0) voutput[nthree][0] = OH2AlX2;
      else if (strcmp(arg[2],"OHAl") == 0) voutput[nthree][0] = OHAl;
      else if (strcmp(arg[2],"OHAlX") == 0) voutput[nthree][0] = OHAlX;
      else if (strcmp(arg[2],"OHAlX2") == 0) voutput[nthree][0] = OHAlX2;
      else if (strcmp(arg[2],"OAl") == 0) voutput[nthree][0] = OAl;
      else if (strcmp(arg[2],"OAlX") == 0) voutput[nthree][0] = OAlX;
      else if (strcmp(arg[2],"OAlX2") == 0) voutput[nthree][0] = OAlX2;

      else error->all(FLERR,"Illegal event command33");

      if (strcmp(arg[3],"O") == 0) vinput[nthree][1] = O;
      else if (strcmp(arg[3],"VAC") == 0) vinput[nthree][1] = VACANCY; 
      else if (strcmp(arg[3],"OH") == 0) vinput[nthree][1] = OH; 
      else if (strcmp(arg[3],"OH2") == 0) vinput[nthree][1] = OH2;
      else if (strcmp(arg[3],"AlX3O") == 0) vinput[nthree][1] = AlX3O;
      else if (strcmp(arg[3],"AlX3OH") == 0) vinput[nthree][1] = AlX3OH;
      else if (strcmp(arg[3],"AlX3OH2") == 0) vinput[nthree][1] = AlX3OH2;
      else if (strcmp(arg[3],"AlX2O") == 0) vinput[nthree][1] = AlX2O;
      else if (strcmp(arg[3],"AlX2OH") == 0) vinput[nthree][1] = AlX2OH;
      else if (strcmp(arg[3],"AlX2OH2") == 0) vinput[nthree][1] = AlX2OH2;
      else if (strcmp(arg[3],"AlXO") == 0) vinput[nthree][1] = AlXO;
      else if (strcmp(arg[3],"AlXOH") == 0) vinput[nthree][1] = AlXOH;
      else if (strcmp(arg[3],"AlXOH2") == 0) vinput[nthree][1] = AlXOH2;
      else if (strcmp(arg[3],"AlO") == 0) vinput[nthree][1] = AlO;
      else if (strcmp(arg[3],"AlOH") == 0) vinput[nthree][1] = AlOH;
      else if (strcmp(arg[3],"AlOH2") == 0) vinput[nthree][1] = AlOH2;
      else if (strcmp(arg[3],"Al") == 0) vinput[nthree][1] = Al;
      else if (strcmp(arg[3],"AlX") == 0) vinput[nthree][1] = AlX;
      else if (strcmp(arg[3],"AlX2") == 0) vinput[nthree][1] = AlX2;
      else if (strcmp(arg[3],"OH2Al") == 0) vinput[nthree][1] = OH2Al;
      else if (strcmp(arg[3],"OH2AlX") == 0) vinput[nthree][1] = OH2AlX;
      else if (strcmp(arg[3],"OH2AlX2") == 0) vinput[nthree][1] = OH2AlX2;
      else if (strcmp(arg[3],"OHAl") == 0) vinput[nthree][1] = OHAl;
      else if (strcmp(arg[3],"OHAlX") == 0) vinput[nthree][1] = OHAlX;
      else if (strcmp(arg[3],"OHAlX2") == 0) vinput[nthree][1] = OHAlX2;
      else if (strcmp(arg[3],"OAl") == 0) vinput[nthree][1] = OAl;
      else if (strcmp(arg[3],"OAlX") == 0) vinput[nthree][1] = OAlX;
      else if (strcmp(arg[3],"OAlX2") == 0) vinput[nthree][1] = OAlX2;

      else error->all(FLERR,"Illegal event command34");

      if (strcmp(arg[4],"O") == 0) voutput[nthree][1] = O;
      else if (strcmp(arg[4],"VAC") == 0) voutput[nthree][1] = VACANCY; 
      else if (strcmp(arg[4],"OH") == 0) voutput[nthree][1] = OH; 
      else if (strcmp(arg[4],"OH2") == 0) voutput[nthree][1] = OH2;
      else if (strcmp(arg[4],"AlX3O") == 0) voutput[nthree][1] = AlX3O;
      else if (strcmp(arg[4],"AlX3OH") == 0) voutput[nthree][1] = AlX3OH;
      else if (strcmp(arg[4],"AlX3OH2") == 0) voutput[nthree][1] = AlX3OH2;
      else if (strcmp(arg[4],"AlX2O") == 0) voutput[nthree][1] = AlX2O;
      else if (strcmp(arg[4],"AlX2OH") == 0) voutput[nthree][1] = AlX2OH;
      else if (strcmp(arg[4],"AlX2OH2") == 0) voutput[nthree][1] = AlX2OH2;
      else if (strcmp(arg[4],"AlXO") == 0) voutput[nthree][1] = AlXO;
      else if (strcmp(arg[4],"AlXOH") == 0) voutput[nthree][1] = AlXOH;
      else if (strcmp(arg[4],"AlXOH2") == 0) voutput[nthree][1] = AlXOH2;
      else if (strcmp(arg[4],"AlO") == 0) voutput[nthree][1] = AlO;
      else if (strcmp(arg[4],"AlOH") == 0) voutput[nthree][1] = AlOH;
      else if (strcmp(arg[4],"AlOH2") == 0) voutput[nthree][1] = AlOH2;
      else if (strcmp(arg[4],"Al") == 0) voutput[nthree][1] = Al;
      else if (strcmp(arg[4],"AlX") == 0) voutput[nthree][1] = AlX;
      else if (strcmp(arg[4],"AlX2") == 0) voutput[nthree][1] = AlX2;
      else if (strcmp(arg[4],"OH2Al") == 0) voutput[nthree][1] = OH2Al;
      else if (strcmp(arg[4],"OH2AlX") == 0) voutput[nthree][1] = OH2AlX;
      else if (strcmp(arg[4],"OH2AlX2") == 0) voutput[nthree][1] = OH2AlX2;
      else if (strcmp(arg[4],"OHAl") == 0) voutput[nthree][1] = OHAl;
      else if (strcmp(arg[4],"OHAlX") == 0) voutput[nthree][1] = OHAlX;
      else if (strcmp(arg[4],"OHAlX2") == 0) voutput[nthree][1] = OHAlX2;
      else if (strcmp(arg[4],"OAl") == 0) voutput[nthree][1] = OAl;
      else if (strcmp(arg[4],"OAlX") == 0) voutput[nthree][1] = OAlX;
      else if (strcmp(arg[4],"OAlX2") == 0) voutput[nthree][1] = OAlX2;

      else error->all(FLERR,"Illegal event command35");

      vA[nthree] = atof(arg[5]);
      vexpon[nthree] = atoi(arg[6]);
      if (vexpon[nthree] != 0.0) error->warning(FLERR,"Illegal vexpon command36");
      vrate[nthree] = atof(arg[7]);
      vcoord[nthree] = atoi(arg[8]);
      vpresson[nthree] = atoi(arg[9]);
      nthree++;

    } else error->all(FLERR,"Illegal event command37");
  } 
  else if (strcmp(command,"pulse_time") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal pulse time");
      T1 = atof(arg[0]);
      T3 = atof(arg[1]);
  }
  else if (strcmp(command,"purge_time") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal purge time");
      T2 = atof(arg[0]);
      T4 = atof(arg[1]);
  }else error->all(FLERR,"Unrecognized command38");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppAldAl2o3::grow_app()
{
  element = iarray[0];
  coord = iarray[1];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppAldAl2o3::init_app()
{
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    firstevent = (int *) memory->smalloc(nlocal*sizeof(int),"app:firstevent");
    //comneigh was defined to avoid double counting of common neighbor in site_propensity
    comneigh = memory->grow(comneigh,12*maxneigh,2,"app/ald:comneigh");
    //comneigh = memory->grow_2d_double_array(comneigh,12*maxneigh,2,"app/ald:comneigh");
    // esites must be large enough for 3 sites and their 1st neighbors
    
    esites = (int *) memory->smalloc(12*maxneigh*sizeof(int),"app:esites");
    //esites = new int[12*maxneigh]; 
  }
  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (coord[i] < -1 || coord[i] > 8) flag = 1;
    if (element[i] < VACANCY) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}
/* ---------------------------------------------------------------------- */

void AppAldAl2o3::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;
  

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;


  if (temperature == 0.0)
    error->all(FLERR,"Temperature cannot be 0.0 for app_ald");
  for (int m = 0; m < none; m++) {
    spropensity[m] = sA[m]*pow(temperature,sexpon[m])*exp(-srate[m]/temperature);
    scount[m] = 0;
  if (spropensity[m] == 0.0) error->warning(FLERR," spropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < ntwo; m++) {
    dpropensity[m] = dA[m]*pow(temperature,dexpon[m])*exp(-drate[m]/temperature);
    dcount[m] = 0;
  if (dpropensity[m] == 0.0) error->warning(FLERR,"dpropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < nthree; m++) {
    vpropensity[m] = vA[m]*pow(temperature,vexpon[m])*exp(-vrate[m]/temperature);
    vcount[m] = 0;
  if (vpropensity[m] == 0.0) error->warning(FLERR,"vpropensity cannot be 0.0 for app_ald");
  }
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppAldAl2o3::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppAldAl2o3::site_propensity(int i)
{
  int j,k,m;


  clear_events(i);

  double proball = 0.0;


  //type I, check species of sites and consider possible events

  // count_coordO was added here to prevent adsorption of metal precursor
  // on the low coordinate oxygen at sublayer

  for (m = 0; m < none; m++) {
    if (element[i] != sinput[m]) continue;
    if ((coord[i] == scoord[m] || scoord[m] == 0) && (spresson[m] == pressureOn || spresson[m] == 0)) {
      add_event(i,1,m,spropensity[m],-1,-1);
      proball += spropensity[m];
    }
  }

  // type II, check species of sites and second neighbor 
  // comneigh variable used to avoid double counting,
  // we consider more than one event between sites, therefore we need 2d comneigh array.

  int nextneib = 1;
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
       for (int kk = 0; kk < numneigh[j]; kk++) {
            k = neighbor[j][kk];
            if (i == k) continue;
            for (m = 0; m < ntwo; m++) {
      if ((element[i] == dinput [m][0] && element[k] == dinput [m][1]) && (dpresson[m] == pressureOn || dpresson[m] == 0) && (coord[i] == dcoord[m] || dcoord[m] == 0)) {
          comevent = 1;
          for (int ii = 0; ii < nextneib; ii++) {
          if ( comneigh[ii][0] == k && comneigh[ii][1] == dpropensity[m]) comevent = 0;
          }
          if (comevent){	
                add_event(i,2,m,dpropensity[m],-1,k);
                proball += dpropensity[m];
            comneigh[nextneib][0] = k;
                comneigh[nextneib][1] = dpropensity[m];
            nextneib++;
        } 
          }
        }
      }
    }
    for (m = 0; m < nextneib; m++) comneigh[m][0] = comneigh[m][1] = 0;

  //type III, check species of sites and first neighbour 


  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
      for (m =0; m < nthree; m++) {
//	if (element[i] == vinput[m][0] && element[j] == vinput[m][1] && (coord[i] == vcoord[m] || vcoord[m] == 0) && (dpresson[m] == pressureOn || dpresson[m] == 0)) { // Zero cn does not allow reaction
  if (element[i] == vinput[m][0] && element[j] == vinput[m][1] && (coord[i] == vcoord[m] ) && (vpresson[m] == pressureOn || vpresson[m] == 0)) {
  add_event(i,3,m,vpropensity[m],j,-1);
  proball += vpropensity[m];
      }
    }
  }


  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppAldAl2o3::site_event(int i, class RandomPark *random)
{
  int j,k,m,n,mm,jj;

  int elcoord = element[i];

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;
  k = events[ievent].kpartner;


  if (rstyle == 1) {
    element[i] = soutput[which];
    scount[which]++;
    } 
  else if (rstyle == 2 && j == -1) {
    element[i] = doutput[which][0];
    element[k] = doutput[which][1];
    dcount[which]++;
    }
  else if (rstyle == 3 && k == -1) {
    element[i] = voutput[which][0];
    element[j] = voutput[which][1];
    vcount[which]++;
    }
  else {printf("Illegal execution event i %d %d %d j %d %d %d k %d %d %d", i, element[i], coord[i], j, element[j], coord[j], k, element[k], coord[k]);
  fflush(stdout); 
  error->all(FLERR,"Illegal execution event"); }

  update_coord(elcoord,i,j,k,which);

  // sequence of ALD, 
  // 1 is metal pulse, 3 purge, 2 oxygen pulse.
  if ((cycle+T1) > time ) {pressureOn = 1;}
  else if ((cycle+T1)<= time && time < (cycle+T1+T2)) {pressureOn = 3;}
  else if ((cycle+T1+T2) <= time && time < (cycle+T1+T2+T3)) {pressureOn = 2;}
  else if ((cycle+T1+T2+T3) <= time && time < (cycle+T1+T2+T3+T4)) {pressureOn = 3;}
  else {cycle += T1+T2+T3+T4; }

  

  int nsites = 0;
  int isite = i2site[i];
  
 

  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  // go from site i to first and second neighbor in all type
  for (n = 0; n < numneigh[i]; n++) {
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
      for (jj = 0; jj< numneigh[m];jj++) {
        mm = neighbor[m][jj];
        isite = i2site[mm];
        if (isite >= 0 && echeck[isite] == 0) {
          propensity[isite] = site_propensity(mm);
          esites[nsites++] = isite;
          echeck[isite] = 1;
        }
     }
  }

  // go from site k to first and second neighbor in type II 
  if (rstyle == 2) {
    for (n = 0; n < numneigh[k]; n++) {
      m = neighbor[k][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
        propensity[isite] = site_propensity(m);
        esites[nsites++] = isite; 
        echeck[isite] = 1;
      }
      for (jj = 0; jj< numneigh[m];jj++) {
        mm = neighbor[m][jj];
        isite = i2site[mm];
        if (isite >= 0 && echeck[isite] == 0) {
          propensity[isite] = site_propensity(mm);
          esites[nsites++] = isite; 
          echeck[isite] = 1;
        }
      }
    }
  }

  // go from site j to first and second neighbor in type III 
  if (rstyle == 3) {
    for (n = 0; n < numneigh[j]; n++) {
      m = neighbor[j][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
        propensity[isite] = site_propensity(m);
        esites[nsites++] = isite;
        echeck[isite] = 1;
      }
      for (jj = 0; jj< numneigh[m];jj++) {
        mm = neighbor[m][jj];
        isite = i2site[mm];
        if (isite >= 0 && echeck[isite] == 0) {
          propensity[isite] = site_propensity(mm);
          esites[nsites++] = isite;
          echeck[isite] = 1;
        }
       }
    }
  }

  solve->update(nsites,esites,propensity);
   // clear echeck array

  for (m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
  
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppAldAl2o3::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppAldAl2o3::add_event(int i, int rstyle, int which, double propensity,
        int jpartner, int kpartner)
{
  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;

  events[freeevent].style = rstyle;
  events[freeevent].which = which;
  events[freeevent].jpartner = jpartner;
  events[freeevent].kpartner = kpartner;
  events[freeevent].propensity = propensity;

  if ( propensity == 0 ) error->all(FLERR,"propensity in add_event wrong app ald");
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   grow list of stored reactions for single and double
------------------------------------------------------------------------- */

void AppAldAl2o3::grow_reactions(int rstyle)
{
  if (rstyle == 1) {
    int n = none + 1;
    srate = (double *) 
      memory->srealloc(srate,n*sizeof(double),"app/ald:srate");
    spropensity = (double *) 
      memory->srealloc(spropensity,n*sizeof(double),"app/ald:spropensity");
    sinput = (int *) 
      memory->srealloc(sinput,n*sizeof(int),"app/ald:sinput");
    soutput = (int *) 
      memory->srealloc(soutput,n*sizeof(int),"app/ald:soutput");
    scount = (int *) 
      memory->srealloc(scount,n*sizeof(int),"app/ald:scount");
    sA = (double *) 
      memory->srealloc(sA,n*sizeof(double),"app/ald:sA");
    sexpon = (int *) 
      memory->srealloc(sexpon,n*sizeof(int),"app/ald:sexpon");
    scoord = (int *) 
      memory->srealloc(scoord,n*sizeof(int),"app/ald:scoord");
    spresson = (int *) 
      memory->srealloc(spresson,n*sizeof(int),"app/ald:spresson");

  } else if (rstyle == 2) {
    int n = ntwo + 1;
    drate = (double *) 
      memory->srealloc(drate,n*sizeof(double),"app/ald:drate");
    dpropensity = (double *) 
      memory->srealloc(dpropensity,n*sizeof(double),"app/ald:dpropensity");
    dinput = memory->grow(dinput,n,2,"app/ald:dinput");
    doutput = memory->grow(doutput,n,2,"app/ald:doutput");
    dcount = (int *) 
      memory->srealloc(dcount,n*sizeof(int),"app/ald:dcount");
    dA = (double *) 
      memory->srealloc(dA,n*sizeof(double),"app/ald:dA");
    dexpon = (int *) 
      memory->srealloc(dexpon,n*sizeof(int),"app/ald:dexpon");
    dcoord = (int *) 
      memory->srealloc(dcoord,n*sizeof(int),"app/ald:dcoord");
    dpresson = (int *) 
      memory->srealloc(dpresson,n*sizeof(int),"app/ald:dpresson");

  } else if (rstyle == 3) {
    int n = nthree + 1;
    vrate = (double *)
      memory->srealloc(vrate,n*sizeof(double),"app/ald:vrate");
    vpropensity = (double *)
      memory->srealloc(vpropensity,n*sizeof(double),"app/ald:vpropensity");
    vinput = memory->grow(vinput,n,2,"app/ald:vinput");
    voutput = memory->grow(voutput,n,2,"app/ald:voutput");
    vcount = (int *)
      memory->srealloc(vcount,n*sizeof(int),"app/ald:vcount");
    vA = (double *)
      memory->srealloc(vA,n*sizeof(double),"app/ald:vA");
    vexpon = (int *)
      memory->srealloc(vexpon,n*sizeof(int),"app/ald:vexpon");
    vcoord = (int *)
      memory->srealloc(vcoord,n*sizeof(int),"app/ald:vcoord");
    vpresson = (int *)
      memory->srealloc(vpresson,n*sizeof(int),"app/ald:vpresson");
  }
}

/* ----------------------------------------------------------------------
   update c.n. for Zn and O, put and remove mask for relative sites
------------------------------------------------------------------------- */
void AppAldAl2o3::update_coord(int elcoord, int i, int j, int k, int which)
{
  // adsorption of TMA yes
  if ((elcoord == O || elcoord == OH || elcoord == OH2) && (element[i] == AlX3O || element[i] == AlX3OH || element[i] == AlX3OH2) && (j == -1)){
    coord[i]=coord[i]+1;
    put_mask(i);
  }
  // desorption of TMA yes
  else if ((elcoord == AlX3O || elcoord == AlX3OH || elcoord == AlX3OH2) && (element[i] == O || element[i] == OH || element[i] == OH2) && (j == -1)){
    coord[i]=coord[i]-1;
    remove_mask(i);
    // count_coordO(i);
  }
  // Ligand elimination from TMA yes
  else if ((elcoord == AlX3OH || elcoord == AlX3OH2) && (element[i] == AlX2O || element[i] == AlX2OH) && (j == -1)){
    remove_mask(i);
    put_mask(i);
  }
  else if ((elcoord == AlX3O || elcoord == AlX3OH || elcoord == AlX3OH2) && (element[i] == AlX2O || element[i] == AlX2OH || element[i] == AlX2OH2) && (element[k] == OH || element[k] == O) && (j == -1)){
    remove_mask(i);
    put_mask(i);
  }
  // Ligand elimination from DMA yes
  else if ((elcoord == AlX2OH || elcoord == AlX2OH2) && (element[i] == AlXO || element[i] == AlXOH) && (j == -1)){
    remove_mask(i);
    put_mask(i);
  }
  else if ((elcoord == AlX2O || elcoord == AlX2OH || elcoord == AlX2OH2) && (element[i] == AlXO || element[i] == AlXOH || element[i] == AlXOH2) && (element[k] == OH || element[k] == O) && (j == -1)){
    remove_mask(i);
    put_mask(i);
  }
  // Ligand elimination from MMA yes
  else if ((elcoord == AlXOH || elcoord == AlXOH2) && (element[i] == AlO || element[i] == AlOH) && (j == -1)){
    remove_mask(i);
  }
  else if ((elcoord == AlXO || elcoord == AlXOH || elcoord == AlXOH2) && (element[i] == AlO || element[i] == AlOH || element[i] == AlOH2) && (element[k] == OH || element[k] == O) && (j == -1)){
    remove_mask(i);
  }
  // TMA dissociation on the surface yes
  else if ((elcoord == AlX3O || elcoord == AlX3OH || elcoord == AlX3OH2) && (element[i] == AlX2O || element[i] == AlX2OH || element[i] == AlX2OH2) && ( element[j] == AlX )){
    remove_mask(i);
    put_mask(i);
    coord[j]=coord[j]+1;
    put_mask(j);
  }
  // DMA dissociation on the surface yes
  else if ((elcoord == AlX2O || elcoord == AlX2OH || elcoord == AlX2OH2) && (element[i] == AlXO || element[i] == AlXOH || element[i] == AlXOH2) && ( element[j] == AlX )){
    remove_mask(i);
    put_mask(i);
    coord[j]=coord[j]+1;
    put_mask(j);
  }
  // MMA dissociation on the surface yes
  else if ((elcoord == AlXO || elcoord == AlXOH || elcoord == AlXOH2) && (element[i] == AlO || element[i] == AlOH || element[i] == AlOH2) && ( element[j] == AlX )){
    remove_mask(i);
    coord[j]=coord[j]+1;
    put_mask(j);
  }
  // water adsorption typeI yes
  else if ((elcoord == Al || elcoord == AlX || elcoord == AlX2) && (element[i] == OH2Al || element[i] == OH2AlX || element[i] == OH2AlX2) && ( j == -1 )) {
    coord[i]=coord[i]+1;
  }
  // water desorption typeI yes
  else if ((elcoord == OH2Al || elcoord == OH2AlX || elcoord == OH2AlX2) && (element[i] == Al || element[i] == AlX || element[i] == AlX2) && ( j == -1 )) {
    coord[i]=coord[i]-1;
  }
  // water desorption typeII no
  else if ((elcoord == OH2) && (element[i] == VACANCY) && ( j == -1 )) {
    count_coord(i,j);
  }  
  // Ligand-exchange reactions during water pulse
  else if ((elcoord == OH2 || elcoord == OH) && (element[i] == OH || element[i] == O) && (element[j] == AlX || element[j] == Al)) {
    if (element[j] == AlX){
      remove_mask(j);
      put_mask(j);
    }
    else {
      remove_mask(j);
      coord[j]=coord[j]-1;
    }
  }
  // MMA with H2O, event I
  else if ((elcoord == OH2AlX || elcoord == OHAlX) && (element[i] == OHAl || element[i] == OAl ) && (j == -1 ) ) { 
    remove_mask(i);
    coord[i]=coord[i]-1;
    }
  // DMA with H2O, event I
  // else if ((elcoord == OH2AlX2 || elcoord == OHAlX2) && (element[i] == OHAlX || element[i] == OAlX ) && (j == -1 ) ) { 
  //   remove_mask(i);
  //   put_mask(i);
  //   }
  // Al densification
  else if ((elcoord == VACANCY) && (element[i] == Al || element[i] == AlX || element[i] == AlX2) && (element[j] == O || element[j] == OH || element[j] == OH2)) {
    if (element[i] == AlX || element[i] == AlX2){
      remove_mask(i, j);
    }
    count_coord(j, i);
    if (element[i] == AlX || element[i] == AlX2){
      put_mask(i);
    }
  }
  // Al reverse densification
  else if ((elcoord == Al || elcoord == AlX || elcoord == AlX2) && (element[i] == VACANCY) && (element[j] >= AlX2O && element[j] <= AlOH2)) {
    if (elcoord == Al){
      count_coord(i, j); // ********* //
    }
    else {
      remove_mask(i, j);
      count_coord(i, j);
      put_mask(j);
    }
  }
  // OH2 densification
  else if ((elcoord == OH2Al || elcoord == OHAl || elcoord == OAl) && (element[i] == Al) && (element[j] == O || element[j] == OH || element[j] == OH2)) {
    count_coord(i, j);
    // count_coordO(j);
  }
  else if ((elcoord == OH2AlX || elcoord == OHAlX || elcoord == OAlX) && (element[i] == AlX) && (element[j] == O || element[j] == OH || element[j] == OH2)) {
    count_coord(i, j);
  }
  else if ((elcoord == OH2AlX2 || elcoord == OHAlX2 || elcoord == OAlX2) && (element[i] == AlX2) && (element[j] == O || element[j] == OH || element[j] == OH2)) {
    count_coord(i, j);
  }
  // OH2 reverse densification
  else if ((elcoord == O || elcoord == OH || elcoord == OH2) && (element[i] == VACANCY) && (element[j] == OAl || element[j] == OHAl || element[j] == OH2Al)) {
    count_coord(i, j);
    // count_coordO(j);
    coord[j] = coord[j] + 1;
  }
  else if ((elcoord == O || elcoord == OH || elcoord == OH2) && (element[i] == VACANCY) && (element[j] == OAlX || element[j] == OHAlX || element[j] == OH2AlX)) {
    count_coord(i, j);
    coord[j] = coord[j] + 1;
  }
  else if ((elcoord == O || elcoord == OH || elcoord == OH2) && (element[i] == VACANCY) && (element[j] == OAlX2 || element[j] == OHAlX2 || element[j] == OH2AlX2)) {
    count_coord(i, j);
    coord[j] = coord[j] + 1;
  }
}
/* ----------------------------------------------------------------------
   put mask for affected sites
------------------------------------------------------------------------- */

void AppAldAl2o3::put_mask(int i)
{
  int isite = i2site[i];
  int nsites = 0;
  esites[nsites++] = isite;
  echeck[isite] = 1;
// Add mask on the second neighbor (oxygen) of the TMA to block adsorption
  if (element[i] == AlX3OH2 || element[i] == AlX3OH || element[i] == AlX3O ){
      for (int n = 0; n < numneigh[i]; n++) {
      int nn = neighbor[i][n];
      isite = i2site[nn];
      if (isite >= 0 && echeck[isite] == 0) { // Cover first neighbour Al site
          coord[isite]=coord[isite]-10;
          esites[nsites++] = isite;
          echeck[isite] = 1;
                        }
      for (int k = 0; k < numneigh[nn]; k++){
        int kk = neighbor[nn][k];
        isite = i2site[kk];
        if (isite >= 0 && echeck[isite] == 0) { // Cover second neighbour O site
          coord[isite]=coord[isite]-10;
          esites[nsites++] = isite;
          echeck[isite] = 1;
        }
      }
      }
    }
// Add mask on the second neighbor (oxygen) of the DMA to block adsorption
  else if (element[i] == AlX2OH2 || element[i] == AlX2OH || element[i] == AlX2O ){
      for (int n = 0; n < numneigh[i]; n++) {
      int nn = neighbor[i][n];
      isite = i2site[nn];
      if (isite >= 0 && echeck[isite] == 0) { // Cover first neighbour Al site
          coord[isite]=coord[isite]-10;
          esites[nsites++] = isite;
          echeck[isite] = 1;
                        }
      for (int k = 0; k < numneigh[nn]; k++){
        int kk = neighbor[nn][k];
        isite = i2site[kk];
        if (isite >= 0 && echeck[isite] == 0) { // Cover second neighbour O site
          coord[isite]=coord[isite]-10;
          esites[nsites++] = isite;
          echeck[isite] = 1;
        }
      }
      }
    }
// Add mask on the second neighbor (oxygen) of the MMA to block adsorption
  else if (element[i] == AlXOH2 || element[i] == AlXOH || element[i] == AlXO ){
      for (int n = 0; n < numneigh[i]; n++) {
      int nn = neighbor[i][n];
      isite = i2site[nn];
      if (isite >= 0 && echeck[isite] == 0) { // Cover first neighbour Al site
          coord[isite]=coord[isite]-10;
          esites[nsites++] = isite;
          echeck[isite] = 1;
                        }
      for (int k = 0; k < numneigh[nn]; k++){
        int kk = neighbor[nn][k];
        isite = i2site[kk];
        if (isite >= 0 && echeck[isite] == 0) { // Cover second neighbour O site
          coord[isite]=coord[isite]-10;
          esites[nsites++] = isite;
          echeck[isite] = 1;
        }
      }
      }
    }
// Add mask to the first neighbor (oxygen) to block adsorption
  else if ( element[i] == AlX || element[i] == AlX2 ){
      for (int n = 0; n < numneigh[i]; n++) {
      int nn = neighbor[i][n];
      isite = i2site[nn];
      if (isite >= 0 && echeck[isite] == 0) { // Cover first neighbour O site
        coord[isite]=coord[isite]-10;
        esites[nsites++] = isite;
        echeck[isite] = 1;
      }
            for (int k = 0; k < numneigh[nn]; k++){
                int kk = neighbor[nn][k];
                isite = i2site[kk];
                if (isite >= 0 && echeck[isite] == 0) { // Cover second neighbour Al site
                    coord[isite] = coord[isite]-10;
                    esites[nsites++] = isite;
                    echeck[isite] = 1;
                }

            }
        }
  }
  for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}

/* ----------------------------------------------------------------------
   remove mask 
------------------------------------------------------------------------- */

void AppAldAl2o3::remove_mask(int i, int j) // j flag for when Al densification
{
  int isite = i2site[i];
  int nsites = 0;
  esites[nsites++] = isite;
  echeck[isite] = 1;

// Remove mask from oxygen sites after TMA desorption, dissociation and ligand elimination
  if ( element[i] == O || element[i] == OH || element[i] == OH2 || element[i] == AlX2O || element[i] == AlX2OH || element[i] == AlX2OH2 ){
      for (int n = 0; n < numneigh[i]; n++) {
      int nn = neighbor[i][n];
      isite = i2site[nn];
      if (isite >= 0 && echeck[isite] == 0) { // Remove first neighbour Al site
                coord[isite]=coord[isite]+10;
                esites[nsites++] = isite;
                echeck[isite] = 1;
            }
      for (int k = 0; k < numneigh[nn]; k++){
        int kk = neighbor[nn][k];
        isite = i2site[kk];
        if (isite >= 0 && echeck[isite] == 0) { // Remove second neighbour O site
          coord[isite]=coord[isite]+10;
          esites[nsites++] = isite;
          echeck[isite] = 1;
        }
        }
        }
  }
// DMA dissociation and ligand elimination
  else if ( element[i] == AlXO || element[i] == AlXOH || element[i] == AlXOH2 ){
      for (int n = 0; n < numneigh[i]; n++) {
      int nn = neighbor[i][n];
      isite = i2site[nn];
      if (isite >= 0 && echeck[isite] == 0) { // Remove first neighbour Al site
                coord[isite]=coord[isite]+10;
                esites[nsites++] = isite;
                echeck[isite] = 1;
            }
      for (int k = 0; k < numneigh[nn]; k++){
        int kk = neighbor[nn][k];
        isite = i2site[kk];
        if (isite >= 0 && echeck[isite] == 0) { // Remove second neighbour O site
          coord[isite]=coord[isite]+10;
          esites[nsites++] = isite;
          echeck[isite] = 1;
        }
        }
        }
  }
// MMA dissociation and ligand elimination
  else if ( element[i] == AlO || element[i] == AlOH || element[i] == AlOH2 ){
      for (int n = 0; n < numneigh[i]; n++) {
      int nn = neighbor[i][n];
      isite = i2site[nn];
      if (isite >= 0 && echeck[isite] == 0) { // Remove first neighbour Al site
                coord[isite]=coord[isite]+10;
                esites[nsites++] = isite;
                echeck[isite] = 1;
            }
      for (int k = 0; k < numneigh[nn]; k++){
        int kk = neighbor[nn][k];
        isite = i2site[kk];
        if (isite >= 0 && echeck[isite] == 0) { // Remove second neighbour O site
          coord[isite]=coord[isite]+10;
          esites[nsites++] = isite;
          echeck[isite] = 1;
        }
        }
        }
  }
// Al densification
  else if ( ( element[i] == AlX || element[i] == AlX2 ) && ( element[j] == O || element[j] == OH || element[j] == OH2 )){ 
      echeck[i2site[i]] = 0;
      for (int n = 0; n < numneigh[j]; n++) {
          int nn = neighbor[j][n];
          isite = i2site[nn];
          if (isite >= 0 && echeck[isite] == 0 ) { // Remove first neighbour Zn site
              coord[isite]=coord[isite]+10;
                esites[nsites++] = isite;
                echeck[isite] = 1;
            }
            for (int k = 0; k < numneigh[nn]; k++){
                int kk = neighbor[nn][k];
                if(kk != j){
                    isite = i2site[kk];
                    if (isite >= 0 && echeck[isite] == 0) { // Remove second neighbour O site
                        if(isite!=j){coord[isite] = coord[isite]+10;}
                        esites[nsites++] = isite;
                        echeck[isite] = 1;
                    }
                }

            }
        }
    }

// Remove mask after second ligand has been removed
  else if ( element[i] == OAl || element[i] == OHAl  ||  element[i] == AlX || element[i] == Al ){
      for (int n = 0; n < numneigh[i]; n++) {
            int nn = neighbor[i][n];
            isite = i2site[nn];
            if (isite >= 0 && echeck[isite] == 0 ) {
                    coord[isite]=coord[isite]+10;
                    esites[nsites++] = isite;
                    echeck[isite] = 1;
            }
            for (int k = 0; k < numneigh[nn]; k++){
                int kk = neighbor[nn][k];
                isite = i2site[kk];
                if (isite >= 0 && echeck[isite] == 0) {
                    coord[isite] = coord[isite]+10;
                    esites[nsites++] = isite;
                    echeck[isite] = 1;
                }

            }
        }
  }


// Remove mask after reverse densification
  else if ( ( element[i] == VACANCY ) && ( element[j] >= AlX2O && element[j] <= AlOH2)){
        for (int n = 0; n < numneigh[i]; n++) {
            int nn = neighbor[i][n];
            isite = i2site[nn];
            if (isite >= 0 && echeck[isite] == 0 ) {
                coord[isite]=coord[isite]+10;
                esites[nsites++] = isite;
                echeck[isite] = 1;
            }
            for (int k = 0; k < numneigh[nn]; k++){
                int kk = neighbor[nn][k];
                isite = i2site[kk];
                if (isite >= 0 && echeck[isite] == 0) {
                    coord[isite] = coord[isite]+10;
                    esites[nsites++] = isite;
                    echeck[isite] = 1;
                }

            }
        }
    }

  for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}

/* ----------------------------------------------------------------------
   count c.n after densification
------------------------------------------------------------------------- */


void AppAldAl2o3::count_coord(int i, int j) // i: Oxygen species, j: Zinc species (does not necessarily hold)
{

// Densification Al
  if (( element[i] == O || element[i] == OH || element[i] == OH2 ) &&  ( element[j] == Al || element[j] == AlX || element[j] == AlX2)  ){
      if (element[j] == AlX || element[j] == AlX2){
        coord[j]=coord[j]+1; // Add one because of X ligand
      }
      for (int s = 0; s < numneigh[j]; s++){
        int nn = neighbor[j][s];
        if (element[nn] >= O && element[nn] <= AlOH2 ) { // Check if neighbouring site is an oxygen site
          coord[j]=coord[j] + 1;
          if (i != nn){ coord[nn]=coord[nn]+1; } // Careful not to change the cn of original site
        }
      }
    }
// Desorption of OH2, event 1 and 3
  else if ( element[i] == VACANCY  && ( j == -1 ) ){
        for (int s = 0; s < numneigh[i]; s++){
            int nn = neighbor[i][s];
            if ( element[nn] >= Al && element[nn] <= OAlX2 ) { // Check if neighbouring site is a Al site
                coord[i]=coord[i] - 1;
                if (i != nn){ coord[nn]=coord[nn] - 1;}
            }
        }
    }
// Reverse densification on ZnX
  else if ( element[i] == VACANCY  && (element[j] >= AlX2O && element[j] <= AlOH2 ) ){
        if ( element[j] == AlO || element[j] == AlOH || element[j] == AlOH2 ){ 
          coord[i]=coord[i] - 0; 
          } 
        else { coord[i]=coord[i] - 1; } // Remove the extra cn from ligand

        for (int s = 0; s < numneigh[i]; s++){
            int nn = neighbor[i][s];
            if ( element[nn] >= O && element[nn] <= AlOH2 ) { // Check if neighbouring site is an oxygen site
                coord[i]=coord[i] - 1;
                if (j != nn){ coord[nn]=coord[nn] - 1;
                }
            }
        }
    }
// Densification of oxygen-species
  else if (( element[i] == Al || element[i] == AlX || element[i] == AlX2) &&  ( element[j] == O || element[j] == OH || element[j] == OH2 )  ){
        for (int s = 0; s < numneigh[j]; s++){
                int nn = neighbor[j][s];
                if ( element[nn] >= Al && element[nn] <= OAlX2) { // Check if neighbouring site is an Al site
                        coord[j]=coord[j] + 1;
                        if (i != nn){ coord[nn]=coord[nn]+1; 
                        }
                }
        }
    }
// Reverse densification of OH2 / OH / O
  else if ( element[i] == VACANCY  && ( element[j] >= OH2Al && element[j] <= OAlX2) ){
        for (int s = 0; s < numneigh[i]; s++){
            int nn = neighbor[i][s];
            if ( element[nn] >= Al && element[nn] <= OAlX2 ) { // Check if neighbouring site is a Al site
                coord[i]=coord[i] - 1;
                if (i != nn){ coord[nn]=coord[nn] - 1;}
            }
        }
    }
}

/* ----------------------------------------------------------------------
   count c.n of oxygen before adsorption
------------------------------------------------------------------------- */
void AppAldAl2o3::count_coordO(int i)
{
    int fullO = 0;
    int emptyO = 0;
    int totalS = 0;

    int isite = i2site[i];
    int nsites = 0;

  for (int m = 0; m < numneigh[i]; m++) {
    int mm = neighbor[i][m];
    for (int s = 0; s < numneigh[mm]; s++) {
      int ss = neighbor[mm][s];
      isite = i2site[ss];
      if (i==ss)  continue;
      if (isite >= 0 && echeck[isite] == 0) {
        if ( element[ss] >= O && element[ss] <= AlOH2 ) {fullO++;}
        else if (element[ss] == VACANCY) {emptyO++;}
              esites[nsites++] = isite;
              echeck[isite] = 1;
            }
    
    }
  }
   totalS = fullO+emptyO;
   if ( float(fullO) > 4*totalS/5 and coord[i] > -20) {coord[i] += -20;} // decrease the coord of the oxygen site to render it inactive for adsorption
   for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}

