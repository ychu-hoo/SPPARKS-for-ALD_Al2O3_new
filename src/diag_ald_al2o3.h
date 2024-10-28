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

#ifdef DIAG_CLASS
DiagStyle(ald/al2o3,DiagAldAl2o3)

#else

#ifndef SPK_DIAG_ALD_Al2o3_H
#define SPK_DIAG_ALD_Al2o3_H

#include "diag.h"

namespace SPPARKS_NS {

class DiagAldAl2o3 : public Diag {
 public:
  DiagAldAl2o3(class SPPARKS *, int, char **);
  ~DiagAldAl2o3();
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:
  class AppAldAl2o3 *appaldal2o3;
  int nlist;
  char **list;
  int *which,*index,*ivector;
  int siteflag;
};

}

#endif
#endif
