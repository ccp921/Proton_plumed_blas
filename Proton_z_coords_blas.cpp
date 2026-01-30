/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* =========================================================================
   Modified from the script Proton.cpp on https://www.acmm.nl/molsim/users/ensing/software/index.html#ptcv
   Change the group allatoms to two groups GROUPA and GROUPB, which contain O atoms and H atoms, respectively
   Add the z coordinates to calculate the r_O(z) for proton in 10.1021/acs.jpcc.9b09715
   The distref is set to (0,0,0), i.e., the position of each atom in GROUPA is used
   ccp, 2025/08
========================================================================= */
#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

#include "tools/Pbc.h"
#include <string>
#include <cmath>

#include <iostream>
#include <iomanip>  // needed to use manipulators with parameters (precision, width)
#include <fstream>
#include <cblas.h>

using namespace std;

namespace PLMD{
namespace vatom{

//+PLUMEDOC VATOM CENTER
/*
Calculate the position of the PROTON CV 

The computed
virtual atom is a proton tracker and be accessed in
an atom list through the label for the PROTON action that creates it.

When running with periodic boundary conditions, it uses the inbuilt plumed PBC




\par Examples

# definition of CVs

c1: PROTON GROUPA=waterO GROUPB=waterH NN=6 MM=12 LAMBDA=100 R_0=1.25 OFFSET=1.8
DUMPATOMS STRIDE=10 FILE=proton.xyz ATOMS=c1

\verbatim
\endverbatim
(See also \ref DISTANCE, \ref MOVINGRESTRAINT and \ref PRINT).

*/
//+ENDPLUMEDOC


class Proton:
  public ActionWithVirtualAtom
{
  std::vector<double> coordnums, weight;
  std::vector< std::vector<Vector> > derivcn;
  std::vector< std::vector<double> > derivw;
  std::vector<Vector> distref;

  int no,nh,nn,mm,j;
  double r0,lambda,offset;
  bool pbc,forcebubble,unbias;
public:
  explicit Proton(const ActionOptions&ao);
  void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Proton,"PROTON")

void Proton::registerKeywords(Keywords& keys){
  ActionWithVirtualAtom::registerKeywords(keys);
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","OFFSET","2","The offset parameter of the weight function. ");
  keys.add("compulsory","LAMBDA","8","The lambda parameter of the weight function. ");
  keys.add("compulsory","NO","64","The number of oxygen atoms ");
  keys.add("compulsory","NH","128","The number of hydrogen atoms.");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("FORCEBUBBLE",false,"ignore the artificial force bubble around imax");
  keys.addFlag("MASS",false,"If set center is mass weighted");
  keys.addFlag("UNBIAS",false,"ignore the calculation of derivatives to accelerate in unbiased simulations");
  // keys.add("atoms","GROUPA","Calculate the distances between all the atoms in GROUPA and all "
  //                             "the atoms in GROUPB. This must be used in conjuction with GROUPB.");
  // keys.add("atoms","GROUPB","Calculate the distances between all the atoms in GROUPA and all "
  //                             "the atoms in GROUPB. This must be used in conjuction with GROUPA.");
  keys.add("atoms","ALLATOMS","Calculate the distances between all the atoms, with the first part of O atoms and second part of H atoms.");
}

Proton::Proton(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao),
  pbc(true),
  forcebubble(false),
  unbias(false)
{

  // group
  vector<AtomNumber> allatoms;
  parseAtomList("ALLATOMS",allatoms);

  // no = groupa.size(); // number of O atoms
  // nh = groupb.size(); // number of H atoms
  
  parse("R_0",r0);
  if(r0<=0.0) error("R_0 should be explicitly specified and positive");
  
  parse("NN",nn);
  parse("MM",mm);
  parse("NO",no);
  parse("NH",nh);
  parse("LAMBDA",lambda);
  parse("OFFSET",offset);
  
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  parseFlag("FORCEBUBBLE",forcebubble);
  parseFlag("UNBIAS",unbias);
  pbc=!nopbc;
  checkRead();
  Vector zero;

  // log.printf("  of groupa of atoms");
  // for(unsigned i=0;i<groupa.size();++i) log.printf(" %d",groupa[i].serial());
  // log.printf("  \n");

  // log.printf("  of groupb of atoms");
  // for(unsigned i=0;i<groupb.size();++i) log.printf(" %d",groupb[i].serial());
  // log.printf("  \n");

  log.printf("  of group of atoms");
  for(unsigned i=0;i<allatoms.size();++i) log.printf(" %d",allatoms[i].serial());
  log.printf("  \n");

  derivcn.resize(no);
  for(unsigned i=0;i<no;i++) derivcn[i].resize(no+nh);

  derivw.resize(no);
  for(unsigned i=0;i<no;i++) derivw[i].resize(no);

  distref.resize(no);
  coordnums.resize(no);
  weight.resize(no);

  log.printf("the offset of the weight function is %f \n", offset);
  log.printf("\n");
 
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  if(nopbc){
    log<<"  PBC will be ignored\n";
  } else {
    log<<"  broken molecules will be rebuilt assuming atoms are in the proper order\n";
  }
  
  if(forcebubble) log.printf(" Forces will be around 3 Ang of imax \n");

  if(unbias) log.printf(" !!!!! Calculations of derivatives will be ignored to accelerate, which can ONLY be used in unbiased simulations !!!!! \n");

  log.printf("  number of H-acceptors: %d\n",no); 
  log.printf("  number of hydrogens  : %d\n",nh); 
  log.printf("  coord.num. r0        : %lf\n",r0); 
  log.printf("  coord.num. nn        : %d\n",nn); 
  log.printf("  coord.num. mm        : %d\n",mm); 
  log.printf("  weighting lambda      : %lf\n",lambda);

  // requestAtoms(groupa);
  requestAtoms(allatoms);
}


void Proton::calculate(){

  Vector pos,temp,poscheck,zero;
  vector<Tensor> deriv(getNumberOfAtoms());
  vector<Tensor> derivcheck(getNumberOfAtoms());
  int io, i;
  double dr=0.0, sum=0.0;
  int imax=0;
  double max=0.0,dfunc=0.0;
  double invr0=1.0/r0;

  /* initialize things to zero */
  for(unsigned i=0;i<getNumberOfAtoms();i++){
    for(unsigned k=0;k<3;k++){
      for(unsigned m=0;m<3;m++){
        deriv[i][k][m]=0.0;
        derivcheck[i][k][m]=0.0;
      }
    }
  }
  for(unsigned j=0;j<3;j++) zero[j] = 0.0;

  sum=0.0;
  for(unsigned io=0;io<no;io++){
    distref[io] = zero;
    coordnums[io] = 0.0;
    weight[io] = 0.0;
    for(unsigned ih=0;ih<no+nh;ih++){
      for(unsigned k=0;k<3;k++){
        derivcn[io][ih][k] = 0.0;
      }
    }
    for(unsigned ih=no;ih<nh+no;ih++){
      Vector distance = pbcDistance(getPosition(io), getPosition(ih));
      dr = distance.modulo();
      const double rdist = dr * invr0;
      double rNdist=Tools::fastpow(rdist,nn-1);
      double rMdist=Tools::fastpow(rdist,mm-1);
      double num = 1.0-rNdist*rdist;
      double iden = 1.0/(1.0-rMdist*rdist);
      double func = num*iden;

      coordnums[io] += func;
      dfunc = ((-nn*rNdist*iden)+(func*(iden*mm)*rMdist));
      dfunc*=invr0;
      dfunc/=dr;
      Vector dd(dfunc*distance);
      derivcn[io][ih] = dd;
      derivcn[io][io] -= dd;
    }
    for(unsigned j=0;j<no;j++){
      derivw[io][j] = 0.0;
    }
    /* weight associated with each oxygen atom */
    weight[io] = exp( lambda * (coordnums[io] - offset) );
    sum += weight[io];
  }
  
  /* find oxygen with highest coordination number */
  max=-999.9;
  for(unsigned i=0;i<no;i++){
    weight[i] /= sum;
    if(weight[i]>max){
      imax=i;
      max=weight[i];
    }
  }
  
  /* make force bubble around imax by setting weights outside to zero */
  if(forcebubble){
    sum=0.0;
    for(unsigned iat=0;iat<no;iat++){
      Vector disderiv = pbcDistance(getPosition(imax), getPosition(iat));
      dr = disderiv.modulo();
      if(dr<3.0){
        weight[iat] = exp( lambda * (coordnums[iat] - offset) );
        sum += weight[iat];
      }else{
        weight[iat]=0.0;
      }
    }
    for(unsigned iat=0;iat<no;iat++) weight[iat] /= sum;
  }  

  if(!unbias){ // calculate the derivatives
    /* derivatives of weight with respect to coordination number */
    std::vector<double> derivw_1d(no*no); // unfold derivw
    std::vector<double> distref_1d_trans(3*no); // unfold distref, including transposation
    pos=getPosition(imax);
    for(unsigned i=0;i<no;i++){
      for(unsigned j=0;j<no;j++){
        if(j!=i){
          derivw[i][j] = -1.0*weight[i]*lambda*weight[j];
          derivw_1d[i*no + j] = -1.0*weight[i]*lambda*weight[j];
        }
      }
      derivw[i][i] = lambda*weight[i]*(1.0-weight[i]);
      derivw_1d[i*no + i] = lambda*weight[i]*(1.0-weight[i]);
      /* sum weighted distances from reference oxygen */
      distref[i] = pbcDistance(getPosition(imax), getPosition(i));
      pos += distref[i] * weight[i];
      for(unsigned m=0;m<3;m++){
        distref_1d_trans[m*no + i] = distref[i][m];
      }
    }
    
    // Temp = distref^T * derivw  (3xno)
    std::vector<double> Temp(3*no, 0.0);
    // call BLAS:  (M=3, N=no, K=no)
    // Temp(no×3) = (distref^T)(3×no) * derivw(no×no)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, no, no,
                1.0, distref_1d_trans.data(), no,
                    derivw_1d.data(), no,
                0.0, Temp.data(), no);

    /* derivatives of weight with respect to coordination number */
    int ntotal = getNumberOfAtoms(); // total number of atoms in group
    for(unsigned k=0;k<3;k++){
      // get ntotal*no slices
      std::vector<double> derivcn_slice_1d(no*ntotal); // unfold derivcn[][][k]
      for(unsigned l=0;l<no;l++){
        for(unsigned i=0;i<ntotal;i++){
          derivcn_slice_1d[l*ntotal + i] = derivcn[l][i][k];
        }
      }
      // Result = Temp * derivcn (3xntotal)
      std::vector<double> Result(3*ntotal, 0.0);
      // call BLAS:  (M=3, N=ntotal, K=no)
      // Result(3×ntotal) = Temp(3×no) * derivcn(no×ntotal)
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                  3, ntotal, no,
                  1.0, Temp.data(), no,
                      derivcn_slice_1d.data(), ntotal,
                  0.0, Result.data(), ntotal);
      // fill the slice in deriv
      for(unsigned m=0;m<3;m++){
        for(unsigned i=0;i<ntotal;i++){
          deriv[i][k][m] = Result[m*ntotal + i];
        }
      }
    }

    for(unsigned k=0;k<3;k++) deriv[imax][k][k]+=1.0;

    for(unsigned j=0;j<no;j++){
      if(j!=imax){ 
        for(unsigned k=0;k<3;k++){
            deriv[j][k][k]+=weight[j];
            deriv[imax][k][k]-=weight[j];
        }
      }
    }

    // // output for comparison with original calculation
    // for(unsigned k=0;k<3;k++){
    //   for(unsigned m=0;m<3;m++){
    //     log << deriv[imax][k][m] << " ";
    //   }
    // }
    // log << "\n";
  }else{ // only calculate the position of proton-attached oxygen atom in unbiased simulations
    /* sum weighted distances from reference oxygen */
    pos=getPosition(imax);
    for(unsigned i=0;i<no;i++){
      distref[i] = pbcDistance(getPosition(imax), getPosition(i));
      pos += distref[i] * weight[i];
    }
  }
  
  /* set final positions and derivatives */
  setPosition(pos);
  //  setMass(mass);
  setAtomsDerivatives(deriv); 

}   /* for calculate */ 
}   /* for namespace vatom */
}   /* for namespace PLMD */
