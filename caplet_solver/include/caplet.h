/*
Created : Jul 28, 2010
Modified: Feb 14, 2013
Author  : Yu-Chung Hsiao
Email   : project.caplet@gmail.com
*/

/*
This file is part of CAPLET.

CAPLET is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CAPLET is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CAPLET.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CAPLET_H_
#define CAPLET_H_

#include <fstream>
#include <iostream>
#include <string>

#include "caplet_blas.h"
#include "caplet_const.h"
#include "caplet_parameter.h"

namespace caplet {

class Caplet
{
 public:  //* enum
  enum MODE
  {
    FAST_GALERKIN,
    DOUBLE_GALERKIN,
    DOUBLE_COLLOCATION
  };

  enum ERROR_REF
  {
    DIAGONAL,
    SELF
  };

 public:  //* functions
  Caplet();
  ~Caplet();

  void extractC(MODE mode = DOUBLE_GALERKIN);

  int getNPanels() const;
  int getNCoefs() const;
  int getSizeP() const;
  int getSizeCoefs() const;
  int getSizeCmat() const;
  const float* const getCmat() const;

  float compareCmatError(const float* cmatRef, ERROR_REF option) const;
  float compareCmatError(const Caplet* const caplet,
                         ERROR_REF option = DIAGONAL) const;
  float compareCmatError(const std::string filename,
                         ERROR_REF option = DIAGONAL) const;

  void clear();
  void loadCapletFile(const std::string filename);
  void loadFastcapFile(const std::string filename);
  void saveCmat(const std::string filename);
  void saveCoefs(const std::string filename);

  void printP(std::ostream& fout = std::cout);
  void printCoefs();
  void printRHS(std::ostream& fout = std::cout);
  void printCmat();
  void printPanel(int index);

  //* Versioned algorithms

 private:  //* variables
  static float epsilon0;
  static float pi;

  bool isInstantiable;
  bool isLoaded;
  bool isSolved;
  MODE mode;

  int nCoefs;
  int nWires;
  int nPanels;
  int* nWirePanels;
  int* nWireCoefs;

  float (*panels)[nDim][nBit];
  int* dirs;
  float* areas;

  int* indexIncrements;

  char* basisTypes;
  int* basisDirs;
  float* basisZs;
  float* basisShifts;

  //* [P] [coefs] = [rhs]
  float* P;      //* symmetric positive definite and diagonally dominant
  float* rhs;    //* right hand side
  float* coefs;  //* solution
  float* Cmat;   //* capacitance matrix

  //* double precision version
  //  [dP] [dcoefs] = [drhs]
  double* dP;
  double* drhs;
  double* dcoefs;
  double* dCmat;

#ifdef CAPLET_TIMER
  double timeStart;
  double timeAfterFilling;
  double timeAfterSolving;

  double fillingTime;
  double solvingTime;
  double totalTime;
#endif

  bool flagMergeProjection1_0;

 private:  //* functions
  void extractCCollocationDouble();
  void extractCGalerkin();
  void extractCGalerkinDouble();

  void generateCollocationPMatrixDouble();
  void generateRHSDouble();

  void generateGalerkinPMatrix();
  void generateGalerkinPMatrixDouble();
  void generateGalerkinPMatrixMPI();
  void generateGalerkinPMatrixDoubleMPI();
  void generateRHS();

  void modifyPanelAspectRatio();
  bool isPanelAspectRatioValid();

  Shape selectShape(int panel);

 private:
  void rotateX2Z(float* coord_ptr[3][4], int panelNo);
  void rotateY2Z(float* coord_ptr[3][4], int panelNo);
  void rotateZ2Z(float* coord_ptr[3][4], int panelNo);
  void mirrorY2X(float* coord_ptr[3][4], int panelNo);
  void mirrorY2Z(float* coord_ptr[3][4], int panelNo);
  void mirrorX2Z(float* coord_ptr[3][4], int panelNo);

  double calCollocationPEntryDouble(int panel1, int panel2);
  float calGalerkinPEntry(int panel1, int panel2);
  double calGalerkinPEntryDouble(int panel1, int panel2);
};

}  // namespace caplet

#endif /* CAPLET_H_ */
