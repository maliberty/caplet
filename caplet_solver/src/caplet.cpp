/*
Created : 2010-07-28
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

#include "caplet.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <string>

#include "caplet_blas.h"
#include "caplet_const.h"
#include "caplet_elem.h"
#include "caplet_gauss.h"
#include "caplet_int.h"
#include "caplet_parameter.h"
#include "caplet_widgets.h"
#include "mpi.h"

namespace caplet {

using namespace std;

//**** Private constants
float Caplet::epsilon0 = 8.8541878176e-12f;
float Caplet::pi = 3.1415926535f;

//**********************************
//*
//* Shape functions
//*
//*
static float flat(const float x, const float w)
{
  return 1;
}

static float arch(const float x, float w)
{
//* Constant arch shapes are good enough for most of the time
#ifdef CAPLET_FLAT_ARCH
  return 0.5;
#endif

#ifdef DEBUG_SHAPE_BOUNDARY_CHECK
  if (x > 1.41e-6 || x < 0 || w < 0) {
    std::cerr << "ERROR: input argument is out of the arch boundary"
              << std::endl;
  }
#endif

  //* Extracted for arch length = 0.4um
  if (w > 10e-6) {
    w = 10e-6;
  }
  return (-1.1065e-06 + (241.98 + 5.1765e+06 * w) * w
          + (350.58 + -2.3076e+08 * x + -5.5044e+07 * w) * x)
         / (-1.0785e-05 + (404.46 + (2.6513e+06 + 2.7367e+11 * w) * w) * w
            + (456.37 + (1.8301e+08 + -3.3833e+14 * x) * x
               + (-1.0181e+07 + 5.4154e+12 * w + 1.0317e+15 * x) * w)
                  * x);
}

static float side(const float x, float w)
{
#ifdef CAPLET_FLAT_SIDE
  return 0.5;
#endif

#ifdef DEBUG_SHAPE_BOUNDARY_CHECK
  if (x > 1.251e-6 || x < 0 || w < 0) {
    std::cerr << "ERROR: input argument is out of the arch boundary"
              << std::endl;
  }
#endif

  //* Extracted for arch length = 0.4um
  if (w > 10e-6) {
    w = 10e-6;
  }
  return (0.00012065 + (237.09 + -8.6825e+06 * w) * w
          + (1580.3 + -1.7957e+09 * x + -2.3054e+07 * w) * x)
         / (0.00029255 + (398.95 + -1.6467e+07 * w) * w
            + (2789.4 + 1.447e+09 * x + 1.0577e+08 * w) * x);
}

//__________________________________________________________
//*
//* Ctor and Dtor
//*

Caplet::Caplet()
    : isLoaded(false), isSolved(false), flagMergeProjection1_0(true)
{
}

Caplet::~Caplet()
{
  clear();
}

void Caplet::clear()
{
  if (isLoaded == true) {
    isLoaded = false;
    isInstantiable = false;

    delete[] nWirePanels;
    delete[] nWireCoefs;

    delete[] panels;

    delete[] dirs;
    delete[] areas;
    delete[] indexIncrements;
    delete[] basisTypes;
    delete[] basisDirs;
    delete[] basisZs;
    delete[] basisShifts;

    if (isSolved == true) {
      isSolved = false;

      //* single precision version
      switch (mode) {
        case FAST_GALERKIN:
          delete[] P;
          delete[] rhs;
          delete[] coefs;

          delete[] Cmat;
          break;
        case DOUBLE_GALERKIN:
        case DOUBLE_COLLOCATION:
          delete[] dP;
          delete[] drhs;
          delete[] dcoefs;

          delete[] dCmat;
          break;
        default:;
      }
    }
  }
  isInstantiable = false;
}

//__________________________________________________________
//*
//* Versioned algorithms
//*

//__________________________________________________________
//*
//* Public functions
//*

void Caplet::loadFastcapFile(const std::string filename)
{
  //* Clean up
  clear();

  ifstream ifile(filename.c_str());
  //* Check if open successfully
  if (!ifile) {
    cerr << "ERROR: Cannot open the file: " << filename << endl;
    return;
  }

  //* Initialize
  nPanels = 0;
  nWires = 0;

  string stringTemp;
  string lineTemp;
  string currentConductorName = "";
  char charTemp;
  list<int> nWirePanelsList;
  int nWirePanelTemp = 0;

  float xCoord[4];
  float yCoord[4];
  float zCoord[4];

  float xmin, xmax, ymin, ymax, zmin, zmax;

  list<float> x1;
  list<float> x2;
  list<float> y1;
  list<float> y2;
  list<float> z1;
  list<float> z2;

  //* skip the first line
  getline(ifile, lineTemp);

  while (getline(ifile, lineTemp)) {
    nWirePanelTemp++;
    nPanels++;

    std::stringstream stringTokenizer(lineTemp);

    stringTokenizer >> charTemp;
    if (charTemp != 'Q') {
      cerr << "ERROR: only 'Q' is allowed in this solver!" << endl;
    }

    stringTokenizer >> stringTemp;

    stringTokenizer >> xCoord[0] >> yCoord[0] >> zCoord[0]

        >> xCoord[1] >> yCoord[1] >> zCoord[1]

        >> xCoord[2] >> yCoord[2] >> zCoord[2]

        >> xCoord[3] >> yCoord[3] >> zCoord[3];

    xmin = xCoord[0];
    xmax = xCoord[0];
    ymin = yCoord[0];
    ymax = yCoord[0];
    zmin = zCoord[0];
    zmax = zCoord[0];

    for (int i = 1; i < 4; i++) {
      if (xmin > xCoord[i]) {
        xmin = xCoord[i];
      }
      if (xmax < xCoord[i]) {
        xmax = xCoord[i];
      }
      if (ymin > yCoord[i]) {
        ymin = yCoord[i];
      }
      if (ymax < yCoord[i]) {
        ymax = yCoord[i];
      }
      if (zmin > zCoord[i]) {
        zmin = zCoord[i];
      }
      if (zmax < zCoord[i]) {
        zmax = zCoord[i];
      }
    }

    x1.push_back(xmin);
    x2.push_back(xmax);
    y1.push_back(ymin);
    y2.push_back(ymax);
    z1.push_back(zmin);
    z2.push_back(zmax);

    if (currentConductorName.compare(stringTemp) != 0) {  // is different
      //* change to another conductor
      currentConductorName = stringTemp;
      nWires++;
      nWirePanelsList.push_back(nWirePanelTemp);
      nWirePanelTemp = 0;
    }
  }
  //* Register the last wire information
  nWirePanelsList.push_back(++nWirePanelTemp);

  //* Collect all information from the .caplet file until here
  //  fill in and allocate each array and field below

  nCoefs = nPanels;

  //* Convert nWirePanelsList to the nWirePanels array
  nWirePanels = new int[nWires];
  nWireCoefs = new int[nWires];
  int i = 0;
  for (std::list<int>::iterator it = ++nWirePanelsList.begin();
       it != nWirePanelsList.end();
       ++it) {
    nWirePanels[i] = *it;
    nWireCoefs[i] = *it;
    i++;
  }

  panels = new float[nPanels][3][4];
  dirs = new int[nPanels];
  areas = new float[nPanels];
  for (int i = 0; i < nPanels; i++) {
    panels[i][X][MIN] = x1.front();
    x1.pop_front();
    panels[i][X][MAX] = x2.front();
    x2.pop_front();
    panels[i][X][LENGTH] = panels[i][X][MAX] - panels[i][X][MIN];
    panels[i][X][CENTER] = (panels[i][X][MAX] + panels[i][X][MIN]) / 2;

    panels[i][Y][MIN] = y1.front();
    y1.pop_front();
    panels[i][Y][MAX] = y2.front();
    y2.pop_front();
    panels[i][Y][LENGTH] = panels[i][Y][MAX] - panels[i][Y][MIN];
    panels[i][Y][CENTER] = (panels[i][Y][MAX] + panels[i][Y][MIN]) / 2;

    panels[i][Z][MIN] = z1.front();
    z1.pop_front();
    panels[i][Z][MAX] = z2.front();
    z2.pop_front();
    panels[i][Z][LENGTH] = panels[i][Z][MAX] - panels[i][Z][MIN];
    panels[i][Z][CENTER] = (panels[i][Z][MAX] + panels[i][Z][MIN]) / 2;

    if (std::abs(panels[i][X][LENGTH]) < zero) {  // in x-dir
      dirs[i] = X;
      areas[i] = panels[i][Y][LENGTH] * panels[i][Z][LENGTH];
    } else if (std::abs(panels[i][Y][LENGTH]) < zero) {  // in y-dir
      dirs[i] = Y;
      areas[i] = panels[i][Z][LENGTH] * panels[i][X][LENGTH];
    } else {  // in z-dir
      dirs[i] = Z;
      areas[i] = panels[i][X][LENGTH] * panels[i][Y][LENGTH];
    }
  }

  //* For a uniform treatment of generateRHS
  indexIncrements = new int[nPanels];
  basisTypes = new char[nPanels];
  for (int i = 0; i < nPanels; i++) {
    indexIncrements[i] = 1;
    basisTypes[i] = 'F';
  }

  //* Allocate dummy array for a uniform clear operation
  basisDirs = new int[nPanels];
  for (int i = 0; i < nPanels; i++) {
    basisDirs[i] = FLAT;
  }
  basisZs = new float[1];
  basisShifts = new float[1];

  //* Set flags
  isSolved = false;
  isLoaded = true;

  ifile.close();

#ifdef DEBUG_ASPECT_RATIO_VALIDITY
  if (isPanelAspectRatioValid() == false) {
    cerr << "ERROR: Panel aspect ratio of FASTCAP file is not valid." << endl;
    exit(1);
  }
#endif
}

void Caplet::loadCapletFile(const std::string filename)
{
  //* Caplet format
  //* Line 1: nWire
  //* Line 2: nShape for each wire
  //* Line 3: nTotalShape
  //* Line n: shape description lines, 12 numbers per line
  //* shapeType, indexIncrement, XL, XU, YL, YU, ZL, ZU, dir, shapeDir,
  // shapeNormalDistance, shapeShift
  //* shapeType: either F for Flat, A for Arch, and S for Side arch
  //* indexIncrement:
  //*   a basis function can consists of multiple shapes. The value here can be
  // 0 or 1.
  //*   1: the beginning of a new basis function
  //*   0: not increment, meaning the current shape is combined with the shape
  // before with value 1
  //* XL, XU, YL, YU, ZL, ZU: the range limit of each direction
  //* dir: normal direction of the rectangle that supports the shape
  //*   0 for x-dir, 1 for y-dir, 2 for z-dir
  //* shapeDir: shape varying direction
  //*   0 for x-dir, 1 for y-dir, 2 for z-dir
  //* shapeNormalDistrance: signed distance between the support rectangle and
  // the neighborhoold rectangle
  //*   positive: decaying in the positive direction
  //*   negative: decaying in the negative direction
  //* shapeShift: shift parameter that indicates the shift of starting point of

  clear();
  isInstantiable = true;

  ifstream ifile(filename.c_str());
  //* Check if open successfully
  if (!ifile) {
    cerr << "ERROR: cannot open the file: " << filename << endl;
    return;
  }

  ifile >> nWires;
  nWirePanels = new int[nWires];
  Cmat = new float[nWires * nWires];
  for (int i = 0; i < nWires; i++) {
    ifile >> nWirePanels[i];
  }
  ifile >> nPanels;
  panels = new float[nPanels][3][4];
  dirs = new int[nPanels];
  areas = new float[nPanels];
  indexIncrements = new int[nPanels];
  basisTypes = new char[nPanels];
  basisDirs = new int[nPanels];
  basisZs = new float[nPanels];
  basisShifts = new float[nPanels];

  nCoefs = 0;

  for (int i = 0; i < nPanels; i++) {
    ifile >> basisTypes[i];
    ifile >> indexIncrements[i];
    nCoefs += indexIncrements[i];

    ifile >> panels[i][X][MIN];
    ifile >> panels[i][X][MAX];
    panels[i][X][LENGTH] = panels[i][X][MAX] - panels[i][X][MIN];
    panels[i][X][CENTER] = (panels[i][X][MAX] + panels[i][X][MIN]) / 2;

    ifile >> panels[i][Y][MIN];
    ifile >> panels[i][Y][MAX];
    panels[i][Y][LENGTH] = panels[i][Y][MAX] - panels[i][Y][MIN];
    panels[i][Y][CENTER] = (panels[i][Y][MAX] + panels[i][Y][MIN]) / 2;

    ifile >> panels[i][Z][MIN];
    ifile >> panels[i][Z][MAX];
    panels[i][Z][LENGTH] = panels[i][Z][MAX] - panels[i][Z][MIN];
    panels[i][Z][CENTER] = (panels[i][Z][MAX] + panels[i][Z][MIN]) / 2;

    ifile >> dirs[i];
    switch (dirs[i]) {
      case X:
        areas[i] = panels[i][Y][LENGTH] * panels[i][Z][LENGTH];
        break;
      case Y:
        areas[i] = panels[i][Z][LENGTH] * panels[i][X][LENGTH];
        break;
      case Z:
        areas[i] = panels[i][X][LENGTH] * panels[i][Y][LENGTH];
    }

    ifile >> basisDirs[i];
    ifile >> basisZs[i];
    ifile >> basisShifts[i];
  }

  // P 	= new float[nCoefs*nCoefs];
  // rhs 	= new float[nCoefs*nWires];
  // coefs = new float[nCoefs*nWires];

  // construct nWireCoefs
  nWireCoefs = new int[nWires];
  int ind = -1;
  for (int i = 0; i < nWires; i++) {
    nWireCoefs[i] = 0;
    for (int j = 0; j < nWirePanels[i]; j++) {
      ind++;
      nWireCoefs[i] += indexIncrements[ind];
    }
  }

  isSolved = false;
  isLoaded = true;

  ifile.close();

#ifdef DEBUG_ASPECT_RATIO_VALIDITY
  if (isPanelAspectRatioValid() == false) {
    cerr << "ERROR: Panel aspect ratio of CAPLET file is not valid." << endl;
    exit(1);
  }
#endif
}

void Caplet::saveCmat(const std::string filename)
{
  if (isSolved) {
    std::ofstream ofile(filename.c_str());
    if (!ofile.is_open()) {
      cerr << "ERROR: cannot write Cmat file: " << filename << endl;
      return;
    }
    //* DOUBLE block was commented before. Double check this.
    switch (mode) {
      case FAST_GALERKIN:
        for (int i = 0; i < nWires; i++) {
          for (int j = 0; j < nWires; j++) {
            ofile << Cmat[i + nWires * j] << " ";
          }
          ofile << endl;
        }
        break;
      case DOUBLE_GALERKIN:
      case DOUBLE_COLLOCATION:
        for (int i = 0; i < nWires; i++) {
          for (int j = 0; j < nWires; j++) {
            ofile << dCmat[i + nWires * j] << " ";
          }
          ofile << endl;
        }
        break;
      default:
        break;
    }

    ofile.close();
  }
}

void Caplet::saveCoefs(const std::string filename)
{
  if (isSolved) {
    ofstream ofile(filename.c_str());
    if (!ofile) {
      cerr << "ERROR: cannot write the file: " << filename << endl;
      return;
    }
    switch (mode) {
      case FAST_GALERKIN:
        for (int i = 0; i < nCoefs; i++) {
          for (int j = 0; j < nWires; j++) {
            ofile << coefs[i + nCoefs * j] * 4 * pi * epsilon0 << " ";
          }
          ofile << endl;
        }
        break;
      case DOUBLE_GALERKIN:
      case DOUBLE_COLLOCATION:
        for (int i = 0; i < nCoefs; i++) {
          for (int j = 0; j < nWires; j++) {
            ofile << dcoefs[i + nCoefs * j] * 4 * pi * epsilon0 << " ";
          }
          ofile << endl;
        }
        break;
      default:;
    }

    ofile.close();
  }
}

void Caplet::extractC(MODE mode)
{
  mode = mode;

#ifdef CAPLET_TIMER
  fillingTime = 0;
  solvingTime = 0;
  totalTime = 0;
#endif

#ifdef CAPLET_INIT_ATAN_LOG
  caplet::atan(1);
  caplet::log(1);
#endif

  //* Init MPI
  MPI::Init();
  int rank = MPI::COMM_WORLD.Get_rank();

  //* Subdivide panels if aspect ratio is too large
  modifyPanelAspectRatio();

  if (rank == 0) {
    std::cout << "Number of conductors        : " << nWires << std::endl;
    std::cout << "Number of basis functions   : " << nCoefs << std::endl;
    std::cout << "Number of basis shapes      : " << nPanels << std::endl;
  }
  for (int iter = 0; iter < N_ITER; iter++) {
    switch (mode) {
      case DOUBLE_GALERKIN:
        extractCGalerkinDouble();
        break;
      case FAST_GALERKIN:
        extractCGalerkin();
        break;
      case DOUBLE_COLLOCATION:
        extractCCollocationDouble();
        break;
    }
  }
  isSolved = true;

  MPI::Finalize();
  if (rank != 0) {
    return;
  }

#ifdef CAPLET_TIMER
  std::cout << "Total extraction time (s)   : " << (totalTime) / N_ITER
            << std::endl;
  std::cout << "  Setup time (s)            : " << (fillingTime) / N_ITER
            << std::endl;
  std::cout << "  Solving time (s)          : " << (solvingTime) / N_ITER
            << std::endl;
#endif

  //* Print Cmat
  std::cout << "Cmat" << std::endl;
  printCmat();
  std::cout << endl;

#ifdef DEBUG_PRINT_P
  std::cout << "Store system matrix in 'pmatrix' and rhs in 'rhs'" << std::endl;
  ofstream pmatrixOut("pmatrix");
  printP(pmatrixOut);
  pmatrixOut.close();
  ofstream rhsOut("rhs");
  printRHS(rhsOut);
  rhsOut.close();
#endif
}

//***************************
//*
//* DOUBLE COOLLOCATION MODE
//*
//*
void Caplet::extractCCollocationDouble()
{
  if (isLoaded == true) {
    //* Allocate system memory for double precision
    if (MPI::COMM_WORLD.Get_rank() == 0) {
      dP = new double[nCoefs * nCoefs];
      drhs = new double[nCoefs * nWires];
      dcoefs = new double[nCoefs * nWires];
      dCmat = new double[nWires * nWires];
    } else {
      dP = new double[1];
      drhs = new double[1];
      dcoefs = new double[1];
      dCmat = new double[1];
      return;
    }

  } else {
    std::cerr << "ERROR: structure file is not yet loaded" << std::endl;
  }
#ifdef CAPLET_TIMER
  timeStart = MPI::Wtime();
#endif

  generateCollocationPMatrixDouble();
  generateRHSDouble();

#ifdef CAPLET_TIMER
  timeAfterFilling = MPI::Wtime();
#endif

  //* Solve the system
  int* ipiv = new int[nCoefs];
  int info;

  dgesv_(&nCoefs, &nWires, dP, &nCoefs, ipiv, dcoefs, &nCoefs, &info);

  delete[] ipiv;

  char transA = 't';
  char transB = 'n';
  double alpha = 4 * pi * epsilon0;
  double beta = 0.0;

  //* Use matrix-matrix product to compute Cmat from coefs
  dgemm_(&transA,
         &transB,
         &nWires,
         &nWires,
         &nCoefs,
         &alpha,
         drhs,
         &nCoefs,
         dcoefs,
         &nCoefs,
         &beta,
         dCmat,
         &nWires);

#ifdef CAPLET_TIMER
  timeAfterSolving = MPI::Wtime();
  fillingTime += timeAfterFilling - timeStart;
  solvingTime += timeAfterSolving - timeAfterFilling;
  totalTime += timeAfterSolving - timeStart;
#endif
}

void Caplet::generateCollocationPMatrixDouble()
{
  //* Set dP to zero
  double zero = 0.0;
  int inc = 1;
  int nC = nCoefs * nCoefs;
  dscal_(&nC, &zero, dP, &inc);

  //* Construct ind_vec from indexIncrements
  int* ind = new int[nPanels];
  ind[0] = 0;
  for (int i = 1; i < nPanels; i++) {
    ind[i] = ind[i - 1] + indexIncrements[i];
  }
  const int nK = nPanels * nPanels;

#ifdef CAPLET_OPENMP
#pragma omp parallel for num_threads(CAPLET_OPENMP_NUM_THREADS)
#endif
  for (int k = 0; k < nK; k++) {
    int j = k / nPanels;
    int i = k - j * nPanels;
    dP[ind[i] + nCoefs * ind[j]] += calGalerkinPEntry(i, j);
  }

  delete[] ind;
}

double Caplet::calCollocationPEntryDouble(int panel1, int panel2)
{
  //* Test point : panel1.center
  //  Integration: panel2

  //* A nDim-element array of pointers pointing to a nBit-element array
  float* coord_ptr_1[3][4];
  float* coord_ptr_2[3][4];

  //* Rotate, mirror, and call proper integrals
  switch (dirs[panel2]) {
    case X:
      rotateX2Z(coord_ptr_1, panel1);
      rotateX2Z(coord_ptr_2, panel2);
      return calColD(coord_ptr_1, coord_ptr_2);
      break;
    case Y:
      rotateY2Z(coord_ptr_1, panel1);
      rotateY2Z(coord_ptr_2, panel2);
      return calColD(coord_ptr_1, coord_ptr_2);
      break;
    case Z:
      return calColD(coord_ptr_1, coord_ptr_2);
      break;
  }
  return 0.0;  // dummy return
}

void Caplet::generateRHSDouble()
{
  int inc = 1;
  int nRHS = nWires * nCoefs;
  double alpha = 0.0f;
  dscal_(&nRHS, &alpha, dcoefs, &inc);

  int rowIndex = -1;
  int columnIndex = 0;
  int coefCounter = 0;
  for (int i = 0; i < nPanels; i++) {
    rowIndex += indexIncrements[i];
    coefCounter += indexIncrements[i];

    if (coefCounter > nWireCoefs[columnIndex]) {
      coefCounter = 1;
      columnIndex++;
    }

    float (*func)(float x, float w) = &arch;
    switch (basisTypes[i]) {
      case 'S':
        func = &side;
      case 'A':
        float b = panels[i][basisDirs[i]][LENGTH];
        float xm = (b) / 2;
        float xr = (b) / 2;
        int dir = (basisDirs[i] + 1) % 3;
        if (dir == dirs[i]) {
          dir = (dir + 1) % 3;
        }
        float w = panels[i][dir][LENGTH];

        //* Gauss quad for computing area integral of shapes
        float val = 0;
        int init_i = 0;
        if ((gauss_n % 2) == 1) {  // odd n
          init_i = 1;
          val += (*gauss::w[gauss_n])[0] * func(xm, w);
        }
        int gauss_n2 = (gauss_n + 1) / 2;
        for (int gi = init_i; gi < gauss_n2; gi++) {
          float dx = xr * (*gauss::p[gauss_n])[gi];
          val += (*gauss::w[gauss_n])[gi]
                 * (func(xm + dx, w) + func(xm - dx, w));
        }
        //* End of Gauss quad

        areas[i] = val * xr * w;
    }
    dcoefs[columnIndex * nCoefs + rowIndex] += areas[i];
  }
  dcopy_(&nRHS, dcoefs, &inc, drhs, &inc);
}

//__________________________________________________________
//*
//* DOUBLE GALERKIN MODE
//*
void Caplet::extractCGalerkinDouble()
{
  int rank = MPI::COMM_WORLD.Get_rank();

  if (isLoaded == true) {
    if (rank == 0) {
      dP = new double[nCoefs * nCoefs];
    } else {
      dP = new double[1];
    }

    drhs = new double[nCoefs * nWires];
    dcoefs = new double[nCoefs * nWires];
    dCmat = new double[nWires * nWires];

  } else {
    std::cerr << "ERROR: structure file is not yet loaded" << std::endl;
  }
#ifdef CAPLET_TIMER
  timeStart = MPI::Wtime();
  ;
#endif

#ifdef CAPLET_MPI
  generateGalerkinPMatrixDoubleMPI();
#else
  generateGalerkinPMatrixDouble();
#endif

  if (MPI::COMM_WORLD.Get_rank() != 0) {
    return;
  }
  //* End of core with non-zero rank

  generateRHSDouble();

#ifdef CAPLET_TIMER
  timeAfterFilling = MPI::Wtime();
  ;
#endif

  //* Solve the system
  int info;
  char uplo = 'u';

  //* Query optimal workspace size
  double* work = new double[1];
  int lwork = -1;
  int* ipiv = new int[nCoefs];
  dsysv_(&uplo,
         &nCoefs,
         &nWires,
         dP,
         &nCoefs,
         ipiv,
         dcoefs,
         &nCoefs,
         work,
         &lwork,
         &info);
  lwork = work[0];
  delete[] work;

  //* Solve system using optimal work length
  work = new double[lwork];
  dsysv_(&uplo,
         &nCoefs,
         &nWires,
         dP,
         &nCoefs,
         ipiv,
         dcoefs,
         &nCoefs,
         work,
         &lwork,
         &info);
  delete[] ipiv;
  delete[] work;

  //* Use matrix-matrix product to compute Cmat from coefs
  char transA = 't';
  char transB = 'n';
  double alpha = 4 * pi * epsilon0;
  double beta = 0.0;
  dgemm_(&transA,
         &transB,
         &nWires,
         &nWires,
         &nCoefs,
         &alpha,
         drhs,
         &nCoefs,
         dcoefs,
         &nCoefs,
         &beta,
         dCmat,
         &nWires);

#ifdef CAPLET_TIMER
  timeAfterSolving = MPI::Wtime();

  fillingTime += timeAfterFilling - timeStart;
  solvingTime += timeAfterSolving - timeAfterFilling;
  totalTime += timeAfterSolving - timeStart;
#endif
}

//__________________________________________________________
//*
//* FAST GALERKIN MODE
//*
void Caplet::extractCGalerkin()
{
  int rank = MPI::COMM_WORLD.Get_rank();

  if (isLoaded == true) {
    if (rank == 0) {
      P = new float[nCoefs * nCoefs];
    } else {
      P = new float[1];
    }

    rhs = new float[nCoefs * nWires];
    coefs = new float[nCoefs * nWires];
    Cmat = new float[nWires * nWires];

  } else {
    std::cerr << "ERROR: structure file is not yet loaded" << std::endl;
  }
#ifdef CAPLET_TIMER
  timeStart = MPI::Wtime();
  ;
#endif

#ifdef CAPLET_MPI
  generateGalerkinPMatrixMPI();
#endif

#ifndef CAPLET_MPI
  generateGalerkinPMatrix();
#endif

  if (MPI::COMM_WORLD.Get_rank() != 0) {
    return;
  }
  generateRHS();

#ifdef CAPLET_TIMER
  timeAfterFilling = MPI::Wtime();
  ;
#endif

  //* Solve the system
  int info;
  char uplo = 'u';

  //* Query optimal workspace size
  float* work = new float[1];
  int lwork = -1;
  int* ipiv = new int[nCoefs];
  ssysv_(&uplo,
         &nCoefs,
         &nWires,
         P,
         &nCoefs,
         ipiv,
         coefs,
         &nCoefs,
         work,
         &lwork,
         &info);
  lwork = work[0];
  delete[] work;

  //* Solve system using optimal work length
  work = new float[lwork];
  ssysv_(&uplo,
         &nCoefs,
         &nWires,
         P,
         &nCoefs,
         ipiv,
         coefs,
         &nCoefs,
         work,
         &lwork,
         &info);
  delete[] ipiv;
  delete[] work;

  //* Use matrix-matrix product to compute Cmat from coefs
  char transA = 't';
  char transB = 'n';
  float alpha = 4 * pi * epsilon0;
  float beta = 0.0;
  sgemm_(&transA,
         &transB,
         &nWires,
         &nWires,
         &nCoefs,
         &alpha,
         rhs,
         &nCoefs,
         coefs,
         &nCoefs,
         &beta,
         Cmat,
         &nWires);

#ifdef CAPLET_TIMER
  timeAfterSolving = MPI::Wtime();

  fillingTime += timeAfterFilling - timeStart;
  solvingTime += timeAfterSolving - timeAfterFilling;
  totalTime += timeAfterSolving - timeStart;
#endif
}

void Caplet::generateGalerkinPMatrix()
{
  float zero = 0.0f;
  int inc = 1;
  int nC = nCoefs * nCoefs;
  sscal_(&nC, &zero, P, &inc);

  //* Construct ind_vec from indexIncrements
  int* ind = new int[nPanels];
  ind[0] = 0;
  for (int i = 1; i < nPanels; i++) {
    ind[i] = ind[i - 1] + indexIncrements[i];
  }
  const int nK = nPanels * (nPanels + 1) / 2;

#ifdef CAPLET_OPENMP
#pragma omp parallel for num_threads(CAPLET_OPENMP_NUM_THREADS)
#endif
  for (int k = 0; k < nK; k++) {
    int j = int((sqrt(double(1 + 8 * k)) - 1) / 2);
    int i = k - j * (j + 1) / 2;

    float result = calGalerkinPEntry(i, j);

    if ((i != j) && (ind[i] == ind[j])) {
      P[ind[i] + nCoefs * ind[j]] += result * 2;
    } else {
      P[ind[i] + nCoefs * ind[j]] += result;
    }
  }

  delete[] ind;
}

void Caplet::generateGalerkinPMatrixDouble()
{
  double zero = 0.0;
  int inc = 1;
  int nC = nCoefs * nCoefs;
  dscal_(&nC, &zero, dP, &inc);

  //* Construct ind_vec from indexIncrements
  int* ind = new int[nPanels];
  ind[0] = 0;
  for (int i = 1; i < nPanels; i++) {
    ind[i] = ind[i - 1] + indexIncrements[i];
  }
  const int nK = nPanels * (nPanels + 1) / 2;

#ifdef CAPLET_OPENMP
#pragma omp parallel for num_threads(CAPLET_OPENMP_NUM_THREADS)
#endif
  for (int k = 0; k < nK; k++) {
    int j = int((sqrt(double(1 + 8 * k)) - 1) / 2);
    int i = k - j * (j + 1) / 2;

    double result = calGalerkinPEntry(i, j);

    if ((i != j) && (ind[i] == ind[j])) {
      dP[ind[i] + nCoefs * ind[j]] += result * 2;
    } else {
      dP[ind[i] + nCoefs * ind[j]] += result;
    }
  }

  delete[] ind;
}

//* Convert lower triangular matrix (column major) index k to subscript i,j
inline void ltind2sub(int k, int& i, int& j)
{
  j = int((std::sqrt(double(1 + 8 * k)) - 1) / 2);
  i = k - j * (j + 1) / 2;
}
void Caplet::generateGalerkinPMatrixMPI()
{
  int rank = MPI::COMM_WORLD.Get_rank();
  int numproc = MPI::COMM_WORLD.Get_size();

  int totalK = nPanels * (nPanels + 1) / 2;
  int nK = totalK / numproc;

  int* startK = new int[numproc];
  int* lastK = new int[numproc];
  int* startC = new int[numproc];
  int* lastC = new int[numproc];

  for (int i = 0; i < numproc; i++) {
    startK[i] = i * nK;
    lastK[i] = (i + 1) * nK - 1;
  }
  lastK[numproc - 1] = totalK - 1;

  //* Construct ind_vec from indexIncrements
  int* ind = new int[nPanels];
  ind[0] = 0;
  for (int i = 1; i < nPanels; i++) {
    ind[i] = ind[i - 1] + indexIncrements[i];
  }

  int maxNC = -1;
  for (int r = 0; r < numproc; r++) {
    int i, startj, lastj;
    ltind2sub(startK[r], i, startj);
    ltind2sub(lastK[r], i, lastj);
    startC[r] = ind[startj];
    lastC[r] = ind[lastj];

    int nC = lastC[r] - startC[r] + 1;
    if (maxNC < nC) {
      maxNC = nC;
    }
  }

  float* ptrP;
  float* tempP;
  if (rank == 0) {
    //* init P
    float zero = 0.0f;
    int inc = 1;
    int nP = nCoefs * nCoefs;
    sscal_(&nP, &zero, P, &inc);
    //* init tempP to cover the largest number of columns among all proccesses.
    int nTempP = maxNC * nCoefs;
    tempP = new float[nTempP];
    ptrP = P;
  } else {
    float zero = 0.0f;
    int inc = 1;
    int nTempP = (lastC[rank] - startC[rank] + 1) * nCoefs;
    tempP = new float[nTempP];
    sscal_(&nTempP, &zero, tempP, &inc);
    ptrP = tempP;
  }

  //* For each k index
  for (int k = startK[rank]; k <= lastK[rank]; k++) {
    //* Convert k index to i,j subscript
    int i, j;
    ltind2sub(k, i, j);
    //* Compute the P entry

    float result = calGalerkinPEntry(i, j);

#ifdef DEBUG_DETECT_NAN_INF_ENTRY
#include <cmath>
    if (isnan(result) || isinf(result)) {
      cout << "Detect nan or inf: (" << i << "," << j << ") = " << result
           << endl;
      cout << "i: " << panels[i][X][MIN] << ", " << panels[i][X][MAX] << ", "
           << panels[i][Y][MIN] << ", " << panels[i][Y][MAX] << ", "
           << panels[i][Z][MIN] << ", " << panels[i][Z][MAX]
           << ". Dir: " << dirs[i] << ", " << basisDirs[i] << endl;
      cout << "j: " << panels[j][X][MIN] << ", " << panels[j][X][MAX] << ", "
           << panels[j][Y][MIN] << ", " << panels[j][Y][MAX] << ", "
           << panels[j][Z][MIN] << ", " << panels[j][Z][MAX]
           << ". Dir: " << dirs[j] << ", " << basisDirs[j] << endl;
    }
#endif

    //* Combine rows or columns in place if consecutive panels
    //  belong to the same basis function
    if ((i != j) && (ind[i] == ind[j])) {
      ptrP[ind[i] + nCoefs * (ind[j] - startC[rank])] += result * 2;
    } else {
      ptrP[ind[i] + nCoefs * (ind[j] - startC[rank])] += result;
    }
  }

  //* Combine sub-matrices to rank0 node
  if (rank == 0) {
    MPI::Status status;
    for (int i = 1; i < numproc; i++) {
      int copylen = (lastC[i] - startC[i] + 1) * nCoefs;
      MPI::COMM_WORLD.Recv(tempP, copylen, MPI::FLOAT, i, 0, status);
      float alpha = 1.0f;
      int inc = 1;
      ptrP = P + (startC[i] * nCoefs);
      saxpy_(&copylen, &alpha, tempP, &inc, ptrP, &inc);
    }
  } else {
    int copylen = (lastC[rank] - startC[rank] + 1) * nCoefs;
    MPI::COMM_WORLD.Send(tempP, copylen, MPI::FLOAT, 0, 0);
  }

  delete[] startK;
  delete[] lastK;
  delete[] startC;
  delete[] lastC;
  delete[] tempP;
  delete[] ind;
}
void Caplet::generateGalerkinPMatrixDoubleMPI()
{
  int rank = MPI::COMM_WORLD.Get_rank();
  int numproc = MPI::COMM_WORLD.Get_size();

  int totalK = nPanels * (nPanels + 1) / 2;
  int nK = totalK / numproc;

  int* startK = new int[numproc];
  int* lastK = new int[numproc];
  int* startC = new int[numproc];
  int* lastC = new int[numproc];

  for (int i = 0; i < numproc; i++) {
    startK[i] = i * nK;
    lastK[i] = (i + 1) * nK - 1;
  }
  lastK[numproc - 1] = totalK - 1;

  //* Construct ind_vec from indexIncrements
  int* ind = new int[nPanels];
  ind[0] = 0;
  for (int i = 1; i < nPanels; i++) {
    ind[i] = ind[i - 1] + indexIncrements[i];
  }

  int maxNC = -1;
  for (int r = 0; r < numproc; r++) {
    int i, startj, lastj;
    ltind2sub(startK[r], i, startj);
    ltind2sub(lastK[r], i, lastj);
    startC[r] = ind[startj];
    lastC[r] = ind[lastj];

    int nC = lastC[r] - startC[r] + 1;
    if (maxNC < nC) {
      maxNC = nC;
    }
  }

  double* ptrP;
  double* tempP;
  if (rank == 0) {
    //* init dP
    double zero = 0.0;
    int inc = 1;
    int nP = nCoefs * nCoefs;
    dscal_(&nP, &zero, dP, &inc);
    //* init tempP to cover the largest number of columns among all proccesses.
    int nTempP = maxNC * nCoefs;
    tempP = new double[nTempP];
    ptrP = dP;
  } else {
    double zero = 0.0;
    int inc = 1;
    int nTempP = (lastC[rank] - startC[rank] + 1) * nCoefs;
    tempP = new double[nTempP];
    dscal_(&nTempP, &zero, tempP, &inc);
    ptrP = tempP;
  }

  //* For each k index
  for (int k = startK[rank]; k <= lastK[rank]; k++) {
    //* Convert k index to i,j subscript
    int i, j;
    ltind2sub(k, i, j);
    //* Compute the P entry

    double result = calGalerkinPEntry(i, j);

#ifdef DEBUG_DETECT_NAN_INF_ENTRY
#include <cmath>
    if (isnan(result) || isinf(result)) {
      cout << "Detect nan or inf: (" << i << "," << j << ") = " << result
           << endl;
      cout << "i: " << panels[i][X][MIN] << ", " << panels[i][X][MAX] << ", "
           << panels[i][Y][MIN] << ", " << panels[i][Y][MAX] << ", "
           << panels[i][Z][MIN] << ", " << panels[i][Z][MAX]
           << ". Dir: " << dirs[i] << ", " << basisDirs[i] << endl;
      cout << "j: " << panels[j][X][MIN] << ", " << panels[j][X][MAX] << ", "
           << panels[j][Y][MIN] << ", " << panels[j][Y][MAX] << ", "
           << panels[j][Z][MIN] << ", " << panels[j][Z][MAX]
           << ". Dir: " << dirs[j] << ", " << basisDirs[j] << endl;
    }
#endif

    //* Combine rows or columns in place if consecutive panels
    //  belong to the same basis function
    if ((i != j) && (ind[i] == ind[j])) {
      ptrP[ind[i] + nCoefs * (ind[j] - startC[rank])] += result * 2;
    } else {
      ptrP[ind[i] + nCoefs * (ind[j] - startC[rank])] += result;
    }
  }

  //* Combine sub-matrices to rank0 node
  if (rank == 0) {
    MPI::Status status;
    for (int i = 1; i < numproc; i++) {
      int copylen = (lastC[i] - startC[i] + 1) * nCoefs;
      MPI::COMM_WORLD.Recv(tempP, copylen, MPI::DOUBLE, i, 0, status);
      double alpha = 1.0f;
      int inc = 1;
      ptrP = dP + (startC[i] * nCoefs);
      daxpy_(&copylen, &alpha, tempP, &inc, ptrP, &inc);
    }
  } else {
    int copylen = (lastC[rank] - startC[rank] + 1) * nCoefs;
    MPI::COMM_WORLD.Send(tempP, copylen, MPI::DOUBLE, 0, 0);
  }

  delete[] startK;
  delete[] lastK;
  delete[] startC;
  delete[] lastC;
  delete[] tempP;
  delete[] ind;
}

float Caplet::calGalerkinPEntry(int panel1, int panel2)
{
  //* A nDim-element array of pointers pointing to a nBit-element array
  float* coord_ptr_1[3][4];
  float* coord_ptr_2[3][4];

  if (basisTypes[panel1] == 'F') {
    int temp = panel1;
    panel1 = panel2;
    panel2 = temp;
  }

  //* Select shapes
  Shape shape1 = selectShape(panel1);
  Shape shape2 = selectShape(panel2);

  //* Rotate, mirror, and call proper integrals
  switch (dirs[panel1]) {
    case X:
      switch (basisDirs[panel1]) {
        case Y:  // Xy
          rotateX2Z(coord_ptr_1, panel1);
          rotateX2Z(coord_ptr_2, panel2);
          switch (dirs[panel2]) {
            case X:  // Xy_X
              switch (basisDirs[panel2]) {
                case Y:  // Xy_Xy
                  return intZXZX(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Z:  // Xy_Xz
                  return intZXZY(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:  // Xy_Xf
                  return intZXZF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
            case Y:  // Xy_Y
              switch (basisDirs[panel2]) {
                case X:  // Xy_Yx
                  return intZXXZ(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Z:  // Xy_Yz
                  return intZXXY(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:  // Xy_Yf
                  return intZXXF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
            case Z:
              switch (basisDirs[panel2]) {
                case X:  // Xy_Zx
                  return intZXYZ(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Y:  // Xy_Zy
                  return intZXYX(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:  // Xy_Zf
                  return intZXYF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
          }
        case Z:  // Xz
          mirrorX2Z(coord_ptr_1, panel1);
          mirrorX2Z(coord_ptr_2, panel2);
          switch (dirs[panel2]) {
            case X:
              switch (basisDirs[panel2]) {
                case Y:  // Xz_Xy
                  return intZXZY(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Z:  // Xz_Xz
                  return intZXZX(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:  // Xz_Xf
                  return intZXZF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
            case Y:
              switch (basisDirs[panel2]) {
                case X:  // Xz_Yx
                  return intZXYZ(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Z:  // Xz_Yz
                  return intZXYX(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:  // Xz_Yf
                  return intZXYF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
              break;
            case Z:
              switch (basisDirs[panel2]) {
                case X:  // Xz_Zx
                  return intZXXZ(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Y:  // Xz_Zy
                  return intZXXY(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:  // Xz_Zf
                  return intZXXF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
          }
        case FLAT:  // Xf-?f
          switch (dirs[panel2]) {
            case X:  // Xf_Xf
              mirrorX2Z(coord_ptr_1, panel1);
              mirrorX2Z(coord_ptr_2, panel2);
              return intZFZF(coord_ptr_1, coord_ptr_2);
            case Y:  // Xf_Yf
              rotateX2Z(coord_ptr_1, panel1);
              rotateX2Z(coord_ptr_2, panel2);
              return intZFXF(coord_ptr_1, coord_ptr_2);
            case Z:  // Xf_Zf
              mirrorX2Z(coord_ptr_1, panel1);
              mirrorX2Z(coord_ptr_2, panel2);
              return intZFXF(coord_ptr_1, coord_ptr_2);
          }
      }
    case Y:
      switch (basisDirs[panel1]) {
        case X:  // Yx
          mirrorY2Z(coord_ptr_1, panel1);
          mirrorY2Z(coord_ptr_2, panel2);
          switch (dirs[panel2]) {
            case X:  // Yx_X
              switch (basisDirs[panel2]) {
                case Y:  // Yx_Xy
                  return intZXXZ(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Z:  // Yx_Xz
                  return intZXXY(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:  // Yx_Xf
                  return intZXXF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
            case Y:  // Yx_Y
              switch (basisDirs[panel2]) {
                case X:  // Yx_Yx
                  return intZXZX(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Z:  // Yx_Yz
                  return intZXZY(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:  // Yx_Yf
                  return intZXZF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
            case Z:  // Yx_Z
              switch (basisDirs[panel2]) {
                case X:  // Yx_Zx
                  return intZXYX(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Y:  // Yx_Zy
                  return intZXYZ(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:  // Yx_Zf
                  return intZXYF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
          }
        case Z:  // Yz
          rotateY2Z(coord_ptr_1, panel1);
          rotateY2Z(coord_ptr_2, panel2);
          switch (dirs[panel2]) {
            case X:  // Yz_X
              switch (basisDirs[panel2]) {
                case Y:  // Yz_Xy
                  return intZXYZ(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Z:  // Yz_Xz
                  return intZXYX(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:  // Yz_Xf
                  return intZXYF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
            case Y:  // Yz_Y
              switch (basisDirs[panel2]) {
                case X:  // Yz_Yx
                  return intZXZY(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Z:  // Yz_Yz
                  return intZXZX(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:  // Yz_Yf
                  return intZXZF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
            case Z:  // Yz_Z
              switch (basisDirs[panel2]) {
                case X:  // Yz_Zx
                  return intZXXY(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Y:  // Yz_Zy
                  return intZXXZ(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:  // Yz_Zf
                  return intZXXF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
          }
        case FLAT:  // Yf_?f
          switch (dirs[panel2]) {
            case X:
              mirrorY2Z(coord_ptr_1, panel1);
              mirrorY2Z(coord_ptr_2, panel2);
              return intZFXF(coord_ptr_1, coord_ptr_2);
            case Y:
              mirrorY2Z(coord_ptr_1, panel1);
              mirrorY2Z(coord_ptr_2, panel2);
              return intZFZF(coord_ptr_1, coord_ptr_2);
            case Z:
              rotateY2Z(coord_ptr_1, panel1);
              rotateY2Z(coord_ptr_2, panel2);
              return intZFXF(coord_ptr_1, coord_ptr_2);
          }
      }
    case Z:  // Z
      switch (basisDirs[panel1]) {
        case X:  // Zx
          rotateZ2Z(coord_ptr_1, panel1);
          rotateZ2Z(coord_ptr_2, panel2);
          switch (dirs[panel2]) {
            case X:  // Zx_X
              switch (basisDirs[panel2]) {
                case Y:  // Zx_Xy
                  return intZXXY(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Z:  // Zx_Xz
                  return intZXXZ(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:  // Zx_Xf
                  return intZXXF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
            case Y:  // Zx_Y
              switch (basisDirs[panel2]) {
                case X:  // Zx_Yx
                  return intZXYX(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Z:  // Zx_Yz
                  return intZXYZ(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:
                  return intZXYF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
            case Z:  // Zx_Z
              switch (basisDirs[panel2]) {
                case X:  // Zx_Zx
                  return intZXZX(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Y:  // Zx_Zy
                  return intZXZY(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:
                  return intZXZF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
          }
        case Y:  // Zy
          mirrorY2X(coord_ptr_1, panel1);
          mirrorY2X(coord_ptr_2, panel2);
          switch (dirs[panel2]) {
            case X:  // Zy_X
              switch (basisDirs[panel2]) {
                case Y:  // Zy_Xy
                  return intZXYX(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Z:  // Zy_Xz
                  return intZXYZ(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:
                  return intZXYF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
            case Y:  // Zy_Y
              switch (basisDirs[panel2]) {
                case X:  // Zy_Yx
                  return intZXXY(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Z:  // Zy_Yz
                  return intZXXZ(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:
                  return intZXXF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
            case Z:  // Zy_Z
              switch (basisDirs[panel2]) {
                case X:  // Zy_Zx
                  return intZXZY(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case Y:  // Zy_Zy
                  return intZXZX(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2,
                                 basisZs[panel2],
                                 basisShifts[panel2],
                                 shape2);
                case FLAT:
                  return intZXZF(coord_ptr_1,
                                 basisZs[panel1],
                                 basisShifts[panel1],
                                 shape1,
                                 coord_ptr_2);
              }
          }
        case FLAT:  // Zf_?f
          switch (dirs[panel2]) {
            case X:
              rotateZ2Z(coord_ptr_1, panel1);
              rotateZ2Z(coord_ptr_2, panel2);
              return intZFXF(coord_ptr_1, coord_ptr_2);
            case Y:
              mirrorY2X(coord_ptr_1, panel1);
              mirrorY2X(coord_ptr_2, panel2);
              return intZFXF(coord_ptr_1, coord_ptr_2);
            case Z:
              rotateZ2Z(coord_ptr_1, panel1);
              rotateZ2Z(coord_ptr_2, panel2);
              return intZFZF(coord_ptr_1, coord_ptr_2);
          }
      }
  }

  return 0.0;  // dummy return
}

void Caplet::generateRHS()
{
  int inc = 1;
  int nRHS = nWires * nCoefs;
  float alpha = 0.0f;
  sscal_(&nRHS, &alpha, coefs, &inc);

  int rowIndex = -1;
  int columnIndex = 0;
  int coefCounter = 0;
  for (int i = 0; i < nPanels; i++) {
    rowIndex += indexIncrements[i];
    coefCounter += indexIncrements[i];

    if (coefCounter > nWireCoefs[columnIndex]) {
      coefCounter = 1;
      columnIndex++;
    }

    float (*func)(float x, float w) = &arch;
    switch (basisTypes[i]) {
      case 'S':
        func = &side;
      case 'A':
        float b = panels[i][basisDirs[i]][LENGTH];
        float xm = (b) / 2;
        float xr = (b) / 2;
        int dir = (basisDirs[i] + 1) % 3;
        if (dir == dirs[i]) {
          dir = (dir + 1) % 3;
        }
        float w = panels[i][dir][LENGTH];

        //* Gauss quad for computing area integral of shapes
        float val = 0;
        int init_i = 0;
        if ((gauss_n % 2) == 1) {  // odd n
          init_i = 1;
          val += (*gauss::w[gauss_n])[0] * func(xm, w);
        }
        int gauss_n2 = (gauss_n + 1) / 2;
        for (int gi = init_i; gi < gauss_n2; gi++) {
          float dx = xr * (*gauss::p[gauss_n])[gi];
          val += (*gauss::w[gauss_n])[gi]
                 * (func(xm + dx, w) + func(xm - dx, w));
        }
        //* End of Gauss quad

        areas[i] = val * xr * w;
    }

    coefs[columnIndex * nCoefs + rowIndex] += areas[i];
  }
  scopy_(&nRHS, coefs, &inc, rhs, &inc);
}

void Caplet::modifyPanelAspectRatio()
{
  //* Aspect ratio of all panels cannot exceed MAX_ASPECT_RAIO
  if (isLoaded == false) {
    std::cerr << "ERROR: not yet load any structure " << std::endl;
    std::cerr << "       Unable to check structure validity" << std::endl;

    std::exit(1);
  }
  const float maxRatio = kMaxAspectRatio;

  std::list<float> panelsTemp[3][2];
  std::list<int> dirsTemp;
  std::list<float> areasTemp;
  std::list<int> indexIncrementsTemp;
  std::list<char> basisTypesTemp;
  std::list<int> basisDirsTemp;
  std::list<float> basisZsTemp;
  std::list<float> basisShiftsTemp;

  int currentWireIndex = 0;
  int currentWirePanelCount = 0;
  int nCurrentWirePanels = nWirePanels[0];
  int nFinalPanels = nPanels;

  for (int i = 0; i < nPanels; i++) {
    currentWirePanelCount++;
    if (currentWirePanelCount > nWirePanels[currentWireIndex]) {
      currentWirePanelCount = 1;
      nWirePanels[currentWireIndex] = nCurrentWirePanels;
      currentWireIndex++;
      nCurrentWirePanels = nWirePanels[currentWireIndex];
    }
    float ba = panels[i][(dirs[i] + 1) % 3][LENGTH]
               / panels[i][(dirs[i] + 2) % 3][LENGTH];
    if ((ba > maxRatio) || (1 / ba > maxRatio)) {
      //* need to split to reduce aspact ratio of a panel
      int nSplit = 2;
      while ((ba / nSplit > maxRatio) || (1 / ba / nSplit > maxRatio)) {
        nSplit++;
      }
      int splitDir = ((ba > 1) ? ((dirs[i] + 1) % 3) : ((dirs[i] + 2) % 3));
      float currentCoord = panels[i][splitDir][0];
      float sublength = panels[i][splitDir][LENGTH] / nSplit;

      for (int j = 0; j < nSplit; j++) {
        panelsTemp[splitDir][0].push_back(currentCoord);
        currentCoord += sublength;
        panelsTemp[splitDir][1].push_back(currentCoord);

        panelsTemp[(splitDir + 1) % 3][0].push_back(
            panels[i][(splitDir + 1) % 3][0]);
        panelsTemp[(splitDir + 1) % 3][1].push_back(
            panels[i][(splitDir + 1) % 3][1]);
        panelsTemp[(splitDir + 2) % 3][0].push_back(
            panels[i][(splitDir + 2) % 3][0]);
        panelsTemp[(splitDir + 2) % 3][1].push_back(
            panels[i][(splitDir + 2) % 3][1]);

        dirsTemp.push_back(dirs[i]);
        areasTemp.push_back(areas[i] / nSplit);
        if (j == 0) {
          //* The first sub-panel after split. Needs index increment.
          indexIncrementsTemp.push_back(indexIncrements[i]);
        } else {
          //* The rest of sub-panel after split. No index increment needed.
          indexIncrementsTemp.push_back(0);
        }
        basisTypesTemp.push_back(basisTypes[i]);
        basisDirsTemp.push_back(basisDirs[i]);
        basisZsTemp.push_back(basisZs[i]);

        if (basisDirs[i] == splitDir) {
          if (basisZs[i] > 0) {  // decaying in the positive direction
            basisShiftsTemp.push_back(basisShifts[i] + j * sublength);
          } else {  // decaying in the negative direction
            basisShiftsTemp.push_back(basisShifts[i]
                                      + (nSplit - 1 - j) * sublength);
          }
        } else {
          basisShiftsTemp.push_back(basisShifts[i]);
        }

        if (j != 0) {
          nFinalPanels++;
          nCurrentWirePanels++;
        }
      }
    } else {
      //* no need to split. simply push back
      panelsTemp[X][0].push_back(panels[i][X][0]);
      panelsTemp[X][1].push_back(panels[i][X][1]);
      panelsTemp[Y][0].push_back(panels[i][Y][0]);
      panelsTemp[Y][1].push_back(panels[i][Y][1]);
      panelsTemp[Z][0].push_back(panels[i][Z][0]);
      panelsTemp[Z][1].push_back(panels[i][Z][1]);

      dirsTemp.push_back(dirs[i]);
      areasTemp.push_back(areas[i]);
      indexIncrementsTemp.push_back(indexIncrements[i]);
      basisTypesTemp.push_back(basisTypes[i]);
      basisDirsTemp.push_back(basisDirs[i]);
      basisZsTemp.push_back(basisZs[i]);
      basisShiftsTemp.push_back(basisShifts[i]);
    }
  }
  nPanels = nFinalPanels;
  nWirePanels[currentWireIndex] = nCurrentWirePanels;

  //* Rebuild panel information based on new panel list with split panels
  delete[] panels;
  delete[] dirs;
  delete[] indexIncrements;
  delete[] basisTypes;
  delete[] basisDirs;
  delete[] basisZs;
  delete[] basisShifts;

  panels = new float[nPanels][3][4];
  dirs = new int[nPanels];
  areas = new float[nPanels];
  indexIncrements = new int[nPanels];
  basisTypes = new char[nPanels];
  basisDirs = new int[nPanels];
  basisZs = new float[nPanels];
  basisShifts = new float[nPanels];

  for (int i = 0; i < nPanels; i++) {
    panels[i][X][0] = panelsTemp[X][0].front();
    panelsTemp[X][0].pop_front();
    panels[i][X][1] = panelsTemp[X][1].front();
    panelsTemp[X][1].pop_front();
    panels[i][Y][0] = panelsTemp[Y][0].front();
    panelsTemp[Y][0].pop_front();
    panels[i][Y][1] = panelsTemp[Y][1].front();
    panelsTemp[Y][1].pop_front();
    panels[i][Z][0] = panelsTemp[Z][0].front();
    panelsTemp[Z][0].pop_front();
    panels[i][Z][1] = panelsTemp[Z][1].front();
    panelsTemp[Z][1].pop_front();

    panels[i][X][LENGTH] = panels[i][X][1] - panels[i][X][0];
    panels[i][Y][LENGTH] = panels[i][Y][1] - panels[i][Y][0];
    panels[i][Z][LENGTH] = panels[i][Z][1] - panels[i][Z][0];

    panels[i][X][CENTER] = (panels[i][X][1] + panels[i][X][0]) / 2;
    panels[i][Y][CENTER] = (panels[i][Y][1] + panels[i][Y][0]) / 2;
    panels[i][Z][CENTER] = (panels[i][Z][1] + panels[i][Z][0]) / 2;

    dirs[i] = dirsTemp.front();
    dirsTemp.pop_front();
    areas[i] = areasTemp.front();
    areasTemp.pop_front();

    indexIncrements[i] = indexIncrementsTemp.front();
    indexIncrementsTemp.pop_front();

    basisTypes[i] = basisTypesTemp.front();
    basisTypesTemp.pop_front();
    basisDirs[i] = basisDirsTemp.front();
    basisDirsTemp.pop_front();
    basisZs[i] = basisZsTemp.front();
    basisZsTemp.pop_front();
    basisShifts[i] = basisShiftsTemp.front();
    basisShiftsTemp.pop_front();
  }
}

Shape Caplet::selectShape(int panel)
{
  switch (basisTypes[panel]) {
    case 'A':
      return &arch;
      break;
    case 'S':
      return &side;
      break;
    default:
      return 0;
  }
}

int Caplet::getNPanels() const
{
  return nPanels;
}

int Caplet::getNCoefs() const
{
  return nCoefs;
}

int Caplet::getSizeP() const
{
  return nCoefs * nCoefs;
}

int Caplet::getSizeCoefs() const
{
  return nCoefs * nWires;
}

int Caplet::getSizeCmat() const
{
  return nWires * nWires;
}

const float* const Caplet::getCmat() const
{
  return Cmat;
}

float Caplet::compareCmatError(const float* const cmatRef,
                               ERROR_REF option) const
{
  using std::abs;

  float* const cmatError = new float[getSizeCmat()];
  int n = nWires;
  float maxError = 0.0f;
  int maxErrorI = 0;
  int maxErrorJ = 0;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      float errRef = 0.0f;
      switch (option) {
        case DIAGONAL:
          errRef = (cmatRef[i + n * i]);
          break;
        case SELF:
          errRef = (cmatRef[i + n * j]);
      }
      float thisError
          = abs((Cmat[i + n * j] - cmatRef[i + n * j]) / errRef * 100);
      cmatError[i + n * j] = thisError;

      if (maxError < thisError) {
        maxError = thisError;
        maxErrorI = i;
        maxErrorJ = j;
      }
    }
  }
  print_matrix(cmatError, n, n, "Cmat Error (%)", 'a');

  std::cout << "Max error = " << std::fixed << maxError << "% at (" << maxErrorI
            << "," << maxErrorJ << ")" << std::endl;

  switch (option) {
    case DIAGONAL:
      std::cout << "    (w.r.t. the row diagonal)" << std::endl;
      break;
    case SELF:
      std::cout << "    (w.r.t. itself)" << std::endl;
  }
  std::cout << std::endl;

  delete[] cmatError;
  return maxError;
}

float Caplet::compareCmatError(const Caplet* const caplet,
                               ERROR_REF option) const
{
  const float* const cmatRef = caplet->getCmat();
  return compareCmatError(cmatRef, option);
}

float Caplet::compareCmatError(const std::string filename,
                               ERROR_REF option) const
{
  const int n = nWires;

  float* cmatRef = new float[n * n];

  std::ifstream ifile(filename.c_str());
  if (!ifile) {
    std::cerr << "ERROR: cannot open the file: " << filename << std::endl;
    std::exit(1);
  }

  std::string lineTemp;
  for (int i = 0; i < n; i++) {
    getline(ifile, lineTemp);
    std::stringstream stringTokenizer(lineTemp);
    for (int j = 0; j < n; j++) {
      stringTokenizer >> cmatRef[i + n * j];
    }
  }

  float maxError = compareCmatError(cmatRef, option);
  delete[] cmatRef;

  return maxError;
}

//__________________________________________________________
//*
//* COORDINATE FUNCTIONS
//*
void Caplet::rotateX2Z(float* coord_ptr[3][4], int panelNo)
{
  /* the argument float* coord_ptr[3][4] means:
   * 1. It is a 3-element array.
   * 2. Each element is a pointer.
   * 3. The pointer points to a 4-element array.
   *
   * *coord_ptr[1] means:
   * 1. Look at the 2nd pointer of the array.
   * 2. We care about the address where the pointer points
   * */
  *coord_ptr[X] = panels[panelNo][Y];
  *coord_ptr[Y] = panels[panelNo][Z];
  *coord_ptr[Z] = panels[panelNo][X];
}

void Caplet::rotateY2Z(float* coord_ptr[3][4], int panelNo)
{
  *coord_ptr[X] = panels[panelNo][Z];
  *coord_ptr[Y] = panels[panelNo][X];
  *coord_ptr[Z] = panels[panelNo][Y];
}

void Caplet::rotateZ2Z(float* coord_ptr[3][4], int panelNo)
{
  *coord_ptr[X] = panels[panelNo][X];
  *coord_ptr[Y] = panels[panelNo][Y];
  *coord_ptr[Z] = panels[panelNo][Z];
}

void Caplet::mirrorY2X(float* coord_ptr[3][4], int panelNo)
{
  *coord_ptr[X] = panels[panelNo][Y];
  *coord_ptr[Y] = panels[panelNo][X];
  *coord_ptr[Z] = panels[panelNo][Z];
}

void Caplet::mirrorY2Z(float* coord_ptr[3][4], int panelNo)
{
  *coord_ptr[X] = panels[panelNo][X];
  *coord_ptr[Y] = panels[panelNo][Z];
  *coord_ptr[Z] = panels[panelNo][Y];
}

void Caplet::mirrorX2Z(float* coord_ptr[3][4], int panelNo)
{
  *coord_ptr[X] = panels[panelNo][Z];
  *coord_ptr[Y] = panels[panelNo][Y];
  *coord_ptr[Z] = panels[panelNo][X];
}

//__________________________________________________________
//*
//* PRINT UTILITIES
//*
void Caplet::printP(std::ostream& out)
{
  if (isLoaded && isSolved) {
    switch (mode) {
      case FAST_GALERKIN:
        print_matrix(P, nCoefs, nCoefs, "", 'u', out);
        break;
      case DOUBLE_GALERKIN:
      case DOUBLE_COLLOCATION:
        print_matrix(dP, nCoefs, nCoefs, "", 'u', out);
        cout << "here" << endl;
        cout << dP[0] << ", " << dP[1] << endl;
        out << "test" << endl;
        break;
      default:;
    }
  }
}

void Caplet::printCoefs()
{
  if (isLoaded && isSolved) {
    switch (mode) {
      case FAST_GALERKIN:
        print_matrix(coefs, nCoefs, nWires, "");
        break;
      case DOUBLE_GALERKIN:
      case DOUBLE_COLLOCATION:
        print_matrix(dcoefs, nCoefs, nWires, "");
        break;
      default:;
    }
  }
}

void Caplet::printRHS(std::ostream& out)
{
  if (isLoaded && isSolved) {
    switch (mode) {
      case FAST_GALERKIN:
        print_matrix(rhs, nCoefs, nWires, "", 'a', out);
        break;
      case DOUBLE_GALERKIN:
      case DOUBLE_COLLOCATION:
        print_matrix(drhs, nCoefs, nWires, "", 'a', out);
        break;
      default:;
    }
  }
}

void Caplet::printCmat()
{
  if (isLoaded && isSolved) {
    switch (mode) {
      case FAST_GALERKIN:
        print_matrix(Cmat, nWires, nWires, "");
        break;
      case DOUBLE_GALERKIN:
      case DOUBLE_COLLOCATION:
        print_matrix(dCmat, nWires, nWires, "");
        break;
      default:;
    }
  }
}

void Caplet::printPanel(int index)
{
  if (isLoaded == false) {
    return;
  }
  if (index >= nPanels) {
    return;
  }
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 4; ++j) {
      cout << setw(12) << panels[index][i][j];
    }
    cout << endl;
  }
}

//__________________________________________________________
//*
//* DEBUG UTILITIES
//*
bool Caplet::isPanelAspectRatioValid()
{
  if (isLoaded == false) {
    std::cerr << "ERROR: not yet load any structure " << std::endl;
    std::cerr << "       Unable to check structure validity" << std::endl;

    return false;
  }
  const float maxRatio = 51;

  for (int i = 0; i < nPanels; i++) {
    float ratio = panels[i][(dirs[i] + 1) % 3][LENGTH]
                  / panels[i][(dirs[i] + 2) % 3][LENGTH];

    if (ratio > maxRatio || 1 / ratio > maxRatio) {
      std::cerr << "ERROR: panel ratio check fails" << std::endl;
      std::cerr << "    Panel: " << i << std::endl;
      std::cerr << "    x:" << panels[i][X][0] << ", " << panels[i][X][1]
                << ", " << panels[i][X][2] << ", " << panels[i][X][3]
                << std::endl;
      std::cerr << "    y:" << panels[i][Y][0] << ", " << panels[i][Y][1]
                << ", " << panels[i][Y][2] << ", " << panels[i][Y][3]
                << std::endl;
      std::cerr << "    z:" << panels[i][Z][0] << ", " << panels[i][Z][1]
                << ", " << panels[i][Z][2] << ", " << panels[i][Z][3]
                << std::endl;

      return false;
    }
  }
  return true;
}

}  // end of namespace caplet
