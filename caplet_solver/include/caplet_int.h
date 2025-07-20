/*
Created : Jul 31, 2010
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

#ifndef CAPLET_INT_H_
#define CAPLET_INT_H_

#include "caplet_const.h"
#include "caplet_debug.h"

#ifdef DEBUG_SHAPE_BOUNDARY_CHECK
#include <iostream>
#endif

namespace caplet {

//* FAST GALERKIN MODE

//* Flat-Flat integrals
//* - integral 1
float intZFZF(float* coord1[nDim][nBit], float* coord2[nDim][nBit]);
//* - integral 2
float intZFXF(float* coord1[nDim][nBit], float* coord2[nDim][nBit]);

//* Linear-Flat integrals
//* - integral 3
float intZXZF(float* coord1[nDim][nBit],
              float bz,
              float bsh,
              float (*shape)(float, float),
              float* coord2[nDim][nBit]);
//* - integral 4
float intZXXF(float* coord1[nDim][nBit],
              float bz,
              float bsh,
              float (*shape)(float, float),
              float* coord2[nDim][nBit]);
//* - integral 5
float intZXYF(float* coord1[nDim][nBit],
              float bz,
              float bsh,
              float (*shape)(float, float),
              float* coord2[nDim][nBit]);

//* Linear-Linear integrals
//* - integral 6
float intZXZX(float* coord1[nDim][nBit],
              float bz1,
              float bsh1,
              float (*shape1)(float, float),
              float* coord2[nDim][nBit],
              float bz2,
              float bsh2,
              float (*shape2)(float, float));
//* - integral 7
float intZXYX(float* coord1[nDim][nBit],
              float bz1,
              float bsh1,
              float (*shape1)(float, float),
              float* coord2[nDim][nBit],
              float bz2,
              float bsh2,
              float (*shape2)(float, float));
//* - integral 8
float intZXZY(float* coord1[nDim][nBit],
              float bz1,
              float bsh1,
              float (*shape1)(float, float),
              float* coord2[nDim][nBit],
              float bz2,
              float bsh2,
              float (*shape2)(float, float));
//* - integral 9
float intZXXZ(float* coord1[nDim][nBit],
              float bz1,
              float bsh1,
              float (*shape1)(float, float),
              float* coord2[nDim][nBit],
              float bz2,
              float bsh2,
              float (*shape2)(float, float));
//* - integral 10
float intZXYZ(float* coord1[nDim][nBit],
              float bz1,
              float bsh1,
              float (*shape1)(float, float),
              float* coord2[nDim][nBit],
              float bz2,
              float bsh2,
              float (*shape2)(float, float));
//* - integral 11
float intZXXY(float* coord1[nDim][nBit],
              float bz1,
              float bsh1,
              float (*shape1)(float, float),
              float* coord2[nDim][nBit],
              float bz2,
              float bsh2,
              float (*shape2)(float, float));

//***********************************************************
//*
//* Double collocation mode
//*
//*
double calColD(float* p1[nDim][nBit], float* p2[nDim][nBit]);

}  // namespace caplet

#endif /* CAPLET_INT_H_ */
