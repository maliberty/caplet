/*
Created : Jan 29, 2013
Modified: Feb 13, 2013
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

#ifndef CAPLET_CONST_H
#define CAPLET_CONST_H

namespace caplet {

enum DIR
{
  X,
  Y,
  Z,
  FLAT
};
enum
{
  MIN,
  MAX,
  LENGTH,
  CENTER
};

constexpr int nDim = 3;
constexpr int nBit = 4;
constexpr float kMaxAspectRatio = 50;

using Shape = float (*)(float, float);

}  // namespace caplet

#endif  // CAPLET_CONST_H
