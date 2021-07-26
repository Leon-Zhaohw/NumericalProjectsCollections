/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <VECTOR.h>

//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    cout << " USAGE: " << argv[0] << " *.vector <optional second vector>" << endl;
    return 0;
  }

  VECTOR firstVector(argv[1]);

  VECTOR::printVertical = false;
  cout << "vector = "<< firstVector << endl;

  if (argc > 2)
  {
    VECTOR secondVector(argv[2]);
    cout << "vector2 = " << secondVector << endl;

    VECTOR diff;
    diff = firstVector - secondVector;
    diff.fabs();

    cout << " diff = " << diff << endl;
    cout << " diff sum = " << diff.sum() << endl;
    cout << " diff mean = " << diff.sum() / diff.size() << endl;
  }

  return 0;
}
