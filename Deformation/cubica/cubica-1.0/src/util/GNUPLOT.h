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
//////////////////////////////////////////////////////////////////////////////////////////
// my gnuplot interface class
//////////////////////////////////////////////////////////////////////////////////////////
#ifndef GNUPLOT_H
#define GNUPLOT_H

#include <vector>
#include <string>
#include <cstdio>
#include <VECTOR.h>
using namespace std;

class GNUPLOT {

public:
  GNUPLOT() : _totalPlots(0) {
  };

  ~GNUPLOT() {
    system("rm gnuplot.data.*");
    system("rm plotscript");
  };

  //////////////////////////////////////////////////////////////////////
  // add a plot
  //////////////////////////////////////////////////////////////////////
  void addPlot(VECTOR& points, string label) {
    cout << " Dumping plot with label " << label.c_str() << endl;
    char buffer[256];
    sprintf(buffer, "%05d", _totalPlots);
    string filename = string("gnuplot.data.") + string(buffer);
    FILE* file = fopen(filename.c_str(), "w");

    if (file == NULL)
    {
      cout << " Failed to open file " << filename << "!" << endl;
      return;
    }

    for (int x = 0; x < points.size(); x++)
      fprintf(file, "%d %f\n", x, points[x]);

    fclose(file);

    _labels.push_back(label);
    _totalPlots++;
  };

  //////////////////////////////////////////////////////////////////////
  // add a plot
  //////////////////////////////////////////////////////////////////////
  void addPlot(float* xPoints, float* yPoints, int size, string label) {
    char buffer[256];
    sprintf(buffer, "%05d", _totalPlots);
    string filename = string("data") + string(buffer);
    FILE* file = fopen(filename.c_str(), "w");
    for (int x = 0; x < size; x++)
      fprintf(file, "%f %f\n", xPoints[x], yPoints[x]);
    fclose(file);
    
    /*
    ofstream fout;
    fout.open(filename.c_str());
    for (int x = 0; x < size; x++)
      fout << xPoints[x] << " " << yPoints[x] << endl;
    fout.close();
    */
    _labels.push_back(label);
    _totalPlots++;
  };

  //////////////////////////////////////////////////////////////////////
  // add a plot - autolabelled version
  //////////////////////////////////////////////////////////////////////
  void addPlot(float* xPoints, float* yPoints, int size) {
    char buffer[256];
    sprintf(buffer, "%05d", _totalPlots);
    string label = string("plot ") + string(buffer);
    addPlot(xPoints, yPoints, size, label);
  };

  //////////////////////////////////////////////////////////////////////
  // output the plots
  //////////////////////////////////////////////////////////////////////
  void outputPlots(const char* filename, const char* title = "")
  {
    cout << " Outputting Gnuplot graph " << title << " to " << filename << endl;
    string plotscript = string("plotscript");
    FILE* file = fopen(plotscript.c_str(), "w");
    fprintf(file, "set terminal postscript solid color\n");
    fprintf(file, "set output '%s'\n", filename);
    fprintf(file, "set title '%s'\n", title);

    // output the actual plots
    for (int x = 0; x < _totalPlots; x++)
    {
      char buffer[256];
      sprintf(buffer, "%05d", x);
      string dataname = string("gnuplot.data.") + string(buffer);
      if (x == 0)
        fprintf(file, "plot ");
      //fout << "'" << dataname.c_str() << "' title '" << _labels[x].c_str() << "' with linespoints";
      fprintf(file, "'%s' title '%s' with linespoints", dataname.c_str(), _labels[x].c_str());
      if (x != _totalPlots - 1)
        fprintf(file, ", \\");
      fprintf(file, "\n");
    }
    fclose(file);

    // plot the actual data
    system("gnuplot plotscript");
  };

  //////////////////////////////////////////////////////////////////////
  // output the plot and view it
  //////////////////////////////////////////////////////////////////////
  void viewPlots(const char* filename)
  {
    outputPlots(filename);
    string command = string("gsview32 ") + filename;
    system(command.c_str());
  }

  //////////////////////////////////////////////////////////////////////
  // output log plots
  //////////////////////////////////////////////////////////////////////
  void outputLogPlots(const char* filename, const char* title = "")
  {
    cout << " Outputting Gnuplot graph " << title << " to " << filename << endl;
    string plotscript = string("plotscript");
    FILE* file = fopen(plotscript.c_str(), "w");
    fprintf(file, "set terminal postscript solid color\n");
    fprintf(file, "set output '%s'\n", filename);
    fprintf(file, "set title '%s'\n", title);
    fprintf(file, "set logscale xy\n");

    // output the actual plots
    for (int x = 0; x < _totalPlots; x++)
    {
      char buffer[256];
      sprintf(buffer, "%05d", x);
      string dataname = string("gnuplot.data.") + string(buffer);
      if (x == 0)
        fprintf(file, "plot ");
      fprintf(file, "'%s' title '%s' with linespoints", dataname.c_str(), _labels[x].c_str());
      if (x != _totalPlots - 1)
        fprintf(file, ", \\");
      fprintf(file, "\n");
    }
    fclose(file);

    // plot the actual data
    system("gnuplot plotscript");
  }

  //////////////////////////////////////////////////////////////////////
  // output log plots
  //////////////////////////////////////////////////////////////////////
  void outputSemiLogYPlots(const char* filename, const char* title = "")
  {
    cout << " Outputting Gnuplot graph " << title << " to " << filename << endl;
    string plotscript = string("plotscript");
    FILE* file = fopen(plotscript.c_str(), "w");
    fprintf(file, "set terminal postscript solid color\n");
    fprintf(file, "set output '%s'\n", filename);
    fprintf(file, "set title '%s'\n", title);
    fprintf(file, "set logscale y\n");

    // output the actual plots
    for (int x = 0; x < _totalPlots; x++)
    {
      char buffer[256];
      sprintf(buffer, "%05d", x);
      string dataname = string("gnuplot.data.") + string(buffer);
      if (x == 0)
        fprintf(file, "plot ");
      fprintf(file, "'%s' title '%s' with linespoints", dataname.c_str(), _labels[x].c_str());
      if (x != _totalPlots - 1)
        fprintf(file, ", \\");
      fprintf(file, "\n");
    }
    fclose(file);

    // plot the actual data
    system("gnuplot plotscript");
  }

  //////////////////////////////////////////////////////////////////////
  // output the log plot and view it
  //////////////////////////////////////////////////////////////////////
  void viewLogPlots(const char* filename)
  {
    outputLogPlots(filename);
    string command = string("gsview32 ") + filename;
    system(command.c_str());
  }

  //////////////////////////////////////////////////////////////////////
  // member variables
  //////////////////////////////////////////////////////////////////////
  vector<string> _labels;
  int _totalPlots;
};

#endif
