#ifndef TEXTUREOBSTACLE_HPP
#define TEXTUREOBSTACLE_HPP

#include "Obstacle.hpp"
#include <string>

class TextureObstacle: public Obstacle {
protected:
  SDL_Surface *height_field;
  Grid grid;
  int n_rows, n_cols; 
  FLOAT cell_size; // meters
  std::string file_texture;
  
  void setBoundaries(uint w);
  void setEquivalentSources(uint w);

public:
  TextureObstacle();
  TextureObstacle(std::string file, uint nr, uint nc, FLOAT cs = 0.1);
  ~TextureObstacle();

  FLOAT getGrid(int i, int j) const;
  void setCellSize(FLOAT cs);

  std::ofstream &exportMitsuba(std::ofstream &file) const;
};

#endif
