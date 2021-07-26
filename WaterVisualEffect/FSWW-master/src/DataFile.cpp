/* 
 * File: DataFile.cpp
 *
 * Copyright (C) 2019  Camille Schreck
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#include "DataFile.hpp"

DataFile *DataFile::DATA(new DataFile());

DataFile::DataFile() {
  open = false;
  write_file = false;
}

void DataFile::init() {
  if (open) {
    close();
  }
  open = true;
  if (write_file) {
    file.open("data/file.dat");
  }
}

void DataFile::close() {
  if (open) {
    open = false;
    if (write_file) {
      file.close();
    }
  }
}
