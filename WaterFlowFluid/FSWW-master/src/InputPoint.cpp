/* 
 * File: InputPoint.cpp
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

#include "InputPoint.hpp"
#include <fstream>
#include <iostream>

InputPoint::InputPoint() {
  pos = VEC2(0, 0);
  window_size = 0;
  dt = 0;
  nb_frequencies = 0;
 frequency_step = 0;

  spectrum_re = NULL;
  spectrum_im = NULL;

  t = 0;
  name = "test";
}

InputPoint::InputPoint(uint nb_samples, FLOAT time_step) {
  pos = VEC2(0, 0);
  window_size = nb_samples;
  dt = time_step;
  sample_rate = 1.0/dt;
  frequency_step = 1/(dt*window_size);
  nb_frequencies = window_size;
  
  for (uint i  = 0; i < window_size; ++i) {
    samples.push_back(0.0);
  }
  spectrum_re = new FLOAT[nb_frequencies];
  spectrum_im = new FLOAT[nb_frequencies];

  t = 0;
  name = "test";
}

InputPoint::~InputPoint() {
   plotSpectrum();
   plotSpectrogram();
   plotSamples();
 
}

void InputPoint::setPos(VEC2 p) {
  pos = p;
}

void InputPoint::setPos(FLOAT x, FLOAT y) {
  pos = VEC2(x, y);
}

VEC2 InputPoint::getPos() const {
  return pos;
}

void InputPoint::update(FLOAT next_sample) {
   samples.push_back(next_sample);
   samples.pop_front();
   if (t == window_size / 2) {
     computeSpectrum();
     t = 0;
   }
  ++t;
}

void InputPoint::computeSpectrum() {
  // FLOAT *in_re = new FLOAT[window_size];
  // FLOAT *in_im = new FLOAT[window_size];

  // std::list<FLOAT>::iterator it;
  // uint i = 0;
  // for (it = samples.begin(); it != samples.end(); ++it, ++i) {
  //   in_re[i] = (*it);
  //   in_im[i] = 0.0;
  // }

  // // TODO
  // // FFT(window_size, false, in_re, in_im, spectrum_re, spectrum_im);
  
  // std::vector<FLOAT> power_spectrum(nb_frequencies);
  // FLOAT epsilon_db = 1e-18;
  // for (uint i = 0; i < nb_frequencies; ++i) {
  //    FLOAT e = powf(spectrum_re[i]/window_size, 2) + powf(spectrum_im[i]/window_size, 2);
    
  //    FLOAT e_db = - 10*log10(e+epsilon_db);
  //    power_spectrum[i] = e;
  // }
  // spectrogram.push_back(power_spectrum);
  // if (spectrogram.size() > 128) {
  //   spectrogram.pop_front();
  // }
  // delete[] in_re;
  // delete[] in_im;
}

COMPLEX InputPoint::amplitude(FLOAT frequency) const {
  uint index = frequency/frequency_step;
  COMPLEX out(spectrum_re[index], spectrum_im[index]);
  return out;
}


void InputPoint::plotSpectrum() const {
  std::stringstream ss;
  ss <<name<<"_spectrum.txt";
  std::string str(ss.str());
  std::ofstream  out_file;
  out_file.open(str.c_str());

  FLOAT epsilon_db = 1e-18; 
  
  for (uint i = 0; i < nb_frequencies; ++i) {
    FLOAT e = powf(spectrum_re[i]/window_size, 2) + powf(spectrum_im[i]/window_size, 2);    
    FLOAT e_db = - 10*log10(e  + epsilon_db);
    out_file << i*frequency_step << " " << e <<"\n";
  }
  out_file.close();
}

void InputPoint::plotSpectrogram() const {
  std::stringstream ss;
  ss <<name<<"_spectrogram.txt";
  std::string str(ss.str());
  std::ofstream  out_file;
  out_file.open(str.c_str());

  std::list<std::vector<FLOAT> >::const_iterator it;
  uint i = 0;
  for (it = spectrogram.begin(); it != spectrogram.end(); ++it, ++i) {
    for (uint j = 0; j < nb_frequencies; ++j) {
      out_file << i*dt<<" "<<j*frequency_step << " " << (*it)[j] <<"\n";
  }
    out_file <<"\n";
      }
    out_file.close();
}

void InputPoint::plotSamples() const {
  std::stringstream ss;
  ss <<name<<"_samples.txt";
  std::string str(ss.str());
  std::ofstream  out_file;
  out_file.open(str.c_str());
  
  std::list<FLOAT>::const_iterator it;
  uint i = 0;
  for (it = samples.begin(); it != samples.end(); ++it, ++i) {
    out_file << i*dt<<" "<<(*it) <<"\n";
  }

  out_file.close();
}

void InputPoint::move(VEC2 trans) {
  pos += trans;
}
