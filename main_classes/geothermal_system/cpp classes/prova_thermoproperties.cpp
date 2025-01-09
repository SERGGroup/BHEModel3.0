//
// Created by Utente on 17/12/2024.
//
#include "REFPROP.h"
#include <iostream>

int main() {
  int ierr = 0;
  int components = 1;
  char hrf[4] = "DEF";
  char herr[255];
  char fluid[10] = "R134a.fld";
  char hfm[8] = "HMX.BNC";

  // Setup REFPROP: load the fluid file
  SETUPdll(&components, fluid, hfm, hrf, &ierr, herr, 10000, 255, 3, 255);
  if (ierr != 0) {
    std::cerr << "Error: " << herr << std::endl;
    return 1;
  }

  return 0;
}