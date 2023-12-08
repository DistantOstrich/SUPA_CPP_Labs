/*
Header for the main file for SUPA C++ Exercise 2
*/

#pragma once

#include "../../Utility.h"
#include "FiniteFunctions.h"

void testDefault(std::vector<double> data);
void testGaussian(std::vector<double> data, double mean = 0, double sigma = 1);
void testCauchyLorentz(std::vector<double> data, double x0 = 0, double gamma = 1);
void testCrystalBall(std::vector<double> data, double mean = 0, double n = 2, double alpha = 1, double sigma = 1);