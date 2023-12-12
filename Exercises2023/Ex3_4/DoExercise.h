/*
Header for the main file for SUPA C++ Exercise 2
*/

#pragma once

#include "../../Utility.h"
#include "FiniteFunctions.h"

// Plot the various functions
void testDefault(std::vector<double> data);
void testGaussian(std::vector<double> data, double mean = 0, double sigma = 1);
void testCauchyLorentz(std::vector<double> data, double x0 = 0, double gamma = 1);
void testCrystalBall(std::vector<double> data, double mean = 0, double n = 2, double alpha = 1, double sigma = 1);

// Find the distribution type and parameters of the mystery data
void findParameters(std::vector<double> data);

// Generate and plot a sample of the Cauchy-Lorentz distribution
void sampleAndPlotCauchyLorentz(std::vector<double> data, double sigma = 1.0, double xMin = -10.0, double xMax = 10.0, double mean = -2.0, double gamma = 0.75);