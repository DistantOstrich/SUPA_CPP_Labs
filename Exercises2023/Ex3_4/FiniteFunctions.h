#pragma once //Replacement for IFNDEF

#include <string>
#include <vector>
#include <numbers>
#include <random>
#include "../../GNUplot/gnuplot-iostream.h"

class FiniteFunction{

public:
  FiniteFunction(); //Empty constructor
  FiniteFunction(double range_min, double range_max, std::string outfile); //Variable constructor
  ~FiniteFunction(); //Destructor
  double rangeMin(); //Low end of the range the function is defined within
  double rangeMax(); //High end of the range the function is defined within
  double integral(int Ndiv = 1000); 
  std::vector< std::pair<double,double> > scanFunction(int Nscan = 1000); //Scan over function to plot it (slight hack needed to plot function in gnuplot)
  void setRangeMin(double RMin);
  void setRangeMax(double RMax);
  void setOutfile(std::string outfile);
  void plotFunction(); //Plot the function using scanFunction
  void plotFunction(int nSteps);
  std::vector<double> sampleFunctionMetropolis(int nSamples = 10000, double sigma = 1.0); // Generate sample data from the distribution using the Metropolis algorithm
  
  //Plot the supplied data points (either provided data or points sampled from function) as a histogram using NBins
  void plotData(std::vector<double> &points, int NBins, bool isdata=true); //NB! use isdata flag to pick between data and sampled distributions
  virtual void printInfo(); //Dump parameter info about the current function (Overridable)
  virtual double callFunction(double x); //Call the function with value x (Overridable)

  //Protected members can be accessed by child classes but not users
protected:
  double m_RMin;
  double m_RMax;
  double m_Integral;
  int m_IntDiv = 0; //Number of division for performing integral
  std::string m_FunctionName;
  std::string m_OutData; //Output filename for data
  std::string m_OutPng; //Output filename for plot
  std::vector< std::pair<double,double> > m_data; //input data points to plot
  std::vector< std::pair<double,double> > m_samples; //Holder for randomly sampled data 
  std::vector< std::pair<double,double> > m_function_scan; //holder for data from scanFunction (slight hack needed to plot function in gnuplot)
  bool m_plotfunction = false; //Flag to determine whether to plot function
  bool m_plotdatapoints = false; //Flag to determine whether to plot input data
  bool m_plotsamplepoints = false; //Flag to determine whether to plot sampled data 
  double integrate(int Ndiv);
  std::vector< std::pair<double, double> > makeHist(std::vector<double> &points, int Nbins); //Helper function to turn data points into histogram with Nbins
  void checkPath(std::string outstring); //Helper function to ensure data and png paths are correct
  void generatePlot(Gnuplot &gp);
  double generateNextValueMetropolis(double xPrev, double sigma, std::uniform_real_distribution<double> zeroToOne, std::mt19937* rng, bool* success); // Generate the next value in the Metropolis algorithm
  
private:
  double invxsquared(double x); //The default functional form
};


// Daughter class for normal distribution, inheriting from FiniteFunction
class Gaussian : public FiniteFunction {

public:
  Gaussian();
  Gaussian(double mean, double sigma, double rangeMin, double rangeMax, std::string outfile);
  ~Gaussian();
  void setMean(double mean);
  void setSigma(double sigma);
  double mean();
  double sigma();
  double callFunction(double x);
  void printInfo();

protected:
  double m_mean;
  double m_sigma;

private:
  double gaussian(double x);
};


// Daughter class for Cauchy-Lorentz distribution, inheriting from FiniteFunction
class CauchyLorentz : public FiniteFunction {

public:
  CauchyLorentz();
  CauchyLorentz(double x0, double gamma, double rangeMin, double rangeMax, std::string outfile);
  ~CauchyLorentz();
  void setX0(double x0);
  void setGamma(double gamma);
  double x0();
  double gamma();
  double callFunction(double x);
  void printInfo();

protected:
  double m_x0;
  double m_gamma;

private:
  double cauchyLorentz(double x);
};


// Daughter class for negative Crystal Ball distribution, inheriting from FiniteFunction
class CrystalBall : public FiniteFunction {

public:
  CrystalBall();
  CrystalBall(double mean, double n, double alpha, double sigma, double rangeMin, double rangeMax, std::string outfile);
  ~CrystalBall();
  void setMean(double mean);
  void setN(double n);
  void setAlpha(double alpha);
  void setSigma(double sigma);
  double mean();
  double n();
  double alpha();
  double sigma();
  double callFunction(double x);
  void printInfo();

protected:
  double m_mean;
  double m_n;
  double m_alpha;
  double m_sigma;
  double m_magAlpha;
  double m_A;
  double m_B;
  double m_C;
  double m_D;
  double m_N;

private:
  double crystalBall(double x);
};