#include <iostream>
#include <string>
#include <vector>
#include "FiniteFunctions.h"
#include <filesystem> //To check extensions in a nice way

#include "gnuplot-iostream.h" //Needed to produce plots (not part of the course) 

using std::filesystem::path;

//Empty constructor
FiniteFunction::FiniteFunction(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  this->checkPath("DefaultFunction");
  m_Integral = NULL;
}

//initialised constructor
FiniteFunction::FiniteFunction(double range_min, double range_max, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
}

//Plots are called in the destructor
//SUPACPP note: They syntax of the plotting code is not part of the course
FiniteFunction::~FiniteFunction(){
  Gnuplot gp; //Set up gnuplot object
  this->generatePlot(gp); //Generate the plot and save it to a png using "outfile" for naming 
}

/*
###################
//Setters
###################
*/ 
void FiniteFunction::setRangeMin(double RMin) {m_RMin = RMin;};
void FiniteFunction::setRangeMax(double RMax) {m_RMax = RMax;};
void FiniteFunction::setOutfile(std::string Outfile) {this->checkPath(Outfile);};

/*
###################
//Getters
###################
*/ 
double FiniteFunction::rangeMin() {return m_RMin;};
double FiniteFunction::rangeMax() {return m_RMax;};

/*
###################
//Function eval
###################
*/ 
double FiniteFunction::invxsquared(double x) {return 1/(1+x*x);};
double FiniteFunction::callFunction(double x) {return this->invxsquared(x);}; //(overridable)

/*
###################
Integration by hand (output needed to normalise function when plotting)
###################
*/ 
double FiniteFunction::integrate(int Ndiv){ //private
  
  // Define/calculate the integration limits and step size (from range min to range max in Ndiv steps)
  double integral = 0.0;
  double range = this->m_RMax - this->m_RMin;
  double stepSize = range / (double)Ndiv;
  double slice = this->m_RMin;

  // If the range is zero or negative, the integral is not defined so return zero
  if (range <= 0) {
    std::cout << "Integration range = " << range << std::endl;
    return 0.0;
  }

  std::cout << "Step size: " << stepSize << std::endl;

  // For each slice in x...
  while (slice <= this->m_RMax) {

    // Calculate the value of the function at this x value
    double valAtSlice = callFunction(slice);

    // The area of this slice is the function value multiplied by the step size
    double thisSliceIntegral = valAtSlice * stepSize;

    // Add this to the total integral
    integral += thisSliceIntegral;

    // Increment the slice by the step size
    slice += stepSize;
  }
  return integral;  
}
double FiniteFunction::integral(int Ndiv) { //public
  if (Ndiv <= 0){
    std::cout << "Invalid number of divisions for integral, setting Ndiv to 1000" <<std::endl;
    Ndiv = 1000;
  }
  if (m_Integral == NULL || Ndiv != m_IntDiv){
    m_IntDiv = Ndiv;
    m_Integral = this->integrate(Ndiv);
    return m_Integral;
  }
  else return m_Integral; //Don't bother re-calculating integral if Ndiv is the same as the last call
}

/*
###################
//Helper functions 
###################
*/
// Generate paths from user defined stem
void FiniteFunction::checkPath(std::string outfile){
 path fp = outfile;
 m_FunctionName = fp.stem(); 
 m_OutData = m_FunctionName+".data";
 m_OutPng = m_FunctionName+".png";
}

//Print (overridable)
void FiniteFunction::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

/*
###################
//Plotting
###################
*/

//Hack because gnuplot-io can't read in custom functions, just scan over function and connect points with a line... 
void FiniteFunction::plotFunction(){
  m_function_scan = this->scanFunction(10000);
  m_plotfunction = true;
}

void FiniteFunction::plotFunction(int nSteps){
  m_function_scan = this->scanFunction(nSteps);
  m_plotfunction = true;
}

//Transform data points into a format gnuplot can use (histogram) and set flag to enable drawing of data to output plot
//set isdata to true (default) to plot data points in black, set to false to plot sample points in blue
void FiniteFunction::plotData(std::vector<double> &points, int Nbins, bool isdata){
  if (isdata){
    m_data = this->makeHist(points,Nbins);
    m_plotdatapoints = true;
  }
  else{
    m_samples = this->makeHist(points,Nbins);
    m_plotsamplepoints = true;
  }
}


/*
  #######################################################################################################
  ## SUPACPP Note:
  ## The three helper functions below are needed to get the correct format for plotting with gnuplot
  ## In theory you shouldn't have to touch them
  ## However it might be helpful to read through them and understand what they are doing
  #######################################################################################################
 */

//Scan over range of function using range/Nscan steps (just a hack so we can plot the function)
std::vector< std::pair<double,double> > FiniteFunction::scanFunction(int Nscan){
  std::vector< std::pair<double,double> > function_scan;
  double step = (m_RMax - m_RMin)/(double)Nscan;
  double x = m_RMin;
  //We use the integral to normalise the function points
  if (m_Integral == NULL) {
    std::cout << "Integral not set, doing it now" << std::endl;
    this->integral(Nscan);
    std::cout << "integral: " << m_Integral << ", calculated using " << Nscan << " divisions" << std::endl;
  }
  //For each scan point push back the x and y values 
  for (int i = 0; i < Nscan; i++){
    function_scan.push_back( std::make_pair(x,this->callFunction(x)/m_Integral));
    x += step;
  }
  return function_scan;
}

//Function to make histogram out of sampled x-values - use for input data and sampling
std::vector< std::pair<double,double> > FiniteFunction::makeHist(std::vector<double> &points, int Nbins){

  std::vector< std::pair<double,double> > histdata; //Plottable output shape: (midpoint,frequency)
  std::vector<int> bins(Nbins,0); //vector of Nbins ints with default value 0 
  int norm = 0;
  for (double point : points){
    //Get bin index (starting from 0) the point falls into using point value, range, and Nbins
    int bindex = static_cast<int>(floor((point-m_RMin)/((m_RMax-m_RMin)/(double)Nbins)));
    if (bindex<0 || bindex>Nbins){
      continue;
    }
    bins[bindex]++; //weight of 1 for each data point
    norm++; //Total number of data points
  }
  double binwidth = (m_RMax-m_RMin)/(double)Nbins;
  for (int i=0; i<Nbins; i++){
    double midpoint = m_RMin + i*binwidth + binwidth/2; //Just put markers at the midpoint rather than drawing bars
    double normdata = bins[i]/((double)norm*binwidth); //Normalise with N = 1/(Ndata*binwidth)
    histdata.push_back(std::make_pair(midpoint,normdata));
  }
  return histdata;
}

//Function which handles generating the gnuplot output, called in destructor
//If an m_plot... flag is set, the we must have filled the related data vector
//SUPACPP note: They syntax of the plotting code is not part of the course
void FiniteFunction::generatePlot(Gnuplot &gp){

  if (m_plotfunction==true && m_plotdatapoints==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotdatapoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
  }
  else if (m_plotfunction==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title 'function'\n";
    gp.send1d(m_function_scan);
  }

  else if (m_plotdatapoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_data);
  }

  else if (m_plotsamplepoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_samples);
  }
}

// -------------------------------------------Gaussian---------------------------------------------------

// Empty constructor
Gaussian::Gaussian() {
  m_mean = 0.0;
  m_sigma = 1.0;
  m_RMin = -5.0;
  m_RMax = 5.0;
  this->checkPath("Gaussian");
  m_Integral = NULL;
}

// Initialised constructor
Gaussian::Gaussian(double mean, double sigma, double rangeMin, double rangeMax, std::string outfile) {
  m_mean = mean;
  m_sigma = sigma;
  m_RMin = rangeMin;
  m_RMax = rangeMax;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
}

// Plots are called in the destructor
// SUPACPP note: They syntax of the plotting code is not part of the course
Gaussian::~Gaussian() {
  Gnuplot gp; // Set up gnuplot object
  this->generatePlot(gp); // Generate the plot and save it to a png using "outfile" for naming 
}

void Gaussian::setMean(double mean) {m_mean = mean;};
void Gaussian::setSigma(double sigma) {m_sigma = sigma;};

double Gaussian::mean() {return m_mean;};
double Gaussian::sigma() {return m_sigma;};

// Print the gaussian parameters
void Gaussian::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "mean: " << m_mean << std::endl;
  std::cout << "sigma:" << m_sigma << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

// Define the equation of a gaussian (normal distribution)
double Gaussian::gaussian(double x) {

  // Calculate the normalisation constant
  double A = 1.0 / (this->m_sigma * 2 * std::numbers::pi);

  // Calculate the exponent
  double exponent = -0.5 * std::pow((x - this->m_mean) / this->m_sigma, 2);

  // Calculate the value of the gaussian at this x
  return A * std::exp(exponent);
}

double Gaussian::callFunction(double x) { return this->gaussian(x); };

// -------------------------------------------Cauchy-Lorentz---------------------------------------------------

// Empty constructor
CauchyLorentz::CauchyLorentz() {
  m_x0 = 0.0;
  m_gamma = 1.0;
  m_RMin = -5.0;
  m_RMax = 5.0;
  this->checkPath("Cauchy-Lorentz");
  m_Integral = NULL;
}

// Initialised constructor
CauchyLorentz::CauchyLorentz(double x0, double gamma, double rangeMin, double rangeMax, std::string outfile) {
  m_x0 = x0;
  m_gamma = gamma;
  m_RMin = rangeMin;
  m_RMax = rangeMax;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
}

// Plots are called in the destructor
// SUPACPP note: They syntax of the plotting code is not part of the course
CauchyLorentz::~CauchyLorentz() {
  Gnuplot gp; // Set up gnuplot object
  this->generatePlot(gp); // Generate the plot and save it to a png using "outfile" for naming 
}

void CauchyLorentz::setX0(double x0) {m_x0 = x0;};
void CauchyLorentz::setGamma(double gamma) {m_gamma = gamma;};

double CauchyLorentz::x0() {return m_x0;};
double CauchyLorentz::gamma() {return m_gamma;};

// Print the gaussian parameters
void CauchyLorentz::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "x0: " << m_x0 << std::endl;
  std::cout << "gamma:" << m_gamma << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

// Define the equation of a Cauchy-Lorentz distribution
double CauchyLorentz::cauchyLorentz(double x) {

  // Calculate pi * gamma
  double piGamma = std::numbers::pi * this->m_gamma;

  // Calculate the squared bit
  double xMinusX0Bit = std::pow((x - this->m_x0) / this->m_gamma, 2);

  // Calculate the denominator
  double denom = piGamma * (1.0 + xMinusX0Bit);

  // Calculate the value of the Cauchy-Lorentz distribution at this x
  return 1.0 / denom;
}

double CauchyLorentz::callFunction(double x) { return this->cauchyLorentz(x); };

// -------------------------------------------Negative Crystal Ball---------------------------------------------------

// Empty constructor
CrystalBall::CrystalBall() {
  m_mean = 0.0;
  m_n = 2.0;
  m_alpha = 1.0;
  m_magAlpha = std::abs(m_alpha);
  m_sigma = 1.0;
  m_RMin = -5.0;
  m_RMax = 5.0;
  this->checkPath("CrystalBall");
  m_Integral = NULL;

  // Calculate exp(-magAlpha^2 / 2)
  double eMagA = std::exp(-0.5 * std::pow(m_magAlpha, 2));

  // Calculate A, B, C, D, and N
  m_A = std::pow(m_n / m_magAlpha, m_n) * eMagA;
  m_B = (m_n / m_magAlpha) - m_magAlpha;
  m_C = (m_n / m_magAlpha) * (1.0 / (m_n - 1.0)) * eMagA;
  m_D = std::sqrt(std::numbers::pi / 2.0) * (1.0 + std::erf(m_magAlpha / std::sqrt(2.0)));
  m_N = 1.0 / (m_sigma * (m_C + m_D));
}

// Initialised constructor
CrystalBall::CrystalBall(double mean, double n, double alpha, double sigma, double rangeMin, double rangeMax, std::string outfile) {
  m_mean = mean;
  m_n = n;
  m_alpha = alpha;
  m_magAlpha = std::abs(m_alpha);
  m_sigma = sigma;
  m_RMin = rangeMin;
  m_RMax = rangeMax;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files

  // Calculate exp(-magAlpha^2 / 2)
  double eMagA = std::exp(-0.5 * std::pow(m_magAlpha, 2));

  // Calculate A, B, C, D, and N
  m_A = std::pow(m_n / m_magAlpha, m_n) * eMagA;
  m_B = (m_n / m_magAlpha) - m_magAlpha;
  m_C = (m_n / m_magAlpha) * (1.0 / (m_n - 1.0)) * eMagA;
  m_D = std::sqrt(std::numbers::pi / 2.0) * (1.0 + std::erf(m_magAlpha / std::sqrt(2.0)));
  m_N = 1.0 / (m_sigma * (m_C + m_D));
}

// Plots are called in the destructor
// SUPACPP note: They syntax of the plotting code is not part of the course
CrystalBall::~CrystalBall() {
  Gnuplot gp; // Set up gnuplot object
  this->generatePlot(gp); // Generate the plot and save it to a png using "outfile" for naming 
}

void CrystalBall::setMean(double mean) {m_mean = mean;};
void CrystalBall::setN(double n) {m_n = n;};
void CrystalBall::setAlpha(double alpha) {m_alpha = alpha;};
void CrystalBall::setSigma(double sigma) {m_sigma = sigma;};

double CrystalBall::mean() {return m_mean;};
double CrystalBall::n() {return m_n;};
double CrystalBall::alpha() {return m_alpha;};
double CrystalBall::sigma() {return m_sigma;};

// Print the gaussian parameters
void CrystalBall::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "mean: " << m_mean << std::endl;
  std::cout << "n:" << m_n << std::endl;
  std::cout << "alpha: " << m_alpha << std::endl;
  std::cout << "sigma:" << m_sigma << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

// Define the equation of a negative crystal ball distribution
double CrystalBall::crystalBall(double x) {

  // Calculate x - mean
  double a = (x - m_mean);

  // Calculate x-mean / sigma
  double b = a / m_sigma;

  // Check which regime we are in, and calculate the value of the distribution accordingly
  if (b > -1.0 * m_alpha) {
    double exponent = (-1.0 * std::pow(a, 2)) / (2.0 * std::pow(m_sigma, 2));
    return std::exp(exponent);
  }
  else {
    double brackets = m_B - a / m_sigma;
    return m_A * std::pow(brackets, -1.0 * m_n);
  }
}

double CrystalBall::callFunction(double x) { return this->crystalBall(x); };