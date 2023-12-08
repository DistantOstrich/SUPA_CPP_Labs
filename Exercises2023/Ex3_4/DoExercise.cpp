/*
Duncan Robb - 8th December 2023
Main file for SUPA C++ Excercise 2
*/

#include "DoExercise.h"

// Test plotting the default function (1 / 1+x^2)
void testDefault(std::vector<double> data) {

    FiniteFunction ff = FiniteFunction(-10, 10, "Defaut");

    ff.printInfo();
    ff.plotFunction();

    ff.plotData(data, 1000);

    std::cout << "Default: Done!" << std::endl << std::endl;
}

// Test Plotting a Gaussian (normal) distribution
void testGaussian(std::vector<double> data, double mean, double sigma) {

    Gaussian g = Gaussian(mean, sigma, -10, 10, "Gaussian");

    g.printInfo();
    g.plotFunction();

    g.plotData(data, 1000);

    std::cout << "Gaussian: Done!" << std::endl << std::endl;
}

// Test Plotting a Cauchy-Lorentz distribution
void testCauchyLorentz(std::vector<double> data, double mean, double gamma) {

    CauchyLorentz cl = CauchyLorentz(mean, gamma, -10, 10, "Cauchy-Lorentz");

    cl.printInfo();
    cl.plotFunction();

    cl.plotData(data, 1000);

    std::cout << "Cauchy-Lorentz: Done!" << std::endl << std::endl;
}

// Test plotting a negative crystal ball distribution
void testCrystalBall(std::vector<double> data, double mean, double n, double alpha, double sigma) {

    CrystalBall cb = CrystalBall(mean, n, alpha, sigma, -10.0, 10.0, "CrystalBall");

    cb.printInfo();
    cb.plotFunction();

    cb.plotData(data, 1000);

    std::cout << "Crystal Ball: Done!" << std::endl << std::endl;
}

void testAllFunctions(std::vector<double> data) {

    // Default has no parameters
    testDefault(data);

    // Gaussian takes mean and sigma as parameters
    testGaussian(data, -2, 1);

    // Cauchy-Lorentz takes mean and gamma as parameters
    testCauchyLorentz(data, -2, 1);

    // Negative Crystal Ball takes mean, n, alpha, and sigma as parameters
    testCrystalBall(data, -2, 2.0, 1.0, 1.0);
}

// Try different params for the Gaussian distribution to match the data
void findGaussianParams(std::vector<double> data) {

    // Have the user enter their chosen parameters
    std::cout << "Enter the mean and sigma values to try, separated by commas" << std::endl << std::endl;
    std::string inVals;
    std::cin >> inVals;

    bool again = true;
    while (again) {

        // Get the parameters as doubles
        std::vector<double> params = splitToDouble(inVals, ',');

        // Test the parameters (the user will need to look at the output graph CauchyLorentz.png)
        testCauchyLorentz(data, params[0], params[1]);

        // Ask if this looks good or they want to try again with different params
        std::cout << "If these params look good, enter 0" << std::endl;
        std::cout << "Otherwise enter new parameter values to try, with the format as before:" << std::endl;
        std::cout << "mean,sigma" << std::endl;
        std::cin >> inVals;

        if (inVals == "0" || !contains(inVals, ",")) {
            again = false;
        }
    }
}

// Try different params for the Cauchy-Lorentz distribution to match the data
void findCauchyLorentzParams(std::vector<double> data) {

    // Have the user enter their chosen parameters
    std::cout << "Enter the mean and gamma values to try, separated by commas" << std::endl << std::endl;
    std::string inVals;
    std::cin >> inVals;

    bool again = true;
    while (again) {

        // Get the parameters as doubles
        std::vector<double> params = splitToDouble(inVals, ',');

        // Test the parameters (the user will need to look at the output graph CauchyLorentz.png)
        testCauchyLorentz(data, params[0], params[1]);

        // Ask if this looks good or they want to try again with different params
        std::cout << "If these params look good, enter 0" << std::endl;
        std::cout << "Otherwise enter new parameter values to try, with the format as before:" << std::endl;
        std::cout << "mean,gamma" << std::endl;
        std::cin >> inVals;

        if (inVals == "0" || !contains(inVals, ",")) {
            again = false;
        }
    }
}

// Try different params for the Crystal Ball distribution to match the data
void findCrystalBallParams(std::vector<double> data) {

    // Have the user enter their chosen parameters
    std::cout << "Enter the mean, n, alpha, and sigma values to try, separated by commas" << std::endl << std::endl;
    std::string inVals;
    std::cin >> inVals;

    bool again = true;
    while (again) {

        // Get the parameters as doubles
        std::vector<double> params = splitToDouble(inVals, ',');

        // Test the parameters (the user will need to look at the output graph CrystalBall.png)
        testCrystalBall(data, params[0], params[1], params[2], params[3]);

        // Ask if this looks good or they want to try again with different params
        std::cout << "If these params look good, enter 0" << std::endl;
        std::cout << "Otherwise enter new parameter values to try, with the format as before:" << std::endl;
        std::cout << "mean,n,alpha,sigma" << std::endl;
        std::cin >> inVals;

        if (inVals == "0" || !contains(inVals, ",")) {
            again = false;
        }
    }
}

int main() {

    // Read in the data file
    std::string dataFile = "Outputs/data/MysteryData11141.txt";
    std::vector<double> data = readInDataSingles(dataFile);

    // Test each of the functions
    //testAllFunctions(data);

    // Try to find the distribution parameters (It is a Cauchy-Lorentz, mean=-2.0, gamma=0.75)
    bool another = true;
    while (another) {
        int distType = coutToInt({"Gaussian", "Cauchy-Lorentz", "CrystalBall"}, 1);
        switch (distType) {
            case 1:
                findGaussianParams(data);
                break;
            case 2:
                findCauchyLorentzParams(data);
                break;
            case 3:
                findCrystalBallParams(data);
                break;
            default:
                break;
        }
        another = coutToInt({"Done", "Try another distribution"}, 0);
    }
}