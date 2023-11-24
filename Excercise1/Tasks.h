/*
Duncan Robb - 24th November 2023
Tasks for SUPA C++ Excercise 1
*/

#pragma once

#include "math.h"
#include "iostream"
#include "fstream"

#include "../Utility.h"

// Print n lines of data, starting from line m
void printNLinesOfData(std::vector<std::vector<float>> data, int m, int n, bool saveResults) {

    // Open the results file, overwriting if it already exists
    std::ofstream of;
    if (saveResults) {
        of.open("PrintNLines.txt");
    }

    // Check that the user has entered sensible limits
    int nLines = data.size();
    if (m < 0) { 
        std::cout << "Error: Enter a positive number" << std::endl;
        return;
    }
    if (n > data.size()) {
        std::cout << "Error: The data file only has " << nLines << " entries" << std::endl;
        return;
    }

    // Print the requested lines
    for (int i = m; i < m + n; i++) {
        std::vector<float> line = data[i];
        std::cout << line[0] << "\t" << line[1] << std::endl;
        if (saveResults) {
            of << line[0] << "\t" << line[1] << std::endl;
        }
    }

    // Close the results file
    if (saveResults) {
        of.close();
    }
}

// Print the first n lines of data
void printNLinesOfData(std::vector<std::vector<float>> data, int n, bool saveResults) {
    printNLinesOfData(data, 0, n, saveResults);
}

// Choose a range of the data file and print it to the terminal and (optionally) a text file
void displayNLines(std::vector<std::vector<float>> data, bool saveResults) {

	std::cout << "To print the first N lines of the data file, enter a single whole number" << std::endl;
    std::cout << "OR" << std::endl;
    std::cout << "To print a range of lines from the data file, enter two whole numbers separated by a comma to indicate the first line and number of lines to be printed" << std::endl;
	std::cout << "Remember the file is zero-indexed, i.e. the first line is line 0" << std::endl;

    std::string inLine;
    std::cin >> inLine;

    // If the user has entered two, comma separated numbers, separated them and print lines m to m + n
	if (contains(inLine, ",")) {
		std::vector<int> limits = splitToInt(inLine, ',');
		printNLinesOfData(data, limits[0], limits[1], saveResults);
	}
    // Otherwise, print the number of lines requested, starting from line 0
	else {
		int limit = std::stoi(inLine);
		printNLinesOfData(data, limit, saveResults);
	}
}

// Calculate the magnitude of each point in the data file
void calculateMagnitude(std::vector<std::vector<float>> data, bool saveResults) {

    // OPen the results file
    std::ofstream of;
    if (saveResults) {
        of.open("CalculateMagnitude.txt");
    }

    // Print column headers
    std::cout << "x\ty\tMagnitude" << std::endl;
    if (saveResults) {
        of << "x\ty\tMagnitude" << std::endl;
    }

    // For each point, calculate the magnitude and print x, y, and magnitude
    int nPoints = data.size();
    for (int i = 0; i < nPoints; i++) {
        float x = data[i][0];
        float y = data[i][1];
        float mag = std::sqrt(std::pow(x, 2) + std::pow(y, 2));
        std::cout << x << "\t" << y << "\t" << mag << std::endl;
        if (saveResults) {
            of << x << "\t" << y << "\t" << mag << std::endl;
        }
    }

    // Close the results file
    if (saveResults) {
        of.close();
    }
}

// Calculate the sum from i=0 to i=N of xi*yi
float calculateSumXY(std::vector<std::vector<float>> data, int N) {
    float sum = 0.0;
    for (int i = 0; i < N; i++) {
        float xy = data[i][0] * data[i][1];
        sum += xy;
    }
    return sum;
}

// Calculate the sum from i=0 to i=N of xi^2
float calculateSumXSquared(std::vector<std::vector<float>> data, int N) {
    float sum = 0.0;
    for (int i = 0; i < N; i++) {
        float xsquared = data[i][0] * data[i][0];
        sum += xsquared;
    }
    return sum;
}

// Calculate the sum from i=0 to i=N of xi
float calculateSumX(std::vector<std::vector<float>> data, int N) {
    float sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += data[i][0];
    }
    return sum;
}

// Calculate the sum from i=0 to i=N of yi
float calculateSumY(std::vector<std::vector<float>> data, int N) {
    float sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += data[i][1];
    }
    return sum;
}

// Calculate the gradient of the line from the given equation
float calculateP(std::vector<std::vector<float>> data) {

    // Calculate the component parts
    int N = data.size();
    float sumX = calculateSumX(data, N);
    float sumY = calculateSumY(data, N);
    float sumX2 = calculateSumXSquared(data, N);
    float sumXY = calculateSumXY(data, N);

    // Calculate the numerator and denominator
    float top = (N * sumXY) - (sumX * sumY);
    float bottom = (N * sumX2) - (sumX * sumX);

    // Return the final value
    return top / bottom;
}

// Calculate the intercept of the line from the given equation
float calculateQ(std::vector<std::vector<float>> data) {

    // Calculate the component parts
    int N = data.size();
    float sumX = calculateSumX(data, N);
    float sumY = calculateSumY(data, N);
    float sumX2 = calculateSumXSquared(data, N);
    float sumXY = calculateSumXY(data, N);

    // Calculate the numerator and denominator
    float top = (sumX2 * sumY) - (sumXY * sumX);
    float bottom = (N * sumX2) - (sumX * sumX);

    // Return the final value
    return top / bottom;
}

// Calculate the y value of the line y = px + q at a given value of x
float calculateY(float x, float p, float q) {
    return (p * x) + q;
}

// Calculate the contribution to the Chi^2 from a single point
float calculateThisChiSquared(float x, float y, float fitY, float errY) {
    float top = std::pow(y - fitY, 2);
    float bottom = std::pow(errY, 2);
    return top / bottom;
}

// Calculate Chi^2 for the fit to the data
float calculateChiSquared(std::vector<std::vector<float>> data, std::vector<std::vector<float>> errors, float p, float q) {

    // Calculate the contribution of each point in turn, and add it to the total
    int N = data.size();
    float sum = 0.0;
    for (int i = 0; i < N; i++) {
        float x = data[i][0];
        float y = data[i][1];
        float fitY = calculateY(x, p, q);
        float errY = errors[i][1];

        sum += calculateThisChiSquared(x, y, fitY, errY);
    }
    return sum;
}

// Fit a line to the data using the given equations, and determine the Chi^2 of the fit
void fitLine(std::vector<std::vector<float>> data, std::vector<std::vector<float>> errors, bool saveResults) {

    // Calculate the equation of the line
    float p = calculateP(data);
    float q = calculateQ(data);

    // Calculate the Chi^2
    float chi2 = calculateChiSquared(data, errors, p, q);

    // Open the results file
    std::ofstream of;
    if (saveResults) {
        of.open("CalculateChiSquared.txt");
    }

    // Print the results
    std::cout << "Fit: y = " << p << " x + " << q << std::endl;
    std::cout << "Chi Squared = " << chi2 << std::endl;
    if (saveResults) {
        of << "Fit: y = " << p << " x + " << q << std::endl;
        of << "Chi Squared = " << chi2 << std::endl;

        of.close();
    }
}

// Calculate x to the power of y (with y rounded to the nearest integer) for a given x and y
float calculateXToTheY(float x, float Y) {

    // Round y to the nearest integer
    int y = roundToInt(Y);


    // If y id positive or zero, calculate x^y by multiplying 1 by x, y times
    if (y >= 0) {
        int i = 0;
        float product = 1.0;
        while (i < y) {
            product *= x;
            i++;
        }
        return product;
    }
    // If y is negative, calculate x^y by multiplying 1 by x, -y times, then take 1 over the product
    else {
        int i = 0;
        float product = 1.0;
        while (i > y) {
            product *= x;
            i--;
        }
        return 1.0 / product;
    }
}

// Calculate x to the power of y for each coordinate (x, y) in the data file
void xToTheY(std::vector<std::vector<float>> data, bool saveResults) {

    // Open the results file
    std::ofstream of;
    if (saveResults) {
        of.open("CalculateXToTheY.txt");
    }

    // Print column headers
    std::cout << "x\ty\tx^y" << std::endl;
    if (saveResults) {
        of << "x\ty\tx^y" << std::endl;
    }

    // Calculate x^y for each point and print the results
    int nPoints = data.size();
    for (int i = 0; i < nPoints; i++) {
        float x = data[i][0];
        float y = data[i][1];
        float xToY = calculateXToTheY(x, y);
        std::cout << x << "\t" << y << "\t" << xToY << std::endl;
        if (saveResults) {
            of << x << "\t" << y << "\t" << xToY << std::endl;
        }
    }

    // Close the results file
    if (saveResults) {
        of.close();
    }
}