#pragma once

#include "math.h"
#include "iostream"
#include "fstream"

#include "../Utility.h"

void printNLinesOfData(std::vector<std::vector<float>> data, int n, int m, bool saveResults) {

    std::ofstream of;
    if (saveResults) {
        of.open("PrintNLines.txt");
    }

    int nLines = data.size();
    if (n < 0) { 
        std::cout << "Error: Enter a positive number" << std::endl;
        return;
    }
    if (m > data.size()) {
        std::cout << "Error: The data file only has " << nLines << " entries" << std::endl;
        return;
    }

    for (int i = n; i < n + m; i++) {
        std::vector<float> line = data[i];
        std::cout << line[0] << "\t" << line[1] << std::endl;
        if (saveResults) {
            of << line[0] << "\t" << line[1] << std::endl;
        }
    }

    if (saveResults) {
        of.close();
    }
}

void printNLinesOfData(std::vector<std::vector<float>> data, int n, bool saveResults) {
    printNLinesOfData(data, 0, n, saveResults);
}

void displayNLines(std::vector<std::vector<float>> data, bool saveResults) {

	std::cout << "To print the first N lines of the data file, enter a single whole number" << std::endl;
    std::cout << "OR" << std::endl;
    std::cout << "To print a range of lines from the data file, enter two whole numbers separated by a comma to indicate the first line and number of lines to be printed" << std::endl;
	std::cout << "Remember the file is zero-indexed, i.e. the first line is line 0" << std::endl;

    std::string inLine;
    std::cin >> inLine;

	if (contains(inLine, ",")) {
		std::vector<int> limits = splitToInt(inLine, ',');
		printNLinesOfData(data, limits[0], limits[1], saveResults);
	}
	else {
		int limit = std::stoi(inLine);
		printNLinesOfData(data, limit, saveResults);
	}
	std::cout << "Done!" << std::endl;
}

void calculateMagnitude(std::vector<std::vector<float>> data, bool saveResults) {

    std::ofstream of;
    if (saveResults) {
        of.open("CalculateMagnitude.txt");
    }

    std::cout << "x\ty\tMagnitude" << std::endl;
    if (saveResults) {
        of << "x\ty\tMagnitude" << std::endl;
    }

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

    if (saveResults) {
        of.close();
    }
}

float calculateSumXY(std::vector<std::vector<float>> data, int N) {
    float sum = 0.0;
    for (int i = 0; i < N; i++) {
        float xy = data[i][0] * data[i][1];
        sum += xy;
    }
    return sum;
}

float calculateSumXSquared(std::vector<std::vector<float>> data, int N) {
    float sum = 0.0;
    for (int i = 0; i < N; i++) {
        float xsquared = data[i][0] * data[i][0];
        sum += xsquared;
    }
    return sum;
}

float calculateSumX(std::vector<std::vector<float>> data, int N) {
    float sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += data[i][0];
    }
    return sum;
}

float calculateSumY(std::vector<std::vector<float>> data, int N) {
    float sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += data[i][1];
    }
    return sum;
}

float calculateP(std::vector<std::vector<float>> data) {

    int N = data.size();
    float sumX = calculateSumX(data, N);
    float sumY = calculateSumY(data, N);
    float sumX2 = calculateSumXSquared(data, N);
    float sumXY = calculateSumXY(data, N);

    float top = (N * sumXY) - (sumX * sumY);
    float bottom = (N * sumX2) - (sumX * sumX);

    return top / bottom;
}

float calculateQ(std::vector<std::vector<float>> data) {

    int N = data.size();
    float sumX = calculateSumX(data, N);
    float sumY = calculateSumY(data, N);
    float sumX2 = calculateSumXSquared(data, N);
    float sumXY = calculateSumXY(data, N);

    float top = (sumX2 * sumY) - (sumXY * sumX);
    float bottom = (N * sumX2) - (sumX * sumX);

    return top / bottom;
}

float calculateY(float x, float p, float q) {
    return (p * x) + q;
}

float calculateThisChiSquared(float x, float y, float fitY, float errY) {

    float top = std::pow(y - fitY, 2);
    float bottom = std::pow(errY, 2);
    return top / bottom;
}

float calculateChiSquared(std::vector<std::vector<float>> data, std::vector<std::vector<float>> errors, float p, float q) {

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

void fitLine(std::vector<std::vector<float>> data, std::vector<std::vector<float>> errors, bool saveResults) {

    float p = calculateP(data);
    float q = calculateQ(data);

    float chi2 = calculateChiSquared(data, errors, p, q);

    std::ofstream of;
    if (saveResults) {
        of.open("CalculateChiSquared.txt");
    }

    std::cout << "Fit: y = " << p << " x + " << q << std::endl;
    std::cout << "Chi Squared = " << chi2 << std::endl;
    if (saveResults) {
        of << "Fit: y = " << p << " x + " << q << std::endl;
        of << "Chi Squared = " << chi2 << std::endl;

        of.close();
    }
}

int roundToInt(float f) {
    return std::roundf(f);
}

float calculateXToTheY(float x, float Y) {

    int y = roundToInt(Y);
    int i = 0;
    float product = 1.0;
    while (i < y) {
        product *= x;
        i++;
    }
    if (y < 0) {
        return 1.0 / product;
    }
    else {
        return product;
    }
}

void xToTheY(std::vector<std::vector<float>> data, bool saveResults) {

    std::ofstream of;
    if (saveResults) {
        of.open("CalculateXToTheY.txt");
    }

    std::cout << "x\ty\tx^y" << std::endl;
    if (saveResults) {
        of << "x\ty\tx^y" << std::endl;
    }

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

    if (saveResults) {
        of.close();
    }
}