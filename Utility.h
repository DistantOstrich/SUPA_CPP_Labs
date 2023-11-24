/*
Duncan Robb - 24th November 2023
Utility functions for SUPA C++ excercises
*/

#pragma once

#include "sstream"
#include "vector"
#include "fstream"
#include "iostream"

// Split a string into a vector of strings using the given delimeter character
std::vector<std::string> splitToString(const std::string& s, char delim) {

	std::stringstream ss(s);
	std::string item;
	std::vector<std::string> elems;
	while (std::getline(ss, item, delim)) {
		if (item.size() > 0) {
			elems.push_back(item);
		}
	}
	return elems;
}

// Split a string into a vector of floats using the given delimeter character
std::vector<float> splitToFloat(const std::string& s, char delim) {

	std::stringstream ss(s);
	std::string item;
	std::vector<float> elems;
	while (std::getline(ss, item, delim)) {
		if (item.size() > 0) {
			float d = std::stod(item);
			elems.push_back(d);
		}
	}
	return elems;
}

// Split a string into a vector of ints using the given delimeter character
std::vector<int> splitToInt(const std::string& s, char delim) {

	std::stringstream ss(s);
	std::string item;
	std::vector<int> elems;
	while (std::getline(ss, item, delim)) {
		if (item.size() > 0) {
			int i = std::stoi(item);
			elems.push_back(i);
		}
	}
	return elems;
}

// Offer a numbered list of options in the terminal, and return the chosen number
int coutToInt(std::initializer_list<std::string> options, int startIndex) {

	for (std::string o : options) {
		std::cout << startIndex << " - " << o << std::endl;
		startIndex++;
	}
	int result;
	std::cin >> result;
	return result;
}

// Offer a numbered list of options in the terminal, and return the chosen number
int coutToInt(std::vector<std::string> options, int startIndex) {

	for (std::string o : options) {
		std::cout << startIndex << " - " << o << std::endl;
		startIndex++;
	}
	int result;
	std::cin >> result;
	return result;
}

// Check if a test string contains a substring
bool contains(std::string testString, std::string segment) {
	return testString.find(segment) != std::string::npos;
}

// Read in a comma-separated text file to a vector of floats
std::vector<std::vector<float>> readInData(std::string filename) {

    std::ifstream ifs(filename);
    std::string line;
	std::vector<std::vector<float>> pairs = std::vector<std::vector<float>>();
	while (std::getline(ifs, line)) {
		if (line[0] == 'x') { continue; }
		pairs.push_back(splitToFloat(line, ','));
	}
    return pairs;
}

// Round a float to the nearest integer
int roundToInt(float f) {
    return std::roundf(f);
}