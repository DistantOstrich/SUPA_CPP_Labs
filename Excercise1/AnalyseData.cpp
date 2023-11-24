/*
Duncan Robb - 24th November 2023
Main file for SUPA C++ Excercise 1
*/

#include "Tasks.h"
#include "../Utility.h"


int main() {
    
	// Read in the coordinates and errors from the text files
    std::vector<std::vector<float>> coordinates = readInData("input2D_float.txt");
	std::vector<std::vector<float>> errors = readInData("error2D_float.txt");

	bool again = true;
	while (again) {
	std::cout << "Which task would you like to perform?" << std::endl;

		// List the available tasks
    	int task = coutToInt({ "Display n lines of the file", "Calculate the magnitude of each data point", "Fit a line to the data", "Calculate x^y for each data point (with y rounded to the nearest integer)", "All of the above" }, 1);
		int saveResults;
		std::cout << "Do you want to save the results of the task in a text file?" << std::endl;
		std::cout << "0 - No" << std::endl << "1 - Yes" << std::endl;
		std::cin >> saveResults;

		// Complete the selected task or tasks
		switch (task) {
		case 1:
			displayNLines(coordinates, saveResults);
			break;
		case 2:
			calculateMagnitude(coordinates, saveResults);
			break;
		case 3:
			fitLine(coordinates, errors, saveResults);
			break;
		case 4:
			xToTheY(coordinates, saveResults);
			break;
		case 5:
			displayNLines(coordinates, saveResults);
			calculateMagnitude(coordinates, saveResults);
			fitLine(coordinates, errors, saveResults);
			xToTheY(coordinates, saveResults);
			break;
		default:
		std::cout << task << " was not one of the options. Do better" << std::endl;
			break;
		}
		std::cout << "Do you want to do another task?" << std::endl;
		std::cout << "0 - No" << std::endl << "1 - Yes" << std::endl;
		std::cin >> again;		
	}
}