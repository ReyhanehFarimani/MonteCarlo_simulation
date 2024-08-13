#include "input.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

Input::Input(const std::string &filename) {
    parseInputFile(filename);
}

void Input::parseInputFile(const std::string &filename) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Could not open input file: " + filename);
    }

    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string key;
        double value;

        if (!(iss >> key >> value)) {
            continue;  // Skip invalid lines
        }

        constants[key] = value;
    }
}

void Input::parseCommandLineArgs(int argc, char *argv[]) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg.find('=') != std::string::npos) {
            std::string key = arg.substr(0, arg.find('='));
            double value = std::stod(arg.substr(arg.find('=') + 1));
            constants[key] = value;
        }
    }
}

double Input::getConstant(const std::string &key) const {
    auto it = constants.find(key);
    if (it != constants.end()) {
        return it->second;
    } else {
        return 0.0;  // Default value if the key is not found
    }
}

void Input::setConstant(const std::string &key, double value) {
    constants[key] = value;
}
