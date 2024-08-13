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
        std::string value;

        if (!(iss >> key >> value)) {
            continue;  // Skip invalid lines
        }

        // Determine if the value is numeric or a string (e.g., a filename)
        try {
            double num_value = std::stod(value);
            constants[key] = num_value;
        } catch (std::invalid_argument&) {
            filenames[key] = value;
        }
    }
}

void Input::parseCommandLineArgs(int argc, char *argv[]) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg.find('=') != std::string::npos) {
            std::string key = arg.substr(0, arg.find('='));
            std::string value = arg.substr(arg.find('=') + 1);

            // Determine if the value is numeric or a string (e.g., a filename)
            try {
                double num_value = std::stod(value);
                constants[key] = num_value;
            } catch (std::invalid_argument&) {
                filenames[key] = value;
            }
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

std::string Input::getFilename(const std::string &key) const {
    auto it = filenames.find(key);
    if (it != filenames.end()) {
        return it->second;
    } else {
        return "";  // Return empty string if the filename is not found
    }
}

void Input::setConstant(const std::string &key, double value) {
    constants[key] = value;
}

void Input::setFilename(const std::string &key, const std::string &value) {
    filenames[key] = value;
}
