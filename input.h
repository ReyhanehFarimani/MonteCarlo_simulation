#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <map>

/**
 * @brief Manages reading and parsing of input files and command-line arguments.
 */
class Input {
private:
    std::map<std::string, double> constants;  ///< Map to store constant values as key-value pairs

public:
    /**
     * @brief Constructs an Input object and parses the input file.
     * @param filename The name of the input file to be parsed.
     */
    Input(const std::string &filename);

    /**
     * @brief Parses the input file and stores constants.
     * @param filename The name of the input file to be parsed.
     */
    void parseInputFile(const std::string &filename);

    /**
     * @brief Parses command-line arguments and updates constants.
     * @param argc Number of command-line arguments.
     * @param argv Array of command-line arguments.
     */
    void parseCommandLineArgs(int argc, char *argv[]);

    /**
     * @brief Gets the value of a constant.
     * @param key The name of the constant to retrieve.
     * @return The value of the constant, or 0 if not found.
     */
    double getConstant(const std::string &key) const;

    /**
     * @brief Sets or updates the value of a constant.
     * @param key The name of the constant to set.
     * @param value The value to assign to the constant.
     */
    void setConstant(const std::string &key, double value);
};

#endif // INPUT_H
