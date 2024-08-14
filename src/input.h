#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <map>

/**
 * @brief Manages reading and parsing of input files and command-line arguments.
 */
class Input {
private:
    std::map<std::string, std::string> filenames;  ///< Map to store filename values as key-value pairs
    std::map<std::string, double> constants;       ///< Map to store constant values as key-value pairs

public:
    /**
     * @brief Constructs an Input object and parses the input file.
     * @param filename The name of the input file to be parsed.
     */
    Input(const std::string &filename);

    /**
     * @brief Parses the input file and stores constants and filenames.
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
     * @brief Gets the value of a filename.
     * @param key The name of the filename to retrieve.
     * @return The filename as a string, or an empty string if not found.
     */
    std::string getFilename(const std::string &key) const;

    /**
     * @brief Sets or updates the value of a constant.
     * @param key The name of the constant to set.
     * @param value The value to assign to the constant.
     */
    void setConstant(const std::string &key, double value);

    /**
     * @brief Sets or updates the value of a filename.
     * @param key The name of the filename to set.
     * @param value The filename as a string.
     */
    void setFilename(const std::string &key, const std::string &value);
};

#endif // INPUT_H
