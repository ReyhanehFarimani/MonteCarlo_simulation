// unit_test_serial/test_input.cpp
#include "catch.hpp"
#include <fstream>
#include "../serial_src/input.h"

static const std::string tempFile = "test_input.txt";

TEST_CASE("parseInputFile reads numeric and filename entries", "[Input]") {
    // Create a temporary input file
    std::ofstream ofs(tempFile);
    ofs << "num1 3.5\n";
    ofs << "file1 data.out\n";
    ofs << "bad_line\n";            // should be skipped
    ofs << "num2 2.71828\n";
    ofs.close();

    // Parse it
    Input inp(tempFile);

    REQUIRE(inp.getConstant("num1") == Approx(3.5));
    REQUIRE(inp.getConstant("num2") == Approx(2.71828));
    REQUIRE(inp.getFilename("file1") == "data.out");

    // Missing keys give defaults
    REQUIRE(inp.getConstant("missing") == Approx(0.0));
    REQUIRE(inp.getFilename("missing") == "");
}

TEST_CASE("parseCommandLineArgs overrides and adds values", "[Input]") {
    // Start from an empty file
    std::ofstream ofs(tempFile); ofs.close();
    Input inp(tempFile);

    // No entries yet
    REQUIRE(inp.getConstant("x") == Approx(0.0));
    REQUIRE(inp.getFilename("y") == "");

    // Simulate argv: program name + three k=v args
    const char* argv[] = { "prog", "x=1.23", "y=out.txt", "z=notanumber" };
    inp.parseCommandLineArgs(4, const_cast<char**>(argv));

    // Numeric and non-numeric both work
    REQUIRE(inp.getConstant("x") == Approx(1.23));
    REQUIRE(inp.getFilename("y") == "out.txt");
    // “z” isn’t a number, so it goes into filenames
    REQUIRE(inp.getFilename("z") == "notanumber");
}

TEST_CASE("setConstant and setFilename store and retrieve", "[Input]") {
    std::ofstream ofs(tempFile); ofs.close();
    Input inp(tempFile);

    inp.setConstant("c", 10.0);
    inp.setFilename("f", "file.dat");

    REQUIRE(inp.getConstant("c") == Approx(10.0));
    REQUIRE(inp.getFilename("f") == "file.dat");
}
