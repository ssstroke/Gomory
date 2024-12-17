#include "Gomory.hpp"

#include <boost/rational.hpp>

#include <iostream>
#include <vector>
#include <string>

int main()
{
    const std::vector<r64> Z_FUNCTION_COEFFICIENTS = {
        -3, 0
    };
    const std::vector<std::vector<r64>> INITIAL_TABLE_COEFFICIENTS = {
        {0, 3, 4},
        {2, 1, 4},
        {-4, 3, 4},
        {-3, 1, 0}
    };
    const std::vector<Sign> INITIAL_TABLE_SIGNS = {
        GREATER_OR_EQUAL,
        GREATER_OR_EQUAL,
        LESS_OR_EQUAL,
        GREATER_OR_EQUAL
    };
    //const std::vector<r64> Z_FUNCTION_COEFFICIENTS = {
    //    7, 9, 0, 0
    //};
    //const std::vector<std::vector<r64>> INITIAL_TABLE_COEFFICIENTS = {
    //    {-1, 3, 1, 0, 6},
    //    {7, 1, 0, 1, 35},
    //};
    //const std::vector<Sign> INITIAL_TABLE_SIGNS = {
    //    LESS_OR_EQUAL,
    //    LESS_OR_EQUAL
    //};

    const auto solution = Gomory(Z_FUNCTION_COEFFICIENTS, INITIAL_TABLE_COEFFICIENTS, INITIAL_TABLE_SIGNS);

    return 0;
}
