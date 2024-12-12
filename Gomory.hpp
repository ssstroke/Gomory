#pragma once

#include <boost/rational.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

using r64 = boost::rational<int64_t>;

enum Sign
{
    LESS,
    LESS_OR_EQUAL,
    EQUAL,
    GREATER_OR_EQUAL,
    GREATER
};

class Gomory
{
public:
    Gomory(const std::vector<r64>& z_function_coefficients, const std::vector<std::vector<r64>>& initial_table_coefficients, const std::vector<Sign>& initial_table_signs)
    {
        // Construct basis.
        this->basis = std::vector<size_t>(initial_table_coefficients.size());
        for (size_t i = 0; i < this->basis.size(); ++i)
        {
            this->basis[i] = z_function_coefficients.size() + i;
        }

        this->z_function_coefficients = z_function_coefficients;

        // Construct grid. It contains:
        // - a_ij coefficients;
        // - b_j coefficients;
        // - simplex_relation values (call it theta values).
        // It's dimensions are:
        // rows: number of inequalities in the initial table;
        // cols: (number of Z function coefficients) +
        //       (number of inequalities in the initial table) +
        //       1 (for b_j coefficients).
        this->GRID_ROWS = initial_table_coefficients.size();
        this->GRID_COLS = initial_table_coefficients.size() + z_function_coefficients.size() + 1;
        this->b_column = GRID_COLS - 1;

        this->grid = std::vector<std::vector<r64>>(GRID_ROWS, std::vector<r64>(GRID_COLS));

        for (size_t row = 0; row < initial_table_coefficients.size(); ++row)
        {
            for (size_t col = 0; col < initial_table_coefficients[0].size() - 1; ++col)
            {
                this->grid[row][col] = initial_table_coefficients[row][col];
            }
        }
        for (size_t row = 0; row < initial_table_coefficients.size(); ++row)
        {
            if (initial_table_signs[row] == GREATER_OR_EQUAL)
            {
                this->grid[row][z_function_coefficients.size() + row] = r64(-1);
            }
            else
            {
                this->grid[row][z_function_coefficients.size() + row] = r64(1);
            }
        }
        for (size_t row = 0; row < initial_table_coefficients.size(); ++row)
        {
            this->grid[row][this->b_column] = initial_table_coefficients[row][initial_table_coefficients[0].size() - 1];
        }
        for (size_t row = 0; row < GRID_ROWS; ++row)
        {
            if (initial_table_signs[row] == GREATER_OR_EQUAL)
            {
                for (size_t col = 0; col < GRID_COLS; ++col)
                {
                    this->grid[row][col] *= -1;
                }
            }
        }

        this->Print();
        std::cout << std::endl;

BNegativeRemoval:
        // Iterate over b_j column and look for smallest negative element.
        // If there is one, find any negative element in the corresponding
        // row, consider it as pivot and perform rectangle transformation.
        // Repeat while there are negative elements in b_j column.
        bool has_negative_b = true;
        while (has_negative_b)
        {
            auto min_b = r64(std::numeric_limits<int64_t>().max(), 1);
            size_t min_b_idx = 0;
            for (size_t row = 0; row < GRID_ROWS; ++row)
            {
                if (this->grid[row][b_column] < min_b)
                {
                    min_b = this->grid[row][b_column];
                    min_b_idx = row;
                }
            }

            if (min_b >= 0)
            {
                has_negative_b = false;
            }
            else
            {
                bool found_negative_a = false;
                for (size_t col = 0; col < GRID_COLS - 1 && found_negative_a == false; ++col)
                {
                    if (this->grid[min_b_idx][col] < 0)
                    {
                        found_negative_a = true;
                        this->Transform(min_b_idx, col);
                        this->basis[min_b_idx] = col;
                    }
                }

                if (found_negative_a == false)
                {
                    std::cout << "Solution does not exist!" << std::endl;
                    return;
                }

                this->Print();
                std::cout << std::endl;
            }
        }

        // Simplex loop.
        bool plan_is_otimal = false;
        while (plan_is_otimal == false)
        {
            // Find min delta.
            r64 delta_min = std::numeric_limits<int64_t>().max();
            size_t delta_min_idx = 0;
            for (size_t col = 0; col < this->GRID_COLS - 1; ++col)
            {
                r64 delta_j = 0;
                for (size_t row = 0; row < this->GRID_ROWS; ++row)
                {
                    delta_j += this->GetCost(this->basis[row]) * this->grid[row][col];
                }
                delta_j -= this->GetCost(col);
                if (delta_j < delta_min)
                {
                    delta_min = delta_j;
                    delta_min_idx = col;
                }

                std::cout << "Delta_" << col + 1 << " = " << delta_j << std::endl;
            }
            std::cout << std::endl;

            // Perform simplex transformation if necessary.
            if (delta_min < 0)
            {
                r64 theta_min = std::numeric_limits<int64_t>().max();
                size_t theta_min_idx = 0;
                for (size_t row = 0; row < this->GRID_ROWS; ++row)
                {
                    if (this->grid[row][delta_min_idx] > 0 &&
                        (this->grid[row][this->b_column] / this->grid[row][delta_min_idx]) < theta_min
                        )
                    {
                        theta_min = this->grid[row][this->b_column] / this->grid[row][delta_min_idx];
                        theta_min_idx = row;
                    }
                }

                this->Transform(theta_min_idx, delta_min_idx);
                this->basis[theta_min_idx] = delta_min_idx;

                this->Print();
                std::cout << std::endl;
            }
            else
            {
                plan_is_otimal = true;
            }
        }

        // If any of b_j elements are non-integer values, proceed.
        size_t max_fractional_b_idx = 0;
        auto max_fractional_b = r64(0);
        for (size_t row = 0; row < this->GRID_ROWS; ++row)
        {
            if (this->grid[row][this->b_column].numerator() != 0 &&
                this->grid[row][this->b_column].numerator() != this->grid[row][this->b_column].denominator() &&
                this->grid[row][this->b_column] - boost::rational_cast<int64_t>(this->grid[row][this->b_column]) > max_fractional_b)
            {
                max_fractional_b_idx = row;
                max_fractional_b = this->grid[row][this->b_column] - boost::rational_cast<int64_t>(this->grid[row][this->b_column]);
            }
        }

        if (max_fractional_b != 0)
        {
            // Add new row.
            this->grid.push_back(std::vector<r64>());
            this->GRID_ROWS += 1;
            this->basis.push_back(this->GRID_ROWS + 1);
            for (size_t col = 0; col < this->GRID_COLS; ++col)
            {
                this->grid[GRID_ROWS - 1].push_back(-(this->grid[max_fractional_b_idx][col] - boost::rational_cast<int64_t>(this->grid[max_fractional_b_idx][col])));
            }
            this->Print();
            goto BNegativeRemoval;
        }

        std::cout << "Solution:" << std::endl;
        for (size_t col = 0; col < this->GRID_COLS - 1; ++col)
        {
            std::cout << "x_" << col + 1 << " = ";

            const auto it = std::find(this->basis.begin(), this->basis.end(), col);
            if (it != this->basis.end())
            {
                std::cout << this->grid[std::distance(this->basis.begin(), it)][this->b_column];
            }
            else
            {
                std::cout << 0;
            }

            std::cout << std::endl;
        }
    }

    const void Print()
    {
        std::cout << "Basis" << "\t";
        for (size_t col = 0; col < this->GRID_COLS - 1; ++col)
        {
            std::cout << "x_" << col + 1 << "\t";
        }
        std::cout << "B" << std::endl;
        for (size_t row = 0; row < this->GRID_ROWS; ++row)
        {
            std::cout << "x_" << this->basis[row] + 1 << "\t";
            for (size_t col = 0; col < this->GRID_COLS; ++col)
            {
                std::cout << this->grid[row][col] << "\t";
            }
            std::cout << std::endl;
        }
    }

private:
    std::vector<size_t> basis;

    std::vector<r64> z_function_coefficients;

    std::vector<std::vector<r64>> grid;
    size_t GRID_ROWS;
    size_t GRID_COLS;
    size_t b_column;
    size_t theta_column;

    const r64 GetCost(const size_t idx)
    {
        if (idx < this->z_function_coefficients.size())
        {
            return this->z_function_coefficients[idx];
        }
        else
        {
            return r64(0);
        }
    }

    void Transform(const size_t pivot_row, const size_t pivot_col)
    {
        const auto OLD_PIVOT = this->grid[pivot_row][pivot_col];

        auto new_grid = std::vector<std::vector<r64>>(this->GRID_ROWS, std::vector<r64>(this->GRID_COLS));

        for (size_t row = 0; row < this->GRID_ROWS; ++row)
        {
            for (size_t col = 0; col < this->GRID_COLS; ++col)
            {
                if (row == pivot_row)
                {
                    new_grid[row][col] = this->grid[row][col] / OLD_PIVOT;
                }
                else
                {
                    new_grid[row][col] = this->grid[row][col] - (this->grid[pivot_row][col] * this->grid[row][pivot_col]) / OLD_PIVOT;
                }
            }
        }

        this->grid = new_grid;
    }
};
