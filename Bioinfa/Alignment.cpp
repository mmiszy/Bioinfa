#include <string>
#include <vector>
#include <algorithm>
#include "Alignment.h"
#include <map>
#include <ostream>
#include <iostream>


std::vector<std::vector<int>> Alignment::get_default_similarity_matrix()
{
    const int SIMILARITY_MATRIX_SIZE = 5;
    std::vector<std::vector<int>> similarity = std::vector<std::vector<int>>(SIMILARITY_MATRIX_SIZE,
                                                                             std::vector<int>(
                                                                                     SIMILARITY_MATRIX_SIZE,
                                                                                     this->WEIGHTS.MISMATCH));

    for (int i = 0; i < SIMILARITY_MATRIX_SIZE; ++i)
    {
        similarity[i][i] = this->WEIGHTS.MATCH;
    }

    return similarity;
}

void Alignment::generate_alignment_grid(std::vector<std::vector<int>> similarity)
{
    unsigned long len_u = this->str_u.length() + 1;
    unsigned long len_w = this->str_w.length() + 1;

    // row is for iterating the first string, col for the second

//    this->grid = std::vector<std::vector<int>>(len_w, std::vector<int>(len_u));
    this->A = std::vector<std::vector<int>>(len_u, std::vector<int>(len_w, std::numeric_limits<int>::min()));
    this->B = std::vector<std::vector<int>>(len_u, std::vector<int>(len_w, std::numeric_limits<int>::min()));
    this->C = std::vector<std::vector<int>>(len_u, std::vector<int>(len_w, std::numeric_limits<int>::min()));
    this->S = std::vector<std::vector<int>>(len_u, std::vector<int>(len_w, std::numeric_limits<int>::min()));

    this->track = std::vector<std::vector<Alignment::DIRECTIONS>>(len_w, std::vector<Alignment::DIRECTIONS>(len_u));

    // A(i,j) = maksymalne dopasowanie przedrostków u[1..i] w[1..j], przy warunku: w[j] połączono z ‘–‘.
    // B(i,j) = maksymalne dopasowanie przedrostków u[1..i] w[1..j], przy warunku: u[i] połączono z ‘–‘.
    // C(i,j) = maksymalne dopasowanie przedrostków u[1..i] w[1..j], przy warunku: u[i] połączono z w[j].
    // S(i,j) = maksymalne dopasowanie u[1..i] z w[1..j].

    S[0][0] = 0;
    this->track[0][0] = Alignment::DIRECTIONS::DIAGONAL;

    for (int i = 1; i < len_u; ++i)
    {
        S[i][0] = B[i][0] = this->penaltyFn(i);
        this->track[0][i] = Alignment::DIRECTIONS::TOP;
    }
    for (int j = 1; j < len_w; ++j)
    {
        S[0][j] = A[0][j] = this->penaltyFn(j);
        this->track[j][0] = Alignment::DIRECTIONS::LEFT;
    }

    for (int j = 1; j < len_w; ++j)
    {
        for (int i = 1; i < len_u; ++i)
        {
            // col = i
            // row = j
            int left, top, diagonal;
            left = top = diagonal = std::numeric_limits<int>::min();
            {
                // A(i,j)=max k∈{0,j–1} { max{ B(i,k), C(i,k) } + p(j–k) }
                for (int k = 0; k <= (j-1); ++k) {
                    left = std::max(left, std::max(B[i][k], C[i][k]) + this->penaltyFn(j - k));
                }
//                int k = 0;
//                left = std::max(left, std::max(B[i][k], C[i][k]) + this->penaltyFn(j - k));
//                k = j - 1;
//                left = std::max(left, std::max(B[i][k], C[i][k]) + this->penaltyFn(j - k));

                A[i][j] = left;
            }
            {
                // B(i,j)=max k∈{0,i–1} { max { A(k,j), C(k,j) } + p(i–k) }
                for (int k = 0; k <= (i-1); ++k) {
                    top = std::max(top, std::max(A[k][j], C[k][j]) + this->penaltyFn(i - k));
                }
//                int k = 0;
//                top = std::max(top, std::max(A[k][j], C[k][j]) + this->penaltyFn(i - k));
//                k = i - 1;
//                top = std::max(top, std::max(A[k][j], C[k][j]) + this->penaltyFn(i - k));
                B[i][j] = top;
            }
            {
                // C(i,j)= S(i–1,j–1)+ s(u[i],w[j])
                char u = this->str_u[i-1];
                char w = this->str_w[j-1];

                diagonal = S[i - 1][j - 1] +
                           similarity[this->char_to_similarity_index[u]][this->char_to_similarity_index[w]];

                C[i][j] = diagonal;
            }
            {
                // S(i,j)=max{A(i,j), B(i,j), C(i,j)}
                int result = std::max({A[i][j], B[i][j], C[i][j]});
                S[i][j] = result;

                if (result == left)
                {
                    this->track[j][i] = Alignment::DIRECTIONS::LEFT;
                }
                else if (result == top)
                {
                    this->track[j][i] = Alignment::DIRECTIONS::TOP;
                }
                else
                {
                    this->track[j][i] = Alignment::DIRECTIONS::DIAGONAL;
                }
            }
        }
    }

//    for (int j = 0; j < len_w; ++j)
//    {
//        for (int i = 0; i < len_u; ++i)
//        {
//            int top, left, diagonal;
//            top = left = diagonal = INT32_MIN;
//
//            if (i != 0)
//            {
//                left = this->grid[j][i - 1] + this->WEIGHTS.INDEL;
//            }
//
//            if (j != 0)
//            {
//                top = this->grid[j - 1][i] + this->WEIGHTS.INDEL;
//            }
//
//            if (i == 0 && j == 0)
//            {
//                diagonal = 0;
//            }
//            else if (i != 0 && j != 0)
//            {
//                char a = this->str_u[i - 1];
//                char b = this->str_w[j - 1];
//
//                diagonal = this->grid[j - 1][i - 1];
//                if (this->char_to_similarity_index.find(a) == this->char_to_similarity_index.end() ||
//                    this->char_to_similarity_index.find(b) == this->char_to_similarity_index.end())
//                {
//                    diagonal += this->WEIGHTS.MISMATCH;
//                }
//                else
//                {
//                    diagonal += similarity[this->char_to_similarity_index[a]][this->char_to_similarity_index[b]];
//                }
//            }
//
//            int result = std::max({top, left, diagonal});
//
//
//            if (result == top)
//            {
//                this->track[j][i] = Alignment::DIRECTIONS::TOP;
//            }
//            else if (result == left)
//            {
//                this->track[j][i] = Alignment::DIRECTIONS::LEFT;
//            }
//            else
//            {
//                this->track[j][i] = Alignment::DIRECTIONS::DIAGONAL;
//            }
//
//            this->grid[j][i] = result;
//        }
//    }
}

Alignment::Alignment(std::string str_a, std::string str_b, STR_WEIGHTS w,
                     std::vector<std::vector<int>> similarity, penaltyFnType penaltyFn)
{
    this->str_u = str_a;
    this->str_w = str_b;
    this->WEIGHTS = w;
    this->char_to_similarity_index = std::map<char, int>({{'A', 0},
                                                          {'G', 1},
                                                          {'C', 2},
                                                          {'T', 3},
                                                          {' ', 4}});

    this->penaltyFn = penaltyFn;

    if (similarity.size() == 0)
    {
        similarity = this->get_default_similarity_matrix();
    }

    this->generate_alignment_grid(similarity);
}

std::pair<std::string, std::string> Alignment::global_alignment()
{
    int len_u = (int) (this->str_u.length() + 1);
    int len_w = (int) (this->str_w.length() + 1);

    int j = len_w - 1;
    int i = len_u - 1;

    std::vector<std::pair<int, int>> path = std::vector<std::pair<int, int>>();

    while (true)
    {
        if (j == 0 && i == 0)
        {
            break;
        }

        path.push_back(std::pair<int, int>(j, i));

        Alignment::DIRECTIONS dir = this->track[j][i];

        if (dir == Alignment::DIRECTIONS::DIAGONAL)
        {
            j -= 1;
            i -= 1;
        }
        else if (dir == Alignment::DIRECTIONS::LEFT)
        {
            j -= 1;
        }
        else if (dir == Alignment::DIRECTIONS::TOP)
        {
            i -= 1;
        }
    }

    std::reverse(path.begin(), path.end());

    std::string top = "";
    std::string side = "";

    for (auto p : path)
    {
        int row2 = p.first;
        int col2 = p.second;

        if (this->track[row2][col2] == Alignment::DIRECTIONS::DIAGONAL)
        {
            top += this->str_u[col2 - 1];
            side += this->str_w[row2 - 1];
        }
        else if (this->track[row2][col2] == Alignment::DIRECTIONS::TOP)
        {
            top += this->str_u[col2 - 1];
            side += Alignment::GAP_SYMBOL;
        }
        else
        {
            top += Alignment::GAP_SYMBOL;
            side += this->str_w[row2 - 1];
        }
    }

    return std::pair<std::string, std::string>(top, side);
}

int Alignment::get_similarity()
{
    return this->S.back().back();
}

void Alignment::debugS()
{
    std::cout << "\nS\n\t";
    unsigned long len_u = this->str_u.length() + 1;
    unsigned long len_w = this->str_w.length() + 1;

    const std::string track_to_symbol[] = {"↑", "←", "↖"};

    std::cout << "      ";

    for (int j = 1; j < len_w; ++j)
    {
        std::cout << "\t\t" << this->str_w[j - 1];
    }
    std::cout << "\n";

    std::cout << "";

    for (int j = 0; j < len_w; ++j)
    {
        std::cout << "\t\t" << j;
    }

    for (int i = 0; i < len_u; ++i)
    {
        std::cout << "\n";
        if (i > 0)
        {
            std::cout << this->str_u[i - 1] << "\t" << i << " ";
        }
        else
        {
            std::cout << "\t" << i << " ";
        }

        for (int j = 0; j < len_w; ++j)
        {
            int idx = (int) this->track[j][i];
            int s = this->S[i][j];
            if (s >= 0)
            {
                std::cout << " ";
            }
            std::cout << s << " " << track_to_symbol[idx] << "\t";
        }
    }
    std::cout << std::endl;
    std::cout << "similarity: " << this->get_similarity() << "\n";
    std::pair<std::string, std::string> result = this->global_alignment();
    std::cout << result.first << "\n";
    std::cout << result.second << "\n";
    std::cout << std::endl;
}

void Alignment::debugC()
{
    std::cout << "\nC\n\t";
    unsigned long len_u = this->str_u.length() + 1;
    unsigned long len_w = this->str_w.length() + 1;

    const std::string track_to_symbol[] = {"↑", "←", "↖"};

    std::cout << "      ";

    for (int j = 1; j < len_w; ++j)
    {
        std::cout << "\t\t" << this->str_w[j - 1];
    }
    std::cout << "\n";

    std::cout << "";

    for (int j = 0; j < len_w; ++j)
    {
        std::cout << "\t\t" << j;
    }

    for (int i = 0; i < len_u; ++i)
    {
        std::cout << "\n";
        if (i > 0)
        {
            std::cout << this->str_u[i - 1] << "\t" << i << " ";
        }
        else
        {
            std::cout << "\t" << i << " ";
        }

        for (int j = 0; j < len_w; ++j)
        {
            int s = this->C[i][j];
            std::cout << s << "\t";
        }
    }
    std::cout << std::endl;
    std::cout << std::endl;
}
