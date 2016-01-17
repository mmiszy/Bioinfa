//
// Created by Michal Miszczyszyn on 16/01/16.
//

#ifndef BIOINFA_ALIGNMENT_H
#define BIOINFA_ALIGNMENT_H

#include <string>
#include <vector>
#include <algorithm>
#include <map>

using penaltyFnType = std::function<int(int)>;

class Alignment
{
public:
    static struct STR_WEIGHTS
    {
    private:
        static const int DEFAULT_MATCH = 1;
        static const int DEFAULT_MISMATCH = -1;
        static const int DEFAULT_INDEL = -1;

    public:
        int MATCH;
        int MISMATCH;
        int INDEL;

        STR_WEIGHTS() : MATCH(STR_WEIGHTS::DEFAULT_MATCH), MISMATCH(STR_WEIGHTS::DEFAULT_MISMATCH),
                        INDEL(STR_WEIGHTS::DEFAULT_INDEL)
        {
        }

        STR_WEIGHTS(int match, int mismatch, int indel) : MATCH(match), MISMATCH(mismatch), INDEL(indel)
        {
        }
    };

    int penalty(int i);

private:
    static enum class DIRECTIONS
    {
        TOP, LEFT, DIAGONAL
    };

    static const char GAP_SYMBOL = '-';

    static const STR_WEIGHTS DEFAULT_WEIGHTS;

    STR_WEIGHTS WEIGHTS;

    std::string str_u;
    std::string str_w;

    std::vector<std::vector<int>> grid;
    std::vector<std::vector<int>> A;
    std::vector<std::vector<int>> B;
    std::vector<std::vector<int>> C;
    std::vector<std::vector<int>> S;
    std::vector<std::vector<Alignment::DIRECTIONS>> track;
    std::map<char, int> char_to_similarity_index;

    std::vector<std::vector<int>> get_default_similarity_matrix();

    void generate_alignment_grid(std::vector<std::vector<int>> similarity);

public:
    Alignment(std::string str_a, std::string str_b, STR_WEIGHTS w = STR_WEIGHTS(),
              std::vector<std::vector<int>> similarity = std::vector<std::vector<int>>(),
              penaltyFnType penaltyFn = [](int j){ return -(j + 2); });

    std::pair<std::string, std::string> global_alignment();

    int get_similarity();

    void debugC();

    void debugS();

    penaltyFnType penaltyFn;
};

#endif //BIOINFA_ALIGNMENT_H
