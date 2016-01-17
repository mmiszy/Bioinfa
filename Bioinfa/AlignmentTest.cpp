#include <string>
#include <stdexcept>
#include <vector>
#include <iostream>
#include "Alignment.h"

class AlignmentTest
{
    void testFailed(std::string testname)
    {
        std::cout << std::endl;
        throw std::runtime_error(testname + " failed");
    }

public:
    void test0()
    {
        std::string string_a = "AGAGTCAATCCATAG";
        std::string string_b = "CAGAGGTCCATCATG";

        std::vector<std::vector<int>> similarity = {
                //A   G   C   T   -
                {+2, +0, +0, +0, -1},
                {+0, +2, +0, +0, -1},
                {+0, +0, +2, +0, -1},
                {+0, +0, +0, +2, -1},
                {-1, -1, -1, -1, +2},
        };

        Alignment alignment = Alignment(string_a, string_b, Alignment::STR_WEIGHTS(1, -1, -1), similarity);

        std::pair<std::string, std::string> result = alignment.global_alignment();

        std::string expected_a = "-AGAG-TCAATCCATAG";
        std::string expected_b = "CAGAGGTCCATC-AT-G";
        std::cout << alignment.get_similarity();

        if (result.first != expected_a || result.second != expected_b)
        {
            alignment.debugS();
            this->testFailed("Test0");
        }
    }

    void test1()
    {
        std::string string_a = "GCATGCA";
        std::string string_b = "GATTACA";

        std::vector<std::vector<int>> similarity = {
                //A   G   C   T   -
                {+2, +0, +0, +0, -1},
                {+0, +2, +0, +0, -1},
                {+0, +0, +2, +0, -1},
                {+0, +0, +0, +2, -1},
                {-1, -1, -1, -1, +2},
        };

        Alignment alignment = Alignment(string_a, string_b, Alignment::STR_WEIGHTS(1, -1, -1), similarity,
                                        [](int j) { return -(j); });

        std::pair<std::string, std::string> result = alignment.global_alignment();

        std::string expected_a = "GCATG-CA";
        std::string expected_b = "G-ATTACA";
        std::cout << alignment.get_similarity();

        if (result.first != expected_a || result.second != expected_b)
        {
            alignment.debugS();
            this->testFailed("Test1");
        }
    }

    void test2()
    {
        std::string string_a = "CAGCCCTAC";
        std::string string_b = "CCTGTACCC";

        std::vector<std::vector<int>> similarity = {
                {+2, +0, +0, +0, -1},
                {+0, +2, +0, +0, -1},
                {+0, +0, +2, +0, -1},
                {+0, +0, +0, +2, -1},
                {-1, -1, -1, -1, +2}
        };

        Alignment alignment = Alignment(string_a, string_b, Alignment::STR_WEIGHTS(1, -1, -1), similarity,
                                        [](int j) { return -(1); });

        std::pair<std::string, std::string> result = alignment.global_alignment();


        std::string expected_a = "CAGCCCTAC--";
        std::string expected_b = "C--CTGTACCC";

        if (result.first != expected_a || result.second != expected_b)
        {
            alignment.debugS();
            this->testFailed("Test2");
        }
    }

    void test3()
    {
        std::string string_a = "AGACTAGTTAC";
        std::string string_b = "CGAGACGT";
        std::vector<std::vector<int>> similarity = {
                {+2, +0, +0, +0, -1},
                {+0, +2, +0, +0, -1},
                {+0, +0, +2, +0, -1},
                {+0, +0, +0, +2, -1},
                {-1, -1, -1, -1, +2}
        };

        Alignment alignment = Alignment(string_a, string_b, Alignment::STR_WEIGHTS(1, -1, -1), similarity,
                                        [](int j) { return 0; });

        std::pair<std::string, std::string> result = alignment.global_alignment();

        std::string expected_a = "--AGACTAGTTAC";
        std::string expected_b = "CGAGAC--GT---";

        if (result.first != expected_a || result.second != expected_b)
        {
            alignment.debugS();
            this->testFailed("Test3");
        }
    }

    void test4()
    {
        std::string string_a = "AGAGTCAATCCATAG";
        std::string string_b = "CAGAGGTCCATCATG";
        std::vector<std::vector<int>> similarity = {
                {+2, +0, +0, +0, -1},
                {+0, +2, +0, +0, -1},
                {+0, +0, +2, +0, -1},
                {+0, +0, +0, +2, -1},
                {-1, -1, -1, -1, +2}
        };

        Alignment alignment = Alignment(string_a, string_b, Alignment::STR_WEIGHTS(1, -1, -1), similarity,
                                        [](int j) { return -(j+2); });

        std::pair<std::string, std::string> result = alignment.global_alignment();

        std::string expected_a = "-AGAG-TCAATCCATAG";
        std::string expected_b = "CAGAGGTCCATC-AT-G";

        if (result.first != expected_a || result.second != expected_b)
        {
            alignment.debugS();
            this->testFailed("Test4");
        }
    }

    void test5()
    {
        std::string string_a = "GTATT";
        std::string string_b = "TTT";

        std::vector<std::vector<int>> similarity = {
                //A   G   C   T   -
                {+1, -1, -1, -1, -1},
                {-1, +1, -1, -1, -1},
                {-1, -1, +1, -1, -1},
                {-1, -1, -1, +1, -1},
                {-1, -1, -1, -1, +1}
        };

        Alignment alignment = Alignment(string_a, string_b, Alignment::STR_WEIGHTS(1, -1, -1), similarity,
                                        [](int j) { return -(j+2); });

        int result = alignment.get_similarity();

        int expected_result = -3;
        if (result != expected_result)
        {
            alignment.debugS();
            this->testFailed("Test5");
        }
    }
};