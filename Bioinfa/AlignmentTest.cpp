#include <string>
#include <stdexcept>
#include <vector>
#include "Alignment.cpp"

class AlignmentTest {
    void testFailed(std::string testname) {
        throw std::runtime_error(testname + " failed");
    }

public:
    void test1() {
        std::string string_a = "GCATGCU";
        std::string string_b = "GATTACA";

        Alignment alignment = Alignment(string_a, string_b);

        std::pair<std::string, std::string> result = alignment.global_alignment();

        std::string expected_a = "GCATG-CU";
        std::string expected_b = "G-ATTACA";

        if (result.first != expected_a || result.second != expected_b) {
            this->testFailed("Test1");
        }
    }

    void test2() {
        std::string string_a = "ACACACTA";
        std::string string_b = "AGCACACA";

        Alignment alignment = Alignment(string_a, string_b, Alignment::STR_WEIGHTS(2, -1, -1));

        std::pair<std::string, std::string> result = alignment.local_alignment();

        std::string expected_a = "A-CACACTA";
        std::string expected_b = "AGCACAC-A";

        if (result.first != expected_a || result.second != expected_b) {
            this->testFailed("Test2");
        }
    }

    void test3() {
        std::string string_a = "AGACTAGTTAC";
        std::string string_b = "CGAGACGT";
        std::vector<std::vector<int>> similarity = {
            {10, -1, -3, -4},
            {-1, 7, -5, -3},
            {-3, -5, 9, 0},
            {-4, -3, 0, 8}
        };

        Alignment alignment = Alignment(string_a, string_b, Alignment::STR_WEIGHTS(1, -1, -1), similarity);

        std::pair<std::string, std::string> result = alignment.global_alignment();

        std::string expected_a = "--AGACTAGTTAC";
        std::string expected_b = "CGAGAC--GT---";

        if (result.first != expected_a || result.second != expected_b) {
            this->testFailed("Test3");
        }
    }

    void test4() {
        std::string string_a = "CATTCG";
        std::string string_b = "TATTAGG";

        Alignment alignment = Alignment(string_a, string_b, Alignment::STR_WEIGHTS(1, -1, -1));

        int result = alignment.get_edit_distance();

        int expected_result = 3;

        if (result != expected_result) {
            this->testFailed("Test4");
        }
    }
};