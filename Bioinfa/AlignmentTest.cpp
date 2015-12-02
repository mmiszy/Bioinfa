#include <string>
#include <stdexcept>
#include "Alignment.cpp"

class AlignmentTest {
    void error(std::string testname) {
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
            this->error("Test1");
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
            this->error("Test2");
        }
    }
};