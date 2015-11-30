#include <string>
#include <stdexcept>
#include "Alignment.cpp"

class AlignmentTest {
public:
	void test1() {
		std::string string_a = "GCATGCU";
		std::string string_b = "GATTACA";

		Alignment alignment = Alignment(string_a, string_b);

		std::pair<std::string, std::string> result = alignment.local_alignment();

		std::string expected_a = "GCATG-CU";
		std::string expected_b = "G-ATTACA";

		if (result.first != expected_a || result.second != expected_b) {
			throw std::runtime_error("Test1 failed");
		}
	}
};