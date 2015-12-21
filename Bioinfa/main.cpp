#include <string>
#include <iostream>
#include <vector>
//#include "AlignmentTest.cpp"
#include "Alignment.cpp"

int main() {
    /*
    AlignmentTest suite = AlignmentTest();
    suite.test1();
    suite.test2();
    suite.test3();
    suite.test4();
    suite.test5();
    */

    const int MATRIX_SIZE = 5;

    std::string sequence1, sequence2;
    std::cin >> sequence1 >> sequence2;

    /*
        similarity and distance matrix input format:
          A G C T
        A
        G
        C
        T
    */
    
    std::vector<std::vector<int>> similarity = std::vector<std::vector<int>>(5, std::vector<int>(5));
    std::vector<std::vector<int>> distance = std::vector<std::vector<int>>(5, std::vector<int>(5));

    for (unsigned int i = 0; i < MATRIX_SIZE; ++i) {
        for (unsigned int j = 0; j < MATRIX_SIZE; ++j) {
            std::cin >> similarity[i][j];
        }
    }

    for (unsigned int i = 0; i < MATRIX_SIZE; ++i) {
        for (unsigned int j = 0; j < MATRIX_SIZE; ++j) {
            std::cin >> distance[i][j];
        }
    }

    Alignment alignment = Alignment(sequence1, sequence2, Alignment::STR_WEIGHTS(1, -1, -1), similarity, distance);
    std::pair<std::string, std::string> global_alignment = alignment.global_alignment();
    std::pair<std::string, std::string> local_alignment = alignment.local_alignment();

    std::cout << "For strings:" << std::endl << sequence1 << std::endl << sequence2 << std::endl << std::endl;
    std::cout << "Best global alignment:" << std::endl << global_alignment.first << std::endl << global_alignment.second << std::endl << std::endl;
    std::cout << "Best local alignment:" << std::endl << local_alignment.first << std::endl << local_alignment.second << std::endl << std::endl;
    std::cout << "Edit distance: " << alignment.get_edit_distance() << std::endl;
    std::cout << "Similarity: " << alignment.get_similarity() << std::endl;

    return 0;
}