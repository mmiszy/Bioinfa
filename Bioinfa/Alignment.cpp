#include <string>
#include <vector>
#include <algorithm>
#include <map>

class Alignment {
public:
    static struct STR_WEIGHTS {
    private:
        static const int DEFAULT_MATCH = 1;
        static const int DEFAULT_MISMATCH = -1;
        static const int DEFAULT_INDEL = -1;

    public:
        int MATCH;
        int MISMATCH;
        int INDEL;
        STR_WEIGHTS() : MATCH(STR_WEIGHTS::DEFAULT_MATCH), MISMATCH(STR_WEIGHTS::DEFAULT_MISMATCH), INDEL(STR_WEIGHTS::DEFAULT_INDEL) {}
        STR_WEIGHTS(int match, int mismatch, int indel) : MATCH(match), MISMATCH(mismatch), INDEL(indel) {}
    };

private:
    static enum class DIRECTIONS { TOP, LEFT, DIAGONAL };

    static const char GAP_SYMBOL = '-';

    static const STR_WEIGHTS DEFAULT_WEIGHTS;

    STR_WEIGHTS WEIGHTS;

    std::string str_a;
    std::string str_b;

    std::vector<std::vector<int>> grid;
    std::vector<std::vector<int>> distanceGrid;
    std::vector<std::vector<Alignment::DIRECTIONS>> track;
    std::map<char, int> char_to_similarity_index;

    std::vector<std::vector<int>> get_default_similarity_matrix() {
        const int SIMILARITY_MATRIX_SIZE = 5;
        std::vector<std::vector<int>> similarity = std::vector<std::vector<int>>(SIMILARITY_MATRIX_SIZE,
            std::vector<int>(SIMILARITY_MATRIX_SIZE, this->WEIGHTS.MISMATCH));

        for (int i = 0; i < SIMILARITY_MATRIX_SIZE; ++i) {
            similarity[i][i] = this->WEIGHTS.MATCH;
        }

        return similarity;
    }

    std::vector<std::vector<int>> get_default_distance_matrix() {
        const int DISTANCE_MATRIX_SIZE = 5;
        const int DISTANCE_MATCH = 0;
        const int DISTANCE_MISMATCH = 1;

        std::vector<std::vector<int>> distance = std::vector<std::vector<int>>(DISTANCE_MATRIX_SIZE,
            std::vector<int>(DISTANCE_MATRIX_SIZE, DISTANCE_MISMATCH));

        for (int i = 0; i < DISTANCE_MATRIX_SIZE; ++i) {
            distance[i][i] = DISTANCE_MATCH;
        }

        return distance;
    }

    void generate_alignment_grid(std::vector<std::vector<int>> similarity) {
        int len1 = this->str_a.length() + 1;
        int len2 = this->str_b.length() + 1;

        //row is for iterating the first string, col for the second

        this->grid = std::vector<std::vector<int>>(len2, std::vector<int>(len1));
        this->track = std::vector<std::vector<Alignment::DIRECTIONS>>(len2, std::vector<Alignment::DIRECTIONS>(len1));

        for (int row = 0; row < len2; ++row) {
            for (int col = 0; col < len1; ++col) {
                int top, left, diagonal;
                top = left = diagonal = INT32_MIN;

                if (col != 0) {
                    left = this->grid[row][col - 1] + this->WEIGHTS.INDEL;
                }

                if (row != 0) {
                    top = this->grid[row - 1][col] + this->WEIGHTS.INDEL;
                }

                if (col == 0 && row == 0) {
                    diagonal = 0;
                }
                else if (col != 0 && row != 0) {
                    char a = this->str_a[col - 1];
                    char b = this->str_b[row - 1];

                    diagonal = this->grid[row - 1][col - 1];
                    if (this->char_to_similarity_index.find(a) == this->char_to_similarity_index.end() ||
                        this->char_to_similarity_index.find(b) == this->char_to_similarity_index.end()) {
                        diagonal += this->WEIGHTS.MISMATCH;
                    }
                    else {
                        diagonal += similarity[this->char_to_similarity_index[a]][this->char_to_similarity_index[b]];
                    }
                }

                int result = std::max({ top, left, diagonal });

                if (result == top) {
                    this->track[row][col] = Alignment::DIRECTIONS::TOP;
                }
                else if (result == left) {
                    this->track[row][col] = Alignment::DIRECTIONS::LEFT;
                }
                else {
                    this->track[row][col] = Alignment::DIRECTIONS::DIAGONAL;
                }

                this->grid[row][col] = result;
            }
        }
    }

    void generate_distance_grid(std::vector<std::vector<int>> distance) {
        int len1 = this->str_a.length() + 1;
        int len2 = this->str_b.length() + 1;

        this->distanceGrid = std::vector<std::vector<int>>(len2, std::vector<int>(len1));

        const int DISTANCE_INDEL = 1;

        for (int row = 0; row < len2; ++row) {
            for (int col = 0; col < len1; ++col) {
                int top, left, diagonal;
                top = left = diagonal = INT32_MAX;

                if (col != 0) {
                    left = this->distanceGrid[row][col - 1] + DISTANCE_INDEL;
                }

                if (row != 0) {
                    top = this->distanceGrid[row - 1][col] + DISTANCE_INDEL;
                }

                if (col == 0 && row == 0) {
                    diagonal = 0;
                }
                else if (col != 0 && row != 0) {
                    char a = this->str_a[col - 1];
                    char b = this->str_b[row - 1];

                    diagonal = this->distanceGrid[row - 1][col - 1];

                    if (this->char_to_similarity_index.find(a) == this->char_to_similarity_index.end() ||
                        this->char_to_similarity_index.find(b) == this->char_to_similarity_index.end()) {
                        diagonal += a == b ? 0 : 1;
                    } else {
                        diagonal += distance[this->char_to_similarity_index[a]][this->char_to_similarity_index[b]];
                    }
                }

                int result = std::min({ top, left, diagonal });
                this->distanceGrid[row][col] = result;
            }
        }
    }

public:
    Alignment(std::string str_a, std::string str_b, STR_WEIGHTS w = STR_WEIGHTS(), 
        std::vector<std::vector<int>> similarity = std::vector<std::vector<int>>(), 
        std::vector<std::vector<int>> distance = std::vector<std::vector<int>>()) {
        this->str_a = str_a;
        this->str_b = str_b;
        this->WEIGHTS = w;
        this->char_to_similarity_index = std::map<char, int>({ {'A', 0}, {'G', 1}, {'C', 2}, {'T', 3}, {' ', 4} });

        if (similarity.size() == 0) {
            similarity = this->get_default_similarity_matrix();
        }

        if (distance.size() == 0) {
            distance = this->get_default_distance_matrix();
        }

        this->generate_alignment_grid(similarity);
        this->generate_distance_grid(distance);
    }

    std::pair<std::string, std::string> global_alignment() {
        int len1 = this->str_a.length() + 1;
        int len2 = this->str_b.length() + 1;

        int row = len2 - 1;
        int col = len1 - 1;

        std::vector<std::pair<int, int>> path = std::vector<std::pair<int, int>>();

        while (true) {
            if (row == 0 && col == 0) {
                break;
            }

            path.push_back(std::pair<int, int>(row, col));

            Alignment::DIRECTIONS dir = this->track[row][col];

            if (dir == Alignment::DIRECTIONS::DIAGONAL) {
                row -= 1;
                col -= 1;
            }
            else if (dir == Alignment::DIRECTIONS::TOP) {
                row -= 1;
            }
            else if (dir == Alignment::DIRECTIONS::LEFT) {
                col -= 1;
            }
        }

        std::reverse(path.begin(), path.end());

        std::string top = "";
        std::string side = "";

        for (auto p : path) {
            int row = p.first;
            int col = p.second;

            if (this->track[row][col] == Alignment::DIRECTIONS::DIAGONAL) {
                top += this->str_a[col - 1];
                side += this->str_b[row - 1];
            }
            else if (this->track[row][col] == Alignment::DIRECTIONS::LEFT) {
                top += this->str_a[col - 1];
                side += Alignment::GAP_SYMBOL;
            }
            else {
                top += Alignment::GAP_SYMBOL;
                side += this->str_b[row - 1];
            }
        }

        return std::pair<std::string, std::string>(top, side);
    }

    std::pair<std::string, std::string> local_alignment() {
        int len1 = this->str_a.length() + 1;
        int len2 = this->str_b.length() + 1;

        std::vector<std::pair<int, int>> path = std::vector<std::pair<int, int>>();
        int row = 0;
        int col = 0;
        int max_elem = this->grid[0][0];

        for (int i = 0; i < this->grid.size(); ++i) {
            for (int j = 0; j < this->grid[i].size(); ++j) {
                if (this->grid[i][j] > max_elem) {
                    max_elem = this->grid[i][j];
                    row = i;
                    col = j;
                }
            }
        }

        while (true) {
            if (row <= 0 || col <= 0) {
                break;
            }
            else if (this->grid[row][col] < 0) {
                break;
            }

            path.push_back(std::pair<int, int>(row, col));

            Alignment::DIRECTIONS dir = this->track[row][col];

            if (dir == Alignment::DIRECTIONS::DIAGONAL) {
                row -= 1;
                col -= 1;
            }
            else if (dir == Alignment::DIRECTIONS::TOP) {
                row -= 1;
            }
            else if (dir == Alignment::DIRECTIONS::LEFT) {
                col -= 1;
            }
        }

        std::reverse(path.begin(), path.end());

        std::string top = "";
        std::string side = "";

        for (auto p : path) {
            int row = p.first;
            int col = p.second;

            if (this->track[row][col] == Alignment::DIRECTIONS::DIAGONAL) {
                top += this->str_a[col - 1];
                side += this->str_b[row - 1];
            }
            else if (this->track[row][col] == Alignment::DIRECTIONS::LEFT) {
                top += this->str_a[col - 1];
                side += Alignment::GAP_SYMBOL;
            }
            else {
                top += Alignment::GAP_SYMBOL;
                side += this->str_b[row - 1];
            }
        }

        return std::pair<std::string, std::string>(top, side);
    }

    int get_edit_distance() {
        return this->distanceGrid.back().back();
    }

    int get_similarity() {
        return this->grid.back().back();
    }
};