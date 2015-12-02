#include <string>
#include <vector>
#include <algorithm>

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
    std::vector<std::vector<Alignment::DIRECTIONS>> track;

    void generate_grid() {
        int len1 = this->str_a.length() + 1;
        int len2 = this->str_b.length() + 1;

        this->grid = std::vector<std::vector<int>>(len1, std::vector<int>(len2));
        this->track = std::vector<std::vector<Alignment::DIRECTIONS>>(len1, std::vector<Alignment::DIRECTIONS>(len2));

        for (int row = 0; row < len1; ++row) {
            for (int col = 0; col < len2; ++col) {
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
                    diagonal = this->grid[row - 1][col - 1];
                    if (this->str_a[col - 1] == this->str_b[row - 1]) {
                        diagonal += this->WEIGHTS.MATCH;
                    }
                    else {
                        diagonal += this->WEIGHTS.MISMATCH;
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

public:
    Alignment(std::string str_a, std::string str_b) : Alignment(str_a, str_b, STR_WEIGHTS()) {}
    Alignment(std::string str_a, std::string str_b, STR_WEIGHTS w) {
        this->str_a = str_a;
        this->str_b = str_b;
        this->WEIGHTS = STR_WEIGHTS(w);
        this->generate_grid();
    }

    std::pair<std::string, std::string> global_alignment() {
        int len1 = this->str_a.length() + 1;
        int len2 = this->str_b.length() + 1;

        int row = len1 - 1;
        int col = len2 - 1;

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

    bool operator () (std::pair<int, int> i, std::pair<int, int> j) {
        return this->grid[i.first][i.second] < this->grid[j.first][j.second];
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
            } else if (this->grid[row][col] < 0) {
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
};