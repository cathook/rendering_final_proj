#ifndef PBRT_TREE_STYLE_ARRAY_H_
#define PBRT_TREE_STYLE_ARRAY_H_


#include <stdint.h>

#include <unordered_map>


class TreeStyleArray2D {
public:
    TreeStyleArray2D(uint32_t num_rows, size_t num_cols) :
            num_rows_(num_rows), num_cols_(num_cols) {}

    int Query(uint32_t row, uint32_t col) const {
        int sum = 0;
        for (uint32_t r = row; r > 0; r -= (r & -r)) {
            for (uint32_t c = col; c > 0; c -= (c & -c)) {
                auto it = table_.find(PackKeys_(r, c));
                sum += (it == table_.end() ? 0 : it->second);
            }
        }
        return sum;
    }

    void Update(uint32_t row, uint32_t col, int delta) {
        for (uint32_t r = row; r <= num_rows_; r += (r & -r)) {
            for (uint32_t c = col; c <= num_cols_; c += (c & -c)) {
                table_[PackKeys_(r, c)] += delta;
            }
        }
    }

private:
    uint64_t PackKeys_(uint32_t row, uint32_t col) const {
        return (static_cast<uint64_t>(row) << 32) + col;
    }

    uint32_t num_rows_;
    uint32_t num_cols_;
    std::unordered_map<uint64_t, int> table_;
};

#endif  // PBRT_TREE_STYLE_ARRAY_H_
