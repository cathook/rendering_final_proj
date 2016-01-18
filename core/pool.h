#ifndef PBRT_CORE_POOL_H_
#define PBRT_CORE_POOL_H_


#include <new>
#include <utility>

#include "pbrt.h"


template <typename T>
class NonDecreasePool {
public:
    template <typename... Args>
    size_t AddEntry(Args&& ...args) {
        size_t index = pool_.size();
        pool_.emplace_back(std::forward<Args>(args)...);
        return index;
    }

    T &GetEntry(size_t index) { return pool_[index]; }

    const T &GetEntry(size_t index) const { return pool_[index]; }

    T &operator[](size_t index) { return GetEntry(index); }

    const T &operator[](size_t index) const { return GetEntry(index); }

private:
    vector<T> pool_;
};


template <typename T>
class Pool {
public:
    ~Pool() {
        vector<bool> destructed(pool_.size(), false);
        for (size_t i = 0; i < recycled_indexes_.size(); ++i) {
            destructed[i] = true;
        }
        for (size_t i = 0; i < pool_.size(); ++i) {
            if (!destructed[i]) {
                pool_[i].Destruct();
            }
        }
    }

    template <typename... Args>
    size_t AddEntry(Args&& ...args) {
        size_t index;
        if (!recycled_indexes_.empty()) {
            index = recycled_indexes_.back();
            recycled_indexes_.pop_back();
        } else {
            index = pool_.size();
            pool_.resize(pool_.size() + 1);
        }
        pool_[index].Construct(std::forward<Args>(args)...);
        return index;
    }

    void RemoveEntry(size_t index) {
        pool_[index].Destruct();
        recycled_indexes_.push_back(index);
    }

    T &GetEntry(size_t index) { return pool_[index].value(); }

    const T &GetEntry(size_t index) const { return pool_[index].value(); }

    T &operator[](size_t index) { return GetEntry(index); }

    const T &operator[](size_t index) const { return GetEntry(index); }

private:
    class Entry {
    public:
        template <typename... Args>
        void Construct(Args&& ...args) {
            new (buf) T(std::forward<Args>(args)...);
        }

        void Destruct() { value().~T(); }

        T &value() { return *reinterpret_cast<T*>(buf); }

        const T &value() const { return *reinterpret_cast<const T*>(buf); }

    private:
        uint32_t buf[sizeof(T)];
    };

    vector<Entry> pool_;
    vector<size_t> recycled_indexes_;
};


#endif  // PBRT_CORE_POOL_H_
