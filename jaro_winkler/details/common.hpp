/* SPDX-License-Identifier: MIT */
/* Copyright © 2022 Max Bachmann */

#pragma once
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace jaro_winkler {
namespace common {

/**
 * @defgroup Common Common
 * Common utilities shared among multiple functions
 * @{
 */

constexpr double result_cutoff(double result, double score_cutoff)
{
    return (result >= score_cutoff) ? result : 0;
}

constexpr double norm_distance(size_t dist, size_t lensum, double score_cutoff = 0)
{
    return result_cutoff(
        (lensum > 0) ? (1.0 - static_cast<double>(dist) / static_cast<double>(lensum)) : 1.0,
        score_cutoff);
}

static inline size_t score_cutoff_to_distance(double score_cutoff, size_t lensum)
{
    return static_cast<size_t>(std::ceil(static_cast<double>(lensum) * (1.0 - score_cutoff)));
}

template <typename T, typename U>
T ceildiv(T a, U divisor)
{
    return (T)(a / divisor) + (T)((a % divisor) != 0);
}

/**
 * @brief Finds the first mismatching pair of elements from two ranges:
 * one defined by [first1, last1) and another defined by [first2,last2).
 * Similar implementation to std::mismatch from C++14
 *
 * @param first1, last1	-	the first range of the elements
 * @param first2, last2	-	the second range of the elements
 *
 * @return std::pair with iterators to the first two non-equal elements.
 */
template <typename InputIt1, typename InputIt2>
std::pair<InputIt1, InputIt2> mismatch(InputIt1 first1, InputIt1 last1, InputIt2 first2,
                                       InputIt2 last2)
{
    while (first1 != last1 && first2 != last2 && *first1 == *first2) {
        ++first1;
        ++first2;
    }
    return std::pair<InputIt1, InputIt2>(first1, first2);
}

/**
 * Removes common prefix of two string views // todo
 */
template <typename InputIt1, typename InputIt2>
size_t remove_common_prefix(InputIt1& first1, InputIt1 last1, InputIt2& first2, InputIt2 last2)
{
    size_t prefix = static_cast<size_t>(
        std::distance(first1, common::mismatch(first1, last1, first2, last2).first));
    first1 += prefix;
    first2 += prefix;
    return prefix;
}

struct BitvectorHashmap {
    struct MapElem {
        uint64_t key = 0;
        uint64_t value = 0;
    };

    BitvectorHashmap() : m_map()
    {}

    template <typename CharT>
    void insert(CharT key, size_t pos)
    {
        insert_mask(key, 1ull << pos);
    }

    template <typename CharT>
    void insert_mask(CharT key, uint64_t mask)
    {
        size_t i = lookup((uint64_t)key);
        m_map[i].key = key;
        m_map[i].value |= mask;
    }

    template <typename CharT>
    uint64_t get(CharT key) const
    {
        return m_map[lookup((uint64_t)key)].value;
    }

private:
    /**
     * lookup key inside the hasmap using a similar collision resolution
     * strategy to CPython and Ruby
     */
    size_t lookup(uint64_t key) const
    {
        size_t i = key % 128;

        if (!m_map[i].value || m_map[i].key == key) {
            return i;
        }

        size_t perturb = key;
        while (true) {
            i = ((i * 5) + perturb + 1) % 128;
            if (!m_map[i].value || m_map[i].key == key) {
                return i;
            }

            perturb >>= 5;
        }
    }

    std::array<MapElem, 128> m_map;
};

struct PatternMatchVector {
    struct MapElem {
        uint64_t key = 0;
        uint64_t value = 0;
    };

    PatternMatchVector() : m_map(), m_extendedAscii()
    {}

    template <typename InputIt1>
    PatternMatchVector(InputIt1 first, InputIt1 last) : m_map(), m_extendedAscii()
    {
        insert(first, last);
    }

    template <typename InputIt1>
    void insert(InputIt1 first, InputIt1 last)
    {
        uint64_t mask = 1;
        for (size_t i = 0; i < std::distance(first, last); ++i) {
            auto key = first[i];
            if (key >= 0 && key <= 255) {
                m_extendedAscii[key] |= mask;
            }
            else {
                m_map.insert_mask(key, mask);
            }
            mask <<= 1;
        }
    }

    template <typename CharT>
    void insert(CharT key, size_t pos)
    {
        uint64_t mask = 1ull << pos;
        if (key >= 0 && key <= 255) {
            m_extendedAscii[key] |= mask;
        }
        else {
            m_map.insert_mask(key, mask);
        }
    }

    template <typename CharT>
    uint64_t get(CharT key) const
    {
        if (key >= 0 && key <= 255) {
            return m_extendedAscii[key];
        }
        else {
            return m_map.get(key);
        }
    }

private:
    BitvectorHashmap m_map;
    std::array<uint64_t, 256> m_extendedAscii;
};

struct BlockPatternMatchVector {
    BlockPatternMatchVector() : m_block_count(0)
    {}

    template <typename InputIt1>
    BlockPatternMatchVector(InputIt1 first, InputIt1 last) : m_block_count(0)
    {
        insert(first, last);
    }

    template <typename CharT>
    void insert(size_t block, CharT key, int pos)
    {
        uint64_t mask = 1ull << pos;

        assert(block < m_block_count);
        if (key >= 0 && key <= 255) {
            m_extendedAscii[key * m_block_count + block] |= mask;
        }
        else {
            m_map[block].insert_mask(key, mask);
        }
    }

    template <typename InputIt1>
    void insert(InputIt1 first, InputIt1 last)
    {
        size_t len = std::distance(first, last);
        m_block_count = ceildiv(len, 64);
        m_map.resize(m_block_count);
        m_extendedAscii.resize(m_block_count * 256);

        for (size_t i = 0; i < len; ++i) {
            size_t block = i / 64;
            size_t pos = i % 64;
            insert(block, first[i], pos);
        }
    }

    template <typename CharT>
    uint64_t get(size_t block, CharT key) const
    {
        assert(block < m_block_count);
        if (key >= 0 && key <= 255) {
            return m_extendedAscii[key * m_block_count + block];
        }
        else {
            return m_map[block].get(key);
        }
    }

private:
    std::vector<BitvectorHashmap> m_map;
    std::vector<uint64_t> m_extendedAscii;
    size_t m_block_count;
};

/**@}*/

} // namespace common
} // namespace jaro_winkler