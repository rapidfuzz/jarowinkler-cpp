/* SPDX-License-Identifier: MIT */
/* Copyright © 2022 Max Bachmann */

#include <jaro_winkler/details/common.hpp>
#include <jaro_winkler/details/intrinsics.hpp>

namespace jaro_winkler {
namespace detail {

struct FlaggedCharsWord {
    uint64_t P_flag;
    uint64_t T_flag;
};

struct FlaggedCharsMultiword {
    std::vector<uint64_t> P_flag;
    std::vector<uint64_t> T_flag;
};

struct SearchBoundMask {
    size_t words = 0;
    size_t empty_words = 0;
    uint64_t last_mask = 0;
    uint64_t first_mask = 0;
};

struct TextPosition {
    TextPosition(size_t Word_, size_t WordPos_) : Word(Word_), WordPos(WordPos_)
    {}
    size_t Word;
    size_t WordPos;
};

static inline double jaro_calculate_similarity(size_t P_len, size_t T_len, size_t CommonChars,
                                               size_t Transpositions)
{
    Transpositions /= 2;
    double Sim = (double)CommonChars / (double)P_len + (double)CommonChars / (double)T_len +
                 (double)(CommonChars - Transpositions) / (double)CommonChars;
    return Sim / 3.0;
}

/**
 * @brief filter matches below score_cutoff based on string lengths
 */
static inline bool jaro_length_filter(size_t P_len, size_t T_len, double score_cutoff)
{
    if (!T_len || !P_len) return false;

    double min_len = (double)std::min(P_len, T_len);
    double Sim = (double)min_len / (double)P_len + (double)min_len / (double)T_len + 1.0;
    Sim /= 3.0;
    return Sim >= score_cutoff;
}

/**
 * @brief filter matches below score_cutoff based on string lengths and common characters
 */
static inline bool jaro_common_char_filter(size_t P_len, size_t T_len, size_t CommonChars,
                                           double score_cutoff)
{
    if (!CommonChars) return false;

    double Sim = (double)CommonChars / (double)P_len + (double)CommonChars / (double)T_len + 1.0;
    Sim /= 3.0;
    return Sim >= score_cutoff;
}

static inline size_t count_common_chars(const FlaggedCharsWord& flagged)
{
    return intrinsics::popcount64(flagged.P_flag);
}

static inline size_t count_common_chars(const FlaggedCharsMultiword& flagged)
{
    size_t CommonChars = 0;
    if (flagged.P_flag.size() < flagged.T_flag.size()) {
        for (uint64_t flag : flagged.P_flag) {
            CommonChars += intrinsics::popcount64(flag);
        }
    }
    else {
        for (uint64_t flag : flagged.T_flag) {
            CommonChars += intrinsics::popcount64(flag);
        }
    }
    return CommonChars;
}

template <typename InputIt1, typename InputIt2>
static inline FlaggedCharsWord
flag_similar_characters_word(const common::PatternMatchVector& PM, InputIt1 P_first,
                             InputIt1 P_last, InputIt2 T_first, InputIt2 T_last, size_t Bound)
{
    using namespace intrinsics;
    size_t P_len = std::distance(P_first, P_last);
    size_t T_len = std::distance(T_first, T_last);
    assert(P_len <= 64);
    assert(T_len <= 64);
    assert(Bound > P_len || P_len - Bound <= T_len);

    FlaggedCharsWord flagged = {0, 0};

    uint64_t BoundMask = bit_mask_lsb<uint64_t>(Bound + 1);

    size_t j = 0;
    for (; j < std::min(Bound, T_len); ++j) {
        uint64_t PM_j = PM.get(T_first[j]) & BoundMask & (~flagged.P_flag);

        flagged.P_flag |= blsi(PM_j);
        flagged.T_flag |= (uint64_t)(PM_j != 0) << j;

        BoundMask = (BoundMask << 1) | 1;
    }

    for (; j < T_len; ++j) {
        uint64_t PM_j = PM.get(T_first[j]) & BoundMask & (~flagged.P_flag);

        flagged.P_flag |= blsi(PM_j);
        flagged.T_flag |= (uint64_t)(PM_j != 0) << j;

        BoundMask <<= 1;
    }

    return flagged;
}

template <typename InputIt1, typename InputIt2>
static inline FlaggedCharsWord
flag_similar_characters_word(const common::BlockPatternMatchVector& PM, InputIt1 P_first,
                             InputIt1 P_last, InputIt2 T_first, InputIt2 T_last, size_t Bound)
{
    using namespace intrinsics;
    size_t P_len = std::distance(P_first, P_last);
    size_t T_len = std::distance(T_first, T_last);
    assert(P_len <= 64);
    assert(T_len <= 64);
    assert(Bound > P_len || P_len - Bound <= T_len);

    FlaggedCharsWord flagged = {0, 0};

    uint64_t BoundMask = bit_mask_lsb<uint64_t>(Bound + 1);

    size_t j = 0;
    for (; j < std::min(Bound, T_len); ++j) {
        uint64_t PM_j = PM.get(0, T_first[j]) & BoundMask & (~flagged.P_flag);

        flagged.P_flag |= blsi(PM_j);
        flagged.T_flag |= (uint64_t)(PM_j != 0) << j;

        BoundMask = (BoundMask << 1) | 1;
    }

    for (; j < T_len; ++j) {
        uint64_t PM_j = PM.get(0, T_first[j]) & BoundMask & (~flagged.P_flag);

        flagged.P_flag |= blsi(PM_j);
        flagged.T_flag |= (uint64_t)(PM_j != 0) << j;

        BoundMask <<= 1;
    }

    return flagged;
}

template <typename CharT>
static inline void flag_similar_characters_step(const common::BlockPatternMatchVector& PM,
                                                CharT T_j, FlaggedCharsMultiword& flagged, size_t j,
                                                SearchBoundMask BoundMask)
{
    using namespace intrinsics;

    size_t j_word = j / 64;
    size_t j_pos = j % 64;
    size_t word = BoundMask.empty_words;
    size_t last_word = word + BoundMask.words;

    if (BoundMask.words == 1) {
        uint64_t PM_j = PM.get(word, T_j) & BoundMask.last_mask & BoundMask.first_mask &
                        (~flagged.P_flag[word]);

        flagged.P_flag[word] |= blsi(PM_j);
        flagged.T_flag[j_word] |= (uint64_t)(PM_j != 0) << j_pos;
        return;
    }

    if (BoundMask.first_mask) {
        uint64_t PM_j = PM.get(word, T_j) & BoundMask.first_mask & (~flagged.P_flag[word]);

        if (PM_j) {
            flagged.P_flag[word] |= blsi(PM_j);
            flagged.T_flag[j_word] |= 1ull << j_pos;
            return;
        }
        word++;
    }

    for (; word < last_word - 1; ++word) {
        uint64_t PM_j = PM.get(word, T_j) & (~flagged.P_flag[word]);

        if (PM_j) {
            flagged.P_flag[word] |= blsi(PM_j);
            flagged.T_flag[j_word] |= 1ull << j_pos;
            return;
        }
    }

    if (BoundMask.last_mask) {
        uint64_t PM_j = PM.get(word, T_j) & BoundMask.last_mask & (~flagged.P_flag[word]);

        flagged.P_flag[word] |= blsi(PM_j);
        flagged.T_flag[j_word] |= (uint64_t)(PM_j != 0) << j_pos;
    }
};

template <typename InputIt1, typename InputIt2>
static inline FlaggedCharsMultiword
flag_similar_characters_block(const common::BlockPatternMatchVector& PM, InputIt1 P_first,
                              InputIt1 P_last, InputIt2 T_first, InputIt2 T_last, size_t Bound)
{
    using namespace intrinsics;
    size_t P_len = std::distance(P_first, P_last);
    size_t T_len = std::distance(T_first, T_last);
    assert(P_len > 64 || T_len > 64);
    assert(Bound > P_len || P_len - Bound <= T_len);
    assert(Bound >= 31);

    const size_t TextWords = common::ceildiv(T_len, 64);
    const size_t PatternWords = common::ceildiv(P_len, 64);
    size_t FullBoundSize = 2 * Bound + 1;

    FlaggedCharsMultiword flagged;
    flagged.T_flag.resize(TextWords);
    flagged.P_flag.resize(PatternWords);

    SearchBoundMask BoundMask;
    size_t start_range = std::min(Bound + 1, P_len);
    BoundMask.words = common::ceildiv(start_range, 64);
    BoundMask.empty_words = 0;
    BoundMask.last_mask = (1ull << (start_range % 64)) - 1;
    BoundMask.first_mask = (uint64_t)-1;

    for (size_t j = 0; j < T_len; ++j) {
        size_t j_word = j / 64;
        size_t j_pos = j % 64;

        flag_similar_characters_step(PM, T_first[j], flagged, j, BoundMask);

        if (j + Bound + 1 < P_len) {
            BoundMask.last_mask = (BoundMask.last_mask << 1) | 1;
            if (j + Bound + 2 < P_len && BoundMask.last_mask == (uint64_t)-1) {
                BoundMask.last_mask = 0;
                BoundMask.words++;
            }
        }

        if (j >= Bound) {
            BoundMask.first_mask <<= 1;
            if (BoundMask.first_mask == 0) {
                BoundMask.first_mask = (uint64_t)-1;
                BoundMask.words--;
                BoundMask.empty_words++;
            }
        }
    }

    return flagged;
}

template <typename InputIt1>
static inline size_t count_transpositions_word(const common::PatternMatchVector& PM,
                                               InputIt1 T_first, InputIt1 T_last,
                                               const FlaggedCharsWord& flagged)
{
    using namespace intrinsics;
    uint64_t P_flag = flagged.P_flag;
    uint64_t T_flag = flagged.T_flag;
    size_t Transpositions = 0;
    while (T_flag) {
        uint64_t PatternFlagMask = blsi(P_flag);

        Transpositions += !(PM.get(T_first[tzcnt(T_flag)]) & PatternFlagMask);

        T_flag = blsr(T_flag);
        P_flag ^= PatternFlagMask;
    }

    return Transpositions;
}

template <typename InputIt1>
static inline size_t count_transpositions_word(const common::BlockPatternMatchVector& PM,
                                               InputIt1 T_first, InputIt1 T_last,
                                               const FlaggedCharsWord& flagged)
{
    using namespace intrinsics;
    uint64_t P_flag = flagged.P_flag;
    uint64_t T_flag = flagged.T_flag;
    size_t Transpositions = 0;
    while (T_flag) {
        uint64_t PatternFlagMask = blsi(P_flag);

        Transpositions += !(PM.get(0, T_first[tzcnt(T_flag)]) & PatternFlagMask);

        T_flag = blsr(T_flag);
        P_flag ^= PatternFlagMask;
    }

    return Transpositions;
}

template <typename InputIt1>
static inline size_t count_transpositions_block(const common::BlockPatternMatchVector& PM,
                                                InputIt1 T_first, InputIt1 T_last,
                                                const FlaggedCharsMultiword& flagged,
                                                size_t FlaggedChars)
{
    using namespace intrinsics;
    size_t TextWord = 0;
    size_t PatternWord = 0;
    uint64_t T_flag = flagged.T_flag[TextWord];
    uint64_t P_flag = flagged.P_flag[PatternWord];

    size_t Transpositions = 0;
    while (FlaggedChars) {
        while (!T_flag) {
            TextWord++;
            T_first += 64;
            T_flag = flagged.T_flag[TextWord];
        }

        while (T_flag) {
            while (!P_flag) {
                PatternWord++;
                P_flag = flagged.P_flag[PatternWord];
            }

            uint64_t PatternFlagMask = blsi(P_flag);

            Transpositions += !(PM.get(PatternWord, T_first[tzcnt(T_flag)]) & PatternFlagMask);

            T_flag = blsr(T_flag);
            P_flag ^= PatternFlagMask;

            FlaggedChars--;
        }
    }

    return Transpositions;
}

template <typename InputIt1, typename InputIt2>
double jaro_similarity(InputIt1 P_first, InputIt1 P_last, InputIt2 T_first, InputIt2 T_last,
                       double score_cutoff)
{
    size_t P_len = std::distance(P_first, P_last);
    size_t T_len = std::distance(T_first, T_last);

    /* filter out based on the length difference between the two strings */
    if (!jaro_length_filter(P_len, T_len, score_cutoff)) {
        return 0.0;
    }

    if (P_len == 1 && T_len == 1) {
        return (double)(P_first[0] == T_first[0]);
    }

    /* since jaro uses a sliding window some parts of T/P might never be in
     * range an can be removed ahead of time
     */
    size_t Bound = 0;
    if (T_len > P_len) {
        Bound = T_len / 2 - 1;
        if (T_len > P_len + Bound) {
            T_last = T_first + P_len + Bound;
        }
    }
    else {
        Bound = P_len / 2 - 1;
        if (P_len > T_len + Bound) {
            P_last = P_first + T_len + Bound;
        }
    }

    /* common prefix never includes Transpositions */
    size_t CommonChars = common::remove_common_prefix(P_first, P_last, T_first, T_last);
    size_t Transpositions = 0;
    size_t P_view_len = std::distance(P_first, P_last);
    size_t T_view_len = std::distance(T_first, T_last);

    if (!P_view_len || !T_view_len) {
        /* already has correct number of common chars and transpositions */
    }
    else if (P_view_len <= 64 && T_view_len <= 64) {
        common::PatternMatchVector PM(P_first, P_last);
        auto flagged = flag_similar_characters_word(PM, P_first, P_last, T_first, T_last, Bound);
        CommonChars += count_common_chars(flagged);

        if (!jaro_common_char_filter(P_len, T_len, CommonChars, score_cutoff)) {
            return 0.0;
        }

        Transpositions = count_transpositions_word(PM, T_first, T_last, flagged);
    }
    else {
        common::BlockPatternMatchVector PM(P_first, P_last);
        auto flagged = flag_similar_characters_block(PM, P_first, P_last, T_first, T_last, Bound);
        size_t FlaggedChars = count_common_chars(flagged);
        CommonChars += FlaggedChars;

        if (!jaro_common_char_filter(P_len, T_len, CommonChars, score_cutoff)) {
            return 0.0;
        }

        Transpositions = count_transpositions_block(PM, T_first, T_last, flagged, FlaggedChars);
    }

    double Sim = jaro_calculate_similarity(P_len, T_len, CommonChars, Transpositions);
    return common::result_cutoff(Sim, score_cutoff);
}

template <typename InputIt1, typename InputIt2>
double jaro_similarity(const common::BlockPatternMatchVector& PM, InputIt1 P_first, InputIt1 P_last,
                       InputIt2 T_first, InputIt2 T_last, double score_cutoff)
{
    size_t P_len = std::distance(P_first, P_last);
    size_t T_len = std::distance(T_first, T_last);

    /* filter out based on the length difference between the two strings */
    if (!jaro_length_filter(P_len, T_len, score_cutoff)) {
        return 0.0;
    }

    if (P_len == 1 && T_len == 1) {
        return (double)(P_first[0] == T_first[0]);
    }

    /* since jaro uses a sliding window some parts of T/P might never be in
     * range an can be removed ahead of time
     */
    size_t Bound = 0;
    if (T_len > P_len) {
        Bound = T_len / 2 - 1;
        if (T_len > P_len + Bound) {
            T_last = T_first + P_len + Bound;
        }
    }
    else {
        Bound = P_len / 2 - 1;
        if (P_len > T_len + Bound) {
            P_last = P_first + T_len + Bound;
        }
    }

    /* common prefix never includes Transpositions */
    size_t CommonChars = 0;
    size_t Transpositions = 0;
    size_t P_view_len = std::distance(P_first, P_last);
    size_t T_view_len = std::distance(T_first, T_last);

    if (!P_view_len || !T_view_len) {
        /* already has correct number of common chars and transpositions */
    }
    else if (P_view_len <= 64 && T_view_len <= 64) {
        auto flagged = flag_similar_characters_word(PM, P_first, P_last, T_first, T_last, Bound);
        CommonChars += count_common_chars(flagged);

        if (!jaro_common_char_filter(P_len, T_len, CommonChars, score_cutoff)) {
            return 0.0;
        }

        Transpositions = count_transpositions_word(PM, T_first, T_last, flagged);
    }
    else {
        auto flagged = flag_similar_characters_block(PM, P_first, P_last, T_first, T_last, Bound);
        size_t FlaggedChars = count_common_chars(flagged);
        CommonChars += FlaggedChars;

        if (!jaro_common_char_filter(P_len, T_len, CommonChars, score_cutoff)) {
            return 0.0;
        }

        Transpositions = count_transpositions_block(PM, T_first, T_last, flagged, FlaggedChars);
    }

    double Sim = jaro_calculate_similarity(P_len, T_len, CommonChars, Transpositions);
    return common::result_cutoff(Sim, score_cutoff);
}

template <typename InputIt1, typename InputIt2>
double jaro_winkler_similarity(InputIt1 P_first, InputIt1 P_last, InputIt2 T_first, InputIt2 T_last,
                               double prefix_weight, double score_cutoff)
{
    size_t P_len = std::distance(P_first, P_last);
    size_t T_len = std::distance(T_first, T_last);
    size_t min_len = std::min(P_len, T_len);
    size_t prefix = 0;
    size_t max_prefix = (min_len >= 4) ? 4 : min_len;

    for (; prefix < max_prefix; ++prefix) {
        if (T_first[prefix] != P_first[prefix]) {
            break;
        }
    }

    double jaro_score_cutoff = score_cutoff;
    if (jaro_score_cutoff > 0.7) {
        double prefix_sim = prefix * prefix_weight;

        if (prefix_sim >= 1.0) {
            jaro_score_cutoff = 0.7;
        }
        else {
            jaro_score_cutoff =
                std::max(0.7, (prefix_sim - jaro_score_cutoff) / (prefix_sim - 1.0));
        }
    }

    double Sim = jaro_similarity(P_first, P_last, T_first, T_last, jaro_score_cutoff);
    if (Sim > 0.7) {
        Sim += prefix * prefix_weight * (1.0 - Sim);
    }

    return common::result_cutoff(Sim, score_cutoff);
}

template <typename InputIt1, typename InputIt2>
double jaro_winkler_similarity(const common::BlockPatternMatchVector& PM, InputIt1 P_first,
                               InputIt1 P_last, InputIt2 T_first, InputIt2 T_last,
                               double prefix_weight, double score_cutoff)
{
    size_t P_len = std::distance(P_first, P_last);
    size_t T_len = std::distance(T_first, T_last);
    size_t min_len = std::min(P_len, T_len);
    size_t prefix = 0;
    size_t max_prefix = (min_len >= 4) ? 4 : min_len;

    for (; prefix < max_prefix; ++prefix) {
        if (T_first[prefix] != P_first[prefix]) {
            break;
        }
    }

    double jaro_score_cutoff = score_cutoff;
    if (jaro_score_cutoff > 0.7) {
        double prefix_sim = prefix * prefix_weight;

        if (prefix_sim >= 1.0) {
            jaro_score_cutoff = 0.7;
        }
        else {
            jaro_score_cutoff =
                std::max(0.7, (prefix_sim - jaro_score_cutoff) / (prefix_sim - 1.0));
        }
    }

    double Sim = jaro_similarity(PM, P_first, P_last, T_first, T_last, jaro_score_cutoff);
    if (Sim > 0.7) {
        Sim += prefix * prefix_weight * (1.0 - Sim);
    }

    return common::result_cutoff(Sim, score_cutoff);
}

} // namespace detail
} // namespace jaro_winkler