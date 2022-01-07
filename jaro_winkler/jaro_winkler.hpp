/* SPDX-License-Identifier: MIT */
/* Copyright © 2022 Max Bachmann */

#pragma once
#include <jaro_winkler/details/common.hpp>
#include <jaro_winkler/details/jaro_impl.hpp>

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace jaro_winkler {

/**
 * @defgroup jaro_winkler jaro_winkler
 * @{
 */

/**
 * @brief Calculates the jaro winkler similarity
 *
 * @tparam Sentence1 This is a string that can be converted to
 * basic_string_view<char_type>
 * @tparam Sentence2 This is a string that can be converted to
 * basic_string_view<char_type>
 *
 * @param s1
 *   string to compare with s2 (for type info check Template parameters above)
 * @param s2
 *   string to compare with s1 (for type info check Template parameters above)
 * @param prefix_weight
 *   Weight used for the common prefix of the two strings.
 *   Has to be between 0 and 0.25. Default is 0.1.
 * @param score_cutoff
 *   Optional argument for a score threshold as a float between 0 and 100.
 *   For ratio < score_cutoff 0 is returned instead. Default is 0,
 *   which deactivates this behaviour.
 *
 * @return jaro winkler similarity between s1 and s2
 *   as a float between 0 and 100
 */
template <typename InputIt1, typename InputIt2>
double jaro_winkler_similarity(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2,
                               double prefix_weight = 0.1, double score_cutoff = 0.0)
{
    if (prefix_weight < 0.0 || prefix_weight > 0.25) {
        throw std::invalid_argument("prefix_weight has to be between 0.0 and 0.25");
    }

    return detail::jaro_winkler_similarity(first1, last1, first2, last2, prefix_weight,
                                           score_cutoff);
}

template <typename InputIt1>
struct CachedJaroWinklerSimilarity {
    CachedJaroWinklerSimilarity(InputIt1 first1_, InputIt1 last1_, double prefix_weight_ = 0.1)
        : first1(first1_), last1(last1_), PM(first1, last1), prefix_weight(prefix_weight_)
    {
        if (prefix_weight < 0.0 || prefix_weight > 0.25) {
            throw std::invalid_argument("prefix_weight has to be between 0.0 and 0.25");
        }
    }

    template <typename InputIt2>
    double ratio(InputIt2 first2, InputIt2 last2, double score_cutoff = 0) const
    {
        return detail::jaro_winkler_similarity(PM, first1, last1, first2, last2, prefix_weight,
                                               score_cutoff);
    }

private:
    InputIt1 first1;
    InputIt1 last1;
    common::BlockPatternMatchVector PM;

    double prefix_weight;
};

/**
 * @brief Calculates the jaro similarity
 *
 * @tparam Sentence1 This is a string that can be converted to
 * basic_string_view<char_type>
 * @tparam Sentence2 This is a string that can be converted to
 * basic_string_view<char_type>
 *
 * @param s1
 *   string to compare with s2 (for type info check Template parameters above)
 * @param s2
 *   string to compare with s1 (for type info check Template parameters above)
 * @param score_cutoff
 *   Optional argument for a score threshold as a float between 0 and 100.
 *   For ratio < score_cutoff 0 is returned instead. Default is 0,
 *   which deactivates this behaviour.
 *
 * @return jaro similarity between s1 and s2
 *   as a float between 0 and 100
 */
template <typename InputIt1, typename InputIt2>
double jaro_similarity(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2,
                       double score_cutoff = 0.0)
{
    return detail::jaro_similarity(first1, last1, first2, last2, score_cutoff);
}

template <typename InputIt1>
struct CachedJaroSimilarity {
    CachedJaroSimilarity(InputIt1 first1_, InputIt1 last1_)
        : first1(first1_), last1(last1_), PM(first1, last1)
    {}

    template <typename InputIt2>
    double ratio(InputIt2 first2, InputIt2 last2, double score_cutoff = 0) const
    {
        return detail::jaro_similarity(PM, first1, last1, first2, last2, score_cutoff);
    }

private:
    InputIt1 first1;
    InputIt1 last1;
    common::BlockPatternMatchVector PM;
};

/**@}*/

} // namespace jaro_winkler
