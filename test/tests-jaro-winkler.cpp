#include <bitset>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <jaro_winkler/jaro_winkler.hpp>

using Catch::Approx;

void validate_bitvector_word(const std::vector<int>& a, std::bitset<64> b)
{
    std::bitset<64> bit_a(0);
    for (size_t i = 0; i < a.size(); ++i)
    {
        bit_a[i] = a[i];
    }

    REQUIRE(bit_a == b);
}

struct FlaggedCharsOriginal {
    std::vector<int> P_flag;
    std::vector<int> T_flag;
    int64_t CommonChars;
};

template <typename CharT1, typename CharT2>
size_t get_jaro_bound(const std::basic_string<CharT1>& P, const std::basic_string<CharT2>& T)
{
    size_t Bound = std::max(P.size(), T.size()) / 2;
    if (Bound > 0) Bound--;
    return Bound;
}

template <typename CharT1, typename CharT2>
static inline FlaggedCharsOriginal flag_similar_characters_original(
    const std::basic_string<CharT1>& P, const std::basic_string<CharT2>& T)
{
    std::vector<int> P_flag(P.size() + 1);
    std::vector<int> T_flag(T.size() + 1);

    size_t Bound = get_jaro_bound(P, T);

    int64_t CommonChars = 0;
    for (size_t i = 0; i < T.size(); i++) {
        size_t lowlim = (i >= Bound) ? i - Bound : 0;
        size_t hilim = (i + Bound <= P.size() - 1) ? (i + Bound) : P.size() - 1;
        for (size_t j = lowlim; j <= hilim; j++) {
            if (!P_flag[j] && (P[j] == T[i])) {
                T_flag[i] = 1;
                P_flag[j] = 1;
                CommonChars++;
                break;
            }
        }
    }

    return {P_flag, T_flag, CommonChars};
}


template <typename CharT1, typename CharT2>
double jaro_similarity_original(const std::basic_string<CharT1>& P, const std::basic_string<CharT2>& T,
                                double score_cutoff)
{
    auto flagged = flag_similar_characters_original(P, T);

    // Count the number of transpositions
    int64_t Transpositions = 0;
    size_t k = 0;
    for (size_t i = 0; i < T.size(); i++) {
        if (flagged.T_flag[i]) {
            size_t j = k;
            for (; j < P.size(); j++) {
                if (flagged.P_flag[j]) {
                    k = j + 1;
                    break;
                }
            }
            if (T[i] != P[j]) {
                Transpositions++;
            }
        }
    }

    double sim = jaro_winkler::detail::jaro_calculate_similarity(
        static_cast<int64_t>(P.size()), static_cast<int64_t>(T.size()),
        flagged.CommonChars, Transpositions);
    return jaro_winkler::common::result_cutoff(sim, score_cutoff);
}

/**
 * @name JaroWinklerFlagCharsTest
 */
TEST_CASE("JaroWinklerTest")
{
    std::array<std::string, 19> names = {
        "james",
        "robert",
        "john",
        "michael",
        "william",
        "david",
        "joseph",
        "thomas",
        "charles",
        "mary",
        "patricia",
        "jennifer",
        "linda",
        "elizabeth",
        "barbara",
        "susan",
        "jessica",
        "sarah",
        "karen"
    };

    SECTION("testFlagCharsWord")
    {
        for (const auto& name1 : names)
        {
            jaro_winkler::common::PatternMatchVector PM(name1.begin(), name1.end());

            for (const auto& name2 : names)
            {
                auto P_first = name1.begin();
                auto P_last = name1.end();
                auto T_first = name2.begin();
                auto T_last = name2.end();

                int Bound = static_cast<int>(jaro_winkler::detail::jaro_bounds(P_first, P_last, T_first, T_last));
                auto flagged_original = flag_similar_characters_original(name1, name2);
                auto flagged_bitparallel = jaro_winkler::detail::flag_similar_characters_word(PM, P_first, P_last, T_first, T_last, Bound);

                INFO("Name1: " << name1 << ", Name2: " << name2);

                validate_bitvector_word(flagged_original.P_flag, flagged_bitparallel.P_flag);
                validate_bitvector_word(flagged_original.T_flag, flagged_bitparallel.T_flag);
                REQUIRE(flagged_original.CommonChars == jaro_winkler::detail::count_common_chars(flagged_bitparallel));
            }    
        }
    }

    // todo write test for blockwise implementation
    SECTION("testFlagCharsBlock")
    {
    }

    SECTION("testFullResult")
    {
        for (const auto& name1 : names)
        {
            for (const auto& name2 : names)
            {
                double Sim_original = jaro_similarity_original(name1, name2, 0);
                double Sim_bitparallel = jaro_winkler::detail::jaro_similarity(
                    name1.begin(), name1.end(), name2.begin(), name2.end(), 0);

                INFO("Name1: " << name1 << ", Name2: " << name2);
                REQUIRE(Sim_original == Approx(Sim_bitparallel));
            }    
        }
    }

    SECTION("testFullResultWithScoreCutoff")
    {
        for (const auto& name1 : names)
        {
            for (const auto& name2 : names)
            {
                double Sim_original = jaro_similarity_original(name1, name2, 0.9);
                double Sim_bitparallel = jaro_winkler::detail::jaro_similarity(
                    name1.begin(), name1.end(), name2.begin(), name2.end(), 0.9);

                INFO("Name1: " << name1 << ", Name2: " << name2);
                REQUIRE(Sim_original == Sim_bitparallel);
            }    
        }
    }

}