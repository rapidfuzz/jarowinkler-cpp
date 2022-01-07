#include <bitset>

#include <catch2/catch_test_macros.hpp>
#include <jaro_winkler/jaro_winkler.hpp>

void validate_bitvector_word(const std::vector<int>& a, std::bitset<64> b)
{
    std::bitset<64> bit_a(0);
    for (int i = 0; i < a.size(); ++i)
    {
        bit_a[i] = a[i];
    }

    REQUIRE(bit_a == b);
}

struct FlaggedCharsOriginal {
    std::vector<int> P_flag;
    std::vector<int> T_flag;
    size_t CommonChars;
};

template <typename CharT1, typename CharT2>
static inline FlaggedCharsOriginal flag_similar_characters_original(jaro_winkler::basic_string_view<CharT1> P,
                                                                    jaro_winkler::basic_string_view<CharT2> T)
{
    std::vector<int> P_flag(P.size() + 1);
    std::vector<int> T_flag(T.size() + 1);

    size_t Bound = std::max(P.size(), T.size()) / 2;
    if (Bound > 0) Bound--;

    size_t CommonChars = 0;
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
double jaro_similarity_original(jaro_winkler::basic_string_view<CharT2> P, jaro_winkler::basic_string_view<CharT1> T,
                                double score_cutoff)
{
    auto flagged = flag_similar_characters_original(P, T);

    // Count the number of transpositions
    size_t Transpositions = 0;
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

    return jaro_winkler::common::result_cutoff(
        jaro_winkler::detail::jaro_calculate_similarity(P, T, flagged.CommonChars, Transpositions), score_cutoff);
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

    std::array<std::string, 1> long_string = {
        "00000000000000000000000000000000000000000000000000000000000000000"
    };

    SECTION("testFlagCharsWord")
    {
        for (const auto& name1 : names)
        {
            jaro_winkler::string_view P(name1);
            jaro_winkler::common::PatternMatchVector PM(P);

            for (const auto& name2 : names)
            {
                jaro_winkler::string_view T(name2);

                auto flagged_original = flag_similar_characters_original(P, T);
                auto flagged_bitparallel = jaro_winkler::detail::flag_similar_characters_word(PM, P, T);

                INFO("Name1: " << name1 << ", Name2: " << name2);

                validate_bitvector_word(flagged_original.P_flag, flagged_bitparallel.P_flag);
                validate_bitvector_word(flagged_original.T_flag, flagged_bitparallel.T_flag);
                REQUIRE(flagged_original.CommonChars == flagged_bitparallel.CommonChars);
            }    
        }
    }

    SECTION("testFlagCharsBlock")
    {
        for (const auto& name1 : names)
        {
            jaro_winkler::string_view P(name1);
            jaro_winkler::common::PatternMatchVector PM(P);

            for (const auto& name2 : names)
            {
                jaro_winkler::string_view T(name2);

                auto flagged_original = flag_similar_characters_original(P, T);
                auto flagged_bitparallel = jaro_winkler::detail::flag_similar_characters_word(PM, P, T);

                INFO("Name1: " << name1 << ", Name2: " << name2);

                validate_bitvector_word(flagged_original.P_flag, flagged_bitparallel.P_flag);
                validate_bitvector_word(flagged_original.T_flag, flagged_bitparallel.T_flag);
                REQUIRE(flagged_original.CommonChars == flagged_bitparallel.CommonChars);
            }    
        }
    }

    SECTION("testFullResult")
    {
        for (const auto& name1 : names)
        {
            jaro_winkler::string_view P(name1);

            for (const auto& name2 : names)
            {
                jaro_winkler::string_view T(name2);

                double Sim_original = jaro_similarity_original(P, T, 0);
                double Sim_bitparallel = jaro_winkler::detail::jaro_similarity_word(P, T, 0);

                INFO("Name1: " << name1 << ", Name2: " << name2);
                REQUIRE(Sim_original == Sim_bitparallel);
            }    
        }
    }

    SECTION("testFullResultWithScoreCutoff")
    {
        for (const auto& name1 : names)
        {
            jaro_winkler::string_view P(name1);

            for (const auto& name2 : names)
            {
                jaro_winkler::string_view T(name2);

                double Sim_original = jaro_similarity_original(P, T, 0.9);
                double Sim_bitparallel = jaro_winkler::detail::jaro_similarity_word(P, T, 0.9);

                INFO("Name1: " << name1 << ", Name2: " << name2);
                REQUIRE(Sim_original == Sim_bitparallel);
            }    
        }
    }

}