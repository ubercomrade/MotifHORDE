#include "sitega_train.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr int kDinucleotideCount = 16;
constexpr int kInvalidSlot = 16;
constexpr int kPrefixSlots = 17;
constexpr int kMaxLpdLen = 6;
constexpr int kMaxMotifLen = 20;
constexpr int kMaxFeatureCount = 120;
constexpr int kMaxPopulation = 500;
constexpr int kDefaultPopulation = 100;
constexpr int kDefaultNumMotifs = 20;
constexpr double kMinPivot = 1e-12;
constexpr double kScoreEpsilon = 1e-12;
constexpr double kFeatureScoreRidge = 1e-9;
constexpr int kDirectedFeatureAttempts = 64;
constexpr int kFeaturePoolTopMultiplier = 12;
constexpr int kPlacementRefineSequenceLimit = 32;
constexpr int kPlacementRefineSamples = 12;

using PrefixRow = std::array<int, kPrefixSlots>;
using PrefixTable = std::vector<PrefixRow>;

struct Rng {
    uint64_t state = 0x9E3779B97F4A7C15ULL;

    explicit Rng(uint64_t seed) {
        reseed(seed);
    }

    void reseed(uint64_t seed) {
        state = seed == 0 ? 0x9E3779B97F4A7C15ULL : seed;
    }

    uint32_t next() {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        return static_cast<uint32_t>(
            (state * 0x2545F4914F6CDD1DULL) >> 32
        );
    }

    int range(int limit) {
        if (limit <= 0) {
            throw std::runtime_error("RNG range limit must be positive");
        }
        return static_cast<int>(next() % static_cast<uint32_t>(limit));
    }

    double unit() {
        return static_cast<double>(next()) /
            static_cast<double>(std::numeric_limits<uint32_t>::max());
    }
};

struct SearchConfig {
    std::string fg_path;
    std::string bg_path;
    std::string out_path;
    std::string log_path;
    std::string fg_file;
    std::string model_name;
    int max_lpd = 0;
    int motif_len = 0;
    int feature_count = 0;
    int olig_bg = 0;
    int max_peak_len = 0;
    int pop_size = 0;
    int num_motifs = 0;
    int generations = 0;
    int mutation_attempts = 0;
    int stale_generations = 0;
    bool verbose = false;
    uint64_t seed = 0;
};

struct EncodedSequence {
    std::string forward;
    std::string reverse;
    std::array<std::vector<int>, 2> dinucs;
    std::array<PrefixTable, 2> prefixes;
    std::vector<double> window_weights;
    std::vector<int> candidate_positions;
    std::vector<int> cumulative_weights;
};

struct EncodedSequences {
    std::vector<EncodedSequence> records;
};

struct BackgroundStats {
    std::vector<std::array<double, kDinucleotideCount>> mean;
    std::vector<std::array<double, kDinucleotideCount>> covariance;
    std::array<double, kDinucleotideCount> expected{};
};

struct Feature {
    int start = 0;
    int end = 0;
    int code = 0;
};

struct FeaturePoolEntry {
    Feature feature;
    double score = 0.0;
    double signed_effect = 0.0;
    double fg_mean = 0.0;
    double bg_mean = 0.0;
    double variance = 0.0;
};

struct Candidate {
    std::vector<Feature> features;
    std::vector<int> positions;
    std::vector<unsigned char> orientations;
    std::vector<double> feature_sum;
    std::vector<double> second_moment;
    int invalid_intervals = 0;
    double kmer_sum = 0.0;
    double mah = 0.0;
    double fpr = 1.0;
    double fit = 0.0;
    uint64_t fingerprint = 0;
    bool stats_valid = false;
};

struct ModelWeights {
    std::vector<double> values;
    double minimum = 0.0;
    double range = 1.0;
};

ModelWeights model_weights_from_stats(
    const Candidate& candidate,
    const EncodedSequences& sequences,
    const BackgroundStats& background
);

int base_code(char nucleotide) {
    switch (nucleotide) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case 'N':
            return -1;
        default:
            return -2;
    }
}

char complement_base(char nucleotide) {
    switch (nucleotide) {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        default:
            return 'N';
    }
}

std::string reverse_complement(const std::string& sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        result.push_back(complement_base(*it));
    }
    return result;
}

std::string dinucleotide_label(int code) {
    static constexpr std::array<char, 4> letters = {'a', 'c', 'g', 't'};
    return {
        letters[static_cast<std::size_t>(code / 4)],
        letters[static_cast<std::size_t>(code % 4)],
    };
}

std::string make_path(
    const char* directory,
    const char* filename,
    const char* fallback_directory = nullptr
) {
    const auto* resolved_directory =
        directory != nullptr ? directory : fallback_directory;
    std::string path = resolved_directory != nullptr ? resolved_directory : "";
    if (filename != nullptr) {
        path += filename;
    }
    return path;
}

std::string model_name_from_file(const std::string& filename) {
    const auto slash = filename.find_last_of("/\\");
    const auto start = slash == std::string::npos ? 0 : slash + 1;
    const auto dot = filename.find('.', start);
    const auto stop = dot == std::string::npos ? filename.size() : dot;
    auto name = filename.substr(start, stop - start);
    return name.empty() ? "sitega" : name;
}

void copy_result_path(char* destination, const std::string& value) {
    std::snprintf(destination, 512, "%s", value.c_str());
}

SearchConfig make_config(const TrainParams& params) {
    SearchConfig config;
    config.fg_file = params.fg_file != nullptr ? params.fg_file : "";
    config.fg_path = make_path(params.fg_path, params.fg_file);
    config.bg_path = make_path(params.bg_path, params.bg_file, params.fg_path);
    config.out_path = params.out_path != nullptr ? params.out_path : "";
    config.log_path = config.out_path;
    config.log_path += params.log_file != nullptr ? params.log_file : "sitega.log";
    config.model_name = model_name_from_file(config.fg_file);
    config.max_lpd = params.max_lpd;
    config.motif_len = params.motif_len;
    config.feature_count = params.size;
    config.olig_bg = params.olig_bg;
    config.max_peak_len = params.max_peak_len;
    config.pop_size = params.pop_size > 0 ? params.pop_size : kDefaultPopulation;
    config.pop_size = std::min(config.pop_size, kMaxPopulation);
    config.num_motifs = params.num_motifs > 0 ? params.num_motifs : kDefaultNumMotifs;
    config.num_motifs = std::max(1, std::min(config.num_motifs, config.pop_size));
    config.generations = params.generations;
    config.mutation_attempts = params.mutation_attempts;
    config.stale_generations = params.stale_generations;
    config.verbose = params.verbose != 0;
    config.seed = params.seed != 0
        ? static_cast<uint64_t>(params.seed)
        : static_cast<uint64_t>(std::time(nullptr));
    return config;
}

void validate_config(const SearchConfig& config) {
    if (config.fg_path.empty()) {
        throw std::runtime_error("foreground FASTA path is empty");
    }
    if (config.bg_path.empty()) {
        throw std::runtime_error("background FASTA path is empty");
    }
    if (config.out_path.empty()) {
        throw std::runtime_error("output path is empty");
    }
    if (config.max_lpd <= 0 || config.max_lpd > kMaxLpdLen) {
        throw std::runtime_error("max_lpd must be in 1..6");
    }
    if (config.motif_len <= 1 || config.motif_len > kMaxMotifLen) {
        throw std::runtime_error("motif_len must be in 2..20");
    }
    if (config.feature_count <= 0 || config.feature_count > kMaxFeatureCount) {
        throw std::runtime_error("size must be in 1..120");
    }
    if (config.feature_count > kDinucleotideCount * (config.motif_len - 1)) {
        throw std::runtime_error("size exceeds available dinucleotide positions");
    }
    if (config.olig_bg <= 0 || config.olig_bg > config.motif_len) {
        throw std::runtime_error("olig_bg must be in 1..motif_len");
    }
    if (config.max_peak_len < config.motif_len) {
        throw std::runtime_error("max_peak_len must be at least motif_len");
    }
    if (config.pop_size <= 0 || config.pop_size > kMaxPopulation) {
        throw std::runtime_error("pop_size must be in 1..500");
    }
    if (config.generations < 0) {
        throw std::runtime_error("generations must be non-negative");
    }
    if (config.mutation_attempts < 0) {
        throw std::runtime_error("mutation_attempts must be non-negative");
    }
    if (config.stale_generations < 0) {
        throw std::runtime_error("stale_generations must be non-negative");
    }
}

std::vector<std::string> read_fasta(
    const std::string& path,
    int motif_len,
    int max_peak_len,
    std::ostream& log
) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open FASTA file: " + path);
    }

    std::vector<std::string> sequences;
    std::string current;
    std::string line;
    int record_index = 0;

    const auto flush = [&]() {
        if (current.empty()) {
            return;
        }

        ++record_index;
        bool valid = true;
        for (char& nucleotide : current) {
            nucleotide = static_cast<char>(
                std::toupper(static_cast<unsigned char>(nucleotide))
            );
            if (base_code(nucleotide) == -2) {
                valid = false;
                break;
            }
        }

        if (!valid) {
            log << "Skipping " << path << " record " << record_index
                << ": unsupported character\n";
        } else if (static_cast<int>(current.size()) < motif_len) {
            log << "Skipping " << path << " record " << record_index
                << ": shorter than motif_len\n";
        } else if (static_cast<int>(current.size()) > max_peak_len) {
            log << "Skipping " << path << " record " << record_index
                << ": longer than max_peak_len\n";
        } else {
            sequences.push_back(current);
        }
        current.clear();
    };

    while (std::getline(input, line)) {
        if (!line.empty() && line.front() == '>') {
            flush();
            continue;
        }
        for (char symbol : line) {
            if (!std::isspace(static_cast<unsigned char>(symbol))) {
                current.push_back(symbol);
            }
        }
    }
    flush();

    if (sequences.empty()) {
        throw std::runtime_error("FASTA file has no usable records: " + path);
    }
    return sequences;
}

std::vector<int> encode_dinucleotides(const std::string& sequence) {
    std::vector<int> encoded;
    if (sequence.size() < 2) {
        return encoded;
    }

    encoded.reserve(sequence.size() - 1);
    int previous = base_code(sequence.front());
    for (std::size_t index = 1; index < sequence.size(); ++index) {
        const int current = base_code(sequence[index]);
        if (previous >= 0 && current >= 0) {
            encoded.push_back(4 * previous + current);
        } else {
            encoded.push_back(-1);
        }
        previous = current;
    }
    return encoded;
}

PrefixTable build_prefix_table(const std::vector<int>& dinucs) {
    PrefixTable prefix(dinucs.size() + 1);
    for (std::size_t index = 0; index < dinucs.size(); ++index) {
        prefix[index + 1] = prefix[index];
        const int code = dinucs[index];
        if (code >= 0) {
            ++prefix[index + 1][static_cast<std::size_t>(code)];
        } else {
            ++prefix[index + 1][kInvalidSlot];
        }
    }
    return prefix;
}

int interval_count(
    const PrefixTable& prefix,
    int code,
    int start,
    int end
) {
    return prefix[static_cast<std::size_t>(end + 1)][static_cast<std::size_t>(code)] -
        prefix[static_cast<std::size_t>(start)][static_cast<std::size_t>(code)];
}

bool interval_has_invalid(const PrefixTable& prefix, int start, int end) {
    return (
        prefix[static_cast<std::size_t>(end + 1)][kInvalidSlot] -
        prefix[static_cast<std::size_t>(start)][kInvalidSlot]
    ) > 0;
}

bool feature_value_for_placement(
    const Feature& feature,
    const EncodedSequence& sequence,
    int position,
    unsigned char orientation,
    double& value
) {
    const auto strand = static_cast<std::size_t>(orientation);
    if (strand >= sequence.prefixes.size()) {
        return false;
    }

    const int start = position + feature.start;
    const int end = position + feature.end;
    const auto& prefix = sequence.prefixes[strand];
    if (start < 0 ||
        end < start ||
        static_cast<std::size_t>(end + 1) >= prefix.size()) {
        return false;
    }
    if (interval_has_invalid(prefix, start, end)) {
        return false;
    }

    value = static_cast<double>(interval_count(prefix, feature.code, start, end)) /
        static_cast<double>(feature.end - feature.start + 1);
    return true;
}

int kmer_code_at(const std::string& sequence, int start, int kmer_len) {
    int code = 0;
    for (int offset = 0; offset < kmer_len; ++offset) {
        const int nucleotide = base_code(sequence[static_cast<std::size_t>(start + offset)]);
        if (nucleotide < 0) {
            return -1;
        }
        code = 4 * code + nucleotide;
    }
    return code;
}

int power4(int exponent) {
    int value = 1;
    for (int index = 0; index < exponent; ++index) {
        value *= 4;
    }
    return value;
}

std::vector<double> kmer_frequencies(
    const std::vector<std::string>& sequences,
    int kmer_len
) {
    const int kmer_count = power4(kmer_len);
    std::vector<double> counts(static_cast<std::size_t>(kmer_count), 1.0);

    for (const auto& sequence : sequences) {
        const std::array<std::string, 2> strands = {
            sequence,
            reverse_complement(sequence),
        };
        for (const auto& strand : strands) {
            const int limit = static_cast<int>(strand.size()) - kmer_len + 1;
            for (int start = 0; start < limit; ++start) {
                const int code = kmer_code_at(strand, start, kmer_len);
                if (code >= 0) {
                    counts[static_cast<std::size_t>(code)] += 1.0;
                }
            }
        }
    }

    const double total = std::accumulate(counts.begin(), counts.end(), 0.0);
    for (double& count : counts) {
        count /= total;
    }
    return counts;
}

std::vector<double> kmer_log_ratios(
    const std::vector<std::string>& foreground,
    const std::vector<std::string>& background,
    int kmer_len
) {
    auto fg_freq = kmer_frequencies(foreground, kmer_len);
    auto bg_freq = kmer_frequencies(background, kmer_len);
    std::vector<double> ratios(fg_freq.size(), 0.0);
    for (std::size_t index = 0; index < ratios.size(); ++index) {
        ratios[index] = std::log10(fg_freq[index]) - std::log10(bg_freq[index]);
    }
    return ratios;
}

std::vector<double> sequence_window_weights(
    const std::string& sequence,
    const std::vector<double>& kmer_ratios,
    int motif_len,
    int kmer_len
) {
    const int window_count = static_cast<int>(sequence.size()) - motif_len + 1;
    const int kmer_windows_per_motif = motif_len - kmer_len + 1;
    std::vector<double> weights(static_cast<std::size_t>(window_count), 0.0);

    for (int window = 0; window < window_count; ++window) {
        double total = 0.0;
        for (int offset = 0; offset < kmer_windows_per_motif; ++offset) {
            const int code = kmer_code_at(sequence, window + offset, kmer_len);
            if (code >= 0) {
                total += kmer_ratios[static_cast<std::size_t>(code)];
            }
        }
        weights[static_cast<std::size_t>(window)] =
            total / static_cast<double>(kmer_windows_per_motif);
    }
    return weights;
}

void select_candidate_positions(EncodedSequence& sequence) {
    const int window_count = static_cast<int>(sequence.window_weights.size());
    if (window_count <= 0) {
        throw std::runtime_error("sequence has no motif windows");
    }

    auto thresholds = sequence.window_weights;
    const int rank = std::max(0, window_count / 3 - 1);
    std::nth_element(
        thresholds.begin(),
        thresholds.begin() + rank,
        thresholds.end(),
        std::greater<double>()
    );
    const double threshold = thresholds[static_cast<std::size_t>(rank)];

    double max_delta = 0.0;
    sequence.candidate_positions.clear();
    for (int position = 0; position < window_count; ++position) {
        const double delta =
            sequence.window_weights[static_cast<std::size_t>(position)] - threshold;
        if (delta >= 0.0) {
            sequence.candidate_positions.push_back(position);
            max_delta = std::max(max_delta, delta);
        }
    }
    if (sequence.candidate_positions.empty()) {
        sequence.candidate_positions.push_back(0);
    }

    const double scale = max_delta > 0.0 ? std::max(1.0, 5.0 / max_delta) : 1.0;
    sequence.cumulative_weights.resize(sequence.candidate_positions.size());
    int cumulative = 0;
    for (std::size_t index = 0; index < sequence.candidate_positions.size(); ++index) {
        const int position = sequence.candidate_positions[index];
        const double delta =
            sequence.window_weights[static_cast<std::size_t>(position)] - threshold;
        cumulative += 1 + static_cast<int>(scale * std::max(0.0, delta));
        sequence.cumulative_weights[index] = cumulative;
    }
}

EncodedSequences encode_sequences(
    const std::vector<std::string>& foreground,
    const std::vector<double>& kmer_ratios,
    int motif_len,
    int kmer_len
) {
    EncodedSequences encoded;
    encoded.records.reserve(foreground.size());
    for (const auto& sequence : foreground) {
        EncodedSequence record;
        record.forward = sequence;
        record.reverse = reverse_complement(sequence);
        record.dinucs[0] = encode_dinucleotides(record.forward);
        record.dinucs[1] = encode_dinucleotides(record.reverse);
        record.prefixes[0] = build_prefix_table(record.dinucs[0]);
        record.prefixes[1] = build_prefix_table(record.dinucs[1]);
        record.window_weights = sequence_window_weights(
            record.forward,
            kmer_ratios,
            motif_len,
            kmer_len
        );
        select_candidate_positions(record);
        encoded.records.push_back(std::move(record));
    }
    return encoded;
}

BackgroundStats build_background_stats(
    const std::vector<std::string>& background,
    int max_lpd
) {
    BackgroundStats stats;
    stats.mean.resize(static_cast<std::size_t>(max_lpd));
    stats.covariance.resize(static_cast<std::size_t>(max_lpd));
    std::vector<std::array<double, kDinucleotideCount>> second(
        static_cast<std::size_t>(max_lpd)
    );
    std::vector<int> interval_count_by_len(static_cast<std::size_t>(max_lpd), 0);
    std::array<double, kDinucleotideCount> expected_counts{};

    for (const auto& sequence : background) {
        const std::array<std::string, 2> strands = {
            sequence,
            reverse_complement(sequence),
        };
        for (const auto& strand : strands) {
            const auto dinucs = encode_dinucleotides(strand);
            const auto prefix = build_prefix_table(dinucs);

            for (int code : dinucs) {
                if (code >= 0) {
                    expected_counts[static_cast<std::size_t>(code)] += 1.0;
                }
            }

            for (int span = 0; span < max_lpd; ++span) {
                const int lpd_len = span + 1;
                const int limit = static_cast<int>(dinucs.size()) - lpd_len + 1;
                for (int start = 0; start < limit; ++start) {
                    const int end = start + lpd_len - 1;
                    if (interval_has_invalid(prefix, start, end)) {
                        continue;
                    }

                    ++interval_count_by_len[static_cast<std::size_t>(span)];
                    for (int code = 0; code < kDinucleotideCount; ++code) {
                        const double frequency = static_cast<double>(
                            interval_count(prefix, code, start, end)
                        ) / static_cast<double>(lpd_len);
                        stats.mean[static_cast<std::size_t>(span)]
                            [static_cast<std::size_t>(code)] += frequency;
                        second[static_cast<std::size_t>(span)]
                            [static_cast<std::size_t>(code)] += frequency * frequency;
                    }
                }
            }
        }
    }

    for (int span = 0; span < max_lpd; ++span) {
        const int count = interval_count_by_len[static_cast<std::size_t>(span)];
        if (count == 0) {
            continue;
        }
        for (int code = 0; code < kDinucleotideCount; ++code) {
            auto& mean = stats.mean[static_cast<std::size_t>(span)]
                [static_cast<std::size_t>(code)];
            auto& covariance = stats.covariance[static_cast<std::size_t>(span)]
                [static_cast<std::size_t>(code)];
            mean /= static_cast<double>(count);
            covariance = second[static_cast<std::size_t>(span)]
                [static_cast<std::size_t>(code)] / static_cast<double>(count);
            covariance -= mean * mean;
            if (covariance < 0.0 && covariance > -1e-12) {
                covariance = 0.0;
            }
        }
    }

    double expected_total =
        std::accumulate(expected_counts.begin(), expected_counts.end(), 0.0);
    if (expected_total <= 0.0) {
        expected_total = 1.0;
    }
    for (int code = 0; code < kDinucleotideCount; ++code) {
        stats.expected[static_cast<std::size_t>(code)] =
            expected_counts[static_cast<std::size_t>(code)] / expected_total;
    }
    return stats;
}

std::vector<FeaturePoolEntry> build_feature_pool(
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config
) {
    std::vector<FeaturePoolEntry> pool;
    const int motif_dinuc_count = config.motif_len - 1;
    pool.reserve(
        static_cast<std::size_t>(
            kDinucleotideCount * motif_dinuc_count * config.max_lpd
        )
    );

    for (int code = 0; code < kDinucleotideCount; ++code) {
        for (int start = 0; start < motif_dinuc_count; ++start) {
            const int max_end = std::min(
                start + config.max_lpd - 1,
                motif_dinuc_count - 1
            );
            for (int end = start; end <= max_end; ++end) {
                const Feature feature{start, end, code};
                double fg_sum = 0.0;
                double fg_second = 0.0;
                int sequence_count = 0;

                for (const auto& sequence : sequences.records) {
                    double sequence_sum = 0.0;
                    int placement_count = 0;
                    for (int position : sequence.candidate_positions) {
                        for (unsigned char orientation = 0;
                             orientation < 2;
                             ++orientation) {
                            double value = 0.0;
                            if (feature_value_for_placement(
                                    feature,
                                    sequence,
                                    position,
                                    orientation,
                                    value
                                )) {
                                sequence_sum += value;
                                ++placement_count;
                            }
                        }
                    }
                    if (placement_count == 0) {
                        continue;
                    }

                    const double sequence_mean =
                        sequence_sum / static_cast<double>(placement_count);
                    fg_sum += sequence_mean;
                    fg_second += sequence_mean * sequence_mean;
                    ++sequence_count;
                }

                if (sequence_count == 0) {
                    continue;
                }

                const double fg_mean = fg_sum / static_cast<double>(sequence_count);
                double fg_variance =
                    fg_second / static_cast<double>(sequence_count) -
                    fg_mean * fg_mean;
                fg_variance = std::max(0.0, fg_variance);

                const int span = end - start;
                const double bg_mean = background.mean[static_cast<std::size_t>(span)]
                    [static_cast<std::size_t>(code)];
                const double bg_variance = std::max(
                    0.0,
                    background.covariance[static_cast<std::size_t>(span)]
                        [static_cast<std::size_t>(code)]
                );
                const double diff = fg_mean - bg_mean;
                const double variance = fg_variance + bg_variance + kFeatureScoreRidge;
                const double score = diff * diff / variance;
                const double signed_effect = diff / variance;
                if (!std::isfinite(score) || !std::isfinite(signed_effect)) {
                    continue;
                }

                pool.push_back({
                    feature,
                    score,
                    signed_effect,
                    fg_mean,
                    bg_mean,
                    variance,
                });
            }
        }
    }

    std::sort(
        pool.begin(),
        pool.end(),
        [](const FeaturePoolEntry& lhs, const FeaturePoolEntry& rhs) {
            if (lhs.score != rhs.score) {
                return lhs.score > rhs.score;
            }
            if (lhs.feature.code != rhs.feature.code) {
                return lhs.feature.code < rhs.feature.code;
            }
            if (lhs.feature.start != rhs.feature.start) {
                return lhs.feature.start < rhs.feature.start;
            }
            return lhs.feature.end < rhs.feature.end;
        }
    );
    return pool;
}

bool feature_less(const Feature& lhs, const Feature& rhs) {
    if (lhs.code != rhs.code) {
        return lhs.code < rhs.code;
    }
    if (lhs.start != rhs.start) {
        return lhs.start < rhs.start;
    }
    return lhs.end < rhs.end;
}

bool same_feature(const Feature& lhs, const Feature& rhs) {
    return lhs.start == rhs.start && lhs.end == rhs.end && lhs.code == rhs.code;
}

void sort_features(Candidate& candidate) {
    std::sort(candidate.features.begin(), candidate.features.end(), feature_less);
}

bool overlaps(const Feature& lhs, const Feature& rhs) {
    return lhs.start <= rhs.end && rhs.start <= lhs.end;
}

bool can_use_feature(
    const std::vector<Feature>& features,
    const Feature& proposed,
    int ignored_index = -1
) {
    if (proposed.start < 0 || proposed.end < proposed.start) {
        return false;
    }
    for (std::size_t index = 0; index < features.size(); ++index) {
        if (static_cast<int>(index) == ignored_index) {
            continue;
        }
        const auto& existing = features[index];
        if (existing.code == proposed.code && overlaps(existing, proposed)) {
            return false;
        }
    }
    return true;
}

uint64_t mix_hash(uint64_t seed, uint64_t value) {
    value += 0x9E3779B97F4A7C15ULL;
    value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9ULL;
    value = (value ^ (value >> 27)) * 0x94D049BB133111EBULL;
    value ^= value >> 31;
    return seed ^ (value + 0x9E3779B97F4A7C15ULL + (seed << 6) + (seed >> 2));
}

uint64_t fingerprint_candidate(const Candidate& candidate) {
    uint64_t hash = 0xCBF29CE484222325ULL;
    for (const auto& feature : candidate.features) {
        hash = mix_hash(hash, static_cast<uint64_t>(feature.start));
        hash = mix_hash(hash, static_cast<uint64_t>(feature.end));
        hash = mix_hash(hash, static_cast<uint64_t>(feature.code));
    }
    hash = mix_hash(hash, 0xF00DCAFEULL);
    for (std::size_t index = 0; index < candidate.positions.size(); ++index) {
        hash = mix_hash(hash, static_cast<uint64_t>(candidate.positions[index]));
        hash = mix_hash(hash, static_cast<uint64_t>(candidate.orientations[index]));
    }
    return hash;
}

bool same_candidate(const Candidate& lhs, const Candidate& rhs) {
    if (lhs.features.size() != rhs.features.size() ||
        lhs.positions.size() != rhs.positions.size()) {
        return false;
    }
    for (std::size_t index = 0; index < lhs.features.size(); ++index) {
        const auto& left = lhs.features[index];
        const auto& right = rhs.features[index];
        if (!same_feature(left, right)) {
            return false;
        }
    }
    return lhs.positions == rhs.positions && lhs.orientations == rhs.orientations;
}

bool duplicate_candidate(
    const Candidate& candidate,
    const std::vector<Candidate>& population,
    int ignored_index
) {
    for (std::size_t index = 0; index < population.size(); ++index) {
        if (static_cast<int>(index) == ignored_index) {
            continue;
        }
        const auto& existing = population[index];
        if (candidate.fingerprint == existing.fingerprint &&
            same_candidate(candidate, existing)) {
            return true;
        }
    }
    return false;
}

int sample_weighted_position(
    const EncodedSequence& sequence,
    Rng& rng
) {
    const int total = sequence.cumulative_weights.back();
    const int value = rng.range(total);
    const auto it = std::upper_bound(
        sequence.cumulative_weights.begin(),
        sequence.cumulative_weights.end(),
        value
    );
    const std::size_t index = static_cast<std::size_t>(
        std::distance(sequence.cumulative_weights.begin(), it)
    );
    return sequence.candidate_positions[index];
}

Candidate make_random_candidate(
    const EncodedSequences& sequences,
    const SearchConfig& config,
    Rng& rng
) {
    Candidate candidate;
    const auto sequence_count = sequences.records.size();
    candidate.positions.resize(sequence_count);
    candidate.orientations.resize(sequence_count);

    for (std::size_t index = 0; index < sequence_count; ++index) {
        candidate.positions[index] = sample_weighted_position(
            sequences.records[index],
            rng
        );
        candidate.orientations[index] = static_cast<unsigned char>(rng.range(2));
    }

    const int motif_dinuc_count = config.motif_len - 1;
    int attempts = 0;
    while (static_cast<int>(candidate.features.size()) < config.feature_count &&
           attempts < 20000) {
        ++attempts;
        Feature feature;
        feature.code = rng.range(kDinucleotideCount);
        const int max_span = std::min(config.max_lpd - 1, motif_dinuc_count - 1);
        const int span = rng.range(max_span + 1);
        feature.start = rng.range(motif_dinuc_count - span);
        feature.end = feature.start + span;
        if (can_use_feature(candidate.features, feature)) {
            candidate.features.push_back(feature);
        }
    }

    if (static_cast<int>(candidate.features.size()) < config.feature_count) {
        for (int code = 0; code < kDinucleotideCount; ++code) {
            for (int start = 0; start < motif_dinuc_count; ++start) {
                Feature feature{start, start, code};
                if (can_use_feature(candidate.features, feature)) {
                    candidate.features.push_back(feature);
                    if (static_cast<int>(candidate.features.size()) ==
                        config.feature_count) {
                        break;
                    }
                }
            }
            if (static_cast<int>(candidate.features.size()) ==
                config.feature_count) {
                break;
            }
        }
    }

    if (static_cast<int>(candidate.features.size()) != config.feature_count) {
        throw std::runtime_error("failed to initialize candidate features");
    }

    sort_features(candidate);
    candidate.fingerprint = fingerprint_candidate(candidate);
    return candidate;
}

double placement_kmer_weight(
    const EncodedSequence& sequence,
    int motif_len,
    int position,
    unsigned char orientation
) {
    int forward_position = position;
    if (orientation != 0) {
        forward_position =
            static_cast<int>(sequence.forward.size()) - motif_len - position;
    }
    if (forward_position < 0 ||
        forward_position >= static_cast<int>(sequence.window_weights.size())) {
        return 0.0;
    }
    return sequence.window_weights[static_cast<std::size_t>(forward_position)];
}

int feature_values_for_placement(
    const Candidate& candidate,
    const EncodedSequence& sequence,
    int position,
    unsigned char orientation,
    std::vector<double>& values
) {
    values.assign(candidate.features.size(), 0.0);
    int invalid_count = 0;

    for (std::size_t index = 0; index < candidate.features.size(); ++index) {
        double value = 0.0;
        if (!feature_value_for_placement(
                candidate.features[index],
                sequence,
                position,
                orientation,
                value
            )) {
            ++invalid_count;
            continue;
        }
        values[index] = value;
    }
    return invalid_count;
}

void reset_stats(Candidate& candidate) {
    const auto feature_count = candidate.features.size();
    candidate.feature_sum.assign(feature_count, 0.0);
    candidate.second_moment.assign(feature_count * feature_count, 0.0);
    candidate.invalid_intervals = 0;
    candidate.kmer_sum = 0.0;
    candidate.stats_valid = true;
}

void add_outer_product(
    std::vector<double>& matrix,
    const std::vector<double>& values,
    double sign
) {
    const std::size_t size = values.size();
    for (std::size_t row = 0; row < size; ++row) {
        const double row_value = values[row];
        if (row_value == 0.0) {
            continue;
        }
        for (std::size_t col = 0; col < size; ++col) {
            const double col_value = values[col];
            if (col_value != 0.0) {
                matrix[row * size + col] += sign * row_value * col_value;
            }
        }
    }
}

bool solve_linear_system(
    std::vector<double> matrix,
    std::vector<double> rhs,
    std::vector<double>& solution,
    int size
) {
    solution.assign(static_cast<std::size_t>(size), 0.0);

    for (int col = 0; col < size; ++col) {
        int pivot = col;
        double best = std::fabs(matrix[static_cast<std::size_t>(col * size + col)]);
        for (int row = col + 1; row < size; ++row) {
            const double value = std::fabs(
                matrix[static_cast<std::size_t>(row * size + col)]
            );
            if (value > best) {
                best = value;
                pivot = row;
            }
        }
        if (best < kMinPivot) {
            return false;
        }
        if (pivot != col) {
            for (int k = col; k < size; ++k) {
                std::swap(
                    matrix[static_cast<std::size_t>(col * size + k)],
                    matrix[static_cast<std::size_t>(pivot * size + k)]
                );
            }
            std::swap(rhs[static_cast<std::size_t>(col)],
                      rhs[static_cast<std::size_t>(pivot)]);
        }

        for (int row = col + 1; row < size; ++row) {
            const double factor =
                matrix[static_cast<std::size_t>(row * size + col)] /
                matrix[static_cast<std::size_t>(col * size + col)];
            if (factor == 0.0) {
                continue;
            }
            matrix[static_cast<std::size_t>(row * size + col)] = 0.0;
            for (int k = col + 1; k < size; ++k) {
                matrix[static_cast<std::size_t>(row * size + k)] -=
                    factor * matrix[static_cast<std::size_t>(col * size + k)];
            }
            rhs[static_cast<std::size_t>(row)] -=
                factor * rhs[static_cast<std::size_t>(col)];
        }
    }

    for (int row = size - 1; row >= 0; --row) {
        double value = rhs[static_cast<std::size_t>(row)];
        for (int col = row + 1; col < size; ++col) {
            value -= matrix[static_cast<std::size_t>(row * size + col)] *
                solution[static_cast<std::size_t>(col)];
        }
        const double diagonal =
            matrix[static_cast<std::size_t>(row * size + row)];
        if (std::fabs(diagonal) < kMinPivot) {
            return false;
        }
        solution[static_cast<std::size_t>(row)] = value / diagonal;
    }
    return true;
}

bool solve_regularized(
    const std::vector<double>& covariance,
    const std::vector<double>& diff,
    std::vector<double>& solution,
    int size
) {
    constexpr std::array<double, 7> ridge_values = {
        0.0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5,
    };
    for (double ridge : ridge_values) {
        auto matrix = covariance;
        if (ridge > 0.0) {
            for (int index = 0; index < size; ++index) {
                matrix[static_cast<std::size_t>(index * size + index)] += ridge;
            }
        }
        if (solve_linear_system(matrix, diff, solution, size)) {
            return true;
        }
    }
    return false;
}

double evaluate_from_stats(
    Candidate& candidate,
    const EncodedSequences& sequences,
    const BackgroundStats& background
) {
    const int feature_count = static_cast<int>(candidate.features.size());
    const int sequence_count = static_cast<int>(sequences.records.size());

    candidate.mah = 0.0;
    candidate.fpr = 1.0;
    candidate.fit = 0.0;

    if (candidate.invalid_intervals > 0) {
        return 0.0;
    }

    std::vector<double> mean(static_cast<std::size_t>(feature_count), 0.0);
    std::vector<double> diff(static_cast<std::size_t>(feature_count), 0.0);
    std::vector<double> covariance(
        static_cast<std::size_t>(feature_count * feature_count),
        0.0
    );

    for (int row = 0; row < feature_count; ++row) {
        mean[static_cast<std::size_t>(row)] =
            candidate.feature_sum[static_cast<std::size_t>(row)] /
            static_cast<double>(sequence_count);
    }

    for (int row = 0; row < feature_count; ++row) {
        for (int col = 0; col < feature_count; ++col) {
            covariance[static_cast<std::size_t>(row * feature_count + col)] =
                candidate.second_moment[
                    static_cast<std::size_t>(row * feature_count + col)
                ] / static_cast<double>(sequence_count) -
                mean[static_cast<std::size_t>(row)] *
                    mean[static_cast<std::size_t>(col)];
        }
    }

    for (int index = 0; index < feature_count; ++index) {
        const auto& feature = candidate.features[static_cast<std::size_t>(index)];
        const int span = feature.end - feature.start;
        const double bg_mean = background.mean[static_cast<std::size_t>(span)]
            [static_cast<std::size_t>(feature.code)];
        diff[static_cast<std::size_t>(index)] =
            mean[static_cast<std::size_t>(index)] - bg_mean;
        covariance[static_cast<std::size_t>(index * feature_count + index)] +=
            background.covariance[static_cast<std::size_t>(span)]
                [static_cast<std::size_t>(feature.code)];
    }

    for (double& value : covariance) {
        value *= 0.5;
    }

    std::vector<double> solution;
    if (!solve_regularized(covariance, diff, solution, feature_count)) {
        return 0.0;
    }

    double mah = 0.0;
    for (int index = 0; index < feature_count; ++index) {
        mah += diff[static_cast<std::size_t>(index)] *
            solution[static_cast<std::size_t>(index)];
    }
    if (!std::isfinite(mah) || mah <= 0.0) {
        return 0.0;
    }

    candidate.mah = mah;
    candidate.fpr = std::pow(
        10.0,
        candidate.kmer_sum / static_cast<double>(sequence_count)
    );
    if (!std::isfinite(candidate.fpr) || candidate.fpr <= 0.0) {
        candidate.fpr = 1.0;
    }
    candidate.fit = candidate.mah * candidate.fpr;
    return candidate.fit;
}

double recompute_candidate_stats(
    Candidate& candidate,
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config
) {
    reset_stats(candidate);
    std::vector<double> values;
    for (std::size_t seq_index = 0; seq_index < sequences.records.size(); ++seq_index) {
        const auto& sequence = sequences.records[seq_index];
        candidate.invalid_intervals += feature_values_for_placement(
            candidate,
            sequence,
            candidate.positions[seq_index],
            candidate.orientations[seq_index],
            values
        );
        for (std::size_t index = 0; index < values.size(); ++index) {
            candidate.feature_sum[index] += values[index];
        }
        add_outer_product(candidate.second_moment, values, 1.0);
        candidate.kmer_sum += placement_kmer_weight(
            sequence,
            config.motif_len,
            candidate.positions[seq_index],
            candidate.orientations[seq_index]
        );
    }
    candidate.fingerprint = fingerprint_candidate(candidate);
    return evaluate_from_stats(candidate, sequences, background);
}

double update_placement_stats(
    Candidate& candidate,
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config,
    int sequence_index,
    int old_position,
    unsigned char old_orientation,
    int new_position,
    unsigned char new_orientation
) {
    if (!candidate.stats_valid) {
        return recompute_candidate_stats(candidate, sequences, background, config);
    }

    const auto& sequence = sequences.records[static_cast<std::size_t>(sequence_index)];
    std::vector<double> old_values;
    std::vector<double> new_values;
    const int old_invalid = feature_values_for_placement(
        candidate,
        sequence,
        old_position,
        old_orientation,
        old_values
    );
    const int new_invalid = feature_values_for_placement(
        candidate,
        sequence,
        new_position,
        new_orientation,
        new_values
    );

    for (std::size_t index = 0; index < old_values.size(); ++index) {
        candidate.feature_sum[index] += new_values[index] - old_values[index];
    }
    add_outer_product(candidate.second_moment, old_values, -1.0);
    add_outer_product(candidate.second_moment, new_values, 1.0);
    candidate.invalid_intervals += new_invalid - old_invalid;
    candidate.kmer_sum += placement_kmer_weight(
        sequence,
        config.motif_len,
        new_position,
        new_orientation
    ) - placement_kmer_weight(
        sequence,
        config.motif_len,
        old_position,
        old_orientation
    );
    candidate.fingerprint = fingerprint_candidate(candidate);
    return evaluate_from_stats(candidate, sequences, background);
}

bool mutate_feature_code(Candidate& candidate, const SearchConfig& config, Rng& rng) {
    (void)config;
    const int feature_index = rng.range(static_cast<int>(candidate.features.size()));
    auto feature = candidate.features[static_cast<std::size_t>(feature_index)];

    for (int attempt = 0; attempt < 32; ++attempt) {
        int new_code = rng.range(kDinucleotideCount - 1);
        if (new_code >= feature.code) {
            ++new_code;
        }
        Feature proposed = feature;
        proposed.code = new_code;
        if (can_use_feature(candidate.features, proposed, feature_index)) {
            candidate.features[static_cast<std::size_t>(feature_index)] = proposed;
            sort_features(candidate);
            candidate.stats_valid = false;
            return true;
        }
    }
    return false;
}

bool mutate_feature_interval(
    Candidate& candidate,
    const SearchConfig& config,
    Rng& rng
) {
    const int feature_index = rng.range(static_cast<int>(candidate.features.size()));
    const auto original = candidate.features[static_cast<std::size_t>(feature_index)];
    const int motif_dinuc_count = config.motif_len - 1;

    for (int attempt = 0; attempt < 64; ++attempt) {
        const int max_span = std::min(config.max_lpd - 1, motif_dinuc_count - 1);
        const int span = rng.range(max_span + 1);
        Feature proposed = original;
        proposed.start = rng.range(motif_dinuc_count - span);
        proposed.end = proposed.start + span;
        if (proposed.start == original.start && proposed.end == original.end) {
            continue;
        }
        if (can_use_feature(candidate.features, proposed, feature_index)) {
            candidate.features[static_cast<std::size_t>(feature_index)] = proposed;
            sort_features(candidate);
            candidate.stats_valid = false;
            return true;
        }
    }
    return false;
}

bool try_feature_mutation(
    Candidate& candidate,
    const std::vector<Candidate>& population,
    int population_index,
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config,
    Rng& rng,
    bool mutate_code
) {
    const double old_fit = candidate.fit;
    const auto old_features = candidate.features;
    const uint64_t old_fingerprint = candidate.fingerprint;

    const bool changed = mutate_code
        ? mutate_feature_code(candidate, config, rng)
        : mutate_feature_interval(candidate, config, rng);
    if (!changed) {
        return false;
    }

    recompute_candidate_stats(candidate, sequences, background, config);
    const bool accept = candidate.fit > old_fit + kScoreEpsilon &&
        !duplicate_candidate(candidate, population, population_index);
    if (accept) {
        return true;
    }

    candidate.features = old_features;
    candidate.fingerprint = old_fingerprint;
    candidate.stats_valid = false;
    recompute_candidate_stats(candidate, sequences, background, config);
    return false;
}

bool try_placement_mutation(
    Candidate& candidate,
    const std::vector<Candidate>& population,
    int population_index,
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config,
    Rng& rng
) {
    const int sequence_index = rng.range(static_cast<int>(sequences.records.size()));
    const auto& sequence = sequences.records[static_cast<std::size_t>(sequence_index)];
    const int old_position = candidate.positions[static_cast<std::size_t>(sequence_index)];
    const auto old_orientation =
        candidate.orientations[static_cast<std::size_t>(sequence_index)];
    int new_position = sample_weighted_position(sequence, rng);
    auto new_orientation = static_cast<unsigned char>(rng.range(2));

    if (new_position == old_position && new_orientation == old_orientation) {
        new_orientation = static_cast<unsigned char>(1 - new_orientation);
    }

    const double old_fit = candidate.fit;
    candidate.positions[static_cast<std::size_t>(sequence_index)] = new_position;
    candidate.orientations[static_cast<std::size_t>(sequence_index)] = new_orientation;
    update_placement_stats(
        candidate,
        sequences,
        background,
        config,
        sequence_index,
        old_position,
        old_orientation,
        new_position,
        new_orientation
    );

    const bool accept = candidate.fit > old_fit + kScoreEpsilon &&
        !duplicate_candidate(candidate, population, population_index);
    if (accept) {
        return true;
    }

    candidate.positions[static_cast<std::size_t>(sequence_index)] = old_position;
    candidate.orientations[static_cast<std::size_t>(sequence_index)] = old_orientation;
    update_placement_stats(
        candidate,
        sequences,
        background,
        config,
        sequence_index,
        new_position,
        new_orientation,
        old_position,
        old_orientation
    );
    return false;
}

double placement_model_score(
    const std::vector<double>& values,
    const ModelWeights& weights
) {
    double score = 0.0;
    for (std::size_t index = 0; index < values.size(); ++index) {
        score += weights.values[index] * values[index];
    }
    return (score - weights.minimum) / weights.range;
}

std::vector<int> refinement_sequence_indices(
    const Candidate& candidate,
    const EncodedSequences& sequences,
    Rng& rng
) {
    const int sequence_count = static_cast<int>(sequences.records.size());
    std::vector<int> indices;
    std::vector<int> optional_indices;
    indices.reserve(static_cast<std::size_t>(sequence_count));
    optional_indices.reserve(static_cast<std::size_t>(sequence_count));

    std::vector<double> values;
    for (int index = 0; index < sequence_count; ++index) {
        const auto& sequence = sequences.records[static_cast<std::size_t>(index)];
        const int invalid_count = feature_values_for_placement(
            candidate,
            sequence,
            candidate.positions[static_cast<std::size_t>(index)],
            candidate.orientations[static_cast<std::size_t>(index)],
            values
        );
        if (invalid_count > 0) {
            indices.push_back(index);
        } else {
            optional_indices.push_back(index);
        }
    }

    const int extra_count = std::min(
        kPlacementRefineSequenceLimit,
        static_cast<int>(optional_indices.size())
    );
    for (int index = 0; index < extra_count; ++index) {
        const int chosen = index + rng.range(
            static_cast<int>(optional_indices.size()) - index
        );
        std::swap(optional_indices[static_cast<std::size_t>(index)],
                  optional_indices[static_cast<std::size_t>(chosen)]);
        indices.push_back(optional_indices[static_cast<std::size_t>(index)]);
    }

    return indices;
}

bool choose_kmer_refined_placement(
    Candidate& candidate,
    const EncodedSequence& sequence,
    const SearchConfig& config,
    int sequence_index,
    Rng& rng
) {
    const auto slot = static_cast<std::size_t>(sequence_index);
    const int old_position = candidate.positions[slot];
    const auto old_orientation = candidate.orientations[slot];
    int best_position = old_position;
    unsigned char best_orientation = old_orientation;
    double best_score = 0.0;
    bool found = false;
    std::vector<double> values;

    const auto consider = [&](int position, unsigned char orientation) {
        if (feature_values_for_placement(
                candidate,
                sequence,
                position,
                orientation,
                values
            ) > 0) {
            return;
        }
        const double score = placement_kmer_weight(
            sequence,
            config.motif_len,
            position,
            orientation
        );
        if (!found || score > best_score) {
            found = true;
            best_score = score;
            best_position = position;
            best_orientation = orientation;
        }
    };

    consider(old_position, old_orientation);
    consider(old_position, static_cast<unsigned char>(1 - old_orientation));
    for (int sample = 0; sample < kPlacementRefineSamples; ++sample) {
        const int position = sample_weighted_position(sequence, rng);
        consider(position, 0);
        consider(position, 1);
    }

    if (!found ||
        (best_position == old_position && best_orientation == old_orientation)) {
        return false;
    }
    candidate.positions[slot] = best_position;
    candidate.orientations[slot] = best_orientation;
    return true;
}

bool choose_model_refined_placement(
    Candidate& candidate,
    const EncodedSequence& sequence,
    const ModelWeights& weights,
    int sequence_index,
    Rng& rng
) {
    const auto slot = static_cast<std::size_t>(sequence_index);
    const int old_position = candidate.positions[slot];
    const auto old_orientation = candidate.orientations[slot];
    int best_position = old_position;
    unsigned char best_orientation = old_orientation;
    double best_score = 0.0;
    bool found = false;
    std::vector<double> values;

    const auto consider = [&](int position, unsigned char orientation) {
        if (feature_values_for_placement(
                candidate,
                sequence,
                position,
                orientation,
                values
            ) > 0) {
            return;
        }
        const double score = placement_model_score(values, weights);
        if (!found || score > best_score) {
            found = true;
            best_score = score;
            best_position = position;
            best_orientation = orientation;
        }
    };

    consider(old_position, old_orientation);
    consider(old_position, static_cast<unsigned char>(1 - old_orientation));
    for (int sample = 0; sample < kPlacementRefineSamples; ++sample) {
        const int position = sample_weighted_position(sequence, rng);
        consider(position, 0);
        consider(position, 1);
    }

    if (!found ||
        (best_position == old_position && best_orientation == old_orientation)) {
        return false;
    }
    candidate.positions[slot] = best_position;
    candidate.orientations[slot] = best_orientation;
    return true;
}

void refine_placements_after_feature_change(
    Candidate& candidate,
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config,
    Rng& rng
) {
    recompute_candidate_stats(candidate, sequences, background, config);
    const auto sequence_indices = refinement_sequence_indices(
        candidate,
        sequences,
        rng
    );

    if (candidate.invalid_intervals > 0) {
        bool repaired = false;
        for (int sequence_index : sequence_indices) {
            repaired = choose_kmer_refined_placement(
                candidate,
                sequences.records[static_cast<std::size_t>(sequence_index)],
                config,
                sequence_index,
                rng
            ) || repaired;
        }
        if (repaired) {
            recompute_candidate_stats(candidate, sequences, background, config);
        }
        if (candidate.invalid_intervals > 0) {
            return;
        }
    }

    const auto weights = model_weights_from_stats(candidate, sequences, background);
    bool refined = false;
    for (int sequence_index : sequence_indices) {
        refined = choose_model_refined_placement(
            candidate,
            sequences.records[static_cast<std::size_t>(sequence_index)],
            weights,
            sequence_index,
            rng
        ) || refined;
    }
    if (refined) {
        recompute_candidate_stats(candidate, sequences, background, config);
    }
}

bool try_directed_feature_replacement(
    Candidate& candidate,
    const std::vector<Candidate>& population,
    int population_index,
    const std::vector<FeaturePoolEntry>& feature_pool,
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config,
    Rng& rng
) {
    if (feature_pool.empty() || candidate.features.empty()) {
        return false;
    }

    const auto pool_limit = std::min(
        feature_pool.size(),
        static_cast<std::size_t>(
            std::max(
                config.feature_count * kFeaturePoolTopMultiplier,
                config.feature_count
            )
        )
    );
    const int feature_index = rng.range(static_cast<int>(candidate.features.size()));
    const double old_fit = candidate.fit;
    const auto old_features = candidate.features;
    const auto old_positions = candidate.positions;
    const auto old_orientations = candidate.orientations;
    const uint64_t old_fingerprint = candidate.fingerprint;

    for (int attempt = 0; attempt < kDirectedFeatureAttempts; ++attempt) {
        const int pool_index = rng.range(static_cast<int>(pool_limit));
        const auto& proposed =
            feature_pool[static_cast<std::size_t>(pool_index)].feature;
        if (same_feature(
                candidate.features[static_cast<std::size_t>(feature_index)],
                proposed
            )) {
            continue;
        }
        if (!can_use_feature(candidate.features, proposed, feature_index)) {
            continue;
        }

        candidate.features[static_cast<std::size_t>(feature_index)] = proposed;
        sort_features(candidate);
        candidate.stats_valid = false;
        refine_placements_after_feature_change(
            candidate,
            sequences,
            background,
            config,
            rng
        );

        const bool accept = candidate.fit > old_fit + kScoreEpsilon &&
            !duplicate_candidate(candidate, population, population_index);
        if (accept) {
            return true;
        }

        candidate.features = old_features;
        candidate.positions = old_positions;
        candidate.orientations = old_orientations;
        candidate.fingerprint = old_fingerprint;
        candidate.stats_valid = false;
        recompute_candidate_stats(candidate, sequences, background, config);
        return false;
    }

    return false;
}

bool valid_feature_set(const std::vector<Feature>& features) {
    for (std::size_t index = 0; index < features.size(); ++index) {
        if (!can_use_feature(features, features[index], static_cast<int>(index))) {
            return false;
        }
    }
    return true;
}

Candidate recombine_candidates(
    const Candidate& first,
    const Candidate& second,
    const SearchConfig& config,
    const EncodedSequences& sequences,
    Rng& rng
) {
    Candidate child = first;
    const int mode = rng.range(2);

    if (mode == 0) {
        for (int attempt = 0; attempt < 16; ++attempt) {
            const int child_index = rng.range(static_cast<int>(child.features.size()));
            const int donor_index = rng.range(static_cast<int>(second.features.size()));
            const auto proposed = second.features[static_cast<std::size_t>(donor_index)];
            if (can_use_feature(child.features, proposed, child_index)) {
                child.features[static_cast<std::size_t>(child_index)] = proposed;
                sort_features(child);
                break;
            }
        }
    } else {
        const int sequence_count = static_cast<int>(sequences.records.size());
        const int swaps = std::max(1, sequence_count / 4);
        for (int swap = 0; swap < swaps; ++swap) {
            const int sequence_index = rng.range(sequence_count);
            child.positions[static_cast<std::size_t>(sequence_index)] =
                second.positions[static_cast<std::size_t>(sequence_index)];
            child.orientations[static_cast<std::size_t>(sequence_index)] =
                second.orientations[static_cast<std::size_t>(sequence_index)];
        }
    }

    if (!valid_feature_set(child.features)) {
        return first;
    }
    child.stats_valid = false;
    child.fingerprint = fingerprint_candidate(child);
    (void)config;
    return child;
}

std::vector<Candidate> initialize_population(
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config,
    Rng& rng,
    std::ostream& log
) {
    std::vector<Candidate> population;
    population.reserve(static_cast<std::size_t>(config.pop_size));
    int attempts = 0;
    while (static_cast<int>(population.size()) < config.pop_size &&
           attempts < config.pop_size * 200) {
        ++attempts;
        auto candidate = make_random_candidate(sequences, config, rng);
        recompute_candidate_stats(candidate, sequences, background, config);
        if (!duplicate_candidate(candidate, population, -1)) {
            population.push_back(std::move(candidate));
        }
    }

    if (population.empty()) {
        throw std::runtime_error("failed to initialize SiteGA population");
    }
    if (static_cast<int>(population.size()) < config.pop_size) {
        log << "Initialized " << population.size() << " unique candidates out of "
            << config.pop_size << " requested\n";
    }
    return population;
}

void sort_population(std::vector<Candidate>& population) {
    std::sort(
        population.begin(),
        population.end(),
        [](const Candidate& lhs, const Candidate& rhs) {
            return lhs.fit > rhs.fit;
        }
    );
}

void run_search(
    std::vector<Candidate>& population,
    const std::vector<FeaturePoolEntry>& feature_pool,
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config,
    Rng& rng,
    std::ostream& log
) {
    sort_population(population);
    const int generations = config.generations > 0
        ? config.generations
        : std::max(6, std::min(24, 8 + config.motif_len / 2));
    const int mutation_attempts = config.mutation_attempts > 0
        ? config.mutation_attempts
        : 3;
    const int stale_generation_limit = config.stale_generations > 0
        ? config.stale_generations
        : 3;
    const int mutation_type_count = feature_pool.empty() ? 3 : 4;
    int stale_generations = 0;
    log << "Search generations=" << generations
        << " mutation_attempts=" << mutation_attempts
        << " stale_generations=" << stale_generation_limit << '\n';

    for (int generation = 0; generation < generations; ++generation) {
        const double best_before = population.front().fit;
        int accepted_mutations = 0;
        int accepted_directed_mutations = 0;
        int accepted_recombinations = 0;

        for (std::size_t index = 0; index < population.size(); ++index) {
            for (int attempt = 0; attempt < mutation_attempts; ++attempt) {
                const int type = rng.range(mutation_type_count);
                bool accepted = false;
                if (type == 0) {
                    accepted = try_feature_mutation(
                        population[index],
                        population,
                        static_cast<int>(index),
                        sequences,
                        background,
                        config,
                        rng,
                        true
                    );
                } else if (type == 1) {
                    accepted = try_feature_mutation(
                        population[index],
                        population,
                        static_cast<int>(index),
                        sequences,
                        background,
                        config,
                        rng,
                        false
                    );
                } else if (type == 2) {
                    accepted = try_placement_mutation(
                        population[index],
                        population,
                        static_cast<int>(index),
                        sequences,
                        background,
                        config,
                        rng
                    );
                } else {
                    accepted = try_directed_feature_replacement(
                        population[index],
                        population,
                        static_cast<int>(index),
                        feature_pool,
                        sequences,
                        background,
                        config,
                        rng
                    );
                }
                if (accepted) {
                    if (type == 3) {
                        ++accepted_directed_mutations;
                    } else {
                        ++accepted_mutations;
                    }
                }
            }
        }

        sort_population(population);

        const int population_size = static_cast<int>(population.size());
        const int recombination_attempts = population_size > 1
            ? std::max(1, population_size / 2)
            : 0;
        const int parent_pool = std::max(1, population_size / 2);
        for (int attempt = 0; attempt < recombination_attempts; ++attempt) {
            const int first_index = rng.range(parent_pool);
            int second_index = rng.range(population_size - 1);
            if (second_index >= first_index) {
                ++second_index;
            }
            const int target_index =
                population_size - 1 - rng.range(parent_pool);
            auto child = recombine_candidates(
                population[static_cast<std::size_t>(first_index)],
                population[static_cast<std::size_t>(second_index)],
                config,
                sequences,
                rng
            );
            recompute_candidate_stats(child, sequences, background, config);
            if (child.fit > population[static_cast<std::size_t>(target_index)].fit +
                    kScoreEpsilon &&
                !duplicate_candidate(child, population, target_index)) {
                population[static_cast<std::size_t>(target_index)] = std::move(child);
                ++accepted_recombinations;
            }
        }

        sort_population(population);
        const double improvement = population.front().fit - best_before;
        log << "Gen " << (generation + 1)
            << " Fit " << std::setprecision(8) << population.front().fit
            << " Mut " << accepted_mutations
            << " Dir " << accepted_directed_mutations
            << " Rec " << accepted_recombinations
            << " Delta " << improvement << '\n';

        if (improvement <= std::max(1e-9, std::fabs(best_before) * 1e-5)) {
            ++stale_generations;
            if (stale_generations >= stale_generation_limit) {
                break;
            }
        } else {
            stale_generations = 0;
        }
    }
}

std::vector<Candidate> run_island_searches(
    const std::vector<FeaturePoolEntry>& feature_pool,
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config,
    Rng& rng,
    std::ostream& log
) {
    std::vector<Candidate> winners;
    winners.reserve(static_cast<std::size_t>(config.num_motifs));

    for (int island = 0; island < config.num_motifs; ++island) {
        log << "Island " << (island + 1) << "/" << config.num_motifs << '\n';

        auto population = initialize_population(
            sequences,
            background,
            config,
            rng,
            log
        );
        run_search(
            population,
            feature_pool,
            sequences,
            background,
            config,
            rng,
            log
        );
        sort_population(population);

        winners.push_back(std::move(population.front()));
        log << "Island " << (island + 1)
            << " winner_fit=" << std::setprecision(8)
            << winners.back().fit << '\n';
    }

    sort_population(winners);
    return winners;
}

ModelWeights model_weights_from_stats(
    const Candidate& candidate,
    const EncodedSequences& sequences,
    const BackgroundStats& background
) {
    ModelWeights weights;
    const int feature_count = static_cast<int>(candidate.features.size());
    weights.values.assign(static_cast<std::size_t>(feature_count), 0.0);
    if (candidate.invalid_intervals > 0 || candidate.mah <= 0.0) {
        return weights;
    }

    const int sequence_count = static_cast<int>(sequences.records.size());
    std::vector<double> mean(static_cast<std::size_t>(feature_count), 0.0);
    std::vector<double> diff(static_cast<std::size_t>(feature_count), 0.0);
    std::vector<double> covariance(
        static_cast<std::size_t>(feature_count * feature_count),
        0.0
    );

    for (int row = 0; row < feature_count; ++row) {
        mean[static_cast<std::size_t>(row)] =
            candidate.feature_sum[static_cast<std::size_t>(row)] /
            static_cast<double>(sequence_count);
    }
    for (int row = 0; row < feature_count; ++row) {
        for (int col = 0; col < feature_count; ++col) {
            covariance[static_cast<std::size_t>(row * feature_count + col)] =
                candidate.second_moment[
                    static_cast<std::size_t>(row * feature_count + col)
                ] / static_cast<double>(sequence_count) -
                mean[static_cast<std::size_t>(row)] *
                    mean[static_cast<std::size_t>(col)];
        }
    }
    for (int index = 0; index < feature_count; ++index) {
        const auto& feature = candidate.features[static_cast<std::size_t>(index)];
        const int span = feature.end - feature.start;
        diff[static_cast<std::size_t>(index)] =
            mean[static_cast<std::size_t>(index)] -
            background.mean[static_cast<std::size_t>(span)]
                [static_cast<std::size_t>(feature.code)];
        covariance[static_cast<std::size_t>(index * feature_count + index)] +=
            background.covariance[static_cast<std::size_t>(span)]
                [static_cast<std::size_t>(feature.code)];
    }
    for (double& value : covariance) {
        value *= 0.5;
    }

    std::vector<double> solution;
    if (!solve_regularized(covariance, diff, solution, feature_count)) {
        return weights;
    }

    double mah = 0.0;
    for (int index = 0; index < feature_count; ++index) {
        mah += diff[static_cast<std::size_t>(index)] *
            solution[static_cast<std::size_t>(index)];
    }
    if (!std::isfinite(mah) || mah <= 0.0) {
        return weights;
    }

    double minimum = 0.0;
    double maximum = 0.0;
    for (int index = 0; index < feature_count; ++index) {
        const double value = solution[static_cast<std::size_t>(index)] / mah;
        weights.values[static_cast<std::size_t>(index)] = value;
        if (value < 0.0) {
            minimum += value;
        } else {
            maximum += value;
        }
    }
    weights.minimum = minimum;
    weights.range = maximum - minimum;
    if (weights.range <= 0.0 || !std::isfinite(weights.range)) {
        weights.range = 1.0;
    }
    return weights;
}

ModelWeights compute_model_weights(
    Candidate& candidate,
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config
) {
    recompute_candidate_stats(candidate, sequences, background, config);
    return model_weights_from_stats(candidate, sequences, background);
}

void write_matrix(
    const std::string& path,
    const Candidate& candidate,
    const ModelWeights& weights,
    const SearchConfig& config
) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("cannot open matrix output: " + path);
    }

    output << config.model_name << '\n';
    output << candidate.features.size() << "\tLPD count\n";
    output << config.motif_len << "\tModel length\n";
    output << std::fixed << std::setprecision(12)
           << weights.minimum << "\tMinimum\n";
    output << std::fixed << std::setprecision(12)
           << weights.range << "\tRazmah\n";

    for (std::size_t index = 0; index < candidate.features.size(); ++index) {
        const auto& feature = candidate.features[index];
        output << feature.start << '\t'
               << feature.end << '\t'
               << std::fixed << std::setprecision(12) << weights.values[index] << '\t'
               << feature.code << '\t'
               << dinucleotide_label(feature.code) << '\n';
    }
}

double output_feature_value(
    const Feature& feature,
    const EncodedSequence& sequence,
    int position,
    unsigned char orientation
) {
    double value = 0.0;
    if (!feature_value_for_placement(
            feature,
            sequence,
            position,
            orientation,
            value
        )) {
        return 0.0;
    }
    return value;
}

double score_window(
    const Candidate& candidate,
    const EncodedSequence& sequence,
    const ModelWeights& weights,
    int position,
    unsigned char orientation
) {
    double score = 0.0;
    for (std::size_t index = 0; index < candidate.features.size(); ++index) {
        score += weights.values[index] * output_feature_value(
            candidate.features[index],
            sequence,
            position,
            orientation
        );
    }
    return (score - weights.minimum) / weights.range;
}

std::string motif_string(
    const EncodedSequence& sequence,
    int position,
    unsigned char orientation,
    int motif_len
) {
    const auto& source = orientation == 0 ? sequence.forward : sequence.reverse;
    std::string result =
        source.substr(static_cast<std::size_t>(position), static_cast<std::size_t>(motif_len));
    std::transform(result.begin(), result.end(), result.begin(), [](char value) {
        return static_cast<char>(std::tolower(static_cast<unsigned char>(value)));
    });
    return result;
}

void write_locations(
    const std::string& path,
    const Candidate& candidate,
    const ModelWeights& weights,
    const EncodedSequences& sequences,
    const SearchConfig& config
) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("cannot open location output: " + path);
    }

    for (std::size_t seq_index = 0; seq_index < sequences.records.size(); ++seq_index) {
        const auto& sequence = sequences.records[seq_index];
        const int window_count =
            static_cast<int>(sequence.forward.size()) - config.motif_len + 1;
        double best_score = -std::numeric_limits<double>::infinity();
        int best_position = 0;
        unsigned char best_orientation = 0;

        for (int position = 0; position < window_count; ++position) {
            for (unsigned char orientation = 0; orientation < 2; ++orientation) {
                const double score = score_window(
                    candidate,
                    sequence,
                    weights,
                    position,
                    orientation
                );
                if (score > best_score) {
                    best_score = score;
                    best_position = position;
                    best_orientation = orientation;
                }
            }
        }

        int start = 0;
        int end = 0;
        if (best_orientation == 0) {
            start = best_position + 1;
            end = best_position + config.motif_len;
        } else {
            start = static_cast<int>(sequence.forward.size()) -
                best_position - config.motif_len + 1;
            end = static_cast<int>(sequence.forward.size()) - best_position;
        }

        output << (seq_index + 1) << '\t'
               << start << '\t'
               << end << '\t'
               << (best_orientation == 0 ? '+' : '-') << '\t'
               << std::fixed << std::setprecision(12) << best_score << '\t'
               << motif_string(sequence, best_position, best_orientation, config.motif_len)
               << '\n';
    }
}

std::string output_base(const SearchConfig& config) {
    std::string base = config.out_path;
    if (!base.empty() && base.back() != '/' && base.back() != '\\') {
        base.push_back('/');
    }
    base += config.model_name;
    return base;
}

void write_outputs(
    std::vector<Candidate>& candidates,
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config,
    TrainResult* result
) {
    const std::filesystem::path output_directory(config.out_path);
    if (!config.out_path.empty()) {
        std::filesystem::create_directories(output_directory);
    }

    const std::string base = output_base(config);
    const int num_motifs = std::min(
        config.num_motifs,
        static_cast<int>(candidates.size())
    );
    for (int index = 0; index < num_motifs; ++index) {
        auto& candidate = candidates[static_cast<std::size_t>(index)];
        auto weights = compute_model_weights(candidate, sequences, background, config);
        const std::string suffix = std::to_string(index + 1);
        write_matrix(base + "_mat" + suffix, candidate, weights, config);
        write_locations(base + "_loc" + suffix, candidate, weights, sequences, config);
    }

    if (result != nullptr) {
        copy_result_path(result->mat_path, base + "_mat1");
        copy_result_path(result->loc_path, base + "_loc1");
        result->best_fit = candidates.empty() ? 0.0 : candidates.front().fit;
        result->status = 0;
    }
}

int train_impl(const TrainParams& params, TrainResult* result) {
    auto config = make_config(params);
    validate_config(config);

    const std::filesystem::path output_directory(config.out_path);
    if (!config.out_path.empty()) {
        std::filesystem::create_directories(output_directory);
    }

    std::ofstream log(config.log_path);
    if (!log) {
        throw std::runtime_error("cannot open log file: " + config.log_path);
    }

    log << "SiteGA update backend\n";
    log << "RNG seed: " << config.seed << '\n';
    log << "Foreground: " << config.fg_path << '\n';
    log << "Background: " << config.bg_path << '\n';
    log << "motif_len=" << config.motif_len
        << " size=" << config.feature_count
        << " max_lpd=" << config.max_lpd
        << " olig_bg=" << config.olig_bg
        << " pop_size_per_island=" << config.pop_size
        << " islands=" << config.num_motifs << '\n';

    const auto foreground = read_fasta(
        config.fg_path,
        config.motif_len,
        config.max_peak_len,
        log
    );
    const auto background_sequences = read_fasta(
        config.bg_path,
        config.motif_len,
        config.max_peak_len,
        log
    );
    if (foreground.empty() || background_sequences.empty()) {
        throw std::runtime_error("foreground and background must be non-empty");
    }

    const auto kmer_ratios = kmer_log_ratios(
        foreground,
        background_sequences,
        config.olig_bg
    );
    auto encoded = encode_sequences(
        foreground,
        kmer_ratios,
        config.motif_len,
        config.olig_bg
    );
    auto background = build_background_stats(
        background_sequences,
        config.max_lpd
    );

    for (const auto& sequence : encoded.records) {
        if (sequence.candidate_positions.empty()) {
            throw std::runtime_error("not enough candidate motif windows");
        }
    }

    auto feature_pool = build_feature_pool(encoded, background, config);
    log << "Feature pool size=" << feature_pool.size()
        << " best_score="
        << (feature_pool.empty() ? 0.0 : feature_pool.front().score)
        << '\n';

    Rng rng(config.seed);
    auto winners = run_island_searches(
        feature_pool,
        encoded,
        background,
        config,
        rng,
        log
    );
    write_outputs(winners, encoded, background, config, result);
    log << "Pipeline completed successfully!\n";
    return 0;
}

}  // namespace

extern "C" int sitega_train(const TrainParams* params, TrainResult* result) {
    if (result != nullptr) {
        std::memset(result, 0, sizeof(*result));
        result->status = 1;
    }
    if (params == nullptr) {
        return 1;
    }

    try {
        return train_impl(*params, result);
    } catch (const std::exception& error) {
        if (result != nullptr) {
            result->status = 1;
        }
        std::fprintf(stderr, "SiteGA training failed: %s\n", error.what());
        return 1;
    }
}
