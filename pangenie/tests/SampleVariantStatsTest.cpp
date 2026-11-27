// Regression test for the per-sample variant statistics layout.
//
// Background: statistics of clustered variants were stored in blocks sized by
// the cluster's allele count and read back with the single variant's allele
// count. When the two counts differ, UK/KC/AK can be read from the wrong
// positions.
//
// This test is self-contained (no test framework, no other project sources)
// so that it builds even while the rest of pangenie/tests is disabled.

#include <iostream>
#include <stdexcept>
#include <vector>
#include "samplevariantstats.hpp"

static int failures = 0;

#define CHECK(cond) \
	do { if (!(cond)) { ++failures; std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond << std::endl; } } while (0)

#define CHECK_THROWS(expr, exception_type) \
	do { bool thrown = false; try { expr; } catch (const exception_type&) { thrown = true; } \
	     if (!thrown) { ++failures; std::cerr << __FILE__ << ":" << __LINE__ << ": expected " #exception_type " from: " #expr << std::endl; } } while (0)

// a value that is different for every (sample, singleton, field)
static unsigned short value(size_t sample, size_t singleton, size_t field) {
	return static_cast<unsigned short>(1000 * (sample + 1) + 100 * (singleton + 1) + field);
}

static void fill(SampleVariantStats& stats, const std::vector<size_t>& alleles) {
	for (size_t s = 0; s < stats.nr_of_samples(); ++s) {
		for (size_t j = 0; j < alleles.size(); ++j) {
			std::vector<unsigned short> counts;
			for (size_t a = 0; a < alleles[j]; ++a) counts.push_back(value(s, j, 2 + a));
			stats.set(s, j, value(s, j, 0), value(s, j, 1), counts);
		}
	}
}

static void check_roundtrip(const std::vector<size_t>& alleles, size_t nr_samples) {
	SampleVariantStats stats;
	CHECK(stats.empty());
	stats.init(nr_samples, alleles);
	CHECK(!stats.empty());
	CHECK(stats.nr_of_samples() == nr_samples);
	CHECK(stats.nr_of_singletons() == alleles.size());
	fill(stats, alleles);
	for (size_t s = 0; s < nr_samples; ++s) {
		for (size_t j = 0; j < alleles.size(); ++j) {
			CHECK(stats.nr_of_alleles(j) == alleles[j]);
			CHECK(stats.unique_kmers(s, j) == value(s, j, 0));
			CHECK(stats.coverage(s, j) == value(s, j, 1));
			for (size_t a = 0; a < alleles[j]; ++a) {
				CHECK(stats.allele_kmer_count(s, j, a) == value(s, j, 2 + a));
			}
		}
	}
}

int main() {
	// a single biallelic variant
	check_roundtrip({2}, 1);
	check_roundtrip({2}, 4);

	// the case that used to fail: a cluster of two biallelic variants. Its
	// haplotypes carry three distinct allele combinations, so the cluster has
	// three alleles while each singleton has two. Every sample and every
	// singleton must read back its own values.
	check_roundtrip({2, 2}, 4);
	check_roundtrip({2, 2, 2}, 8);

	// singletons of one cluster with different numbers of alleles
	check_roundtrip({2, 3}, 3);
	check_roundtrip({4, 2, 3}, 5);

	// the second singleton of the first sample must not alias the first
	// singleton of the second sample (what the old index arithmetic did)
	{
		SampleVariantStats stats;
		stats.init(2, {2, 2});
		fill(stats, {2, 2});
		CHECK(stats.unique_kmers(0, 1) != stats.unique_kmers(1, 0));
		CHECK(stats.unique_kmers(0, 1) == value(0, 1, 0));
		CHECK(stats.unique_kmers(1, 0) == value(1, 0, 0));
	}

	// values are zero until set
	{
		SampleVariantStats stats;
		stats.init(2, {2, 3});
		CHECK(stats.unique_kmers(1, 1) == 0);
		CHECK(stats.allele_kmer_count(1, 1, 2) == 0);
	}

	// out-of-range accesses are reported, never read
	{
		SampleVariantStats stats;
		CHECK_THROWS(stats.unique_kmers(0, 0), std::out_of_range);
		stats.init(2, {2, 3});
		CHECK_THROWS(stats.unique_kmers(2, 0), std::out_of_range);
		CHECK_THROWS(stats.unique_kmers(0, 2), std::out_of_range);
		CHECK_THROWS(stats.coverage(0, 2), std::out_of_range);
		CHECK_THROWS(stats.allele_kmer_count(0, 0, 2), std::out_of_range);
		CHECK_THROWS(stats.allele_kmer_count(0, 1, 3), std::out_of_range);
		CHECK_THROWS(stats.nr_of_alleles(2), std::out_of_range);
		CHECK_THROWS(stats.set(0, 0, 1, 1, {1, 2, 3}), std::invalid_argument);
		CHECK_THROWS(stats.set(0, 2, 1, 1, {1, 2}), std::out_of_range);
	}

	if (failures == 0) {
		std::cout << "SampleVariantStatsTest: all checks passed" << std::endl;
		return 0;
	}
	std::cerr << "SampleVariantStatsTest: " << failures << " check(s) failed" << std::endl;
	return 1;
}
