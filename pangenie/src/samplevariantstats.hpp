#ifndef SAMPLEVARIANTSTATS_HPP
#define SAMPLEVARIANTSTATS_HPP

#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <vector>

/**
 * Per-sample k-mer statistics of one (possibly clustered) variant.
 *
 * Variants closer than k-1 bases are genotyped as one cluster and written to
 * the VCF as separate records ("singletons"). For every sample and every
 * singleton this class stores:
 *   - the number of unique k-mers (UK in the VCF),
 *   - the local k-mer coverage (KC),
 *   - one k-mer count per allele of that singleton (AK).
 *
 * Singletons of one cluster can have different numbers of alleles. Each
 * singleton therefore has its own block width, and the position of a block is
 * always computed here, never by the caller. All accesses are bounds-checked.
 **/
class SampleVariantStats {
public:
	SampleVariantStats() : nr_samples(0), block_total(0) {}

	/** true until init() has been called. **/
	bool empty() const { return this->offsets.empty(); }

	/**
	 * @param nr_samples number of samples
	 * @param alleles_per_singleton number of alleles of each singleton, in VCF order
	 **/
	void init(size_t nr_samples, const std::vector<size_t>& alleles_per_singleton) {
		this->nr_samples = nr_samples;
		this->nr_alleles = alleles_per_singleton;
		this->offsets.assign(alleles_per_singleton.size(), 0);
		this->block_total = 0;
		for (size_t j = 0; j < alleles_per_singleton.size(); ++j) {
			this->offsets[j] = this->block_total;
			this->block_total += HEADER + alleles_per_singleton[j];
		}
		this->data.assign(this->block_total * nr_samples, 0);
	}

	size_t nr_of_samples() const { return this->nr_samples; }
	size_t nr_of_singletons() const { return this->offsets.size(); }
	size_t nr_of_alleles(size_t singleton) const { check_singleton(singleton); return this->nr_alleles[singleton]; }

	/** store the statistics of one sample for one singleton. allele_counts must hold one value per allele. **/
	void set(size_t sample, size_t singleton, unsigned short unique_kmers, unsigned short coverage, const std::vector<unsigned short>& allele_counts) {
		size_t base = index(sample, singleton);
		if (allele_counts.size() != this->nr_alleles[singleton]) {
			std::ostringstream oss;
			oss << "SampleVariantStats::set: got " << allele_counts.size() << " allele counts for a singleton with " << this->nr_alleles[singleton] << " alleles.";
			throw std::invalid_argument(oss.str());
		}
		this->data[base + UNIQUE_KMERS] = unique_kmers;
		this->data[base + COVERAGE] = coverage;
		for (size_t a = 0; a < allele_counts.size(); ++a) {
			this->data[base + HEADER + a] = allele_counts[a];
		}
	}

	unsigned short unique_kmers(size_t sample, size_t singleton) const { return this->data[index(sample, singleton) + UNIQUE_KMERS]; }
	unsigned short coverage(size_t sample, size_t singleton) const { return this->data[index(sample, singleton) + COVERAGE]; }
	unsigned short allele_kmer_count(size_t sample, size_t singleton, size_t allele) const {
		size_t base = index(sample, singleton);
		if (allele >= this->nr_alleles[singleton]) {
			std::ostringstream oss;
			oss << "SampleVariantStats::allele_kmer_count: allele " << allele << " requested, singleton has " << this->nr_alleles[singleton] << ".";
			throw std::out_of_range(oss.str());
		}
		return this->data[base + HEADER + allele];
	}

private:
	static const size_t UNIQUE_KMERS = 0;
	static const size_t COVERAGE = 1;
	static const size_t HEADER = 2;

	void check_singleton(size_t singleton) const {
		if (singleton >= this->offsets.size()) {
			std::ostringstream oss;
			oss << "SampleVariantStats: singleton " << singleton << " requested, variant has " << this->offsets.size() << ".";
			throw std::out_of_range(oss.str());
		}
	}

	/** position of the block of (sample, singleton) in data. **/
	size_t index(size_t sample, size_t singleton) const {
		check_singleton(singleton);
		if (sample >= this->nr_samples) {
			std::ostringstream oss;
			oss << "SampleVariantStats: sample " << sample << " requested, " << this->nr_samples << " stored.";
			throw std::out_of_range(oss.str());
		}
		return sample * this->block_total + this->offsets[singleton];
	}

	size_t nr_samples;
	// sum of the block widths of all singletons, i.e. the width of one sample
	size_t block_total;
	// number of alleles of each singleton
	std::vector<size_t> nr_alleles;
	// start of each singleton's block within one sample
	std::vector<size_t> offsets;
	std::vector<unsigned short> data;
};

#endif // SAMPLEVARIANTSTATS_HPP
