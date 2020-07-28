#include "fresslib.hpp"
#include "nthash.hpp"
#include "./kmc_api/kmc_file.h"

#include <chrono>

std::map<uint32_t, std::size_t> compute_histogram(std::string kmc_name)
{
	CKMCFile kmcdb;
	if (!kmcdb.OpenForListing(kmc_name)) {
		throw std::runtime_error("Unable to open the database\n");
	}
	
	unsigned int _kmer_length;
	unsigned int _mode;
	unsigned int _counter_size;
	unsigned int _lut_prefix_length;
	unsigned int _signature_len;
	unsigned int _min_count;
	unsigned long long _max_count;
	unsigned long long _total_kmers;
	kmcdb.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

	std::map<uint32_t, std::size_t> histo;
	CKmerAPI kmer(_kmer_length);
	uint32_t counter = 0;
	while(kmcdb.ReadNextKmer(kmer, counter))
	{
		++histo[counter];
	}

	kmcdb.Close();
	//std::vector<std::pair<uint32_t, std::size_t>> sorted_histo;
	//std::copy(histo.cbegin(), histo.cend(), std::back_inserter(sorted_histo));
	//std::sort(sorted_histo.begin(), sorted_histo.end(), [](auto& left, auto& right) { return left.second > right.second; });
	return histo;
}

std::map<uint32_t, std::size_t> load_histogram(std::string histo_name)
{
	std::map<uint32_t, std::size_t> histo;
	std::string line;
	std::stringstream ss;
	std::ifstream histo_file(histo_name);
	while(std::getline(histo_file, line))
	{
		auto pair = str2set<uint64_t>(line, '\t');
		if(pair.size() != 2) throw std::logic_error("[Error] The specified file is not a histogram");
		histo[static_cast<uint32_t>(pair[0])] = pair[1];
	}
	return histo;
}

std::vector<std::pair<uint32_t, std::size_t>> sort_histogram(const std::map<uint32_t, std::size_t>& histo)
{
	std::vector<std::pair<uint32_t, std::size_t>> sorted_columns(histo.size());
	for(auto hbucket : histo) sorted_columns.push_back(hbucket);
	std::sort(sorted_columns.begin(), sorted_columns.end(), [](auto &left, auto &right) { return left.second > right.second; });
	return sorted_columns;
}

void fill_sketch_small(std::string kmc_filename, std::size_t nrows, std::size_t ncolumns, std::vector<std::string>& str_combinations, std::vector<uint32_t>& sketch, uint32_t heavy_element)
{
	CKMCFile kmcdb;
	if (!kmcdb.OpenForListing(kmc_filename)) {
		throw std::runtime_error("Unable to open the database\n");
	}

	unsigned int _kmer_length;
	unsigned int _mode;
	unsigned int _counter_size;
	unsigned int _lut_prefix_length;
	unsigned int _signature_len;
	unsigned int _min_count;
	unsigned long long _max_count;
	unsigned long long _total_kmers;
	kmcdb.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

	fprintf(stderr, "Starting filling the sketch of size %lu x %lu = %lu\n", nrows, ncolumns, nrows * ncolumns);

	CKmerAPI kmer(_kmer_length);
	uint32_t counter = 0;
	char str_kmer[_kmer_length + 1];
	uint64_t hashes[nrows];

	
	std::map<std::string, uint32_t> set_index;
	set_index.emplace("", 0);
	std::vector<std::string> dsc;
	{//Begin
	std::vector<std::vector<uint32_t>> combinations;
	combinations.push_back(std::vector<uint32_t>(0));
	while(kmcdb.ReadNextKmer(kmer, counter))
	{
		if (heavy_element != counter)
		{
			kmer.to_string(str_kmer);
			NTM64(str_kmer, _kmer_length, nrows, hashes);
			for(std::size_t i = 0; i < nrows; ++i)
			{
				std::size_t bucket_index = hashes[i] % ncolumns + i * ncolumns;
				const auto& vset = combinations[sketch[bucket_index]];
				if(std::find(vset.cbegin(), vset.cend(), counter) == vset.cend())
				{
					auto key = vset;
					key.insert(std::upper_bound(key.begin(), key.end(), counter), counter);
					auto skey = set2str(key);
					const auto invkey_it = set_index.find(skey);
					if(invkey_it != set_index.cend())
					{
						sketch[bucket_index] = invkey_it->second;
					}
					else
					{
						combinations.push_back(key);
						set_index.emplace(skey, combinations.size() - 1);
						sketch[bucket_index] = combinations.size() - 1;
					}
				}
			}
		}
	}
	kmcdb.Close();
	for(const auto& combo : combinations) dsc.push_back(set2str(combo));
	}//End

	//There might be some unused combinations in the combinations vector, te following code is used to remove them
	std::set<std::string> str_set_combos;
	for(auto idx : sketch) str_set_combos.insert(dsc[idx]);
	str_combinations = std::vector<std::string>();
	std::copy(str_set_combos.cbegin(), str_set_combos.cend(), std::back_inserter(str_combinations));
	for(uint32_t i = 0; i < str_combinations.size(); ++i) set_index[str_combinations[i]] = i;
	for(auto& idx : sketch) idx = set_index[dsc[idx]];
}

void check_sketch(std::string kmc_filename, std::size_t nrows, std::size_t ncolumns, const std::vector<uint32_t>& setmap, const std::vector<std::vector<uint32_t>>& frequency_sets)
{
	CKMCFile kmcdb;
	if (!kmcdb.OpenForListing(kmc_filename)) {
		throw std::runtime_error("Unable to open the database\n");
	}

	unsigned int _kmer_length;
	unsigned int _mode;
	unsigned int _counter_size;
	unsigned int _lut_prefix_length;
	unsigned int _signature_len;
	unsigned int _min_count;
	unsigned long long _max_count;
	unsigned long long _total_kmers;
	kmcdb.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

	CKmerAPI kmer(_kmer_length);
	uint32_t counter = 0;
	char str_kmer[_kmer_length + 1];
	uint64_t hashes[nrows];
	bucket_t intersection;
	std::size_t ncolls = 0;
	//using namespace std::chrono;
	while(kmcdb.ReadNextKmer(kmer, counter))
	{
		kmer.to_string(str_kmer);
		NTM64(str_kmer, _kmer_length, nrows, hashes);
		for(std::size_t i = 0; i < nrows; ++i)
		{
			//auto start = high_resolution_clock::now();
			std::size_t bucket_index = hashes[i] % ncolumns + i * ncolumns;
			if(i == 0) intersection = frequency_sets[setmap[bucket_index]];
			else {
				bucket_t dummy;
				const bucket_t& current = frequency_sets[setmap[bucket_index]];
				std::set_intersection(intersection.cbegin(), intersection.cend(), current.cbegin(), current.cend(), std::back_inserter(dummy));
				intersection = dummy;
			}
			//std::cerr << "time for one query: " << std::chrono::duration_cast<nanoseconds>(system_clock::now() - start).count() << "\n";
		}
		//bool wrong_low_hitter = intersection.size() == 0 and counter != 1;
		//bool wrong_value = (intersection.size() == 1) and (counter != intersection[0]);
		//bool unsolved_collisions = intersection.size() > 1;
		if(intersection.size() > 0 and counter != intersection.back()) //FIXME use the min(histo[intersection]) for selecting the right probability
		{
			++ncolls;
			if(std::abs(static_cast<long long>(counter) - intersection.back()) > 3)
				std::cout << str_kmer << ", " << counter << ", " << intersection << "\n";
		}
	}
	kmcdb.Close();
	std::cout << std::endl;
	std::cerr << "Total number of collisions: " << ncolls << std::endl;
}
