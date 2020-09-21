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
	std::vector<std::pair<uint32_t, std::size_t>> sorted_columns;
	for(auto hbucket : histo) sorted_columns.push_back(hbucket);
	std::sort(sorted_columns.begin(), sorted_columns.end(), [](auto &left, auto &right) { return left.second > right.second; });
	return sorted_columns;
}

std::unordered_map<uint32_t, uint32_t> create_inv_index(const std::vector<std::pair<uint32_t, std::size_t>>& sorted_histogram)//FIXME use an Elias-Fano inverted index for query speed
{
	std::unordered_map<uint32_t, uint32_t> toRet;
	for(uint32_t i = 0; i < sorted_histogram.size(); ++i)
	{
		toRet[sorted_histogram[i].first] = i;
	}
	return toRet;
}

void fill_sketch_small(std::string kmc_filename, std::size_t nrows, std::size_t ncolumns, uint32_t heavy_element, std::vector<std::string>& str_combinations, std::vector<uint32_t>& sketch)
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

	//There might be some unused combinations in the combinations vector, the following code is used to remove them
	std::set<std::string> str_set_combos;
	for(auto idx : sketch) str_set_combos.insert(dsc[idx]);
	str_combinations = std::vector<std::string>();
	std::copy(str_set_combos.cbegin(), str_set_combos.cend(), std::back_inserter(str_combinations));
	for(uint32_t i = 0; i < str_combinations.size(); ++i) set_index[str_combinations[i]] = i;
	for(auto& idx : sketch) idx = set_index[dsc[idx]];
}

void check_sketch(std::string kmc_filename, std::size_t nrows, std::size_t ncolumns, const std::vector<uint32_t>& setmap, const std::vector<std::vector<uint32_t>>& frequency_sets, const std::unordered_map<uint32_t, uint32_t>& inverted_index)
{
	using namespace std::chrono;
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
	std::vector<const bucket_t*> sets(nrows);
	bucket_t intersection, dummy;

	std::size_t ncolls = 0;
	std::size_t delta_sum = 0;
	std::size_t nqueries = 0;
	std::size_t bucket_index;

	//std::size_t total_hash_time = 0;
	//std::size_t total_cycle_time = 0;
	//std::size_t total_getidx_time = 0;
	//std::size_t total_intersect_time = 0;
	//std::size_t total_time = 0;
	while(kmcdb.ReadNextKmer(kmer, counter))
	{
		kmer.to_string(str_kmer);
		//auto start = high_resolution_clock::now();
		NTM64(str_kmer, _kmer_length, nrows, hashes);
		//total_hash_time += duration_cast<nanoseconds>(high_resolution_clock::now() - start).count();
		//auto start2 = high_resolution_clock::now();
		for(std::size_t i = 0; i < nrows; ++i)
		{
			//auto getidx_start = high_resolution_clock::now();
			bucket_index = hashes[i] % ncolumns + i * ncolumns;
			//total_getidx_time += duration_cast<nanoseconds>(high_resolution_clock::now() - getidx_start).count();
			//auto intersect_start = high_resolution_clock::now();
			if(i == 0) intersection = frequency_sets[setmap[bucket_index]];
			else {
				dummy.clear();
				const bucket_t& current = frequency_sets[setmap[bucket_index]];
				std::set_intersection(intersection.cbegin(), intersection.cend(), current.cbegin(), current.cend(), std::back_inserter(dummy));
				std::swap(intersection, dummy);
			}
			//total_intersect_time += duration_cast<nanoseconds>(high_resolution_clock::now() - intersect_start).count(); 
		}
		//total_cycle_time += duration_cast<nanoseconds>(high_resolution_clock::now() - start2).count();
		//total_time += duration_cast<nanoseconds>(high_resolution_clock::now() - start).count();
		++nqueries;
		if(intersection.size() > 0 and counter != intersection.back())
		{
			++ncolls;
			dummy.clear();
			for(auto elem : intersection) dummy.push_back(inverted_index.at(elem));
			auto smallest_itr = std::max_element(dummy.cbegin(), dummy.cend());
			std::size_t idx = std::distance(dummy.cbegin(), smallest_itr);
			uint32_t qval = intersection.at(idx);
			auto delta = std::abs(static_cast<long long>(counter) - qval);
			delta_sum += delta;	
			if(delta > 0) std::cout << str_kmer << ", " << counter << ", " << intersection << ", " << qval << ", " << dummy << ", " << idx << "\n";
		}
	}
	kmcdb.Close();
	std::cout << std::endl;
	std::cerr << "Total number of collisions: " << ncolls << "\n";
	std::cerr << "L1 norm of the errors: " << delta_sum << "\n";
	//std::cerr << "Mean time to build the hash vector: " << total_hash_time / nqueries << " nanoseconds\n";
	//std::cerr << "Mean time to run the outer for loop: " << total_cycle_time / nqueries << " nanoseconds\n";
	//std::cerr << "Mean time to get set index: " << total_getidx_time / (nrows * nqueries) << " nanoseconds\n";
	//std::cerr << "Mean time to compute set intersection: " << total_intersect_time / (nrows * nqueries) << " nanoseconds\n";
	//std::cerr << "Mean time to retrieve a frequency: " << total_time / nqueries << " nanoseconds" << std::endl;
}

//bool wrong_low_hitter = intersection.size() == 0 and counter != 1;
//bool wrong_value = (intersection.size() == 1) and (counter != intersection[0]);
//bool unsolved_collisions = intersection.size() > 1;
