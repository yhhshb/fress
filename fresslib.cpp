#include "fresslib.hpp"
#include "nthash.hpp"
#include "./kmc_api/kmc_file.h"

#include <set>
#include <map>
#include <chrono>

hist_t compute_histogram(std::string kmc_name)
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
	hist_t toRet;
	std::copy(histo.cbegin(), histo.cend(), std::back_inserter(toRet));
	return toRet;
}

void store_histogram(std::string histo_name, const hist_t& histo)
{
	std::ofstream hf(histo_name);
	for(auto it = histo.cbegin(); it != histo.cend(); ++it)
	{
		hf << it->first << "\t" << it->second << "\n";
	}
	hf.close();
}

hist_t load_histogram(std::string histo_name)
{
	hist_t histo;
	std::string line;
	std::ifstream histo_file(histo_name);
	while(std::getline(histo_file, line))
	{
		auto pair = str2set<uint64_t>(line, '\t');
		if(pair.size() != 2) throw std::logic_error("[Error] The specified file is not a histogram");
		histo.emplace_back(static_cast<uint32_t>(pair[0]), pair[1]);
	}
	return histo;
}

hist_t sort_histogram(const hist_t& histo)
{
	hist_t sorted_columns = histo;
	std::sort(sorted_columns.begin(), sorted_columns.end(), [](auto &left, auto &right) { return left.second > right.second; });
	return sorted_columns;
}

std::unordered_map<uint32_t, uint32_t> create_inv_index(const hist_t& sorted_histogram)
{
	std::unordered_map<uint32_t, uint32_t> toRet;
	for(uint32_t i = 0; i < sorted_histogram.size(); ++i) toRet[sorted_histogram[i].first] = i;
	return toRet;
}

void store_cmb(std::string comb_name, const comb_t& combinations)
{
	std::ofstream combo(comb_name);
	for(auto s : combinations) combo << s << "\n";
	combo.close();
}

std::vector<std::vector<uint32_t>> load_cmb_for_query(std::string comb_name)
{
	std::ifstream skdump(comb_name);
	if(not skdump.is_open()) throw std::runtime_error("Unable to open sketch combination file");
	std::vector<std::vector<uint32_t>> frequency_sets;
	std::string line;
	bool first = true;
	while(std::getline(skdump, line))
	{
		if(first or line != "") frequency_sets.push_back(str2set<uint32_t>(line, ','));
		if(first) first = false;//The first line of the file could be empty if the heavies element is implicit.
	}
	skdump.close();
	return frequency_sets;
}

void store_setmap(std::string setmap_name, uint64_t nrows, uint64_t ncolumns, const sketch_t& setmap)
{
	std::ofstream skdump(setmap_name, std::ios::binary);
	skdump.write(reinterpret_cast<char*>(&nrows), sizeof(decltype(nrows)));
	skdump.write(reinterpret_cast<char*>(&ncolumns), sizeof(decltype(ncolumns)));
	skdump.write(reinterpret_cast<char*>(const_cast<sketch_t::value_type*>(setmap.data())), setmap.size() * sizeof(sketch_t::value_type));
	skdump.close();
}

sketch_t load_setmap(std::string setmap_name, uint64_t& nrows, uint64_t& ncolumns, bool all)
{
	std::ifstream skdump(setmap_name, std::ios::binary);
	if(not skdump.is_open()) throw std::runtime_error("Unable to open sketch index file");
	skdump.read(reinterpret_cast<char*>(&nrows), sizeof(decltype(nrows)));
	skdump.read(reinterpret_cast<char*>(&ncolumns), sizeof(decltype(ncolumns)));
	sketch_t setmap(nrows * ncolumns);
	if(all) skdump.read(reinterpret_cast<char*>(setmap.data()), setmap.size() * sizeof(decltype(setmap)::value_type));
	skdump.close();
	return setmap;
}

double estimate_error(const hist_t& sorted_histo, uint64_t nrows, uint64_t ncolumns, std::vector<double>& buffer)
{
	double error, rs;
	std::size_t L = sorted_histo.size();
	if(buffer.size() < L) buffer.resize(L);
	for(std::size_t i = 0; i < L; ++i) buffer[i] = 1.0-std::pow(1.0-1.0/ncolumns, sorted_histo.at(i).second);
	error = 0;
	for(std::size_t i = 0; i < L; ++i)
	{
		rs = 0;
		for(std::size_t j = i+1; j < L; ++j)
		{
			rs += std::pow(buffer[j], nrows) * std::abs(static_cast<long long>(sorted_histo.at(j).first) - sorted_histo.at(i).first);
		}
		error += sorted_histo.at(i).second * rs;
	}
	return error;
}

void optimise_r_b(const hist_t& sorted_histo, double target_error, uint64_t& nrows, uint64_t& ncolumns)
{
	std::size_t L1_norm = 0;
	for(auto p : sorted_histo) L1_norm += p.first * p.second;
	const double threshold = L1_norm * target_error;
	std::cerr << "L1 norm of " << L1_norm << " -> " << "threshold = " << threshold << "\n";

	uint64_t r = nrows;
	uint64_t  b = ncolumns;
	bool cr = false, cb = false; 
	if(b == 0) b = sorted_histo[1].second * 1.443;
	else cb = true;
	if(r == 0) r = 1;
	else cr = true;
	std::vector<double> buffer(sorted_histo.size());
	double error = estimate_error(sorted_histo, r, b, buffer);
	std::cerr << "(" << r << ", " << b << ") -> " << error << "\n";
	while(error > threshold and not cr)
	{
		++r;
		error = estimate_error(sorted_histo, r, b, buffer);
		std::cerr << "(" << r << ", " << b << ") -> " << error << "\n";
	}
	std::size_t rb = r * b;
	bool decr = false;
	while(error < threshold and r > 1 and not cr and not cb)
	{
		decr = true;
		--r;
		b = rb / r;
		error = estimate_error(sorted_histo, r, b, buffer);
		std::cerr << "(" << r << ", " << b << ") -> " << error << "\n";
	}
	if(decr) 
	{
		++r;
		b = rb/r;
	}
	if(cr and cb and error > threshold) 
	{
		std::cerr << "Unable to achieve an error below " << threshold << " with (r, b) = (" << r << ", " << "b)\n";
		std::cerr << "The error for the given parameters is " << error << "\n";
	} 
	else 
	{
		nrows = r;
		ncolumns = b;
		std::cerr << "(r, b) = (" << nrows << ", " << ncolumns << ")\n";
	}
}

void fill_sketch_small(std::string kmc_filename, uint64_t nrows, uint64_t ncolumns, uint32_t heavy_element, comb_t& str_combinations, sketch_t& sketch)
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

	
	std::map<std::string, uint32_t> set_index;//string set of combinations used for quick search 
	set_index.emplace("", 0);
	std::vector<std::string> dsc;
	{//Begin
	std::vector<std::vector<uint32_t>> combinations;//exactly the same as set_index but used for quick modification
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

std::vector<std::string> check_sketch(std::string kmc_filename, uint64_t nrows, uint64_t ncolumns, const sketch_t& setmap, const std::vector<std::vector<uint32_t>>& frequency_sets, const std::unordered_map<uint32_t, uint32_t>& inverted_index)
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
	std::size_t delta_max = 0;
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
		if(intersection.size() > 0)
		{
			if((intersection.size() == 1 and intersection[0] != counter) or intersection.size() > 1)
			{
				++ncolls;
				dummy.clear();
				for(auto elem : intersection) dummy.push_back(inverted_index.at(elem));
				auto smallest_itr = std::max_element(dummy.cbegin(), dummy.cend());
				std::size_t idx = std::distance(dummy.cbegin(), smallest_itr);
				uint32_t qval = intersection.at(idx);
				auto delta = static_cast<std::size_t>(std::abs(static_cast<long long>(counter) - qval));
				delta_sum += delta;
				if(delta_max < delta) delta_max = delta;
			}
		}
	}
	kmcdb.Close();
	std::vector<std::string> toRet(4);
	toRet[0] = std::to_string(ncolls);
	toRet[1] = std::to_string(delta_sum);
	toRet[2] = std::to_string(static_cast<double>(delta_sum)/ncolls);
	toRet[3] = std::to_string(delta_max);
	std::cerr << "Total number of collisions: " << toRet[0] << "\n";
	std::cerr << "L1 sum of deltas: " << toRet[1] << "\n";
	std::cerr << "Average delta: " << toRet[2] << "\n";
	std::cerr << "MAX delta: " << toRet[3] << std::endl;
	//std::cerr << "Mean time to build the hash vector: " << total_hash_time / nqueries << " nanoseconds\n";
	//std::cerr << "Mean time to run the outer for loop: " << total_cycle_time / nqueries << " nanoseconds\n";
	//std::cerr << "Mean time to get set index: " << total_getidx_time / (nrows * nqueries) << " nanoseconds\n";
	//std::cerr << "Mean time to compute set intersection: " << total_intersect_time / (nrows * nqueries) << " nanoseconds\n";
	//std::cerr << "Mean time to retrieve a frequency: " << total_time / nqueries << " nanoseconds" << std::endl;
	return toRet;
}
