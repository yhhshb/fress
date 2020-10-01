#include <cmath>
#include <sstream>
#include <numeric>
#include <exception>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>

typedef std::vector<std::pair<uint32_t, std::size_t>> hist_t;
typedef std::vector<std::string> comb_t;
typedef std::vector<uint32_t> bucket_t;
typedef std::vector<uint32_t> sketch_t;

template <class SetType>
std::string set2str(const SetType& aSet) 
{
	return std::accumulate(
			std::cbegin(aSet), 
			std::cend(aSet), 
			std::string{}, 
			[](const std::string& a, typename SetType::value_type b) { return a.empty() ? std::to_string(b) : (a + ',' + std::to_string(b)); }
		);
};

template <typename T>
std::vector<T> str2set(const std::string& setstr, char sep)
{
	std::vector<T> bucket;
	std::stringstream ss(setstr);
	T val;
	while(ss >> val)
	{
		bucket.push_back(val);
		if(ss.peek() == sep) ss.ignore();
	}
	return bucket;
}

hist_t compute_histogram(std::string kmc_name);
void store_histogram(std::string histo_name, const hist_t& histo);
hist_t load_histogram(std::string histo_name);
hist_t sort_histogram(const hist_t& histo);
std::unordered_map<uint32_t, uint32_t> create_inv_index(const hist_t& sorted_histogram);

void store_cmb(std::string comb_name, const comb_t& combinations);
std::vector<std::vector<uint32_t>> load_cmb_for_query(std::string comb_name);
void store_setmap(std::string setmap_name, uint64_t nrows, uint64_t ncolumns, const sketch_t& setmap);
sketch_t load_setmap(std::string setmap_name, uint64_t& nrows, uint64_t& ncolumns, bool all);

double estimate_error(const hist_t& sorted_hist, uint64_t nrows, uint64_t ncolumns);
std::size_t optimise_r_b(const hist_t& sorted_histo, double target_error, uint64_t& nrows, uint64_t& ncolumns);
void fill_sketch_small(std::string kmc_name, uint64_t nrows, uint64_t ncolumns, uint32_t heavy_element, comb_t& combinations, sketch_t& setmap);
std::vector<std::string> check_sketch(std::string kmc_name, uint64_t nrows, uint64_t ncolumns, const sketch_t& setmap, const std::vector<std::vector<uint32_t>>& frequency_sets, const std::unordered_map<uint32_t, uint32_t>& inverted_index);

void fill_cms_sketch(std::string kmc_filename, uint64_t nrows, uint64_t ncolumns, std::vector<uint32_t>& cms);
std::vector<std::string> check_cm_sketch(std::string kmc_filename, uint64_t nrows, uint64_t ncolumns, const std::vector<uint32_t>& cms);
