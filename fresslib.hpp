#include <cmath>
#include <sstream>
#include <numeric>
#include <exception>
#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>
#include "prettyprint.hpp"

typedef std::vector<uint32_t> bucket_t;
typedef std::vector<bucket_t> sketch_t;

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

std::map<uint32_t, std::size_t> compute_histogram(std::string kmc_name);
std::map<uint32_t, std::size_t> load_histogram(std::string histo_name);
std::vector<std::pair<uint32_t, std::size_t>> sort_histogram(const std::map<uint32_t, std::size_t>& histo);
std::unordered_map<uint32_t, uint32_t> create_inv_index(const std::vector<std::pair<uint32_t, std::size_t>>& sorted_histogram);
void fill_sketch_small(std::string kmc_name, std::size_t nrows, std::size_t ncolumns, uint32_t heavy_element, std::vector<std::string>& combinations, std::vector<uint32_t>& sketch);
void check_sketch(std::string kmc_name, std::size_t nrows, std::size_t ncolumns, const std::vector<uint32_t>& setmap, const std::vector<std::vector<uint32_t>>& frequency_sets, const std::unordered_map<uint32_t, uint32_t>& inverted_index);
