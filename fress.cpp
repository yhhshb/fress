#include <cstdio>
#include <cmath>
#include <cstring>

#include <numeric>
#include <sstream>
#include <fstream>
#include <exception>

#include <map>
#include <unordered_map>

#include "prettyprint.hpp"
#include "nthash.hpp"
#include "./kmc_api/kmc_file.h"

extern "C" {
#include "ketopt.h"
}

typedef std::vector<uint32_t> bucket_t;
typedef std::vector<bucket_t> sketch_t;

void print_subcommands()
{
	fprintf(stderr, "histogram\tCompute the k-mer frequency profile from a kmc database\n");
	fprintf(stderr, "sense\tBuild a probabilistic map of the frequencies using a compressed sensing framework\n");
	fprintf(stderr, "check\tCheck the ability of the map to assign frequencies\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Use the subcommand to see its documentation\n");
}

void print_histogram_help()
{
	fprintf(stderr, "histogram options:\n");
	fprintf(stderr, "i\tinput kmc database (without extentions)\n");
	fprintf(stderr, "o\toutput file. A two-column tsv file where the first column contains the frequencies and the second column the total number of k-mers having that frequency\n");
	fprintf(stderr, "h\tshows this help\n");
}

void print_sense_help()
{
	fprintf(stderr, "sense options:\n");
	fprintf(stderr, "i\tinput kmc database (without extentions)\n");
	fprintf(stderr, "o\toutput probabilistic map containing the frequencies\n");
	fprintf(stderr, "r\tnumber of independent bucket rows\n");
	fprintf(stderr, "h\tshows this help\n");
}

void print_check_help()
{
	fprintf(stderr, "check options:\n");
	fprintf(stderr, "i\tinput kmc database (without extentions)\n");
	fprintf(stderr, "d\tinput map built from the input kmc database (without extentions)\n");
	fprintf(stderr, "h\tshows this help\n");
	fprintf(stderr, "\nThe output is on stdout and are all k-mers for which the predicted frequency is wrong or there are multiple frequencies\n");
	fprintf(stderr, "Each k-mer will have its true frequency and the (wrong) intersection of frequencies\n");
}

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
	return histo;
}

int histogram_main(int argc, char* argv[])
{
	static ko_longopt_t longopts[] = {
		{NULL, 0, 0}
	};

	ketopt_t opt = KETOPT_INIT;
	int c;
	std::string kmc_filename, output_filename;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:h:", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'o') {
			output_filename = opt.arg;
		} else if (c == 'h') {
			print_histogram_help();
			return 0;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_histogram_help();
			return 1;
		}
	}
	auto histo = compute_histogram(kmc_filename);
	std::ofstream hf(output_filename);
	for(auto it = histo.cbegin(); it != histo.cend(); ++it)
	{
		hf << it->first << "\t" << it->second << "\n";
	}
	hf.close();
	return 0;
}

std::vector<std::pair<uint32_t, std::size_t>> analyze_histogram(const std::map<uint32_t, std::size_t>& histo)
{
	std::vector<std::pair<uint32_t, std::size_t> > sorted_columns(histo.size());
	for(auto hbucket : histo)
	{
		sorted_columns.push_back(hbucket);
	}
	std::sort(sorted_columns.begin(), sorted_columns.end(), [](auto &left, auto &right) {
    		return left.second > right.second;
	});
	return sorted_columns;
}

std::string set2str(const bucket_t& aSet) 
{
	return std::accumulate(
			std::cbegin(aSet), 
			std::cend(aSet), 
			std::string{}, 
			[](const std::string& a, uint32_t b) { return a.empty() ? std::to_string(b) : (a + ',' + std::to_string(b)); }
		);
};

std::vector<uint32_t> str2set(const std::string& setstr, char sep)
{
	std::vector<uint32_t> bucket;
	std::stringstream ss(setstr);
	uint32_t val;
	while(ss >> val)
	{
		bucket.push_back(val);
		if(ss.peek() == sep) ss.ignore();
	}
	return bucket;
}

int sense_main(int argc, char* argv[])
{
	static ko_longopt_t longopts[] = {
		{NULL, 0, 0}
	};

	ketopt_t opt = KETOPT_INIT;
	int c;
	std::string kmc_filename, output_filename;
	std::size_t nrows = 0;
	std::size_t ncolumns = 0;
	double delta = 0.01;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:b:r:d:h", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'o') {
			output_filename = opt.arg;
		} else if (c == 'r') {
			nrows = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'b') {
			ncolumns = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'd') {
			delta = std::stod(opt.arg);
		} else if (c == 'h') {
			print_sense_help();
			return 0;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_sense_help();
			return 1;
		}
	}

	if(kmc_filename == "" or output_filename == "") throw std::runtime_error("-i and -o are mandatory arguments");
	if(ncolumns == 0) throw std::runtime_error("-b is a mandatory argument");
	if (nrows == 0) throw std::runtime_error("-r is a mandatory argument");
	//nrows = log(sch.size() / delta);
	//auto sch = analyze_histogram(compute_histogram(kmc_filename));
	//ncolumns = static_cast<std::size_t>(sch[1].second); 
	
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
	sketch_t sketch(nrows * ncolumns);

	CKmerAPI kmer(_kmer_length);
	uint32_t counter = 0;
	char str_kmer[_kmer_length + 1];
	uint64_t hashes[nrows];

	while(kmcdb.ReadNextKmer(kmer, counter))
	{
		//fprintf(stderr, "found a low-hitter\n");
		kmer.to_string(str_kmer);
		NTM64(str_kmer, _kmer_length, nrows, hashes);
		//fprintf(stderr, "hashes computed\n");
		for(std::size_t i = 0; i < nrows; ++i)
		{
			std::size_t bucket_index = hashes[i] % ncolumns + i * ncolumns;
			//fprintf(stderr, "row %lu -> bucket %lu\n", i, bucket_index);
			auto& bucket = sketch[bucket_index];
			if(std::find(bucket.cbegin(), bucket.cend(), counter) == bucket.cend()) bucket.push_back(counter); //.insert(counter);
			//fprintf(stderr, "insertion done\n");
		}
	}
	kmcdb.Close();
	
	/*
	counter=0;
	for(const auto& bucket : sketch)
	{
		std::cout << bucket << "\n";
		++counter;
		if(counter % ncolumns == 0) std::cout << "\n\n\n\n";
	}
	*/

	std::vector<std::string> sorted_combinations;
	std::vector<uint32_t> to_be_stored;

	{//Begin
	std::unordered_map<std::string, uint32_t> combinations;
	for(auto& bucket : sketch) 
	{
		std::sort(bucket.begin(), bucket.end());
		combinations[set2str(bucket)] = 0;
	}

	for(auto it = combinations.cbegin(); it != combinations.cend(); ++it) sorted_combinations.push_back(it->first);
	std::sort(sorted_combinations.begin(), sorted_combinations.end());
	for(std::size_t value = 0; value < sorted_combinations.size(); ++value) combinations[sorted_combinations[value]] = value;

	for(const auto& bucket : sketch) to_be_stored.push_back(combinations.at(set2str(bucket)));
	}//End

	std::ofstream combo(output_filename + ".cmb.txt");
	for(auto s : sorted_combinations) combo << s << "\n";
	combo.close();

	std::ofstream skdump(output_filename + ".bin", std::ios::binary);
	skdump.write(reinterpret_cast<char*>(&nrows), sizeof(decltype(nrows)));
	skdump.write(reinterpret_cast<char*>(&ncolumns), sizeof(decltype(ncolumns)));
	skdump.write(reinterpret_cast<char*>(to_be_stored.data()), to_be_stored.size() * sizeof(decltype(to_be_stored)::value_type));
	skdump.close();

	return 0;
}

int check_main(int argc, char* argv[])
{
	static ko_longopt_t longopts[] = {
		{NULL, 0, 0}
	};

	ketopt_t opt = KETOPT_INIT;
	int c;
	std::string kmc_filename, map_filename;
	while((c = ketopt(&opt, argc, argv, 1, "i:d:h", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'd') {
			map_filename = opt.arg;
		} else if (c == 'h') {
			print_check_help();
			return 0;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_check_help();
			return 1;
		}
	}

	if(kmc_filename == "" or map_filename == "") throw std::runtime_error("-i and -d are mandatory arguments");

	std::size_t nrows, ncolumns;
	std::ifstream skdump(map_filename + ".bin", std::ios::binary);
	skdump.read(reinterpret_cast<char*>(&nrows), sizeof(decltype(nrows)));
	skdump.read(reinterpret_cast<char*>(&ncolumns), sizeof(decltype(ncolumns)));
	std::vector<uint32_t> setmap(nrows * ncolumns);
	skdump.read(reinterpret_cast<char*>(setmap.data()), setmap.size() * sizeof(decltype(setmap)::value_type));
	skdump.close();
	skdump.open(map_filename + ".cmb.txt");
	std::vector<std::vector<uint32_t>> frequency_sets;
	std::string line;
	while(std::getline(skdump, line))
	{
		if(line != "")
		{
			//auto dummy = str2set(line, ',');
			//std::cerr << line << " -> " << dummy << "\n";
			frequency_sets.push_back(str2set(line, ','));
		}
	}
	skdump.close();

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
	while(kmcdb.ReadNextKmer(kmer, counter))
	{
		kmer.to_string(str_kmer);
		NTM64(str_kmer, _kmer_length, nrows, hashes);
		for(std::size_t i = 0; i < nrows; ++i)
		{
			std::size_t bucket_index = hashes[i] % ncolumns + i * ncolumns;
			if(i == 0) intersection = frequency_sets[setmap[bucket_index]];
			else {
				bucket_t dummy;
				const bucket_t& current = frequency_sets[setmap[bucket_index]];
				std::set_intersection(intersection.cbegin(), intersection.cend(), current.cbegin(), current.cend(), std::back_inserter(dummy));
				intersection = dummy;
			}
		}
		
		//bool wrong_low_hitter = intersection.size() == 0 and counter != 1;
		//bool wrong_value = intersection.size() == 1 and counter != intersection[0];
		//bool unsolved_collisions = intersection.size() > 1;
		if(intersection.size() > 0 and counter != intersection.back())
		{
			std::cout << str_kmer << ", " << counter << ", " << intersection << "\n";
		}
	}
	std::cout << std::endl;
	kmcdb.Close();

	return 0;
}

int main(int argc, char* argv[])
{
	ketopt_t om = KETOPT_INIT;
	int c;
	while((c = ketopt(&om, argc, argv, 0, "x", 0)) >= 0) {} //parse main options and find subcommand
	if (om.ind == argc) {
		fprintf(stderr, "[Error] Subcommand unavailable\n\n");
		print_subcommands();
		return 1;
	}

	if (std::strcmp(argv[om.ind], "histogram") == 0) {
		return histogram_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "sense") == 0) {
		return sense_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "check") == 0) {
		return check_main(argc - om.ind, &argv[om.ind]);
	} else {
		fprintf(stderr, "Missing subcommand\n\n");
		print_subcommands();
	}
	return 0;
}
