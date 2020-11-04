#include <cstring>
#include "fresslib.hpp"

#include "BooPHF.hpp"
#include "KMCRangeWrapper.hpp"

#include "prettyprint.hpp"

extern "C" {
#include "ketopt.h"
}

void print_subcommands()
{
	fprintf(stderr, "histogram\tCompute the k-mer frequency profile from a kmc database\n");
	fprintf(stderr, "sense\tBuild a probabilistic map of the frequencies using a compressed sensing framework\n");
	fprintf(stderr, "check\tCheck the ability of the map to assign frequencies\n");
	fprintf(stderr, "cms\t build a count-min sketch dimensioned as a set-min sketch\n");
	fprintf(stderr, "cmschk\t check a count-min sketch\n");
	fprintf(stderr, "bbhash\t build a BBHash MPHF + external frequency array\n");//No need for checking because error = 0
	fprintf(stderr, "info\t get r, b and other useful information about one sketch\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Use the subcommand to see its documentation\n");
}

void print_histogram_help()
{
	fprintf(stderr, "histogram options:\n");
	fprintf(stderr, "\t-i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "\t-o\toutput file. A two-column tsv file where the first column contains the frequencies and the second column the total number of k-mers having that frequency\n");
	fprintf(stderr, "\t-h\tshows this help\n");
}

void print_sense_help()
{
	fprintf(stderr, "sense options:\n");
	fprintf(stderr, "\t-i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "\t-o\toutput probabilistic map containing the frequencies\n");
	fprintf(stderr, "\t-e\tepsilon approximation of the L1 sum of errors [0.01]\n");
	fprintf(stderr, "\t-s\tinput histogram generated from the input kmc database using the <histogram> subcommand. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-r\tnumber of independent bucket rows. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-b\tnumber of columns. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-h\tshows this help\n");
}

void print_check_help()
{
	fprintf(stderr, "check options:\n");
	fprintf(stderr, "\t-i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "\t-d\tinput map built from the input kmc database (without extensions)\n");
	fprintf(stderr, "\t-h\tshows this help\n");
	fprintf(stderr, "\nHuman-readable output on stderr, script-friendly output on stdout\n");
}

void print_info_help()
{
	fprintf(stderr, "info options:\n");
	fprintf(stderr, "\t-d\tinput map built from a kmc database (without extensions)\n");
	fprintf(stderr, "\t-h\tshows this help\n");
	fprintf(stderr, "\nOutput the parameters used for construction and some other informations\n");
}

void print_mms_help()
{
	fprintf(stderr, "mms options:\n");
	fprintf(stderr, "\t-i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "\t-o\toutput probabilistic map containing the frequencies\n");
	fprintf(stderr, "\t-e\tepsilon approximation ??? of the L1 sum of errors [0.01]\n");
	fprintf(stderr, "\t-s\tinput histogram generated from the input kmc database using the <histogram> subcommand. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-r\tnumber of independent bucket rows. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-b\tnumber of columns. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-h\tshows this help\n");
}

void print_cms_help()
{
	fprintf(stderr, "cms options:\n");
	fprintf(stderr, "\t-i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "\t-o\toutput count-min sketch binary file\n");
	fprintf(stderr, "\t-e\tepsilon approximation of the L1 sum of errors [0.01]\n");
	fprintf(stderr, "\t-s\tinput histogram generated from the input kmc database using the <histogram> subcommand. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-r\tnumber of independent bucket rows. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-b\tnumber of columns. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-g\tcounter value to be ignored during construction\n");
	fprintf(stderr, "\t-h\tshows this help\n");

}

void print_cmschk_help()
{
	fprintf(stderr, "check options:\n");
	fprintf(stderr, "\t-i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "\t-d\tinput count-min sketch built from the input kmc database\n");
	fprintf(stderr, "\t-g\tdefault value to be restored when encountering a minimum of 0\n");
	fprintf(stderr, "\t-h\tshows this help\n");
	fprintf(stderr, "\nHuman-readable output on stderr, script-friendly output on stdout\n");
}

void print_bbhash_help()
{
	fprintf(stderr, "mphf options:\n");
	fprintf(stderr, "\t-i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "\t-o\toutput name for mphf-related structures (mphf + arrayof frequencies)\n");
	fprintf(stderr, "\t-h\tshow this help\n");
	fprintf(stderr, "\nscript-friendly output on stdout\n");
}

int histogram_main(int argc, char* argv[])
{
	std::string kmc_filename, output_filename;

	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:h", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'o') {
			output_filename = opt.arg;
		} else if (c == 'h') {
			print_histogram_help();
			return EXIT_SUCCESS;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_histogram_help();
			return EXIT_FAILURE;
		}
	}
	auto histo = compute_histogram(kmc_filename);
	store_histogram(output_filename, histo);
	return EXIT_SUCCESS;
}

int sense_main(int argc, char* argv[])
{
	std::string kmc_filename, output_filename;
	std::vector<std::pair<uint32_t, std::size_t>> sorted_hist;
	uint64_t nrows = 0;
	uint64_t ncolumns = 0;
	double epsilon = 0.01;
	
	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:s:r:b:e:h", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'o') {
			output_filename = opt.arg;
		} else if (c == 's') {
			sorted_hist = sort_histogram(load_histogram(opt.arg));
		} else if (c == 'r') {
			nrows = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'b') {
			ncolumns = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'e') {
			epsilon = std::stod(opt.arg);
		} else if (c == 'h') {
			print_sense_help();
			return EXIT_SUCCESS;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_sense_help();
			return EXIT_FAILURE;
		}
	}

	if(kmc_filename == "" or output_filename == "") throw std::runtime_error("-i and -o are mandatory arguments");
	if(sorted_hist.size() == 0) sorted_hist = sort_histogram(compute_histogram(kmc_filename));
	std::size_t L1_norm = optimise_r_b(sorted_hist, epsilon, nrows, ncolumns);
	
	fprintf(stderr, "Starting filling the sketch of size %lu x %lu = %lu\n", nrows, ncolumns, nrows * ncolumns);

	std::vector<std::string> combinations;
	std::vector<uint32_t> sketch(nrows * ncolumns, 0);
	fill_sketch_small(kmc_filename, nrows, ncolumns, sorted_hist.size() != 0 ? sorted_hist[0].first : std::numeric_limits<uint32_t>::max(), combinations, sketch);
	
	store_histogram(output_filename + ".shist.txt", sorted_hist);
	store_cmb(output_filename + ".cmb.txt", combinations);
	store_setmap(output_filename + ".bin", nrows, ncolumns, sketch);
	fprintf(stdout, "%lu %lu", L1_norm, nrows * ncolumns);//script-friendly output
	return EXIT_SUCCESS;
}

int check_main(int argc, char* argv[])
{
	std::string kmc_filename, map_filename;
	uint64_t nrows, ncolumns;
	uint8_t additional_opts = 0;
	std::size_t merged = 0;
	double freq = 0.0;

	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "i:d:g:f:h", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'd') {
			map_filename = opt.arg;
		} else if (c == 'h') {
			print_check_help();
			return EXIT_SUCCESS;
		} else if (c == 'g') {
			merged = std::stoul(opt.arg, nullptr, 0);
			++additional_opts;
		} else if (c == 'f') {
			freq = std::stod(opt.arg);
			++additional_opts;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_check_help();
			return EXIT_FAILURE;
		}
	}

	if(kmc_filename == "" or map_filename == "") throw std::runtime_error("-i and -d are mandatory arguments");
	if(additional_opts == 1) throw std::runtime_error("-g and -f must be used together at the same time");
	auto sorted_histogram = sort_histogram(load_histogram(map_filename + ".shist.txt"));
	std::unordered_map<uint32_t, uint32_t> invidx = create_inv_index(sorted_histogram);
	auto frequency_sets = load_cmb_for_query(map_filename + ".cmb.txt");
	auto setmap = load_setmap(map_filename + ".bin", nrows, ncolumns, true);
	std::vector<std::string> rvals;
	if(additional_opts == 0) rvals = check_sketch(kmc_filename, nrows, ncolumns, setmap, frequency_sets, invidx);
	else 
	{
		std::vector<uint32_t> mcols(merged);
		for(std::size_t i = 0; i < merged; ++i) mcols[i] = sorted_histogram.at(i).first;
		rvals = check_sketch_merge(kmc_filename, nrows, ncolumns, setmap, frequency_sets, invidx, mcols, freq);
	}
	fprintf(stdout, "%s %s %s %s %s", rvals[0].c_str(), rvals[1].c_str(), rvals[2].c_str(), rvals[3].c_str(), rvals[4].c_str());//script-friendly output	
	return EXIT_SUCCESS;
}

int info_main(int argc, char* argv[])
{
	std::string map_filename;
	uint64_t nrows, ncolumns;
	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "d:h", longopts)) >= 0)
	{
		if (c == 'd') {
			map_filename = opt.arg;
		} else if (c == 'h') {
			print_info_help();
			return EXIT_SUCCESS;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_info_help();
			return EXIT_FAILURE;
		}
	}

	if(map_filename == "") throw std::runtime_error("-d is a mandatory argument");
	auto sorted_hist = sort_histogram(load_histogram(map_filename + ".shist.txt"));
	auto frequency_sets = load_cmb_for_query(map_filename + ".cmb.txt");
	auto setmap = load_setmap(map_filename + ".bin", nrows, ncolumns, false);
	std::cerr << "r = " << nrows << " | b = " << ncolumns << "\n";
	std::cerr << "total dimension = " << nrows * ncolumns << "\n";
	fprintf(stdout, "%lu %lu %lu", nrows, ncolumns, nrows*ncolumns);
	return EXIT_SUCCESS;
}

int mms_main(int argc, char* argv[])
{
	std::string kmc_filename, output_filename;
	std::vector<std::pair<uint32_t, std::size_t>> sorted_hist;
	uint64_t nrows = 0;
	uint64_t ncolumns = 0;
	double epsilon = 0.01;
	uint32_t ignored = std::numeric_limits<uint32_t>::max();
	
	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:s:r:b:e:g:h", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'o') {
			output_filename = opt.arg;
		} else if (c == 's') {
			sorted_hist = sort_histogram(load_histogram(opt.arg));
		} else if (c == 'r') {
			nrows = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'b') {
			ncolumns = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'e') {
			epsilon = std::stod(opt.arg);
		} else if (c == 'g') {
			ignored = static_cast<uint32_t>(std::stoul(opt.arg, nullptr, 10));
		} else if (c == 'h') {
			print_cms_help();
			return EXIT_SUCCESS;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_cms_help();
			return EXIT_FAILURE;
		}
	}

	if(kmc_filename == "" or output_filename == "") throw std::runtime_error("-i and -o are mandatory arguments");
	if(sorted_hist.size() == 0) sorted_hist = sort_histogram(compute_histogram(kmc_filename));
	std::size_t L1_norm = optimise_r_b(sorted_hist, epsilon, nrows, ncolumns);
	
	fprintf(stderr, "Starting filling the sketch of size %lu x %lu = %lu\n", nrows, ncolumns, nrows * ncolumns);
	std::vector<uint32_t> sketch(nrows * ncolumns, 0);
	std::unordered_map<uint32_t, uint32_t> invidx = create_inv_index(sorted_hist);
	fill_mms_sketch(kmc_filename, nrows, ncolumns, ignored, invidx, sketch);
	
	store_setmap(output_filename + ".mms", nrows, ncolumns, sketch);
	fprintf(stdout, "%lu %lu %u", L1_norm, nrows * ncolumns, *std::max_element(std::cbegin(sketch), std::cend(sketch)));//script-friendly output
	return EXIT_SUCCESS;
}

int mmschk_main(int argc, char* argv[])
{
	std::string kmc_filename, cms_filename;
	uint64_t nrows, ncolumns;
	uint32_t ignored = std::numeric_limits<uint32_t>::max();

	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "i:d:g:h", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'd') {
			cms_filename = opt.arg;
		} else if (c == 'g') {
			ignored = static_cast<uint32_t>(std::stoul(opt.arg, nullptr, 10));
		} else if (c == 'h') {
			print_cmschk_help();
			return EXIT_SUCCESS;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_cmschk_help();
			return EXIT_FAILURE;
		}
	}

	if(kmc_filename == "" or cms_filename == "") throw std::runtime_error("-i and -d are mandatory arguments");
	auto sorted_histogram = sort_histogram(load_histogram(cms_filename + ".shist.txt"));
	std::unordered_map<uint32_t, uint32_t> invidx = create_inv_index(sorted_histogram);
	auto setmap = load_setmap(cms_filename + ".mms", nrows, ncolumns, true);	
	auto rvals = check_mm_sketch(kmc_filename, nrows, ncolumns, setmap, invidx);
	fprintf(stdout, "%s %s %s %s %s", rvals[0].c_str(), rvals[1].c_str(), rvals[2].c_str(), rvals[3].c_str(), rvals[4].c_str());//script-friendly output
	return EXIT_SUCCESS;
}

int cms_main(int argc, char* argv[])
{
	std::string kmc_filename, output_filename;
	std::vector<std::pair<uint32_t, std::size_t>> sorted_hist;
	uint64_t nrows = 0;
	uint64_t ncolumns = 0;
	double epsilon = 0.01;
	uint32_t ignored = std::numeric_limits<uint32_t>::max();
	
	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:s:r:b:e:g:h", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'o') {
			output_filename = opt.arg;
		} else if (c == 's') {
			sorted_hist = sort_histogram(load_histogram(opt.arg));
		} else if (c == 'r') {
			nrows = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'b') {
			ncolumns = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'e') {
			epsilon = std::stod(opt.arg);
		} else if (c == 'g') {
			ignored = static_cast<uint32_t>(std::stoul(opt.arg, nullptr, 10));
		} else if (c == 'h') {
			print_cms_help();
			return EXIT_SUCCESS;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_cms_help();
			return EXIT_FAILURE;
		}
	}

	if(kmc_filename == "" or output_filename == "") throw std::runtime_error("-i and -o are mandatory arguments");
	if(sorted_hist.size() == 0) sorted_hist = sort_histogram(compute_histogram(kmc_filename));
	std::size_t L1_norm = optimise_r_b(sorted_hist, epsilon, nrows, ncolumns);
	
	fprintf(stderr, "Starting filling the sketch of size %lu x %lu = %lu\n", nrows, ncolumns, nrows * ncolumns);
	std::vector<uint32_t> sketch(nrows * ncolumns, 0);//just to be explicit about the 0
	fill_cms_sketch(kmc_filename, nrows, ncolumns, ignored, sketch);
	
	store_setmap(output_filename + ".cms", nrows, ncolumns, sketch);
	fprintf(stdout, "%lu %lu %u", L1_norm, nrows * ncolumns, *std::max_element(std::cbegin(sketch), std::cend(sketch)));//script-friendly output
	return EXIT_SUCCESS;
}

int cmschk_main(int argc, char* argv[])
{
	std::string kmc_filename, cms_filename;
	uint64_t nrows, ncolumns;
	uint32_t ignored = std::numeric_limits<uint32_t>::max();

	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "i:d:g:h", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'd') {
			cms_filename = opt.arg;
		} else if (c == 'g') {
			ignored = static_cast<uint32_t>(std::stoul(opt.arg, nullptr, 10));
		} else if (c == 'h') {
			print_cmschk_help();
			return EXIT_SUCCESS;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_cmschk_help();
			return EXIT_FAILURE;
		}
	}

	if(kmc_filename == "" or cms_filename == "") throw std::runtime_error("-i and -d are mandatory arguments");
	auto setmap = load_setmap(cms_filename + ".cms", nrows, ncolumns, true);
	auto rvals = check_cm_sketch(kmc_filename, nrows, ncolumns, ignored, setmap);
	fprintf(stdout, "%s %s %s %s %s", rvals[0].c_str(), rvals[1].c_str(), rvals[2].c_str(), rvals[3].c_str(), rvals[4].c_str());//script-friendly output
	return EXIT_SUCCESS;
}

int bbhash_main(int argc, char* argv[])
{
	typedef boomphf::SingleHashFunctor<uint64_t>  hasher_t;
	std::string kmc_filename, output_filename;
	
	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:s:r:b:f:e:h", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'o') {
			output_filename = opt.arg;
		} else if (c == 'h') {
			print_bbhash_help();
			return EXIT_SUCCESS;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_bbhash_help();
			return EXIT_FAILURE;
		}
	}

	if(kmc_filename == "" or output_filename == "") throw std::runtime_error("-i and -o are mandatory arguments");
	
	CKMCFile kmcdb;
	if (!kmcdb.OpenForListing(kmc_filename)) throw std::runtime_error("Unable to open the database\n");
	if (kmcdb.KmerLength() > 32) throw std::length_error("Maximum value of k = 32");
	uint8_t shift = static_cast<uint8_t>(2*kmcdb.KmerLength());
	uint64_t mask = (shift == 64 ? 0ULL : (1ULL<<shift)) - 1;
	//fprintf(stderr, "k-mer length = %u, mask = %lu\n", kmcdb.KmerLength(), mask);
	
	KMCRangeWrapper kmcrw(kmcdb, mask);
	boomphf::mphf<uint64_t, hasher_t> *bbhash = new boomphf::mphf<uint64_t, hasher_t>(kmcdb.KmerCount(), kmcrw, 1, 1.0);

	//saving bbhash
	std::ofstream mphfout(output_filename + ".bbh", std::ios::binary);
	bbhash->save(mphfout);
	mphfout.close();

	//creating external array
	std::vector<uint32_t> payload(kmcdb.KmerCount());
	
	CKmerAPI kmer(kmcdb.KmerLength());
	char str_kmer[kmcdb.KmerLength() + 1];
	uint32_t counter = 0;
	uint32_t max_counter = 0;
	kmcdb.RestartListing();
	while(kmcdb.ReadNextKmer(kmer, counter))
	{
		kmer.to_string(str_kmer);
		uint64_t packed = pack_kmer(str_kmer, kmcdb.KmerLength(), mask);
		//fprintf(stderr, "packed = %lu\n", packed);
		uint64_t index = bbhash->lookup(packed); 
		//fprintf(stderr, "index = %lu\n", index);
		payload[index] = counter;
		if(max_counter < counter) max_counter = counter;
	}
	store_setmap(output_filename + ".pld", payload.size(), 1, payload);
	fprintf(stdout, "%u %lu", max_counter, payload.size());//script-friendly output
	kmcdb.Close();
	delete bbhash;
	return EXIT_SUCCESS;
}


int main(int argc, char* argv[])
{
	ketopt_t om = KETOPT_INIT;
	int c;
	while((c = ketopt(&om, argc, argv, 0, "h", 0)) >= 0) //parse main options and find subcommand
	{
		if(c == 'h') 
		{
			print_subcommands();
			return EXIT_SUCCESS;
		}
	}
	if (om.ind == argc) {
		fprintf(stderr, "[Error] Subcommand unavailable\n\n");
		print_subcommands();
		return EXIT_FAILURE;
	}

	if (std::strcmp(argv[om.ind], "histogram") == 0) {
		return histogram_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "sense") == 0) {
		return sense_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "check") == 0) {
		return check_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "info") == 0) {
		return info_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "mms") == 0) {
		return mms_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "mmschk") == 0) {
		return mmschk_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "cms") == 0) {
		return cms_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "cmschk") == 0) {
		return cmschk_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "bbhash") == 0) {
		return bbhash_main(argc - om.ind, &argv[om.ind]);
	} else {
		fprintf(stderr, "Missing subcommand\n\n");
		print_subcommands();
	}
	return EXIT_SUCCESS;
}
