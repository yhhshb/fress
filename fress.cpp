#include <cstring>
#include "fresslib.hpp"

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
	fprintf(stderr, "mphf\t build a BBHash MPHF + external frequency array\n");//No need for checking because error = 0
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

void print_cms_help()
{
	fprintf(stderr, "cms options:\n");
	fprintf(stderr, "\t-i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "\t-o\toutput count-min sketch binary file\n");
	fprintf(stderr, "\t-e\tepsilon approximation of the L1 sum of errors [0.01]\n");
	fprintf(stderr, "\t-s\tinput histogram generated from the input kmc database using the <histogram> subcommand. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-r\tnumber of independent bucket rows. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-b\tnumber of columns. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-h\tshows this help\n");

}

void print_cmschk_help()
{
	fprintf(stderr, "check options:\n");
	fprintf(stderr, "\t-i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "\t-d\tinput count-min sketch built from the input kmc database\n");
	fprintf(stderr, "\t-h\tshows this help\n");
	fprintf(stderr, "\nHuman-readable output on stderr, script-friendly output on stdout\n");
}

int histogram_main(int argc, char* argv[])
{
	std::string kmc_filename, output_filename;

	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:h:", longopts)) >= 0)
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
	//double f = 0.5;
	double epsilon = 0.01;
	
	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:s:r:b:f:e:h", longopts)) >= 0)
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
		//} else if (c == 'f') {//Expected fraction of empty buckets for the second heaviest element
		//	f = std::stod(opt.arg);//0.5 is the optimal value that minimizes the space for all p
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
	//if(ncolumns == 0) ncolumns = - static_cast<double>(sorted_hist[1].second) / std::log(f);
	//double temp = std::log(p) / std::log(1-f);
	//if(nrows == 0) nrows = std::max(static_cast<std::size_t>(1), static_cast<std::size_t>(std::ceil(temp)));
	std::size_t L1_norm = optimise_r_b(sorted_hist, epsilon, nrows, ncolumns);
	
	fprintf(stderr, "Starting filling the sketch of size %lu x %lu = %lu\n", nrows, ncolumns, nrows * ncolumns);

	std::vector<std::string> combinations;
	std::vector<uint32_t> sketch(nrows * ncolumns);
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

	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "i:d:g:h", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'd') {
			map_filename = opt.arg;
		} else if (c == 'h') {
			print_check_help();
			return EXIT_SUCCESS;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_check_help();
			return EXIT_FAILURE;
		}
	}

	if(kmc_filename == "" or map_filename == "") throw std::runtime_error("-i and -d are mandatory arguments");
	std::unordered_map<uint32_t, uint32_t> invidx = create_inv_index(sort_histogram(load_histogram(map_filename + ".shist.txt")));
	auto frequency_sets = load_cmb_for_query(map_filename + ".cmb.txt");
	auto setmap = load_setmap(map_filename + ".bin", nrows, ncolumns, true);
	auto rvals = check_sketch(kmc_filename, nrows, ncolumns, setmap, frequency_sets, invidx);
	fprintf(stdout, "%s %s %s %s", rvals[0].c_str(), rvals[1].c_str(), rvals[2].c_str(), rvals[3].c_str());//script-friendly output
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
	return EXIT_SUCCESS;
}

int cms_main(int argc, char* argv[])
{
	std::string kmc_filename, output_filename;
	std::vector<std::pair<uint32_t, std::size_t>> sorted_hist;
	uint64_t nrows = 0;
	uint64_t ncolumns = 0;
	//double f = 0.5;
	double epsilon = 0.01;
	
	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:s:r:b:f:e:h", longopts)) >= 0)
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
		//} else if (c == 'f') {//Expected fraction of empty buckets for the second heaviest element
		//	f = std::stod(opt.arg);//0.5 is the optimal value that minimizes the space for all p
		} else if (c == 'e') {
			epsilon = std::stod(opt.arg);
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
	std::vector<uint32_t> sketch(nrows * ncolumns);
	fill_cms_sketch(kmc_filename, nrows, ncolumns, sketch);
	
	store_setmap(output_filename + ".bin", nrows, ncolumns, sketch);
	fprintf(stdout, "%lu %lu", L1_norm, nrows * ncolumns);//script-friendly output
	return EXIT_SUCCESS;
}

int cmschk_main(int argc, char* argv[])
{
	std::string kmc_filename, cms_filename;
	uint64_t nrows, ncolumns;

	static ko_longopt_t longopts[] = {{NULL, 0, 0}};
	ketopt_t opt = KETOPT_INIT;
	int c;
	while((c = ketopt(&opt, argc, argv, 1, "i:d:g:h", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'd') {
			cms_filename = opt.arg;
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
	auto setmap = load_setmap(cms_filename + ".bin", nrows, ncolumns, true);
	auto rvals = check_cm_sketch(kmc_filename, nrows, ncolumns, setmap);
	fprintf(stdout, "%s %s %s %s", rvals[0].c_str(), rvals[1].c_str(), rvals[2].c_str(), rvals[3].c_str());//script-friendly output
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
	} else if (std::strcmp(argv[om.ind], "cms") == 0) {
		return cms_main(argc - om.ind, &argv[om.ind]);
	} else if (std::strcmp(argv[om.ind], "cmschk") == 0) {
		return cmschk_main(argc - om.ind, &argv[om.ind]);
	} else {
		fprintf(stderr, "Missing subcommand\n\n");
		print_subcommands();
	}
	return EXIT_SUCCESS;
}
