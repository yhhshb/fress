#include "fresslib.hpp"

extern "C" {
#include "ketopt.h"
}

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
	fprintf(stderr, "\t-i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "\t-o\toutput file. A two-column tsv file where the first column contains the frequencies and the second column the total number of k-mers having that frequency\n");
	fprintf(stderr, "\t-h\tshows this help\n");
}

void print_sense_help()
{
	fprintf(stderr, "sense options:\n");
	fprintf(stderr, "\t-i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "\t-o\toutput probabilistic map containing the frequencies\n");
	fprintf(stderr, "\t-s\tinput histogram generated from the input kmc database using the <histogram> subcommand\n");
	fprintf(stderr, "\t-r\tnumber of independent bucket rows. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-b\tnumber of columns. If not specified it is computed from the histogram\n");
	fprintf(stderr, "\t-p\tprobability of collision of two elements in the intersection [0.01]\n");
	fprintf(stderr, "\t-h\tshows this help\n");
}

void print_check_help()
{
	fprintf(stderr, "check options:\n");
	fprintf(stderr, "\t-i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "\t-d\tinput map built from the input kmc database (without extensions)\n");
	fprintf(stderr, "\t-h\tshows this help\n");
	fprintf(stderr, "\nThe output is on stdout and are all k-mers for which the predicted frequency is wrong or there are multiple frequencies\n");
	fprintf(stderr, "Each k-mer will have its true frequency and the (wrong) intersection of frequencies\n");
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
			return EXIT_SUCCESS;
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_histogram_help();
			return EXIT_FAILURE;
		}
	}
	auto histo = compute_histogram(kmc_filename);
	std::ofstream hf(output_filename);
	for(auto it = histo.cbegin(); it != histo.cend(); ++it)
	{
		hf << it->first << "\t" << it->second << "\n";
	}
	hf.close();
	return EXIT_SUCCESS;
}

int sense_main(int argc, char* argv[])
{
	static ko_longopt_t longopts[] = {
		{NULL, 0, 0}
	};

	ketopt_t opt = KETOPT_INIT;
	int c;
	std::string kmc_filename, output_filename;
	std::vector<std::pair<uint32_t, std::size_t>> sch;
	std::size_t nrows = 0;
	std::size_t ncolumns = 0;
	double f = 0.5;
	double p = 0.01;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:s:r:b:f:p:h", longopts)) >= 0)
	{
		if (c == 'i') {
			kmc_filename = opt.arg;
		} else if (c == 'o') {
			output_filename = opt.arg;
		} else if (c == 's') {
			sch = sort_histogram(load_histogram(opt.arg));
		} else if (c == 'r') {
			nrows = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'b') {
			ncolumns = std::stoul(opt.arg, nullptr, 10);
		} else if (c == 'f') {//Expected fraction of empty buckets for the second heaviest element
			f = std::stod(opt.arg);//0.5 is the optimal value that minimizes the space for all p
		} else if (c == 'p') {
			p = std::stod(opt.arg);
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
	//if(nrows == 0 or ncolumns == 0)
	{
		if(sch.size() == 0) sch = sort_histogram(compute_histogram(kmc_filename));
		if(ncolumns == 0) ncolumns = - static_cast<double>(sch[1].second) / std::log(f);
		double temp = std::log(p) / std::log(1-f);
		if(nrows == 0) nrows = std::max(static_cast<std::size_t>(1), static_cast<std::size_t>(std::ceil(temp)));
	}
	
	fprintf(stderr, "Starting filling the sketch of size %lu x %lu = %lu\n", nrows, ncolumns, nrows * ncolumns);

	std::vector<std::string> combinations;
	std::vector<uint32_t> sketch(nrows * ncolumns);
	fill_sketch_small(kmc_filename, nrows, ncolumns, sch.size() != 0 ? sch[0].first : std::numeric_limits<uint32_t>::max(), combinations, sketch);

	std::ofstream combo(output_filename + ".cmb.txt");
	for(auto s : combinations) combo << s << "\n";
	combo.close();

	std::ofstream skdump(output_filename + ".bin", std::ios::binary);
	skdump.write(reinterpret_cast<char*>(&nrows), sizeof(decltype(nrows)));
	skdump.write(reinterpret_cast<char*>(&ncolumns), sizeof(decltype(ncolumns)));
	skdump.write(reinterpret_cast<char*>(sketch.data()), sketch.size() * sizeof(decltype(sketch)::value_type));
	skdump.close();
	
	std::ofstream hf(output_filename + ".shist.txt");
	for(auto it = sch.cbegin(); it != sch.cend(); ++it)
	{
		hf << it->first << "\t" << it->second << "\n";
	}
	hf.close();


	return EXIT_SUCCESS;
}

int check_main(int argc, char* argv[])
{
	static ko_longopt_t longopts[] = {
		{NULL, 0, 0}
	};

	ketopt_t opt = KETOPT_INIT;
	int c;
	std::string kmc_filename, map_filename;
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

	std::size_t nrows, ncolumns;
	std::ifstream skdump(map_filename + ".bin", std::ios::binary);
	if(not skdump.is_open()) throw std::runtime_error("Unable to open sketch index file");
	skdump.read(reinterpret_cast<char*>(&nrows), sizeof(decltype(nrows)));
	skdump.read(reinterpret_cast<char*>(&ncolumns), sizeof(decltype(ncolumns)));
	std::vector<uint32_t> setmap(nrows * ncolumns);
	skdump.read(reinterpret_cast<char*>(setmap.data()), setmap.size() * sizeof(decltype(setmap)::value_type));
	skdump.close();

	skdump.open(map_filename + ".cmb.txt");
	if(not skdump.is_open()) throw std::runtime_error("Unable to open sketch combination file");
	std::vector<std::vector<uint32_t>> frequency_sets;
	std::string line;
	c = 0;
	while(std::getline(skdump, line))
	{
		if(c == 0 or line != "") frequency_sets.push_back(str2set<uint32_t>(line, ','));
		if(c == 0) c = 1;//The first line of the file could be empty if the heavies element is implicit.
	}
	skdump.close();

	std::unordered_map<uint32_t, uint32_t> invidx = create_inv_index(sort_histogram(load_histogram(map_filename + ".shist.txt")));
	check_sketch(kmc_filename, nrows, ncolumns, setmap, frequency_sets, invidx);
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
	} else {
		fprintf(stderr, "Missing subcommand\n\n");
		print_subcommands();
	}
	return EXIT_SUCCESS;
}
