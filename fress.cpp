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
	fprintf(stderr, "i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "o\toutput file. A two-column tsv file where the first column contains the frequencies and the second column the total number of k-mers having that frequency\n");
	fprintf(stderr, "h\tshows this help\n");
}

void print_sense_help()
{
	fprintf(stderr, "sense options:\n");
	fprintf(stderr, "i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "o\toutput probabilistic map containing the frequencies\n");
	fprintf(stderr, "s\tinput histogram generated from the input kmc database using the <histogram> subcommand\n");
	fprintf(stderr, "r\tnumber of independent bucket rows. If not specified it is computed from the histogram\n");
	fprintf(stderr, "b\tnumber of columns. If not specified it is computed from the histogram\n");
	fprintf(stderr, "d\tfailing probability [0.01]. Ignored if the number of rows has been specified\n");
	fprintf(stderr, "h\tshows this help\n");
}

void print_check_help()
{
	fprintf(stderr, "check options:\n");
	fprintf(stderr, "i\tinput kmc database (without extensions)\n");
	fprintf(stderr, "d\tinput map built from the input kmc database (without extensions)\n");
	fprintf(stderr, "h\tshows this help\n");
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
	double delta = 0.01;
	while((c = ketopt(&opt, argc, argv, 1, "i:o:s:r:b:d:h", longopts)) >= 0)
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
		} else if (c == 'd') {
			delta = std::stod(opt.arg);
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
	if(nrows == 0 or ncolumns == 0)
	{
		if(sch.size() == 0) sch = sort_histogram(compute_histogram(kmc_filename));
		if(nrows == 0) nrows = log(sch.size() / delta); //FIXME find good dimensioning of the rows
		if(ncolumns == 0) ncolumns = static_cast<std::size_t>(sch[1].second); //FIXME is this the best number of columns?
	}
	//FIXME directly use an enumeration of the sets. Keep a hash table to store the sets and insert a new set when it is not found. 
	sketch_t sketch(nrows * ncolumns);
	uint32_t heavy_element_to_exclude = 0;
	if (sch.size() != 0) heavy_element_to_exclude = sch[0].first; //If we know the histogram then skip the heavy hitter to reduce space in the buckets
	fill_sketch(kmc_filename, nrows, ncolumns, sketch, heavy_element_to_exclude);

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
	while((c = ketopt(&opt, argc, argv, 1, "i:d:h", longopts)) >= 0)
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
	skdump.read(reinterpret_cast<char*>(&nrows), sizeof(decltype(nrows)));
	skdump.read(reinterpret_cast<char*>(&ncolumns), sizeof(decltype(ncolumns)));
	std::vector<uint32_t> setmap(nrows * ncolumns);
	skdump.read(reinterpret_cast<char*>(setmap.data()), setmap.size() * sizeof(decltype(setmap)::value_type));
	skdump.close();

	skdump.open(map_filename + ".cmb.txt");
	std::vector<std::vector<uint32_t>> frequency_sets;
	std::string line;

	c = 0;
	while(std::getline(skdump, line))
	{
		if(c == 0 or line != "") frequency_sets.push_back(str2set<uint32_t>(line, ','));
		if(c == 0) c = 1;//The first line of the file could be empty if the heavies element is implicit.
	}
	skdump.close();
	check_sketch(kmc_filename, nrows, ncolumns, setmap, frequency_sets);
	return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
	ketopt_t om = KETOPT_INIT;
	int c;
	while((c = ketopt(&om, argc, argv, 0, "x", 0)) >= 0) {} //parse main options and find subcommand
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
