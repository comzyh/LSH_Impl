#include <cstdio>
#include <string>
#include <fstream>
#include "vendor/cmdline/cmdline.h"  // argparser
#include "lsh_impl/lsh_impl.hpp"

using namespace std;

int main(int argc, char **argv)
{
	// CommandLine parser
	cmdline::parser parser;
	parser.add<string>("base", 'b', "Base data path.", true);
	parser.add<string>("query", 'q', "Query data path.", true);
	parser.parse_check(argc, argv);

	//LoadData
	//http://corpus-texmex.irisa.fr/
	size_t filesize, dimension;

	lsh_impl::Matrix<float> base = lsh_impl::Matrix<float>::load_from_file(parser.get<string>("base").c_str());
	lsh_impl::Matrix<float> query = lsh_impl::Matrix<float>::load_from_file(parser.get<string>("query").c_str());

	
	return 0;
}
