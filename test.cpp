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
	parser.add<string>("train", '\0', "Train data path.", true);
	parser.add<string>("test", '\0', "Test data path.", true);
	parser.parse_check(argc, argv);

	//LoadData
	//http://corpus-texmex.irisa.fr/
	size_t filesize, dimension;

	lsh_impl::Matrix<float> train = lsh_impl::Matrix<float>::load_from_file(parser.get<string>("train").c_str());
	for (size_t i = 0; i < train.col; i ++) {
		cout << train[0][i] << " ";
	}

}
