#include <dirent.h>

#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/unordered_map.hpp>

#include <boost/iostreams/filter/counter.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <sstream>

using namespace std;
using namespace boost::iostreams;
namespace io = boost::iostreams;

class USM {
public:
  USM() {}

public:
  void load_data(int file_count, char *pdb_dir, char **filenames);
  void calculate_pairwise_distances();
  void calculate_pairwise_distances(std::set<std::string> lines);
  double calculate_pairwise_distances(int &x, int &y, std::vector< std::string > &vector1, std::vector< std::string > &vector2);  
  double compute_distance(int x, int y, int xy, int yx);
  
private:
  boost::unordered_map< std::string, std::vector<std::string> > data_map;
  boost::unordered_map< std::string, int > complexity_map;
};
//
