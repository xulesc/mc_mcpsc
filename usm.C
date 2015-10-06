#include "usm.hpp"
#include <sys/timeb.h>

using namespace std;
using namespace boost::iostreams;
namespace io = boost::iostreams;

//////////////////////////////////////////////////////////
template<typename T>
inline int compressIt(std::vector<T> s){
    std::stringstream uncompressed, compressed;
    for (typename std::vector<T>::iterator it = s.begin();
         it != s.end(); it++)
        uncompressed << *it;

    io::filtering_streambuf<io::input> o;
    o.push(io::gzip_compressor());
    o.push(uncompressed);
    io::copy(o, compressed);

    return compressed.str().length();
}

void USM::load_data(int file_count, char *pdb_dir, char **filenames) {
	cout << "loading data for USM" << endl;
	int i;
	for(i = 0; i < file_count; i++) {
		string filepath = std::string(pdb_dir) + std::string(filenames[i]);
		ifstream in(filepath.c_str(), ios_base::in | ios_base::binary);
		std::vector<std::string> cm_lines;
		std::copy(std::istream_iterator<std::string>(in), std::istream_iterator<std::string>(), std::back_inserter(cm_lines));
		data_map[filenames[i]] = cm_lines;
		complexity_map[filenames[i]] = compressIt(cm_lines);
	}
	cout << "done loading data for USM" << endl;
}

int iComp(string name, std::set<std::string> lines) {
	for(auto line : lines)
		if(name.find(line) == 0)
			return 1;
	return 0;
}

void USM::calculate_pairwise_distances(std::set<std::string> lines) {
        cout << "calculating pairwise distances" << endl;
        struct timeb t1, t2;
        for(auto entry1: data_map) {
		if(iComp(entry1.first,lines) == 0)
			continue;
                int x = complexity_map[entry1.first];
                std::vector< std::string > vector1 = entry1.second;
                for(auto entry2: data_map) {
			if(iComp(entry2.first,lines)==1)
				continue;
                        int y = complexity_map[entry2.first];
                        std::vector<std::string> vector2 = entry2.second;
                        ftime(&t1);
                        double d = calculate_pairwise_distances(x, y, vector1, vector2);
                        ftime(&t2);
                        cout << "Algo: 3, Protein1: " << entry1.first << ", Protein2: " << entry2.first 
                             << ", len1: " << vector1.size() << ", len2: " << vector2.size() << ", USM: " << d << ", TM1: -1, TM2: -1" 
                             << ", id: -1, sec_job_start: "<< t1.time << ", sec_result_send: " << t2.time << ", processing_time_msec: "
                             << ( (1000 * t2.time + t2.millitm) - (1000 * t1.time + t1.millitm) )
                             << ", data_collect_time_msec: -1, idle_time_msec: -1" << endl;
                }
        }
        cout << "done calculating pairwise distances" << endl;
}

void USM::calculate_pairwise_distances() {
	cout << "calculating pairwise distances" << endl;
	struct timeb t1, t2;
	for(auto entry1: data_map) {
		int x = complexity_map[entry1.first];
		std::vector< std::string > vector1 = entry1.second;
		for(auto entry2: data_map) {
			int y = complexity_map[entry2.first];
			std::vector<std::string> vector2 = entry2.second;
			ftime(&t1);
			double d = calculate_pairwise_distances(x, y, vector1, vector2);
			ftime(&t2);
			cout << "Algo: 3, Protein1: " << entry1.first << ", Protein2: " << entry2.first 
			     << ", len1: " << vector1.size() << ", len2: " << vector2.size() << ", USM: " << d << ", TM1: -1, TM2: -1"
			     << ", id: -1, sec_job_start: "<< t1.time << ", sec_result_send: " << t2.time << ", processing_time_msec: " 
			     << ( (1000 * t2.time + t2.millitm) - (1000 * t1.time + t1.millitm) )
			     << ", data_collect_time_msec: -1, idle_time_msec: -1" << endl;
		}
	}
	cout << "done calculating pairwise distances" << endl;
}

double USM::calculate_pairwise_distances(int &x, int &y, std::vector< std::string > &vector1, std::vector< std::string > &vector2) {
	std::vector< std::string > v1;
	v1.insert(v1.end(), vector1.begin(), vector1.end());
	v1.insert(v1.end(), vector2.begin(), vector2.end());	
	int xy = compressIt(v1);

	v1.clear();
        v1.insert(v1.end(), vector2.begin(), vector2.end());
	v1.insert(v1.end(), vector1.begin(), vector1.end());
	int yx = compressIt(vector2);

	return compute_distance(x, y, xy, yx);
}

double USM::compute_distance(int x, int y, int xy, int yx) {
	double maxNum = ((yx - y) > (xy - x)) ? yx - y : xy - x;
	double maxDen = (x > y) ? x : y;

        return (maxNum/maxDen);
}
//////////////////////////////////////////////////////////
