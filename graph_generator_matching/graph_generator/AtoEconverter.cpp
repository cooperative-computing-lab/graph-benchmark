//
// Created by Page, Brian Andrew on 2019-02-04.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

int main(int argc, char** argv){
    std::ifstream infile(argv[1]);
    if (!infile){
        std::cout << "FAILED TO OPEN FILE!" << std::endl;
        exit(1);
    }
    std::string line;

    int i = 0, j = 0, uSize, vsize, edgecount;
    double tempData;
    bool previousLineCommented = false;
    while (std::getline(infile, line)) {
        if (line[0] != '%') {
            if (previousLineCommented == true) {
                //some stuff to get row,col and size data, etc.
                previousLineCommented = false;

                size_t pos = 0;
                std::string token;
                i = 0;
                while ((pos = line.find(',')) != std::string::npos) {
                    token = line.substr(0, pos);
                    line.erase(0, pos + 1);

                    if (i == 0) {
                        uSize = std::stoi(token);
                    } else {
                        vsize = std::stoi(token);
                    }

                    i++;
                }
                edgecount = std::stoi(line);
		
		cout << "%" << endl;
                cout << uSize << " " << vsize << " " << edgecount << endl;
            } else {
                std::vector<int> adjacencies;
                size_t pos = 0;
                std::string token;
                int uVert, temp;

		pos = line.find(',');
		if (pos != std::string::npos){
                    token = line.substr(0, pos);
                    line.erase(0, pos + 1);
                    uVert = ::atoi(token.c_str());
                    if (uVert < 0) uVert *= -1;
		}
                while ((pos = line.find(',')) != std::string::npos) {
                    token = line.substr(0, pos);
                    line.erase(0, pos + 1);
                    int temp = ::atoi(token.c_str());
                    if (temp < 0) temp *= -1;
                    //cout << j << " " << temp << " " << 1 << endl;
                    adjacencies.push_back(temp);
                }

                temp = ::atoi(line.c_str());
                if (temp < 0) temp *= -1;
                //cout << j << " " << temp << " " << 1 << endl;
                adjacencies.push_back(temp);
                sort(adjacencies.begin(), adjacencies.end());

                for (int k = 0; k < adjacencies.size(); ++k){
                    cout << uVert << " " << adjacencies[k] << " " << 1 << endl;
                }

                ++j;
            }
        } else {
            previousLineCommented = true;
        }
    }
}
