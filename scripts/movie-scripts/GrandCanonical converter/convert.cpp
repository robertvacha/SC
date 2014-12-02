#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

/*
 *  Converts a "movie" to "movieFix"
 *
 *  Used for grandCanonical simulations
 *
 *  unifies numbers of particles in all snapshots by adding particles to 0.0
 */

int main()
{
    string line;
    ifstream in("movie");
    ofstream out("movieFix");
    int num;
    int largest = 0;
    int i=0;

    if (in.is_open() && out.is_open())
    {
        while ( getline (in,line) ) {
            istringstream ( line ) >> num;
            if(num > largest)
                largest = num;
        }

        in.clear();
        in.seekg (0, in.beg);
        while ( getline (in,line) ) {

            if(line.size() < 10) { // not pretty
                if(i==0 ) {
                    out << largest << "\n";
                    continue;
                }
                for(;i<=largest; i++) {
                    out << " 0.00000001e+00  0.00000001e+00  0.00000001e+00"
                        << "    0.00000001e+00  0.00000001e+00  0.00000001e+00"
                        << "    0.00000001e+00  0.00000001e+00  0.00000001e+00 0\n";
                }
                i=0;
                out << largest << "\n";
            } else {
                out << line << "\n";
                i++;
            }
        }
        in.close();
        out.close();
    } else cout << "Unable to open file";

    cout << "Particles: " << largest << endl;

    return 0;
}

