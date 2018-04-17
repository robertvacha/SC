/** @file mesh.h*/

#ifndef MESH_H
#define MESH_H

#include "../structures/Conf.h"


/**
 * @brief Mesh for hole order parameter
 */
class Mesh
{
public:

    int dim[2]; ///< \@brief Mesh dimensions
    int* data;  ///< \@brief Mesh data
    int* tmp;   ///< \@brief tempporary list for hole search
    double alpha_init;

private:
    bool verbose;
    int _hist[200] = {0};

public:

    Mesh(bool verbose) : verbose(verbose) {}
    ~Mesh() {
        if(verbose) {
            ofstream myfile;
            string file= "pore_distrib";
            string num = std::to_string(mcout.rank);
            file.insert(0, num);

            myfile.open (file, std::fstream::app);

            if( myfile.is_open() ) { // histogram of holes

                for(int i=0; i<200; i++){
                    myfile << _hist[i] << " ";
                }
                myfile << endl;
                myfile.close();
            }
        }
    }

    /*..............................................................................*/
    /*........................HOLE IN MESH-MEMBRANE ORDER PARAM.....................*/

    /**
     * @brief meshInit
     * @param meshsize
     * @param npart
     * @param wlmtype
     * @param box
     * @param particleStore
     * @return
     */
    int meshInit(double meshsize, long npart, int wlmtype, Vector box, std::vector<Particle >* pvec);

    /**
     * @brief meshFill filling the mesh
     * @param npart
     * @param wlmtype
     * @param particleStore
     */
    void meshFill(long npart, int wlmtype, std::vector<Particle >* pvec);

    /**
     * @brief mesh_addpart add particle on coordinates posx posy to mesh return 0 if it was placed on empty spot
     * @param posx
     * @param posy
     * @param mesh
     * @param dim
     * @return
     */
    int addPart(double posx, double posy);

    /**
     * @brief mesh_removepart remove particle on coordinates posx posy from mesh and return 0 if there is a empty spot now
     * @param posx
     * @param posy
     * @param mesh
     * @param dim
     * @return
     */
    int removePart(double posx, double posy);

    /**
     * @brief mesh_square
     * @param x
     * @param y
     * @param dim
     */
    void meshSquare(int x, int y, int square[9]);

    /**
     * @brief mesh_neighbors
     * @param pos
     * @param dim
     * @param neighbors
     */
    void meshNeighbors(int pos, int neighbors[4]);

    /**
     * @brief mesh_findholes returns the number of holes and a list of mesh points belonging to each of them
     * @return
     */
    int findHoles();

    /**
     * @brief mesh_findholesdistrib prints out distribution of holes and returns the number of holes and a list of mesh points belonging to each of them
     * @return
     */
    int findHolesDistrib();

    /**
     * @brief mesh_print
     */
    void print();

    /**
     * @brief mesh_cpy
     * @param source
     */
    void operator=(Mesh& source);

    /**
     * @brief mesh_end
     * @return
     */
    int end();

};

#endif // MESH_H
