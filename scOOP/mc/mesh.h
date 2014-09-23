/** @file mesh.h*/

#ifndef MESH_H
#define MESH_H

#include <string.h>
#include "../structures/structures.h"
#include "../structures/Conf.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Mesh for hole order parameter
 */
class Mesh
{
public:

    int dim[2]; ///< \@brief Mesh dimensions
    int* data;  ///< \@brief Mesh data
    int* tmp;   ///< \@brief tempporary list for hole search

    Mesh() {}

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
    int meshInit(double meshsize, long npart, int wlmtype, Vector box, std::vector<Particle >* particleStore);

    /**
     * @brief meshFill filling the mesh
     * @param npart
     * @param wlmtype
     * @param particleStore
     */
    void meshFill(long npart, int wlmtype, std::vector<Particle >* particleStore);

    /**
     * @brief mesh_addpart add particle on coordinates posx posy to mesh return 0 if it was placed on empty spot
     * @param posx
     * @param posy
     * @param mesh
     * @param dim
     * @return
     */
    static int addPart(double posx, double posy, int **mesh, int dim[2]);

    /**
     * @brief mesh_removepart remove particle on coordinates posx posy from mesh and return 0 if there is a empty spot now
     * @param posx
     * @param posy
     * @param mesh
     * @param dim
     * @return
     */
    static int removePart(double posx, double posy, int **mesh, int dim[2]);

    /**
     * @brief mesh_square
     * @param x
     * @param y
     * @param dim
     */
    static void meshSquare(int x, int y, int dim[2], int (*square)[9]);

    /**
     * @brief mesh_neighbors
     * @param pos
     * @param dim
     * @param neighbors
     */
    static void meshNeighbors(int pos, int dim[2], int neighbors[4]);

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
