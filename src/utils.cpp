//
// Created by Asunción Gómez on 18-01-24.
//
#include<bits/stdc++.h>

/**
 * Given a position in the bitvector and the dimension of the quadtree, return the coordinates (x_1,x_2,..., x_d) in the quadtree.
 * @param leafPos the index of the leaf position. For example, the leaf in 0010 will have the 2 index.
 * @param k_d the arity of the quadtree (number of children per node).
 * @param height the height of the quadtree
 * @param gridSide the grid side of the quadtree
 * @param gridDimension of the grid
 * @return the point coordinates of the leaf in the quadtree.
 * @example getPointCoordinates(30, 8, 4, 3, 2) = [6,3]
 */
uint64_t* getPointCoordinates(uint64_t leafPos, uint16_t k_d, uint16_t height, uint64_t gridSide, uint64_t gridDimension){
    // TODO:
}

/**
 * See paper
 * @param bv
 * @param k_d
 * @return an array of k_d dimension with the point coordinates in the grid according to the bit path.
 */
uint64_t* getPointCoordFromBitvectorPath(std::vector<bool> bv, uint16_t k_d){
    uint64_t* pointCoord = new uint64_t[k_d];
    for(uint16_t i = 0 ; i < k_d; i++){
        // TODO: create bitmask with 1 each i+k_d position
        uint64_t bitmask = 0;
        /*if(j == i) {
            bitmask = bitmask | (1 << j);
        }*/
        pointCoord[i] = 0;
    }
    // TODO: retrieve each i, i+k_d, i+2k_d, ... bit from the bitvector and convert it to an integer
    // TODO: pass bits to an integer
    return pointCoord;

}


// function to convert decimal to binary
void decToBinary(int n)
{
    // array to store binary number
    int binaryNum[32];

    // counter for binary array
    int i = 0;
    while (n > 0) {

        // storing remainder in binary array
        binaryNum[i] = n % 2;
        n = n / 2;
        i++;
    }

    // printing binary array in reverse order
    for (int j = i - 1; j >= 0; j--)
        cout << binaryNum[j];
}