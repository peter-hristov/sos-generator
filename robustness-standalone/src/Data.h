#pragma once

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "./CGALTypedefs.h"

class Data
{
  public:
    Data() {}

    // Min/max range coordinates
    float minF, maxF;
    float minG, maxG;

    // Min/max domain coordinates
    float minX, maxX;
    float minY, maxY;
    float minZ, maxZ;


    std::string longnameF, longnameG, units;

    std::vector<std::vector<size_t>> tetrahedra;
    std::vector<std::vector<float>> vertexDomainCoordinates;

    std::vector<float> vertexCoordinatesF;
    std::vector<float> vertexCoordinatesG;


    void readData(const std::string&, const float&);
    void readDataVTU(const std::string&, const float&);

    void printMesh();
    void sortVertices();
    // Compute the min/max F, G and X, Y, Z coordinates
    void computeMinMaxRangeDomainCoordinates();


    // Reeb stuff relted content
    Arrangement_2 arr;
};
