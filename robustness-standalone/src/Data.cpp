#include "Data.h"
#include "./Timer.h"

#include <filesystem>
#include <cassert>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <queue>

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <random>
#include <ranges>

#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>

#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

using namespace std;

double randomPerturbation(double epsilon) {
    static std::mt19937 gen(std::random_device{}()); // Seeded generator
    std::uniform_real_distribution<double> dist(-epsilon, epsilon);
    return dist(gen);
}




void
Data::computeMinMaxRangeDomainCoordinates()
{
    // Compute the min/max domain coordinates
    this->minX = this->vertexDomainCoordinates[0][0];
    this->maxX = this->vertexDomainCoordinates[0][0];

    this->minY = this->vertexDomainCoordinates[0][1];
    this->maxY = this->vertexDomainCoordinates[0][1];

    this->minZ = this->vertexDomainCoordinates[0][2];
    this->maxZ = this->vertexDomainCoordinates[0][2];

    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        this->minX = std::min(this->minX, this->vertexDomainCoordinates[i][0]);
        this->maxX = std::max(this->maxX, this->vertexDomainCoordinates[i][0]);

        this->minY = std::min(this->minY, this->vertexDomainCoordinates[i][1]);
        this->maxY = std::max(this->maxY, this->vertexDomainCoordinates[i][1]);

        this->minZ = std::min(this->minZ, this->vertexDomainCoordinates[i][2]);
        this->maxZ = std::max(this->maxZ, this->vertexDomainCoordinates[i][2]);
    }


    // Compute the min/max range coordinates
    this->minF = this->vertexCoordinatesF[0];
    this->maxF = this->vertexCoordinatesF[0];

    this->minG = this->vertexCoordinatesG[0];
    this->maxG = this->vertexCoordinatesG[0];

    for (int i = 0 ; i < this->vertexCoordinatesF.size() ; i++)
    {
        this->minF = std::min(this->minF, this->vertexCoordinatesF[i]);
        this->maxF = std::max(this->maxF, this->vertexCoordinatesF[i]);

        this->minG = std::min(this->minG, this->vertexCoordinatesG[i]);
        this->maxG = std::max(this->maxG, this->vertexCoordinatesG[i]);
    }

    // Add some padding to the range coordinates for better visibility
    //this->minF -= .2;
    //this->maxF += .2;
    //this->minG -= .2;
    //this->maxG += .2;

    this->minF -= 0.1 * (this->maxF - this->minF);
    this->maxF += 0.1 * (this->maxF - this->minF);
    this->minG -= 0.1 * (this->maxG - this->minG);
    this->maxG += 0.1 * (this->maxG - this->minG);
}


void Data::printMesh()
{
    // Print vertex domain coordinates
    cout << "Vertex domain coordinates: " << endl;
    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        printf("%d - (%f, %f, %f)\n", i, this->vertexDomainCoordinates[i][0], this->vertexDomainCoordinates[i][1], this->vertexDomainCoordinates[i][2]);
    }
    // Print vertex range coordinates
    cout << "Vertex range coordinates: " << endl;
    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        printf("%d - (%f, %f)\n", i, this->vertexCoordinatesF[i], this->vertexCoordinatesG[i]);
    }

    // Print tetrahedra
    cout << "Tetrahedra: " << endl;
    for (int i = 0 ; i < this->tetrahedra.size() ; i++)
    {
        printf("%d - (%ld, %ld, %ld, %ld)\n", i, this->tetrahedra[i][0], this->tetrahedra[i][1], this->tetrahedra[i][2], this->tetrahedra[i][3]);
    }
}
void Data::readDataVTU(const string &filename, const float &perturbationEpsilon)
{
    // Read the VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();


    // Set deault names for the range axis
    this->longnameF = "f";
    this->longnameG = "g";

    int numVertices = mesh->GetPoints()->GetNumberOfPoints(); 
    int numTets = mesh->GetNumberOfCells();

    // Initialize all the data arrays
    this->vertexCoordinatesF = std::vector<float>(numVertices, 0);
    this->vertexCoordinatesG = std::vector<float>(numVertices, 0);
    this->tetrahedra = std::vector<std::vector<size_t>>(numTets, {0, 0, 0, 0});
    this->vertexDomainCoordinates = std::vector<std::vector<float>>(numVertices, {0, 0, 0});



    // Print vertices
    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
    //std::cout << "Vertices:\n";
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
        double p[3];
        points->GetPoint(i, p);
        //std::cout << "Vertex " << i << ": (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";

        this->vertexDomainCoordinates[i][0] = p[0];
        this->vertexDomainCoordinates[i][1] = p[1];
        this->vertexDomainCoordinates[i][2] = p[2];
    }

    // Print tetrahedra
    //std::cout << "\nTetrahedra:\n";
    for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); i++) {
        vtkCell* cell = mesh->GetCell(i);
        if (cell->GetNumberOfPoints() == 4) { // Tetrahedron check
            //std::cout << "Tetrahedron " << i << ": ";
            for (vtkIdType j = 0; j < 4; j++) {
                //std::cout << cell->GetPointId(j) << " ";
                this->tetrahedra[i][j] = cell->GetPointId(j);
            }
            //std::cout << "\n";
        }
    }

    // Print vertex data arrays
    //std::cout << "\nVertex Data Arrays:\n";
    vtkPointData* pointData = mesh->GetPointData();

    assert(pointData->GetNumberOfArrays() >= 2);

    vtkDataArray* fDataArray = pointData->GetArray(1);
    vtkDataArray* gDataArray = pointData->GetArray(0);

    assert(fDataArray->GetNumberOfTuples() == numVertices);
    assert(gDataArray->GetNumberOfTuples() == numVertices);

    for (vtkIdType i = 0; i < fDataArray->GetNumberOfTuples(); i++) 
    {
        this->vertexCoordinatesF[i] = fDataArray->GetTuple1(i);
    }

    for (vtkIdType i = 0; i < gDataArray->GetNumberOfTuples(); i++) 
    {
        this->vertexCoordinatesG[i] = gDataArray->GetTuple1(i);
    }

    // Add numerical perturbation, taken from Tierny 2027 Jacobi Fiber Surfaces
    for (vtkIdType i = 0; i < gDataArray->GetNumberOfTuples(); i++) 
    {
        //const float e = randomPerturbation(1e-8);
        //const float iFloat = static_cast<float>(i);

        // Tierny 2017 - Jacobi Fiber Surfaces Version - not that useful
        //this->vertexCoordinatesF[i] += iFloat * e;
        //this->vertexCoordinatesG[i] += iFloat * e * iFloat * e;

        //this->vertexCoordinatesF[i] += randomPerturbation(perturbationEpsilon);
        //this->vertexCoordinatesG[i] += randomPerturbation(perturbationEpsilon);
    }

    // Now we can sort
    this->sortVertices();

    // Compute the ranges for a bounding box in the range.
    this->computeMinMaxRangeDomainCoordinates();
}

void
Data::readData(const string &filename, const float &perturbationEpsilon)
{
    // Set deault names for the range axis
    this->longnameF = "f";
    this->longnameG = "g";

    // Open data file
    std::ifstream dataFile (filename);
    if (false == dataFile.is_open()) { throw "Could not open data file."; }

    // Read in data in a string and skip the comments
    string rawStringData;
    string myline;
    while (dataFile) {
        std::getline (dataFile, myline);
        if (myline[0] == '#')
        {
            //std::cout << myline << '\n';
        }
        else
        {
            rawStringData += " " + myline;
        }
    }

    // Set up the inputstream from the string
    std::istringstream dataStream(rawStringData);

    // Read in the number of vertices and tets
    int numVertices, numTets;
    dataStream >> numVertices >> numTets;

    // Initialize all the data arrays
    this->vertexCoordinatesF = std::vector<float>(numVertices, 0);
    this->vertexCoordinatesG = std::vector<float>(numVertices, 0);
    this->tetrahedra = std::vector<std::vector<size_t>>(numTets, {0, 0, 0, 0});
    this->vertexDomainCoordinates = std::vector<std::vector<float>>(numVertices, {0, 0, 0});

    // Read in the domain coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> this->vertexDomainCoordinates[i][0];
        dataStream >> this->vertexDomainCoordinates[i][1];
        dataStream >> this->vertexDomainCoordinates[i][2];
    }

    // Read in the range coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> this->vertexCoordinatesF[i];
        dataStream >> this->vertexCoordinatesG[i];
    }
    
    // Read in the tetrahedron configuration
    for  (int i = 0 ; i < numTets ; i++)
    {
        dataStream >> this->tetrahedra[i][0];
        dataStream >> this->tetrahedra[i][1];
        dataStream >> this->tetrahedra[i][2];
        dataStream >> this->tetrahedra[i][3];
    }

    this->computeMinMaxRangeDomainCoordinates();

    this->sortVertices();
}

void Data::sortVertices()
{
    vector<pair<CartesianPoint, int>> pointWithIndices(this->vertexCoordinatesG.size());

    for (int i = 0 ; i < pointWithIndices.size() ; i++)
    {
        pointWithIndices[i] = {CartesianPoint(this->vertexCoordinatesF[i], this->vertexCoordinatesG[i]), i};
    }

    //cout << "Before sorting" << endl;
    //this->printMesh();

    std::sort(pointWithIndices.begin(), pointWithIndices.end(),
            [](const pair<CartesianPoint, int> &a, const pair<CartesianPoint, int>& b) {
            return CGAL::compare_xy(a.first, b.first) == CGAL::SMALLER;
            });

    //cout << "Printing points... " << endl;
    //for (int i = 0 ; i < pointWithIndices.size() ; i++)
    //{
        //printf("%d - (%f, %f) - %d\n", i, pointWithIndices[i].first.x(), pointWithIndices[i].first.y(), pointWithIndices[i].second);
    //}


    // Set up the inverse index search
    vector<int> meshIDtoSortIndex(pointWithIndices.size());
    for (int i = 0 ; i < pointWithIndices.size() ; i++)
    {
        meshIDtoSortIndex[pointWithIndices[i].second] = i;
    }

    //cout << "Swapping..." << endl;
    //for (int i = 0 ; i < pointWithIndices.size() ; i++)
    //{
        //printf("%d -> %d\n", i, meshIDtoSortIndex[i]);

    //}

    //
    // Now we can swap things around
    //

    // Set up copies of the originals for the swap, otherwise editin in place causes errors
    std::vector<std::vector<float>> vertexDomainCoordinatesOriginal = this->vertexDomainCoordinates;
    std::vector<float> vertexCoordinatesFOriginal = this->vertexCoordinatesF;
    std::vector<float> vertexCoordinatesGOriginal = this->vertexCoordinatesG;

    // Swap tet indices
    for (int i = 0 ; i < this->tetrahedra.size() ; i++)
    {
        for (int j = 0 ; j < this->tetrahedra[i].size() ; j++)
        {
            this->tetrahedra[i][j] = meshIDtoSortIndex[this->tetrahedra[i][j]];
        }
    }

    // Swap domain positions
    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        this->vertexDomainCoordinates[i] = vertexDomainCoordinatesOriginal[pointWithIndices[i].second];
    }
    
    // Swap range positions
    for (int i = 0 ; i < this->vertexCoordinatesF.size() ; i++)
    {
        this->vertexCoordinatesF[i] = vertexCoordinatesFOriginal[pointWithIndices[i].second];
        this->vertexCoordinatesG[i] = vertexCoordinatesGOriginal[pointWithIndices[i].second];
    }
    //cout << "After sorting..." << endl;
    //this->printMesh();

}
