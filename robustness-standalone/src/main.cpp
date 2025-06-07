#include <cstddef>
#include <set>
#include <filesystem>
#include <unordered_map>
#include <omp.h>

#include "./CGALTypedefs.h"
#include "./Timer.h"
#include "./Data.h"
#include "./utility/CLI11.hpp"



struct LineLess {
    bool operator()(const Line_2& l1, const Line_2& l2) const {
        auto a1 = l1.a();
        auto b1 = l1.b();
        auto c1 = l1.c();

        auto a2 = l2.a();
        auto b2 = l2.b();
        auto c2 = l2.c();

        if (a1 < a2) return true;
        if (a2 < a1) return false;
        if (b1 < b2) return true;
        if (b2 < b1) return false;
        return c1 < c2;
    }
};


using namespace std;
namespace fs = std::filesystem;

size_t computeDegenerateIntersectionsBruteForce(const vector<Segment_2> &segments)
{
    auto segments_intersect = [](const Segment_2& s1, const Segment_2& s2) 
    {
        // If the bounding boxes don't intersect, skip
        if (false == CGAL::do_overlap(s1.bbox(), s2.bbox()))
        {
            return false;
        }

        // If they intersect at their endpoints, skip
        if (s1.source() == s2.source() || s1.source() == s2.target() || s1.target() == s2.source() || s1.target() == s2.target())
        {
            return false;
        }

        return true;
    };

    Timer::start();
    int concurrentSegments = 0;

    for (std::size_t i = 0; i < segments.size(); ++i) {
        for (std::size_t j = i + 1; j < segments.size(); ++j) {

            if (!segments_intersect(segments[i], segments[j]))
                continue;  // skip if bounding boxes don't intersect

            auto i1 = CGAL::intersection(segments[i], segments[j]);

            // If there is no intersection just continue
            if (!i1) { continue; }

            for (std::size_t k = j + 1; k < segments.size(); ++k) {

                if (!segments_intersect(segments[j], segments[k]) ||
                        !segments_intersect(segments[i], segments[k]))
                    continue;

                auto i2 = CGAL::intersection(segments[j], segments[k]);

                // If there is no intersection just continue
                if (!i2) { continue; }

                auto i3 = CGAL::intersection(segments[i], segments[k]);

                // If there is no intersection just continue
                if (!i3) { continue; }

                auto p1 = std::get_if<Point_2>(&*i1);
                auto p2 = std::get_if<Point_2>(&*i2);
                auto p3 = std::get_if<Point_2>(&*i3);

                if (p1 && p2 && p3 && *p1 == *p2 && *p2 == *p3) {
                    concurrentSegments++;
                }
            }
        }
    }

    return concurrentSegments;
}


pair<size_t, size_t> computeCollienarPoints(const vector<Point_2> &points)
{
    int collinearPoints = 0;
    std::set<Line_2, LineLess> unique_lines; 
    for (int i = 0; i < points.size(); ++i) 
    {
        for (int j = i + 1; j < points.size(); ++j) 
        {
            Line_2 line(points[i], points[j]);

            if (unique_lines.contains(line))
            {
                collinearPoints++;
            }
            else
            {
                unique_lines.insert(line);
            }
        }
    }

    return {collinearPoints, unique_lines.size()};
}

pair<size_t, size_t> computeDegenerateIntersections(const vector<Segment_2> &segments)
{
    auto can_segments_intersect = [](const Segment_2& s1, const Segment_2& s2) 
    {
        if (!CGAL::do_overlap(s1.bbox(), s2.bbox()))
            return false;

        if (s1.source() == s2.source() || s1.source() == s2.target() ||
            s1.target() == s2.source() || s1.target() == s2.target())
            return false;

        return true; 
    };

    const std::size_t n = segments.size();

    // Thread-local sets
    std::vector<std::set<Point_2>> local_seen;
    std::vector<std::set<Point_2>> local_duplicates;

    const int max_threads = omp_get_max_threads();
    local_seen.resize(max_threads);
    local_duplicates.resize(max_threads);

    #pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < n; ++i) {
        int tid = omp_get_thread_num();
        for (std::size_t j = i + 1; j < n; ++j) {
            if (!can_segments_intersect(segments[i], segments[j]))
                continue;

            const auto i1 = CGAL::intersection(segments[i], segments[j]);
            if (!i1) continue;

            if (const auto p = std::get_if<Point_2>(&(*i1))) {
                Point_2 pt = *p;

                if (local_seen[tid].contains(pt)) {
                    local_duplicates[tid].insert(pt);
                } else {
                    local_seen[tid].insert(pt);
                }
            }
        }
    }

    // Merge thread-local sets into global ones
    std::set<Point_2> seen_all;
    std::set<Point_2> degenerate_all;

    for (int t = 0; t < max_threads; ++t) {

        for (const auto &pt : local_seen[t]) {
            if (seen_all.contains(pt)) {
                degenerate_all.insert(pt);
            } else {
                seen_all.insert(pt);
            }
        }

        for (const auto &pt : local_duplicates[t]) {
            degenerate_all.insert(pt);
        }
    }

    return {seen_all.size(), degenerate_all.size()};
}

size_t computeRepeatedVertices(const vector<Point_2> &points)
{
    std::set<Point_2> uniqueVertices;
    for (std::size_t i = 0; i < points.size(); ++i) 
    {
        uniqueVertices.insert(points[i]);
    }
    return points.size() - uniqueVertices.size();
}

pair<vector<Point_2>, vector<Segment_2>> computePointsAndIndices(Data *data)
{
    std::vector<Point_2> arrangementPoints;
    std::map<Point_2, int> arrangementPointsIdices;

    arrangementPoints.resize(data->vertexCoordinatesF.size());
    for (int i = 0 ; i < data->vertexCoordinatesF.size() ; i++)
    {
        const float u = data->vertexCoordinatesF[i];
        const float v = data->vertexCoordinatesG[i];
        const Point_2 point(u, v);

        arrangementPoints[i] = point;
        arrangementPointsIdices[point] = i;
    };

    std::map<std::set<size_t>, bool> uniqueEdges;
    for (int i = 0 ; i < data->tetrahedra.size() ; i++)
    {
        // Bottom face
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][0], data->tetrahedra[i][1]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][1], data->tetrahedra[i][2]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][2], data->tetrahedra[i][0]})] = true;

        // Connect bottom face to top vertex
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][0], data->tetrahedra[i][3]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][1], data->tetrahedra[i][3]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][2], data->tetrahedra[i][3]})] = true;
    }

    std::vector<Segment_2> segments;
    for (const auto& edge : uniqueEdges) 
    {
        // Put in a vector for easy access
        std::vector<int> edgeVector(edge.first.begin(), edge.first.end());
        assert(edgeVector.size() == 2);

        segments.push_back(Segment_2(arrangementPoints[edgeVector[0]], arrangementPoints[edgeVector[1]]));
        //std::cout << "Adding edge " << edgeVector[0] << " - " << edgeVector[1] << std::endl;
    }

    return {arrangementPoints, segments};
}


int main(int argc, char* argv[])
{

    // Parse the command line arguments
    CLI::App cliApp("Fiber Visualiser");

    string filename;
    cliApp.add_option("--file, -f", filename, "Input data filename. Has to be either .txt of .vti.")->required();

    CLI11_PARSE(cliApp, argc, argv);

    fs::path filePath(filename);
    
    if (!fs::exists(filePath)) 
    {
        std::cerr << "Error: File does not exist: " << filename << std::endl;
        return 0;
    }

    // Set up the data class
    Data* data = new Data();

    std::string extension = filePath.extension().string();
    if (extension == ".vtu") 
    {
        data->readDataVTU(filename, 0.0);
    } 
    else if (extension == ".txt") 
    {
        data->readData(filename, 0.0);
    } 
    else 
    {
        std::cerr << "Error: Unsupported file type: " << extension << std::endl;
    }


    const auto [points, segments] = computePointsAndIndices(data);

    printf("There are %ld points and %ld segments.\n", points.size(), segments.size());
    //return 0;


    Timer::start();
    cout << "Computing repeated vertices..." << endl;
    size_t repeatedVertices = computeRepeatedVertices(points);
    Timer::stop("Computed repeated vertices             :");

    Timer::start();
    cout << "Computing degenerate intersections..." << endl;
    auto [intersectionPoints, degenerateIntesectionPoints] = computeDegenerateIntersections(segments);
    Timer::stop("Computed degenerate intersections      :");






    Timer::start();
    cout << "Computing collinear points..." << endl;
    auto [collinearPoints, uniqueLines] = computeDegenerateIntersections(segments);
    Timer::stop("Computed collinear points              :");











    //Timer::start();
    //size_t concurrentSegmentsBruteForce = computeDegenerateIntersectionsBruteForce(segments);
    //Timer::stop("Computed concurrent segments           :");


    //Timer::start();
    //Arrangement_2 arr;
    //CGAL::insert(arr, segments.begin(), segments.end());
    //Timer::stop("Computed arrangement                   :");

    cout << endl;
    printf("There are %ld points and %ld segments and %ld intersection points.\n", points.size(), segments.size(), intersectionPoints);
    //printf("The arrangement has %ld vertices, %ld edges and %ld faces.\n\n", arr.number_of_vertices(), arr.number_of_edges(), arr.number_of_faces());
    printf("There are %ld repeated vertices out of %ld.\n", repeatedVertices, points.size());
    printf("There are %ld collinear points and this many unique lines %ld.\n", collinearPoints, uniqueLines);
    printf("There are %ld degenerate intersections out of %ld.\n", degenerateIntesectionPoints, intersectionPoints);
    //printf("There are %ld concurrent segments.\n", concurrentSegmentsBruteForce);

    std::cout << "Press Enter to continue...";
    std::cin.get();

    return 0;
}
