#include <cstddef>
#include <set>
#include <filesystem>
#include <unordered_map>
#include <random>
#include <omp.h>
#include <boost/multiprecision/cpp_int.hpp>

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


bool canSegmentsIntersect(const Segment_2& s1, const Segment_2& s2) 
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
}

bool doSegmentsIntersect(const Segment_2& si, const Segment_2& sj, const Segment_2& sk)
{
    if (!canSegmentsIntersect(si, sj) || !canSegmentsIntersect(sj, sk) || !canSegmentsIntersect(si, sk))
    {
        return false;
    }

    auto i1 = CGAL::intersection(si, sj);

    // If there is no intersection just continue
    if (!i1) 
    { 
        return false; 
    }

    auto i2 = CGAL::intersection(sj, sk);

    // If there is no intersection just continue
    if (!i2) 
    { 
        return false; 
    }

    auto i3 = CGAL::intersection(si, sk);

    if (!i3) 
    { 
        return false; 
    }

    auto p1 = std::get_if<Point_2>(&*i1);
    auto p2 = std::get_if<Point_2>(&*i2);
    auto p3 = std::get_if<Point_2>(&*i3);

    if (p1 && p2 && p3 && *p1 == *p2 && *p2 == *p3) 
    {
        return true;
    }

    return false;

}

long long computeConcurrentSegmentsRandom(const std::vector<Segment_2>& segments, const size_t numberOfSamples)
{
    long long  n = segments.size();
    if (n < 3) return 0;

    long long  concurrentCount = 0;

    // Parallel region
    #pragma omp parallel
    {
        // Each thread gets its own RNG seeded differently
        std::mt19937 rng;
        {
            auto seed = std::chrono::steady_clock::now().time_since_epoch().count();
            // Use thread num to diversify seeds
            rng.seed(seed + omp_get_thread_num());
        }
        std::uniform_int_distribution<long long > dist(0, n - 1);

        long long  localCount = 0;

        #pragma omp for
        for (long long  sample = 0; sample < numberOfSamples; ++sample)
        {
            // Sample unordered tripplet
            long long  i, j, k;
            do { i = dist(rng); } while (false);
            do { j = dist(rng); } while (j == i);
            do { k = dist(rng); } while (k == i || k == j);

            if (j < i) std::swap(i, j);
            if (k < j) std::swap(j, k);
            if (j < i) std::swap(i, j);

            const Segment_2 &si = segments[i];
            const Segment_2 &sj = segments[j];
            const Segment_2 &sk = segments[k];


            if (doSegmentsIntersect(si, sj, sk))
            {
                ++localCount;
            }
        }

        // Reduce local counts into global count atomically
        #pragma omp atomic
        concurrentCount += localCount;
    }

    return concurrentCount;
}

long long  computeCollinearPointsRandom(const std::vector<Point_2> &points, const size_t numberOfSamples)
{
    long long  n = points.size();
    if (n < 3) return 0;

    long long  collinearCount = 0;

    // Parallel region
    #pragma omp parallel
    {
        // Each thread gets its own RNG seeded differently
        std::mt19937 rng;
        {
            auto seed = std::chrono::steady_clock::now().time_since_epoch().count();
            // Use thread num to diversify seeds
            rng.seed(seed + omp_get_thread_num());
        }
        std::uniform_int_distribution<long long > dist(0, n - 1);

        long long  localCount = 0;

        #pragma omp for
        for (long long  sample = 0; sample < numberOfSamples; ++sample)
        {
            // Sample unordered tripplet
            long long  i, j, k;
            do { i = dist(rng); } while (false);
            do { j = dist(rng); } while (j == i);
            do { k = dist(rng); } while (k == i || k == j);

            if (j < i) std::swap(i, j);
            if (k < j) std::swap(j, k);
            if (j < i) std::swap(i, j);

            const Point_2 &A = points[i];
            const Point_2 &B = points[j];
            const Point_2 &C = points[k];

            if (CGAL::orientation(A, B, C) == CGAL::COLLINEAR)
            {
                ++localCount;
            }
        }

        // Reduce local counts into global count atomically
        #pragma omp atomic
        collinearCount += localCount;
    }

    return collinearCount;
}












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

pair<size_t, size_t> computeCollinearPointsBruteForce(const vector<Point_2>& points)
{
    set<Point_2> collinearPoints;
    const size_t n = points.size();

    for (size_t i = 0; i < n; ++i)
    {
        cout << i << " / " << n << endl;
        for (size_t j = i + 1; j < n; ++j)
        {
            for (size_t k = j + 1; k < n; ++k)
            {
                if (CGAL::collinear(points[i], points[j], points[k]))
                {
                    collinearPoints.insert(points[i]);
                    collinearPoints.insert(points[j]);
                    collinearPoints.insert(points[k]);
                }
            }
        }
    }


    return {collinearPoints.size(), -1};
}


pair<size_t, size_t> computeCollinearPoints(const vector<Point_2> &points)
{
    //int collinearPoints = 0;
    std::set<Line_2, LineLess> unique_lines; 
    std::set<Point_2> collinearPoints;

    for (int i = 0; i < points.size(); ++i) 
    {
        for (int j = i + 1; j < points.size(); ++j) 
        {
            Line_2 line(points[i], points[j]);

            if (unique_lines.contains(line))
            {
                //collinearPoints++;
                collinearPoints.insert(points[i]);
                collinearPoints.insert(points[j]);
            }
            else
            {
                unique_lines.insert(line);
            }
        }
    }

    return {collinearPoints.size(), unique_lines.size()};
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


size_t computeRepeatedVerticesBruteForce(const vector<Point_2> &points)
{
    int duplicatePairs = 0;
    for (std::size_t i = 0; i < points.size(); ++i) 
    {
        cout << i << " / " << points.size() << endl;
        for (std::size_t j = i + 1; j < points.size(); ++j) 
        {
            if (points[i] == points[j])
            {
                duplicatePairs++;
            }

        }
    }

    return duplicatePairs;
}


pair<vector<Point_2>, vector<Segment_2>> computePointsAndIndices(Data *data)
{
    std::vector<Point_2> arrangementPoints;

    arrangementPoints.resize(data->vertexCoordinatesF.size());
    for (int i = 0 ; i < data->vertexCoordinatesF.size() ; i++)
    {
        const float u = data->vertexCoordinatesF[i];
        const float v = data->vertexCoordinatesG[i];
        const Point_2 point(u, v);

        arrangementPoints[i] = point;
    };


    int duplicatePairs = 0;
    for (std::size_t i = 0; i < arrangementPoints.size(); ++i) 
    {
        cout << i << " / " << arrangementPoints.size() << endl;
        for (std::size_t j = i + 1; j < arrangementPoints.size(); ++j) 
        {
            if (arrangementPoints[i].x() == arrangementPoints[j].x() && arrangementPoints[i].y() == arrangementPoints[j].y())
            {
                duplicatePairs++;
            }

        }
    }

    cout << "There are " << duplicatePairs << " duplicate pairs.\n";


    return {};






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


    bool runRepeatedPoints = false;
    cliApp.add_flag("--repeated-points, -r", runRepeatedPoints, "Only compute the number of repeated vertices.");

    bool runCollinearPoints = false;
    cliApp.add_flag("--collinear-points, -p", runCollinearPoints, "Only compute the number of collinear points.");

    bool runConcurrentSegments = false;
    cliApp.add_flag("--concurrent-segments, -s", runConcurrentSegments, "Only compute the number of concurrent segments.");

    bool runAll = false;
    cliApp.add_flag("--all, -a", runAll, "Compute everything.");

    long long numberOfSamples = -1;
    cliApp.add_option("--number-Samples, -n", numberOfSamples, "Compute everything.");

    CLI11_PARSE(cliApp, argc, argv);

    fs::path filePath(filename);
    
    if (!fs::exists(filePath)) 
    {
        std::cerr << "Error: File does not exist: " << filename << std::endl;
        return 0;
    }

    // Set up the data class
    Data* data = new Data();

    printf("Reading data...\n");
    Timer::start();
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
    Timer::stop("Read data                              :");

    //cout << "Tets = " << data->tetrahedra.size() << endl;
    //return 0;

    printf("Computing points and segments...\n");
    Timer::start();
    const auto [points, segments] = computePointsAndIndices(data);
    Timer::stop("Computed points and segments             :");

    // Triplles of over a miliion will overflow long long 
    const boost::multiprecision::uint128_t numberOfPointTripplets = points.size() * (points.size() - 1) * (points.size() - 2) / 6;
    const boost::multiprecision::uint128_t numberOfSegmentTripplets = segments.size() * (segments.size() - 1) * (segments.size() - 2) / 6;

    printf("There are %ld points and %ld segments in dataset %s.\n", points.size(), segments.size(), filename.c_str());
    //return 0;

    if (true == runAll || true == runRepeatedPoints)
    {
        printf("Computing duplicated pints...\n");
        Timer::start();
        cout << "Computing repeated vertices..." << endl;
        //size_t repeatedVertices = computeRepeatedVertices(points);
        size_t repeatedVertices = computeRepeatedVerticesBruteForce(points);
        Timer::stop("Computed repeated vertices             :");
        printf("------------------------------------------------------------------------There are %ld repeated vertices out of %ld.\n", repeatedVertices, points.size());
    }

    if (true == runAll || true == runCollinearPoints)
    {
        printf("Computing collinear pints...\n");
        if (numberOfSamples != -1)
        {
            Timer::start();
            long long  collinearPoints =  computeCollinearPointsRandom(points, numberOfSamples);
            Timer::stop("Computed collinear points              :");
            printf("------------------------------------------------------------------------There are %lld collinear points from %lld samples and %s tripplets.\n", collinearPoints, numberOfSamples, numberOfPointTripplets.str().c_str());
        }
        else
        {
            Timer::start();
            cout << "Computing collinear points..." << endl;
            auto [collinearPoints, uniqueLines] = computeCollinearPoints(points);
            Timer::stop("Computed collinear points              :");
            printf("------------------------------------------------------------------------There are %ld collinear points and this many unique lines %ld.\n", collinearPoints, uniqueLines);

        }
    }

    if (true == runAll || true == runConcurrentSegments)
    {
        printf("Computing concurrent segments...\n");
        if (numberOfSamples != -1)
        {
            Timer::start();
            long long  concurrentSegments =  computeConcurrentSegmentsRandom(segments, numberOfSamples);
            Timer::stop("Computed collinear points              :");
            printf("------------------------------------------------------------------------There are %lld concurrent segments from %lld samples and %s tripplets.\n", concurrentSegments, numberOfSamples, numberOfSegmentTripplets.str().c_str());
        }
        else
        {
            Timer::start();
            cout << "Computing degenerate intersections..." << endl;
            auto [intersectionPoints, degenerateIntesectionPoints] = computeDegenerateIntersections(segments);
            Timer::stop("Computed degenerate intersections      :");
            printf("------------------------------------------------------------------------There are %ld degenerate intersections out of %ld.\n", degenerateIntesectionPoints, intersectionPoints);
        }
    }

    return 0;
}
