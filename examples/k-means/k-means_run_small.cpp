/*******************************************************************************
 * examples/k-means/k-means_run.cpp
 *
 * Part of Project Thrill - http://project-thrill.org
 *
 * Copyright (C) 2015 Matthias Stumpp <mstumpp@gmail.com>
 * Copyright (C) 2016 Timo Bingmann <tb@panthema.net>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <thrill/api/gather.hpp>
#include <thrill/api/generate.hpp>
#include <thrill/api/read_lines.hpp>
#include <thrill/api/zip_with_index.hpp>
#include <thrill/common/logger.hpp>
#include <thrill/common/string.hpp>
#include <tlx/cmdline_parser.hpp>
#include <thrill/api/all_gather.hpp>
#include <thrill/api/cache.hpp>
#include <thrill/api/collapse.hpp>
#include <thrill/api/reduce_by_key.hpp>
#include <thrill/api/sample.hpp>
#include <thrill/api/sum.hpp>
#include <thrill/api/size.hpp>
#include <thrill/api/zip.hpp>
#include <thrill/common/vector.hpp>

#include <cereal/types/vector.hpp>
#include <thrill/data/serialization_cereal.hpp>

#include <limits>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace thrill; // NOLINT

template <size_t D>
using Point = thrill::common::Vector<D, double>;

template <typename Point>
using PointClusterId = std::pair<Point, size_t>;

//! A point which contains "count" accumulated vectors.
template <typename Point>
struct CentroidAccumulated {
    Point  p;
    size_t count;

    template <typename Archive>
    void serialize(Archive& archive) {
        archive(p, count);
    }
};

//! Assignment of a point to a cluster, which is the input to
template <typename Point>
struct ClosestCentroid {
    size_t                     cluster_id;
    CentroidAccumulated<Point> center;

    template <typename Archive>
    void serialize(Archive& archive) {
        archive(cluster_id, center);
    }
};

//! Model returned by KMeans algorithm containing results.
template <typename Point>
class KMeansModel
{
public:
    KMeansModel(size_t dimensions, size_t num_clusters, size_t iterations,
                const std::vector<Point>& centroids)
        : dimensions_(dimensions), num_clusters_(num_clusters),
          iterations_(iterations), centroids_(centroids)
    { }

    //! \name Accessors
    //! \{

    //! Returns dimensions_
    size_t dimensions() const { return dimensions_; }

    //! Returns number of clusters
    size_t num_clusters() const { return num_clusters_; }

    //! Returns iterations_
    size_t iterations() const { return iterations_; }

    //! Returns centroids_
    const std::vector<Point>& centroids() const { return centroids_; }

    //! \}

    //! \name Classification
    //! \{

    //! Calculate closest cluster to point
    size_t Classify(const Point& p) const {
        double min_dist = p.DistanceSquare(centroids_[0]);
        size_t closest_id = 0;
        for (size_t i = 1; i < centroids_.size(); ++i) {
            double dist = p.DistanceSquare(centroids_[i]);
            if (dist < min_dist) {
                min_dist = dist;
                closest_id = i;
            }
        }
        return closest_id;
    }

    //! Calculate closest cluster to all points, returns DIA containing only the
    //! cluster ids.
    template <typename PointDIA>
    auto Classify(const PointDIA& points) const {
        return points
               .Map([this](const Point& p) { return Classify(p); });
    }

    //! Calculate closest cluster to all points, returns DIA contains pairs of
    //! points and their cluster id.
    template <typename PointDIA>
    auto ClassifyPairs(const PointDIA& points) const {
        return points
               .Map([this](const Point& p) {
                        return PointClusterId<Point>(p, Classify(p));
                    });
    }

    //! Calculate the k-means cost: the squared distance to the nearest center.
    double ComputeCost(const Point& p) const {
        double min_dist = p.DistanceSquare(centroids_[0]);
        for (size_t i = 1; i < centroids_.size(); ++i) {
            double dist = p.DistanceSquare(centroids_[i]);
            if (dist < min_dist) {
                min_dist = dist;
            }
        }
        return min_dist;
    }

    //! Calculate the overall k-means cost: the sum of squared distances to
    //! their nearest center.
    template <typename PointDIA>
    double ComputeCost(const PointDIA& points) const {
        return points
               .Map([this](const Point& p) { return ComputeCost(p); })
               .Sum();
    }

    //! \}

private:
    //! dimensions of space
    size_t dimensions_;

    //! number of clusters
    size_t num_clusters_;

    //! number of iterations
    size_t iterations_;

    //! computed centroids in cluster id order
    std::vector<Point> centroids_;
};

//! Calculate k-Means using Lloyd's Algorithm. The DIA centroids is both an
//! input and an output parameter. The method returns a std::pair<Point2D,
//! size_t> = Point2DClusterId into the centroids for each input point.
template <typename Point, typename InStack>
auto KMeans(const DIA<Point, InStack>& input_points, size_t dimensions, size_t num_clusters, size_t iterations, double epsilon = 0.0) {

    thrill::api::Context& ctx = input_points.context();

    auto points = input_points.Cache();

    bool break_condition = false;

    using ClosestCentroid = ClosestCentroid<Point>;
    using CentroidAccumulated = CentroidAccumulated<Point>;

    std::vector<Point> local_centroids =
        points.Keep().Sample(num_clusters).AllGather();

    //auto time_start = std::chrono::high_resolution_clock::now();
    for (size_t iter = 0; iter < iterations && !break_condition; ++iter) {

        std::vector<Point> old_centroids = local_centroids;

        // calculate the closest centroid for each point
        auto closest = points.Keep().Map(
            [&local_centroids](const Point& p) {
                assert(local_centroids.size());
                double min_dist = p.DistanceSquare(local_centroids[0]);
                size_t closest_id = 0;

                for (size_t i = 1; i < local_centroids.size(); ++i) {
                    double dist = p.DistanceSquare(local_centroids[i]);
                    if (dist < min_dist) {
                        min_dist = dist;
                        closest_id = i;
                    }
                }
                return ClosestCentroid {
                    closest_id, CentroidAccumulated { p, 1 }
                };
            });

        // Calculate new centroids as the mean of all points associated with it.
        auto centroids =
            closest
            .ReduceByKey(
                [](const ClosestCentroid& cc) { return cc.cluster_id; },
                [](const ClosestCentroid& a, const ClosestCentroid& b) {
                    return ClosestCentroid {
                        a.cluster_id,
                        CentroidAccumulated { a.center.p + b.center.p,
                                              a.center.count + b.center.count }
                    };
                })
            .Map([](const ClosestCentroid& cc) {
                     return CentroidAccumulated {
                         cc.center.p / static_cast<double>(cc.center.count),
                         cc.cluster_id
                     };
                 })
            .Collapse();

        // collect centroids again, and put back into cluster order
        for (const CentroidAccumulated& uc : centroids.AllGather()) {
            local_centroids[uc.count] = uc.p;
        }

        // Check whether centroid positions changed significantly, if yes do
        // another iteration. only check if epsilon > 0, otherwise we run a
        // fixed number of iterations.
        if (epsilon > 0) {
            break_condition = true;

            for (size_t i = 0; i < local_centroids.size(); ++i) {
                if (local_centroids[i].Distance(old_centroids[i]) > epsilon) {
                    break_condition = false;
                    break;
                }
            }
        }
        //auto time_now = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> diff = time_now - time_start;
        //if (ctx.my_rank() == 0) {
        //    std::cout << "step " << iter << ", time: " << diff.count() << " s" << std::endl;
        //}
    }

    return KMeansModel<Point>(
        dimensions, num_clusters, iterations,
        local_centroids);
}

template <typename Point>
static void RunKMeansFile(
    thrill::Context& ctx,
    size_t dimensions, size_t num_clusters, size_t iterations, double eps,
    const std::vector<std::string>& input_paths) {

    common::StatsTimerStart timer;
    auto points =
        ReadLines(ctx, input_paths).Map(
            [dimensions](const std::string& input) {
                // parse "<pt> <pt> <pt> ..." lines
                Point p = Point::Make(dimensions);
                char* endptr = const_cast<char*>(input.c_str());
                // jump first item
                while (*endptr != ' ') ++endptr;
                for (size_t i = 0; i < dimensions; ++i) {
                    while (*endptr == ' ') ++endptr;
                    p.x[i] = std::strtod(endptr, &endptr);
                    if (!endptr || (*endptr != ' ' && i != dimensions - 1)) {
                        die("Could not parse point coordinates: " << input);
                    }
                }
                while (*endptr == ' ') ++endptr;
                if (!endptr || *endptr != 0) {
                    die("Could not parse point coordinates: " << input);
                }
                return p;
            });

    auto result = KMeans(points.Keep(), dimensions, num_clusters, iterations, eps);

    double cost = result.ComputeCost(points.Keep());
    if (ctx.my_rank() == 0)
        std::cout << "k-means cost: " << cost << std::endl;

    ctx.net.Barrier();
    timer.Stop();
    //if (ctx.my_rank() == 0) {
    //    auto traffic = ctx.net_manager().Traffic();
    //    std::cout << "RESULT"
    //         << " benchmark=k-means"
    //         << " bisecting=" << bisecting
    //         << " dimensions=" << dimensions
    //         << " num_clusters=" << num_clusters
    //         << " iterations=" << iterations
    //         << " eps=" << eps
    //         << " cost=" << cost
    //         << " time=" << timer
    //         << " traffic=" << traffic.total()
    //         << " hosts=" << ctx.num_hosts();
    //}
}

int main(int argc, char* argv[]) {

    tlx::CmdlineParser clp;

    size_t iterations = 10;
    clp.add_size_t('n', "iterations", iterations,
                   "iterations, default: 10");

    size_t dimensions = 64;
    clp.add_param_size_t("dim", dimensions,
                         "dimensions of points 2-10, default: 2");

    size_t num_clusters;
    clp.add_param_size_t("clusters", num_clusters, "Number of clusters");

    double epsilon = 0;
    clp.add_double('e', "epsilon", epsilon,
                   "centroid position delta for break condition, default: 0");

    std::vector<std::string> input_paths;
    clp.add_param_stringlist("input", input_paths,
                             "input file pattern(s)");

    if (!clp.process(argc, argv)) {
        return -1;
    }

    clp.print_result();

    auto start_func =
        [&](thrill::Context& ctx) {
            ctx.enable_consume();
            RunKMeansFile<Point<64> >(
                ctx, dimensions, num_clusters, iterations,
                epsilon, input_paths);
        };

    return thrill::Run(start_func);
}

/******************************************************************************/
