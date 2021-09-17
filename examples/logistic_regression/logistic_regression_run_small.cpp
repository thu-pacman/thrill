/*******************************************************************************
 * examples/logistic_regression/logistic_regression.cpp
 *
 * Part of Project Thrill - http://project-thrill.org
 *
 * Copyright (C) 2016 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <thrill/api/cache.hpp>
#include <thrill/api/dia.hpp>
#include <thrill/api/generate.hpp>
#include <thrill/api/print.hpp>
#include <thrill/api/read_lines.hpp>
#include <thrill/common/logger.hpp>
#include <thrill/common/string.hpp>
#include <tlx/cmdline_parser.hpp>
#include <thrill/api/size.hpp>
#include <thrill/api/sum.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <array>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace thrill;                        // NOLINT

// Dimensions of the data
constexpr size_t dim = 784;
using T = double;

using Element = std::array<T, dim>;
using DataObject = std::pair<bool, Element>;

#define LOGM LOGC(debug && ctx.my_rank() == 0)

template <typename T>
inline T sigmoid(const T& x) { // can't make it constexpr, exp isn't one
    return 1.0 / (1.0 + exp(-x));
}

template <typename T, size_t dim>
T calc_norm(const std::array<T, dim>& weights,
            const std::array<T, dim>& new_weights) {
    T sum = 0.;
    for (size_t i = 0; i < dim; ++i) {
        T diff = weights[i] - new_weights[i];
        sum += (diff * diff);
    }
    return std::sqrt(sum);
}

template <typename T, size_t dim>
auto gradient(const bool& y, const std::array<T, dim>& x,
              const std::array<T, dim>& w) {
    std::array<T, dim> grad;
    T dot_product = std::inner_product(w.begin(), w.end(), x.begin(), T { 0.0 });
    T s = sigmoid(dot_product) - y;
    for (size_t i = 0; i < dim; ++i) {
        grad[i] = s * x[i];
    }
    return grad;
}

template <typename Input>
static auto ReadInputFile(api::Context& ctx, const Input& input_path) {
    return ReadLines(ctx, input_path)
           .Map([](const std::string& line) {
                    // parse "value,dim_1,dim_2,...,dim_n" lines
                    char* endptr;
                    DataObject obj;
                    // yikes C stuff, TODO template
                    obj.first = common::from_cstr<T>(line.c_str(), &endptr);
                    die_unless(endptr && *endptr == ' ' &&
                               "Could not parse input line");

                    for (size_t i = 0; i < dim; ++i) {
                        T value = common::from_cstr<T>(endptr + 1, &endptr);
                        die_unless(endptr &&
                                   ((i + 1 <= dim && *endptr == ',') ||
                                    (i + 1 == dim && *endptr == 0)) &&
                                   "Could not parse input line");
                        obj.second[i] = value;
                    }
                    return obj;
                })
           .Cache();
}

template <typename T, size_t dim, typename InStack,
          typename Element = std::array<T, dim> >
auto logit_train(const DIA<std::pair<bool, Element>, InStack>& data,
                 size_t max_iterations, double gamma = 0.002,
                 double epsilon = 0.0001) {
    thrill::api::Context& ctx = data.context();

    // weights, initialized to zero
    Element weights, new_weights;
    for (int i = 0; i < dim; ++i)
        weights[i] = 0.0;
    size_t iter = 0;
    T norm = 0.0;

    //auto time_start = std::chrono::high_resolution_clock::now();
    while (iter < max_iterations) {
        Element grad =
            data.Keep()
            .Map([&weights](const std::pair<bool, Element>& elem) -> Element {
                     return gradient(elem.first, elem.second, weights);
                 })
            .Sum([](const Element& a, const Element& b) -> Element {
                     Element result;
                     std::transform(a.begin(), a.end(), b.begin(),
                                    result.begin(), std::plus<T>());
                     return result;
                 });

        std::transform(weights.begin(), weights.end(), grad.begin(),
                       new_weights.begin(),
                       [&gamma](const T& a, const T& b) -> T
                       { return a - gamma * b; });

        //norm = calc_norm(new_weights, weights);
        weights = new_weights;

        //auto time_now = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> diff = time_now - time_start;
        //if (ctx.my_rank() == 0) {
        //    std::cout << "step " << iter << ", time: " << diff.count() << " s" << std::endl;
        //}
        iter++;
        //if (norm < epsilon) break;
    }

    return std::make_tuple(weights, norm, iter);
}

template <typename InputDIA>
auto TrainLogit(api::Context& ctx,
                const InputDIA& input_dia,
                size_t max_iterations, double gamma, double epsilon) {

    Element weights;
    double norm;
    size_t iterations;
    std::tie(weights, norm, iterations) =
        logit_train<T, dim>(input_dia.Keep(), max_iterations, gamma, epsilon);

    if (ctx.my_rank() == 0) {
        std::cout << "Iterations: " << iterations << std::endl;
        std::cout << "Norm: " << norm << std::endl;
        std::cout << "Final weights (model):" << std::endl;
    }
    return weights;
}

int main(int argc, char* argv[]) {
    tlx::CmdlineParser clp;

    std::string training_path;
    std::vector<std::string> test_path;
    clp.add_param_string("input", training_path, "training file pattern(s)");
    //clp.add_param_stringlist("test", test_path, "test file pattern(s)");

    size_t max_iterations = 1000;
    clp.add_size_t('n', "iterations", max_iterations,
                   "Maximum number of iterations, default: 1000");

    double gamma = 0.002, epsilon = 0.0001;
    clp.add_double('g', "gamma", gamma, "Gamma, default: 0.002");
    clp.add_double('e', "epsilon", epsilon, "Epsilon, default: 0.0001");

    if (!clp.process(argc, argv)) {
        return -1;
    }

    clp.print_result();

    return api::Run(
        [&](api::Context& ctx) {
            ctx.enable_consume();

            common::StatsTimerStart timer;

            Element weights;

            weights = TrainLogit(ctx, ReadInputFile(ctx, training_path),
                                 max_iterations, gamma, epsilon);

            timer.Stop();
            //if (ctx.my_rank() == 0) {
            //    std::cout << "whole time: " << timer << std::endl;
            //}
        });
}

/******************************************************************************/
