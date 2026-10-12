/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2013 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// test for merging two empty histogram observables

#include <alps/alea.h>
#include <gtest/gtest.h>

TEST(Histogram, MergingEmptyObservablesPreservesRangeAndZeroCounts)
{
    alps::ObservableSet first, second;
    first << alps::IntHistogramObservable("histogram", 0, 10);
    second << alps::IntHistogramObservable("histogram", 0, 10);
    first << second;
    const alps::IntHistogramObsevaluator result(first["histogram"]);
    EXPECT_EQ(result.count(), 0u);
    ASSERT_EQ(result.size(), 10u);
    EXPECT_EQ(result.min(), 0);
    EXPECT_EQ(result.max(), 10);
    EXPECT_EQ(result.stepsize(), 1);
    for (std::size_t i = 0; i < result.size(); ++i) EXPECT_EQ(result[i], 0u);
}

#include <alps/alea/histogramdata.h>
#include <limits>
#include <sstream>

TEST(Histogram, IntegralStepsCoverTheWholeHalfOpenRange)
{
    for (int step : {1, 2, 3, 5, 6}) {
        alps::IntHistogramObservable histogram("h", -2, 3, step);
        ASSERT_EQ(histogram.size(), (5u + step - 1) / step);
        for (int x = -3; x <= 3; ++x) histogram << x;
        EXPECT_EQ(histogram.count(), 5u);
        unsigned count = 0;
        for (auto value : histogram) count += value;
        EXPECT_EQ(count, 5u);
        EXPECT_GT(histogram[histogram.size() - 1], 0u);
    }
    alps::IntHistogramObservable wide("wide", std::numeric_limits<int>::min(),
        std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
    ASSERT_EQ(wide.size(), 3u);
    wide << std::numeric_limits<int>::max() - 1;
    EXPECT_EQ(wide[2], 1u);
    EXPECT_THROW(alps::IntHistogramObservable("h", 0, 5, 0), std::invalid_argument);
}

TEST(Histogram, XmlEntryCountMustMatchTheDeclaredSize)
{
    const std::string entry = "<ENTRY><COUNT>7</COUNT><VALUE>3</VALUE></ENTRY>";
    for (const auto& xml : {
            "<HISTOGRAM nvalues=\"1\">" + entry + entry + "</HISTOGRAM>",
            "<HISTOGRAM nvalues=\"0\">" + entry + "</HISTOGRAM>",
            "<HISTOGRAM nvalues=\"2\">" + entry + "</HISTOGRAM>"}) {
        std::istringstream in(xml);
        EXPECT_THROW(alps::HistogramObservableData<int> data(in, alps::parse_tag(in)), std::runtime_error);
    }
    std::istringstream in("<HISTOGRAM nvalues=\"2\">" + entry + entry + "</HISTOGRAM>");
    alps::HistogramObservableData<int> data(in, alps::parse_tag(in));
    EXPECT_EQ(data[0], 3u);
    EXPECT_EQ(data[1], 3u);
    EXPECT_EQ(data.count(), 7u); // COUNT is the total, repeated in each ENTRY.
    std::istringstream empty("<HISTOGRAM nvalues=\"0\"></HISTOGRAM>");
    EXPECT_NO_THROW(alps::HistogramObservableData<int> data(empty, alps::parse_tag(empty)));
}
