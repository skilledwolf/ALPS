// SPDX-License-Identifier: MIT
// Unit tests for alps/scheduler/seeds.h against the reference table that the
// pyalps implementation of task seeds is also checked against.
#include <alps/parameter.h>
#include <alps/scheduler/seeds.h>
#include <boost/lexical_cast.hpp>
#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace alps::scheduler;

std::vector<std::vector<std::string> > rows(const std::string& layer) {
    std::ifstream table(ALPS_SEED_VECTORS);
    EXPECT_TRUE(table.good()) << ALPS_SEED_VECTORS;
    std::vector<std::vector<std::string> > out;
    std::string line;
    while (std::getline(table, line)) {
        std::istringstream fields(line);
        std::vector<std::string> row;
        std::string field;
        while (fields >> field)
            row.push_back(field);
        if (!row.empty() && row[0] == layer)
            out.push_back(std::vector<std::string>(row.begin() + 1, row.end()));
    }
    EXPECT_FALSE(out.empty()) << layer;
    return out;
}

alps::Parameters with(const std::string& name, const std::string& value) {
    alps::Parameters p;
    if (value != "-")
        p[name] = value;
    return p;
}

TEST(SchedulerSeeds, DerivesTaskSeedsLikePyalps) {
    for (const auto& row : rows("task")) {
        const auto baseseed = static_cast<uint64_t>(boost::lexical_cast<long long>(row[0]));
        EXPECT_EQ(task_seed(baseseed, boost::lexical_cast<uint64_t>(row[1]),
                            boost::lexical_cast<uint64_t>(row[2])),
                  boost::lexical_cast<uint32_t>(row[3]))
            << row[0] << " " << row[1] << " " << row[2];
    }
}

TEST(SchedulerSeeds, AdvancesRunSeeds) {
    for (const auto& row : rows("run")) {
        const alps::StringValue seed(row[0]);
        const auto cpus = boost::lexical_cast<uint32_t>(row[1]);
        if (row[2] == "error")
            EXPECT_THROW(next_run_seed(seed, cpus), boost::bad_lexical_cast) << row[0];
        else
            EXPECT_EQ(next_run_seed(seed, cpus), boost::lexical_cast<uint32_t>(row[2])) << row[0];
    }
}

TEST(SchedulerSeeds, SeedsWorkerEngines) {
    for (const auto& row : rows("worker")) {
        const auto p = with("SEED", row[0]);
        if (row[1] == "error")
            EXPECT_THROW(worker_seed(p), boost::bad_lexical_cast) << row[0];
        else
            EXPECT_EQ(worker_seed(p), boost::lexical_cast<uint32_t>(row[1])) << row[0];
    }
}

TEST(SchedulerSeeds, SeedsTheDisorderGenerator) {
    for (const auto& row : rows("disorder"))
        EXPECT_EQ(disorder_seed(with("DISORDERSEED", row[0])),
                  boost::lexical_cast<unsigned int>(row[1]))
            << row[0];
}

TEST(SchedulerSeeds, LeavesParametersUnchanged) {
    const alps::Parameters empty;
    EXPECT_THROW(worker_seed(empty), boost::bad_lexical_cast);
    EXPECT_EQ(disorder_seed(empty), 0u);
    EXPECT_TRUE(empty.empty());
}

// pyalps divides by zero from 2^30 tasks on; the C++ version says why.
TEST(SchedulerSeeds, RejectsTaskCountsWithoutRoomForTheIndex) {
    EXPECT_NO_THROW(task_seed(1, 0, (uint64_t(1) << 30) - 1));
    EXPECT_THROW(task_seed(1, 0, uint64_t(1) << 30), std::invalid_argument);
}

} // namespace
