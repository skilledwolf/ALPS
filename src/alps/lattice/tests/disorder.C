// SPDX-License-Identifier: MIT
#include <alps/lattice/disorder.h>
#include <gtest/gtest.h>
#include <sstream>

TEST(LatticeDisorder, ChangedEdgePreservesBothEndpoints)
{
    for (const std::string endpoints : {
        "<SOURCE cell=\"0\" vertex=\"2\"/><TARGET cell=\"3\" vertex=\"7\"/>",
        "<SOURCE cell=\"0\" vertex=\"2\"></SOURCE><TARGET cell=\"3\" vertex=\"7\"></TARGET>"}) {
        std::istringstream in("<EDGE type=\"4\">" + endpoints + "</EDGE>");
        alps::detail::EdgeReference edge(alps::parse_tag(in), in);
        EXPECT_EQ(edge.source().vertex(), 2u);
        EXPECT_EQ(edge.target().vertex(), 7u);
        EXPECT_EQ(edge.source().cell_offset()[0], 0);
        EXPECT_EQ(edge.target().cell_offset()[0], 3);
        EXPECT_EQ(edge.new_type(), 4);
        std::istringstream descriptor_input("<CHANGED><EDGE type=\"4\">" + endpoints + "</EDGE></CHANGED></LATTICEGRAPH>");
        auto tag = alps::parse_tag(descriptor_input);
        alps::InhomogeneityDescriptor descriptor(tag, descriptor_input);
        std::ostringstream text;
        alps::oxstream out(text);
        descriptor.write_xml(out);
        EXPECT_NE(text.str().find("vertex=\"2\""), std::string::npos);
        EXPECT_NE(text.str().find("vertex=\"7\""), std::string::npos);
    }
}
