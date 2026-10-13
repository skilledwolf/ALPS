// SPDX-License-Identifier: MIT
// Unit tests for alps::parse_tag, the pull parser that reads job and task
// files and that XMLParser uses underneath, for processing instructions:
// the XML declaration and xml-stylesheet lines at the top of a file.
#include <alps/parser/parser.h>
#include <alps/parser/xmlhandler.h>
#include <alps/parser/xmlparser.h>
#include <gtest/gtest.h>

#include <sstream>
#include <string>
#include <vector>

namespace {

const std::string declaration = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
const std::string stylesheet = "<?xml-stylesheet type=\"text/xsl\" href=\"ALPS.xsl\"?>\n";
const std::string body = "<SIMULATION>\n<PARAMETERS/>\n</SIMULATION>\n";

// The names of all tags in a text, skipping content.
std::vector<std::string> tags(const std::string& text, bool skip_comments) {
    std::istringstream in(text);
    std::vector<std::string> names;
    char c;
    while (in >> c) {
        in.putback(c);
        if (c == '<')
            names.push_back(alps::parse_tag(in, skip_comments).name);
        else
            alps::parse_content(in);
    }
    return names;
}

const std::vector<std::string> elements = {"SIMULATION", "PARAMETERS", "/SIMULATION"};

TEST(XmlParseTag, ReadsTheAttributesOfTheDeclaration) {
    std::istringstream in(declaration + body);
    alps::XMLTag tag = alps::parse_tag(in, false);
    EXPECT_EQ(tag.type, alps::XMLTag::PROCESSING);
    EXPECT_EQ(tag.name, "xml");
    EXPECT_EQ(tag.attributes["version"], "1.0");
    EXPECT_EQ(tag.attributes["encoding"], "UTF-8");
    EXPECT_EQ(alps::parse_tag(in, false).name, "SIMULATION");
}

// A processing instruction ends at its own "?>", so the declaration alone,
// or followed by a comment, does not hide the rest of the file.
TEST(XmlParseTag, SkipsADeclarationWithoutAStylesheet) {
    EXPECT_EQ(tags(declaration + body, true), elements);
    EXPECT_EQ(tags(declaration + "<!-- a comment -->\n" + body, true), elements);
    EXPECT_EQ(tags(stylesheet + body, true), elements);
    EXPECT_EQ(tags(declaration + stylesheet + body, true), elements);
}

TEST(XmlParseTag, ReportsEachProcessingInstruction) {
    EXPECT_EQ(tags(declaration + stylesheet + body, false),
              (std::vector<std::string>{"xml", "xml-stylesheet", "SIMULATION", "PARAMETERS",
                                        "/SIMULATION"}));
}

class recording_handler : public alps::XMLHandlerBase {
public:
    recording_handler() : alps::XMLHandlerBase("record") {}
    void start_element(const std::string& name, const alps::XMLAttributes& attributes,
                       alps::xml::tag_type type) {
        events.push_back(name);
        if (type == alps::xml::stylesheet)
            style = attributes["href"];
    }
    void end_element(const std::string&, alps::xml::tag_type) {}
    void text(const std::string&) {}
    std::vector<std::string> events;
    std::string style;
};

TEST(XmlParseTag, XmlParserSeesTheElementsAndTheStylesheet) {
    recording_handler only_declaration;
    std::istringstream first(declaration + body);
    alps::XMLParser(only_declaration).parse(first);
    EXPECT_EQ(only_declaration.events,
              (std::vector<std::string>{"xml", "SIMULATION", "PARAMETERS"}));

    recording_handler with_stylesheet;
    std::istringstream second(declaration + stylesheet + body);
    alps::XMLParser(with_stylesheet).parse(second);
    EXPECT_EQ(with_stylesheet.events,
              (std::vector<std::string>{"xml", "xml-stylesheet", "SIMULATION", "PARAMETERS"}));
    EXPECT_EQ(with_stylesheet.style, "ALPS.xsl");
}

} // namespace
