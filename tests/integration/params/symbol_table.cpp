// SPDX-License-Identifier: MIT
#include <alps/expression/symbol_table.h>
#include <stdexcept>
int main() {
    // The lattice and model libraries evaluate text symbols; typed
    // parameters reach them through this one conversion.
    alps::params p;
    p["LATTICE"]="chain lattice"; p["L"]=10; p["T"]=2.25;
    const alps::SymbolTable symbols(p);
    if(std::string(symbols["LATTICE"])!="chain lattice" ||
       std::string(symbols["L"])!="10" || std::string(symbols["T"])!="2.25")
        throw std::runtime_error("typed parameters to symbol table conversion failed");
}
