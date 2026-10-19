// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// XML lattice and model libraries with n-site terms (<NSITE>, <NSITETERM>;
// ALPSim/ALPS PR 223) for the catalog importer's tests: groups whose roles are
// out of site order or wrap around a periodic ladder, fermion products on
// two, three and four sites, two fermion operators on one role, a constant,
// spin ring exchange and truncated bosons.
#pragma once

#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>
#include <unistd.h>

namespace nsite_test {

inline char const* const lattices = R"xml(<LATTICES>
<LATTICE name="chain lattice" dimension="1"><BASIS><VECTOR>1</VECTOR></BASIS></LATTICE>
<UNITCELL name="ladder cell" dimension="1" vertices="2">
  <VERTEX id="1"/><VERTEX id="2"/>
  <EDGE type="1"><SOURCE vertex="1" offset="0"/><TARGET vertex="2" offset="0"/></EDGE>
  <EDGE type="0"><SOURCE vertex="1" offset="0"/><TARGET vertex="1" offset="1"/></EDGE>
  <EDGE type="0"><SOURCE vertex="2" offset="0"/><TARGET vertex="2" offset="1"/></EDGE>
  <NSITE type="0">
    <SITE vertex="1" offset="0"/><SITE vertex="2" offset="0"/>
    <SITE vertex="2" offset="1"/><SITE vertex="1" offset="1"/>
  </NSITE>
</UNITCELL>
<LATTICEGRAPH name="open ladder">
  <FINITELATTICE><LATTICE ref="chain lattice"/><EXTENT dimension="1" size="L"/>
    <BOUNDARY type="open"/></FINITELATTICE>
  <UNITCELL ref="ladder cell"/>
</LATTICEGRAPH>
<LATTICEGRAPH name="periodic ladder">
  <FINITELATTICE><LATTICE ref="chain lattice"/><EXTENT dimension="1" size="L"/>
    <BOUNDARY type="periodic"/></FINITELATTICE>
  <UNITCELL ref="ladder cell"/>
</LATTICEGRAPH>
<GRAPH name="scrambled" vertices="5">
  <EDGE source="1" target="2"/><EDGE source="2" target="3"/><EDGE source="3" target="4"/><EDGE source="4" target="5"/>
  <NSITE type="0"><SITE id="4"/><SITE id="1"/><SITE id="3"/></NSITE>
  <NSITE type="0"><SITE id="2"/><SITE id="5"/><SITE id="1"/></NSITE>
  <NSITE type="1"><SITE id="5"/><SITE id="3"/><SITE id="2"/><SITE id="1"/></NSITE>
</GRAPH>
</LATTICES>
)xml";

inline char const* const models = R"xml(<MODELS>
<SITEBASIS name="fermion">
  <QUANTUMNUMBER name="Nup" min="0" max="1" type="fermionic"/>
  <QUANTUMNUMBER name="Ndown" min="0" max="1" type="fermionic"/>
  <OPERATOR name="c_down" matrixelement="1"><CHANGE quantumnumber="Ndown" change="-1"/></OPERATOR>
  <OPERATOR name="cdag_down" matrixelement="1"><CHANGE quantumnumber="Ndown" change="1"/></OPERATOR>
  <OPERATOR name="c_up" matrixelement="1"><CHANGE quantumnumber="Nup" change="-1"/></OPERATOR>
  <OPERATOR name="cdag_up" matrixelement="1"><CHANGE quantumnumber="Nup" change="1"/></OPERATOR>
  <OPERATOR name="n" matrixelement="Nup+Ndown"/>
  <OPERATOR name="n_up" matrixelement="Nup"/>
  <OPERATOR name="n_down" matrixelement="Ndown"/>
</SITEBASIS>
<BASIS name="fermion"><SITEBASIS ref="fermion"/>
  <CONSTRAINT quantumnumber="Nup" value="Nup_total"/><CONSTRAINT quantumnumber="Ndown" value="Ndown_total"/></BASIS>
<SITEBASIS name="spin half">
  <QUANTUMNUMBER name="Sz" min="-1/2" max="1/2"/>
  <OPERATOR name="Splus" matrixelement="1"><CHANGE quantumnumber="Sz" change="1"/></OPERATOR>
  <OPERATOR name="Sminus" matrixelement="1"><CHANGE quantumnumber="Sz" change="-1"/></OPERATOR>
  <OPERATOR name="Sz" matrixelement="Sz"/>
</SITEBASIS>
<BASIS name="spin half"><SITEBASIS ref="spin half"/><CONSTRAINT quantumnumber="Sz" value="Sz_total"/></BASIS>
<SITEBASIS name="boson"><PARAMETER name="Nmax" default="2"/>
  <QUANTUMNUMBER name="N" min="0" max="Nmax"/>
  <OPERATOR name="bdag" matrixelement="sqrt(N+1)"><CHANGE quantumnumber="N" change="1"/></OPERATOR>
  <OPERATOR name="b" matrixelement="sqrt(N)"><CHANGE quantumnumber="N" change="-1"/></OPERATOR>
  <OPERATOR name="n" matrixelement="N"/>
</SITEBASIS>
<BASIS name="boson"><SITEBASIS ref="boson"/><CONSTRAINT quantumnumber="N" value="N_total"/></BASIS>
<HAMILTONIAN name="pair hopping">
  <PARAMETER name="t" default="1"/><PARAMETER name="U" default="4"/><PARAMETER name="K" default="0.7"/>
  <PARAMETER name="J" default="0.3"/><PARAMETER name="V" default="0.45"/><PARAMETER name="C" default="0.25"/>
  <BASIS ref="fermion"/>
  <SITETERM site="i">U*n_up(i)*n_down(i)</SITETERM>
  <BONDTERM source="i" target="j">-t*(cdag_up(i)*c_up(j)+cdag_up(j)*c_up(i)+cdag_down(i)*c_down(j)+cdag_down(j)*c_down(i))</BONDTERM>
  <NSITETERM sites="i j k l">
    -K*( (cdag_up(l)*cdag_down(k) - cdag_down(l)*cdag_up(k)) * (c_down(j)*c_up(i) - c_up(j)*c_down(i))
       + (cdag_up(i)*cdag_down(j) - cdag_down(i)*cdag_up(j)) * (c_down(k)*c_up(l) - c_up(k)*c_down(l)) )
    + J*(cdag_up(i)*c_down(i)*cdag_down(k)*c_up(k) + cdag_down(i)*c_up(i)*cdag_up(k)*c_down(k))
    + V*n_up(j)*(cdag_down(i)*c_down(l) + cdag_down(l)*c_down(i)) + C
  </NSITETERM>
</HAMILTONIAN>
<HAMILTONIAN name="three-site fermions">
  <PARAMETER name="t" default="1"/><PARAMETER name="U" default="2"/><PARAMETER name="V" default="0.6"/>
  <PARAMETER name="W" default="0.35"/>
  <BASIS ref="fermion"/>
  <SITETERM site="i">U*n_up(i)*n_down(i)</SITETERM>
  <BONDTERM source="i" target="j">-t*(cdag_up(i)*c_up(j)+cdag_up(j)*c_up(i)+cdag_down(i)*c_down(j)+cdag_down(j)*c_down(i))</BONDTERM>
  <NSITETERM type="0" sites="a b c">
    V*(cdag_up(a)*c_up(b) + cdag_up(b)*c_up(a))*n_down(c) + W*(cdag_down(c)*cdag_up(a)*c_up(b)*c_down(c) + cdag_down(c)*cdag_up(b)*c_up(a)*c_down(c))
  </NSITETERM>
  <NSITETERM type="1" sites="a b c d">
    W*(cdag_up(a)*cdag_down(c)*c_down(d)*c_up(b) + cdag_up(b)*cdag_down(d)*c_down(c)*c_up(a))
  </NSITETERM>
</HAMILTONIAN>
<HAMILTONIAN name="ring exchange">
  <PARAMETER name="J" default="1"/><PARAMETER name="K" default="0.8"/><PARAMETER name="Jz" default="0.5"/>
  <BASIS ref="spin half"/>
  <BONDTERM source="i" target="j">J*(Sz(i)*Sz(j)+(Splus(i)*Sminus(j)+Sminus(i)*Splus(j))/2)</BONDTERM>
  <NSITETERM sites="i j k l">
    K*(Splus(i)*Sminus(j)*Splus(k)*Sminus(l) + Sminus(i)*Splus(j)*Sminus(k)*Splus(l)) + Jz*Sz(i)*Sz(j)*Sz(k)*Sz(l)
  </NSITETERM>
</HAMILTONIAN>
<HAMILTONIAN name="boson ring hop">
  <PARAMETER name="t" default="1"/><PARAMETER name="U" default="1.5"/><PARAMETER name="K" default="0.4"/>
  <BASIS ref="boson"/>
  <SITETERM site="i">U/2*n(i)*(n(i)-1)</SITETERM>
  <BONDTERM source="i" target="j">-t*(bdag(i)*b(j)+bdag(j)*b(i))</BONDTERM>
  <NSITETERM sites="i j k l">K*(bdag(i)*bdag(k)*b(j)*b(l) + bdag(j)*bdag(l)*b(i)*b(k)) + 0.1</NSITETERM>
</HAMILTONIAN>
</MODELS>
)xml";

// The libraries, written to a temporary directory for the object's lifetime.
struct Libraries {
    std::filesystem::path dir;
    Libraries() {
        static int count = 0;
        dir = std::filesystem::temp_directory_path()
            / ("alps_nsite_libraries_" + std::to_string(::getpid()) + "_" + std::to_string(count++));
        std::filesystem::create_directories(dir);
        std::ofstream(dir / "lattices.xml") << lattices;
        std::ofstream(dir / "models.xml") << models;
    }
    ~Libraries() {
        std::error_code ignored;
        std::filesystem::remove_all(dir, ignored);
    }
};

} // namespace nsite_test
