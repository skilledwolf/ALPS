#!/bin/sh
#  Copyright Matthias Troyer 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

for directory in */ alpsize/*/ code/*/ dmft/*/ dmrg/*/ ed/*/ hybridization/[0-9]*/ intro/*/ mc/*/
do
  [ -d "$directory" ] || continue
  rm -f "$directory"/parm*xml "$directory"/*clone* "$directory"/*.xsl \
    "$directory"/*.run* "$directory"/*.dat "$directory"/*.h5* \
    "$directory"/*out.xml "$directory"/*.mp4 "$directory"/*in.xml "$directory"/selfenergy*
  rm -rf "$directory"/*.chkp
done
rm -f parm*xml *clone* *.xsl *.run* *.dat *.h5* *out.xml *in.xml
