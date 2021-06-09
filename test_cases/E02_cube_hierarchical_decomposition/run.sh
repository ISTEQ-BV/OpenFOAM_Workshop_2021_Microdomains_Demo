#!/bin/bash

set -e

# II. Cell renumbering

pushd original_case

# 1. Multilevel decomposition: scotch + scotch
decomposePar -cellDist -decomposeParDict system/decomposeParDict_multilevel -fileHandler collated -force

# 2. Reconstruct cellDist field to be used for the next step
#reconstructPar -fields '(cellDist)' -withZero

# 3. Generate new cell ordering based on cellDist
# This creates file 'constant/cell_order.txt'
# This also creates 'constant/nodes.txt' for the final decomposition
# And fields 'node' and 'microdomain'
# Number of nodes must match number of nodes used in decomposition
../../build/src/cell_ordering/cell_ordering -nodes 20 -fileHandler collated

# Reorder cells based on generated ordering
# This writes renumbered mesh and fields into directory '1'
renumberMesh -dict system/renumberMeshDict -fileHandler collated

popd

# III. Prepare new case

mkdir -p new_case
cp -r original_case/system new_case
cp -r original_case/constant new_case
rm -rf new_case/0
cp -r original_case/1 new_case/0
rm -r new_case/constant/polyMesh
mv new_case/0/polyMesh new_case/constant/

pushd new_case

decomposePar -cellDist -decomposeParDict system/decomposeParDict_manual -force
#decomposePar -cellDist -decomposeParDict system/decomposeParDict_manual -fileHandler collated -force

popd
