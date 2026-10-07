# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Ping Nang Ma
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************
from pyalps.dataset import DataSet
from pyalps.hlist import HList
from pyalps.tools import groupSets


def test_hlist():
    
    hl = HList([[1,2,3],[4,5]])
    
    print(hl)
    assert list(hl) == [1,2,3,4,5]
    # [[1,2,3],[4,5]]
    
    # !!! Testing linear access
    
    print(hl[0])
    assert hl[0] == 1
    # 1
    
    print(hl[0:2])
    assert hl[0:2] == [1, 2]
    # [1, 2]
    
    # !!! Testing 'recursive' access
    
    print(hl[0,0])
    assert hl[0,0] == 1
    # 1
    print(hl[1,1])
    assert hl[1,1] == 5
    # 5
    
    # !!! Linear assignment
    hl[0] = 27
    print(hl[0])
    assert hl[0] == 27
    print(hl[0,0])
    assert hl[0,0] == 27
    # 27
    
    hl[1,1] = 13
    print(hl[1,1])
    assert hl[1,1] == 13 
    print(hl[4])
    assert hl[4] == 13
    # 13



def test_group_sets_regroups_each_level():
    sets = []
    for a in (1, 2):
        for b in (1, 2):
            sets.append(DataSet())
            sets[-1].props.update(a=a, b=b)
    by_a = groupSets(sets, ['a'])
    assert [[s.props['a'] for s in group] for group in by_a] == [[1, 1], [2, 2]]
    nested = groupSets(by_a, ['b'])
    assert [[[(s.props['a'], s.props['b']) for s in group] for group in outer] for outer in nested] == \
        [[[(1, 1)], [(1, 2)]], [[(2, 1)], [(2, 2)]]]


if __name__ == '__main__':
    test_hlist()
