from pybustools.pybus import merge_iterators
import pytest


def test_merge():
    """
    check that values are grouped correctly
    """
    I1 = zip(['A', 'B', 'C'], [1, 2, 3])
    I2 = zip(['A', 'B', 'D'], [10, 20, 30])

    M = merge_iterators({'I1': I1, 'I2': I2})

    mlist = list(M)

    assert mlist[0] == ('A', {'I1': 1, 'I2': 10})
    assert mlist[1] == ('B', {'I1': 2, 'I2': 20})
    assert mlist[2] == ('C', {'I1': 3})
    assert mlist[3] == ('D', {'I2': 30})


def test_merge_fail_unsorted():
    """
    should raise exeption if the two iterators are not sorted
    """

    I1 = zip(['B', 'A', 'C'], [1, 2, 3])
    I2 = zip(['A', 'B', 'D'], [10, 20, 30])

    with pytest.raises(AssertionError):
        list(merge_iterators({'I1': I1, 'I2': I2}))  # list( ) to pull the items of the generator!
