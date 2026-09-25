# SPDX-License-Identifier: MIT
"""Plot renderers must keep each dataset's values and uncertainties together."""
import re

import numpy as np
import pytest

from pyalps.dataset import DataSet
from pyalps.floatwitherror import FloatWithError
from pyalps.plot_core import makeGracePlot, makeGnuplotPlot


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("scatter", [False, True])
@pytest.mark.parametrize("labelled", [False, True])
def test_export_mixed_uncertainties(reverse, scatter, labelled):
    cases = [
        ("xy", [1., 2.], "1:2", ""),
        ("xydy", [1., 2., .7], "1:2:3", "y"),
        ("xydx", [1., 2., .1], "1:2:3", "x"),
        ("xydxdy", [1., 2., .1, .7], "1:2:3:4", "xy"),
    ]
    if reverse:
        cases.reverse()
    data = []
    for kind, row, using, errors in cases:
        props = {"line": "scatter" if scatter else "line"}
        if labelled:
            props["label"] = kind
        data.append(DataSet(
            [FloatWithError(1., .1) if "x" in errors else 1.],
            [FloatWithError(2., .7) if "y" in errors else 2.], props))
    # Empty datasets must not leave dangling plot entries or change error flags.
    data.insert(1, DataSet())
    grace = makeGracePlot([data])
    blocks = re.findall(r"@type (\w+)\n([^&]+)&", grace)
    assert [kind for kind, _ in blocks] == [case[0] for case in cases]
    for (_, block), (_, row, _, _) in zip(blocks, cases):
        assert [float(value) for value in block.split()] == row

    gnuplot = makeGnuplotPlot([data])
    command = next(line for line in gnuplot.splitlines() if line.startswith("plot"))
    entries = command[4:].split(",")
    assert len(entries) == len(cases)
    rows = re.findall(r"# X Y[^\n]*\n([^\n]+)\nend", gnuplot)
    assert len(rows) == len(cases)
    for entry, row, (kind, expected, using, errors) in zip(entries, rows, cases):
        assert [float(value) for value in row.split()] == expected
        title = f'title "{kind}"' if labelled else "notitle"
        style = f" w {errors}error{'bars' if scatter else 'lines'}" if errors else ""
        assert entry.strip() == f'"-" using {using}{style} {title}'


@pytest.mark.parametrize("representation", ["scalar", "numpy_scalar", "array", "list", "numeric_strings", "error_scalar", "error_array", "error_list"])
def test_plot_representations(representation):
    values = {
        "scalar": 2., "numpy_scalar": np.array(2.), "array": np.array([2.]), "list": [2.],
        "numeric_strings": ["2.0"],
        "error_scalar": FloatWithError(2., .7),
        "error_array": FloatWithError(np.array([2.]), np.array([.7])),
        "error_list": [FloatWithError(2., .7)],
    }
    data = DataSet(np.array(1.), values[representation])
    expected = [1., 2., .7] if representation.startswith("error") else [1., 2.]
    for export in (makeGracePlot, makeGnuplotPlot):
        text = export([data])
        row = next(line for line in text.splitlines() if line.startswith("1.0\t"))
        assert [float(value) for value in row.split()] == expected

    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import pyalps.plot as plot
    for core in (False, True):
        figure, axes = plot.plt.subplots()
        try:
            if core:
                renderer = plot.MplXYPlot_core()
                renderer.plt = {"data": [data]}
                renderer.draw_lines()
            else:
                plot.plot(data)
            np.testing.assert_array_equal(axes.lines[0].get_xdata(), [1.])
            np.testing.assert_array_equal(axes.lines[0].get_ydata(), [2.])
            bars = axes.containers[0].lines[2]
            assert len(bars) == (1 if len(expected) == 3 else 0)
            if bars:
                np.testing.assert_allclose(bars[0].get_segments()[0], [[1., 1.3], [1., 2.7]])
        finally:
            plot.plt.close(figure)


@pytest.mark.parametrize("export", [makeGracePlot, makeGnuplotPlot])
def test_empty_and_mismatched_exports(export):
    for data in ([], [DataSet()]):
        text = export(data)
        assert '@type ' not in text and '\nplot' not in text
    with pytest.raises(ValueError, match="same length"):
        export([DataSet([1., 2.], [3.])])
