import sys

import kindel

import pandas as pd

import plotly.offline as py
import plotly.graph_objs as go


def coverage(weights):
    return [sum(weight.values()) for weight in weights]


def plot(coverage, l_clip_starts, r_clip_starts):
    x_axis = list(range(9100))
    t0 = go.Scattergl(
        x = x_axis,
        y = coverage,
        mode = 'lines',
        name = 'coverage')
    t1 = go.Scattergl(
        x = x_axis,
        y = l_clip_starts,
        mode = 'markers',
        name = 'left clip starts')
    t2 = go.Scattergl(
        x = x_axis,
        y = r_clip_starts,
        mode = 'markers',
        name = 'right clip starts')
    layout = go.Layout(
        xaxis=dict(
            type='linear',
            autorange=True),
        yaxis=dict(
            type='linear',
            autorange=True))

    data = [t0, t1, t2]
    fig = go.Figure(data=data, layout=layout)
    py.plot(fig, filename='softclips.html')


if __name__ == '__main__':
    print(sys.argv[1])
    aln = list(kindel.parse_bam(sys.argv[1]).items())[0][1]
    cov = coverage(aln.weights)
    cov_smoothed = pd.Series(cov).rolling(window=10).mean().tolist()
    plot(cov_smoothed, aln.clip_starts, aln.clip_ends)
