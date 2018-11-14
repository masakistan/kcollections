#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

num_kmers = np.array(
    [
        #1000000,
        #10000000,
        100000000,
        500000000,
        1000000000,
        1500000000,
        2000000000,
        2400000000
    ],
    dtype = np.float
)
xticklabels = num_kmers / 1000000000
# in megabytes
kset_mem = [
    #0.03018,
    #0.06381,
    0.56733,
    2.47606,
    4.491792,
    6.510888,
    8.531176,
    10.171508
]
bft_mem = [
    #0.02172,
    #0.08460,
    0.77856,
    3.61322,
    6.963256,
    10.189268,
    13.347128,
    15.884996
]
set_mem = [
    #0.10582,
    #0.90696,
    11.980,
    48.540,
    97.07,
    191.61,
    194.14,
    220.060
]

# in seconds
kset_insert_time = [
    #3,
    #34,
    408,
    2338,
    4567,
    6733,
    9133,
    11163
]
bft_insert_time =  [
    #3,
    #29,
    351,
    1895,
    4301,
    6977,
    9930,
    12518
]
set_insert_time =  [
    #1,
    #13,
    140,
    790,
    1981,
    3462,
    4416,
    5426
]

kset_query_time = [
    #3,
    #33,
    377,
    2193,
    4720,
    7326,
    9870,
    11812
]
set_query_time = [
    #1,
    #11,
    122,
    746,
    1824,
    3300,
    4559,
    5583
]

fig, ax = plt.subplots()

index = np.arange(len(num_kmers))
bar_width = 0.3

ax.bar(index, kset_mem, bar_width, label='kset')
ax.bar(index + bar_width, bft_mem, bar_width, label='BFT')
ax.bar(index + (bar_width * 2), set_mem, bar_width, label='Python set')
ax.set(
    xlabel='Number of Kmers (Billions)',
    ylabel='Memory Used (GB)',
    title='Memory Usage of kset'
)
ax.set_xticks(index + bar_width)
ax.set_xticklabels(xticklabels)
ax.legend()

fig.savefig('../memory_fig.png')

fig, ax = plt.subplots()

ax.bar(index, kset_insert_time, bar_width, label='kset')
ax.bar(index + bar_width, bft_insert_time, bar_width, label='BFT')
ax.bar(index + (bar_width * 2), set_insert_time, bar_width, label='Python set')
ax.set(
    xlabel='Number of Kmers (Billions)',
    ylabel='Time (seconds)',
    title='Insertion Time of kset'
)
ax.set_xticks(index + bar_width)
ax.set_xticklabels(xticklabels)
ax.legend()

fig.savefig('../insert_fig.png')

fig, ax = plt.subplots()

ax.bar(index, kset_query_time, bar_width, label='kset')
ax.bar(index + bar_width, set_query_time, bar_width, label='Python set')
ax.set(
    xlabel='Number of Kmers (Billions)',
    ylabel='Time (seconds)',
    title='Query Time of kset'
)
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(xticklabels)
ax.legend()

fig.savefig('../query_fig.png')

plt.show()
