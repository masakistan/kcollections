#!/usr/bin/env python

import matplotlib.pyplot as plt

num_kmers = [1000000, 10000000, 100000000, 500000000, 2400000000]

# in megabytes
kset_mem = [0.02532, 0.06374, 0.560, 2.420, 10.080]
set_mem = [0.10582, 0.90696, 11.980, 48.540, 220.060]

# in seconds
kset_insert_time = [3, 36, 465, 2473, 12298]
set_insert_time = [1, 13, 140, 790, 5426]

kset_query_time = [3, 33, 377, 2193, 60720]
set_query_time = [1, 11, 122, 746, 7334]

fig, ax = plt.subplots()

ax.plot(num_kmers, kset_mem, label='kset')
ax.plot(num_kmers, set_mem, label='Python set')
ax.set(xlabel='Number of Kmers', ylabel='Memory Used (GB)',
       title='Memory Usage of kset')
ax.legend()

fig.savefig('../memory_fig.png')

fig, ax = plt.subplots()

ax.plot(num_kmers, kset_insert_time, label='kset')
ax.plot(num_kmers, set_insert_time, label='Python set')
ax.set(xlabel='Number of Kmers', ylabel='Time (seconds)',
       title='Insertion Time of kset')
ax.legend()

fig.savefig('../insert_fig.png')

fig, ax = plt.subplots()

ax.plot(num_kmers, kset_query_time, label='kset')
ax.plot(num_kmers, set_query_time, label='Python set')
ax.set(xlabel='Number of Kmers', ylabel='Time (seconds)',
       title='Query Time of kset')
ax.legend()

fig.savefig('../query_fig.png')
