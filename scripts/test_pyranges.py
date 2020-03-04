
#%%
import pyranges as pr
from pyranges import PyRanges

import pandas as pd

from io import StringIO

#%%
f1 = """Chromosome Start End Score Strand
chr1 4 7 23.8 +
chr1 6 11 0.13 -
chr2 0 14 42.42 +"""


df1 = pd.read_csv(StringIO(f1), sep="\s+")

gr1 = PyRanges(df1)

#%%
f2 = """Chromosome Start End Score Strand
chr2 14 18 15.3 +
"""


df2 = pd.read_csv(StringIO(f2), sep="\s+")

gr2 = PyRanges(df2)

#%%
gr1.overlap(gr2)


#%%
f3 = """Chromosome Start End Score Strand
chr1 4 7 23.8 +
chr1 6 11 0.13 -
chr2 0 14 42.42 +
chr2 14 18 15.3 +"""


df3 = pd.read_csv(StringIO(f3), sep="\s+")

gr3 = PyRanges(df3)

#%%
gr3.merge(count=True)

#%%
#test genome_bounds function
f4 = """Chromosome Start End Score Strand
chr1 4 7 23.8 +
chr1 6 11 0.13 -
chr2 0 14 42.42 +
chr2 14 190000000 15.3 +"""


df4 = pd.read_csv(StringIO(f4), sep="\s+")

gr4 = PyRanges(df4)



#%%
import assembly

chrom_sizes = assembly.build('mm10', 1)._chromsizes

#chromsizes to pyranges
chroms = []
start = []
end = []
for k, v in chrom_sizes.items():
    chroms.append(k)
    start.append(0)
    end.append(v)

chromsize_gr = pr.PyRanges(chromosomes=chroms, starts=start, ends=end)

#%%
pr.gf.genome_bounds(gr4, chromsize_gr, clip=True)

#%%
