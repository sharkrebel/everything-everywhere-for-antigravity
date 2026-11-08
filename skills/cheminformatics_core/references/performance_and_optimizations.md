# High-Performance RDKit & Optimization Best Practices

This guide details computational optimizations, C++ multithreading,
constant-memory streaming, and algorithmic speedups to scale 2D cheminformatics
workflows across millions of molecules.

--------------------------------------------------------------------------------

## 1. Constant-Memory Streaming for Multi-Gigabyte Files

Loading multi-gigabyte SDF files into a standard `SDMolSupplier` consumes
excessive RAM. `Chem.ForwardSDMolSupplier` streams molecules lazily with
near-zero memory footprint, supporting `.sdf` and `.sdf.gz` on the fly:

```python
import gzip
from rdkit import Chem

# Lazy streaming reader from gzip-compressed file (Memory: <50 MB regardless of file size)
with gzip.open("huge_library.sdf.gz", "rb") as gz_in:
    suppl = Chem.ForwardSDMolSupplier(gz_in, sanitize=True)
    for mol in suppl:
        if mol is None:
            continue
        process_molecule(mol)
```

For high-throughput multi-threaded parsing of uncompressed files on disk:

```python
# C++ multithreaded reader (spawns internal reader threads)
suppl = Chem.MultithreadedSDMolSupplier("large_library.sdf", numWriterThreads=4)
for mol in suppl:
    if mol is not None:
        process_molecule(mol)
```

--------------------------------------------------------------------------------

## 2. Two-Tier Substructure Search Optimization

Substructure graph isomorphism (VF2 algorithm) is NP-complete. Screening
millions of molecules directly with `mol.HasSubstructMatch(query)` is slow.
Using a **Two-Tier Pattern Fingerprint Prescreen** eliminates >95% of
non-matching molecules in vectorized C++ before graph matching:

```python
from rdkit import Chem, DataStructs

query = Chem.MolFromSmarts("c1ccccc1[#8]")
q_pat_fp = Chem.PatternFingerprint(query)

# 1. Fast binary bitmask check (microseconds per molecule)
matching_mols = []
for mol, target_pat_fp in zip(candidate_mols, target_pat_fps):
    if DataStructs.AllProbeBitsMatch(q_pat_fp, target_pat_fp):
        # 2. Exact VF2 graph isomorphism only on surviving candidates
        if mol.HasSubstructMatch(query):
            matching_mols.append(mol)
```

--------------------------------------------------------------------------------

## 3. Vectorized Bulk Similarity

Never use nested Python `for` loops to compare fingerprints. Use
`DataStructs.BulkTanimotoSimilarity`:

```python
from rdkit import DataStructs

# C++ SIMD accelerated bulk similarity (Millions of comparisons/sec)
similarities = DataStructs.BulkTanimotoSimilarity(query_fp, target_fps_list)
```

--------------------------------------------------------------------------------

## 4. 1D Lower-Triangular Distance Lists for Clustering

An $N \times N$ numpy matrix for $N = 50,000$ requires 20 GB of RAM.
`rdkit.ML.Cluster.Butina` takes a 1D lower-triangular list of length $N(N-1)/2$,
consuming a fraction of the memory:

```python
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

n = len(fps)
dist_1d = []
for i in range(1, n):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
    dist_1d.extend([1.0 - s for s in sims])

clusters = Butina.ClusterData(dist_1d, n, distThresh=0.4, isDistData=True)
```
