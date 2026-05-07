# 🧬 HelixForge — De-Novo Genome Assembler

A complete genome assembler combining a **C++ DAA engine** with a **Python/Streamlit frontend** and an **interactive 3D De Bruijn graph visualizer**.

---

## 🏗️ Architecture

<img width="1536" height="1024" alt="image" src="https://github.com/user-attachments/assets/89b597f8-f48d-4d7c-8682-ceb411318230" />


---

## 🚀 Quick Start

### Prerequisites
- **GCC 11+** (`g++ --version`)
- **Python 3.9+** (`python3 --version`)

### One-command launch
```bash
chmod +x run.sh && ./run.sh
```

This script:
1. Compiles `assembler_standalone.cpp` → `./assembler`
2. Installs Python requirements
3. Starts Streamlit at `http://localhost:8501`

### Manual steps
```bash
# Compile
g++ -O2 -std=c++17 -o assembler assembler_standalone.cpp

# Install Python deps
pip install -r requirements.txt

# Run frontend
streamlit run app.py

# Or run assembler directly
./assembler input.fastq 19 output/
```

---

## 📁 Output Files

| File | Contents |
|------|---------|
| `genome.fasta` | Assembled DNA sequence, 60-char wrapped |
| `stats.txt` | Read count, k-mer count, full complexity + timing breakdown |
| `graph_data.json` | Up to 500 nodes / 2000 edges for 3D visualization |
| `repeats.txt` | Repeat region analysis (suffix array + LCP) |

### stats.txt format
```
=== HelixForge Assembly Statistics ===

Input
  Reads processed   : 9707
  Unique k-mer hashes: 1973878
  k-mer size (k)    : 19
  Solid threshold   : 3x

Graph
  Nodes (V)         : 73833
  Edges (E)         : 75145

Assembly
  Method selected   : Hierholzer (Eulerian)
  Dijkstra length   : 43 bp
  Hierholzer length : 1143 bp
  Final length      : 1143 bp
  GC content        : 38.15 %
  N50               : 1143 bp

Repeat Analysis (Suffix Array + LCP)
  Min repeat length : 19 bp
  Repeat regions    : 0

Time Complexity (Theoretical)
  Rolling Hash      : O(N)           one slide per character
  Freq Counter      : O(N)           exact unordered_map count
  Graph Build       : O(V + E)       adjacency list
  Dijkstra          : O((V+E) log V) min-heap relaxation
  Hierholzer        : O(E)            each edge visited once
  DP Correction     : O(N)           single pass
  Suffix Array      : O(n log² n)    prefix doubling
  LCP Array         : O(n)           Kasai's algorithm
  Repeat Finding    : O(n)           LCP scan
  Overall           : O(N + (V+E) log V + n log² n)

Execution Time (Measured)
  Total             : 2408.22 ms
  Hashing           : 881.69 ms
  Graph Build       : 1249.92 ms
  Dijkstra          : 42.66 ms
  Hierholzer        : 34.64 ms
  DP Correction     : 0.00 ms
  Suffix Array+LCP  : 0.15 ms
```

---

## 🧩 DAA Module Details

### 1. Rolling Hash (Rabin-Karp) — `O(1)` per k-mer
Uses Mersenne prime modulus `2⁶¹ − 1` and `__uint128_t` multiplication to avoid overflow. The hash slides by removing the contribution of the leftmost base and adding the new rightmost base in constant time.

**Why it matters:** Without rolling hash, extracting N k-mers of length k costs O(N·k). Rolling hash reduces this to O(N) total.

### 2. Exact Frequency Counter — `O(N)` total
Two-pass approach using `unordered_map<uint64_t, int>`:

- **Pass 1:** Count every forward k-mer hash **and** every reverse-complement k-mer hash in the frequency map. This ensures k-mers that only appear on the reverse strand still meet the solid threshold.
- **Pass 2:** Build the De Bruijn graph using **forward-only** edges for solid k-mers (frequency ≥ 3). Reverse-complement edges are deliberately excluded — adding them makes every edge symmetric, which creates artificial Eulerian circuits spanning the entire graph instead of one gene.

**Solid threshold = 3×:** Oxford Nanopore sequencing has ~10–15% per-base error rate. A threshold of 2× still admits many error k-mers; 3× reliably removes them.

**Why not Bloom filter?** A Bloom filter can only detect whether a k-mer was seen "at least once" or "at least twice" — it cannot store exact counts needed to distinguish noise (1–2×) from real signal (3×+) on noisy ONT data.

### 3. De Bruijn Graph — `O(V + E)`
Each unique (k-1)-mer becomes a node. Each solid k-mer creates a **directed edge** from its prefix `[0..k-2]` to its suffix `[1..k-1]`. Built as an adjacency list with an edge frequency map tracking how many reads support each transition.

**Key insight:** If a genome has an Eulerian path through its De Bruijn graph, that path spells out the assembled sequence. Euler's theorem guarantees this when every node has balanced in/out degrees.

**No canonical k-mers in graph:** Taking `min(kmer, rev_comp(kmer))` randomly reverses edge direction — consecutive k-mers in a read no longer form a connected chain. Forward-only edges preserve read-order continuity.

### 4. Multi-Start Greedy Assembly — `O(E)`
Sorts all nodes by total outgoing edge frequency (highest-coverage first), then runs a greedy walk from each of the top-80 seed nodes, always following the highest-frequency unused edge at each step. Returns the longest sequence found across all starts.

**Why multi-start?** Real sequencing graphs are fragmented — not a single connected Eulerian graph. Single-start greedy gets stuck at the first dead end. Multi-start explores the whole graph, stitching together the longest possible chain.

### 5. Dijkstra's Coverage-Weighted Assembly — `O((V+E) log V)`
Assigns edge weight `= (max_freq + 1) − freq(u→v)` so high-coverage edges get low cost. Standard Dijkstra with a min-heap finds the path from source to sink that follows the most-supported transitions.

**Why it matters:** In a fragmented graph, the Dijkstra path may be shorter than Greedy but represents the highest-confidence route through the data.

### 6. Hierholzer's Eulerian Path — `O(E)`
Finds a path that visits every edge exactly once using a stack-based iterative DFS. When the De Bruijn graph is truly Eulerian (all nodes balanced), this produces the longest possible assembly in linear time.

**Sanity check:** If Hierholzer's output is more than 10× the Greedy output, it detected a full-graph circuit (not a gene path) and is discarded. The longest of the remaining methods is used instead.

**Why it matters:** Hamiltonian path (visiting every node) is NP-hard. Eulerian path (visiting every edge) is O(E). De Bruijn graphs convert genome assembly into the tractable problem.

### 7. DP Error Correction — `O(N)`
Single-pass correction: if a base disagrees with both its immediate neighbours (which agree with each other), it is corrected to match them. Fixes isolated substitution errors introduced by sequencing.

### 8. Suffix Array (Prefix-Doubling) — `O(n log² n)`
Builds a sorted array of all suffix starting positions using O(log n) doubling rounds, each sorting by the pair `(rank[i], rank[i+gap])`.

### 9. LCP Array (Kasai's Algorithm) — `O(n)`
Computes the Longest Common Prefix between adjacent suffixes in linear time using the invariant that LCP drops by at most 1 when moving from suffix i to i+1.

### 10. Repeat Finder — `O(n)`
Scans the LCP array for runs of values ≥ min_len. Each such run represents a group of suffixes sharing a long common prefix — a repeat region in the assembly.

---

## ⚠️ Large File Handling

The assembler is designed for large FASTQ files:
- FASTQ is read **line-by-line** (streaming, O(1) memory per read)
- Frequency counter pre-reserves 8M hash buckets to minimize rehashing
- Graph is built **on the fly** during the second pass
- JSON output is **capped at 500 nodes / 2000 edges** for visualization (full graph used for assembly)

---

## 📐 Complexity Summary

| Module | Theoretical | Notes |
|--------|------------|-------|
| FASTQ read | O(N) | N = total characters |
| Rolling hash | O(N) | O(1) per k-mer step |
| Freq counter | O(N) | Exact unordered_map count |
| Graph build | O(V + E) | Adjacency list, forward edges only |
| Multi-start greedy | O(E) | Top-80 seed nodes |
| Dijkstra | O((V+E) log V) | Min-heap, coverage-weighted |
| Hierholzer | O(E) | Each edge visited once |
| DP correction | O(N) | Single pass |
| Suffix array | O(n log² n) | Prefix doubling |
| LCP array | O(n) | Kasai's algorithm |
| Repeat finder | O(n) | LCP scan |
| **Overall** | **O(N + (V+E) log V + n log² n)** | Dominated by Dijkstra + SA |
