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
./assembler input.fastq 21 output/
```

---

## 📁 Output Files

| File | Contents |
|------|---------|
| `genome.fasta` | Assembled DNA sequence, 60-char wrapped |
| `stats.txt` | Read count, k-mer count, full complexity + timing breakdown |
| `graph_data.json` | Up to 500 nodes / 2000 edges for 3D visualization |

### stats.txt format
```
Total Reads Processed: 1000
Total k-mers: 95000
Graph Nodes (V): 4500
Graph Edges (E): 12300
Assembly Length: 5023 bp
k-mer size: 21

Time Complexity (Theoretical):
  - Rolling Hash:        O(N)       - one hash update per character
  - Bloom Filter:        O(N)       - constant-time insert/query per k-mer
  - Graph Construction:  O(V + E)   - adjacency list insertion
  - Traversal:           O(E)       - Hierholzer visits each edge once
  - DP Correction:       O(n*m)     - LCS over sliding windows
  Overall:               O(N + V + E)

Execution Time (Measured):
  Total Time:            142.3 ms
  Hashing Time:          98.5 ms
  Graph Build Time:      31.2 ms
  Traversal Time:        9.8 ms
  DP Time:               2.8 ms
```

---

## 🧩 DAA Module Details

### 1. Rolling Hash (Rabin-Karp) — `O(1)` per k-mer
Uses Mersenne prime modulus `2⁶¹ − 1` and `__uint128_t` multiplication to avoid overflow. The hash slides by removing the contribution of the leftmost base and adding the new rightmost base in constant time.

**Why it matters:** Without rolling hash, extracting N k-mers of length k costs O(N·k). Rolling hash reduces this to O(N) total.

### 2. Bloom Filter — `O(1)` insert/query
Three hash functions over a 64-million-bit array (~8 MB). Two-pass counting: first pass marks seen k-mers, second pass marks k-mers seen twice. Only k-mers present twice are added to the graph, eliminating most sequencing errors (which appear exactly once).

**Why it matters:** A hash map of all k-mers would require O(N·k) memory. The Bloom filter uses a fixed 8 MB regardless of N.

### 3. De Bruijn Graph — `O(V + E)`
Each unique (k-1)-mer becomes a node. Each k-mer becomes a directed edge from its prefix to its suffix. Built as an adjacency list (`unordered_map<string, vector<string>>`). Node imbalance (out_degree − in_degree) identifies Eulerian path endpoints.

### 4. Hierholzer's Algorithm — `O(E)`
Finds an Eulerian path (visiting every edge exactly once) using a stack-based DFS. Correct for graphs where at most two nodes have imbalanced degree. Falls back to greedy DFS for disconnected or complex graphs.

**Why it matters:** Hamiltonian path (visiting every node) is NP-hard. Eulerian path (visiting every edge) is linear. De Bruijn graphs convert genome assembly into the tractable problem.

### 5. DP Error Correction — `O(n·m)`
Slides a window over the assembled sequence and computes LCS between adjacent windows to detect inconsistency (LCS < 80% of window size signals a potential assembly artifact). Production systems would re-assemble or pull a corrected k-mer; this implementation preserves valid assemblies.

---

## 🌐 3D Visualization

- `networkx.DiGraph` builds the graph structure from `graph_data.json`
- `nx.spring_layout(G, dim=3)` computes 3D positions
- `plotly.graph_objects.Scatter3d` renders edges as thin lines and nodes as coloured spheres
- Node colour = degree (dark teal = high connectivity = repeat regions)
- Interactive: rotate, zoom, pan, hover for k-mer label

---

## ⚠️ Large File Handling

The assembler is designed for 100 MB+ files:
- FASTQ is read **line-by-line** (streaming, O(1) memory per read)
- Bloom filter uses **fixed 8 MB** regardless of input size
- Graph is built **on the fly** during the second pass
- JSON output is **capped at 500 nodes / 2000 edges** for visualization (full graph used for assembly)

---

## 📐 Complexity Summary

| Module | Theoretical | Notes |
|--------|------------|-------|
| FASTQ read | O(N) | N = total characters |
| Rolling hash | O(N) | O(1) per k-mer step |
| Bloom filter | O(N) | Constant-time hash |
| Graph build | O(V + E) | Adjacency list |
| Hierholzer | O(E) | Each edge visited once |
| DP correction | O(n·m) | n, m = window size (≤50) |
| **Overall** | **O(N + V + E)** | Linear in input + graph |
