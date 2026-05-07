/*
 * HelixForge — Genome Assembler
 * DAA Modules:
 *   Rolling Hash          — Rabin-Karp O(1) per k-mer slide
 *   Bloom Filter          — probabilistic solid-kmer filter, O(1) insert/query
 *   Exact Freq Counter    — two-pass unordered_map solid-kmer filter, O(N)
 *   De Bruijn Graph       — adjacency list with edge-frequency tracking
 *   Hierholzer            — Eulerian path traversal, O(E)
 *   Dijkstra              — coverage-weighted shortest path, O((V+E) log V)
 *   Multi-Start Greedy    — coverage-ranked greedy walk, handles fragmented graphs
 *   DP Error Correction   — sliding window majority vote, O(N)
 *   Suffix Array          — prefix-doubling O(n log² n)
 *   LCP Array             — Kasai's algorithm, O(n)
 *   Repeat Finder         — LCP-scan for repeated substrings
 *   Genome Metrics        — GC content, N50
 *
 * Compile:
 *   g++ -O2 -std=c++17 -o assembler assembler_standalone.cpp
 *   (Windows/MinGW AV lock: g++ -O2 -std=c++17 -static -o assembler assembler_standalone.cpp)
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <queue>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <chrono>
#include <filesystem>

using namespace std;
using namespace chrono;

/* =========================================================
   TIMING
   ========================================================= */
struct Timer {
    time_point<high_resolution_clock> t0;
    void start() { t0 = high_resolution_clock::now(); }
    double ms() const {
        return duration<double, milli>(high_resolution_clock::now() - t0).count();
    }
};

/* =========================================================
   REVERSE COMPLEMENT — canonical k-mer representation
   Used in Pass 1 (frequency counting) so k-mers that only
   appear on the reverse strand still meet the solid threshold.
   NOT used to add edges to the graph (would create symmetric
   back-edges → artificial full-graph Eulerian circuits).
   ========================================================= */
string rev_comp(const string& s) {
    string r; r.reserve(s.size());
    for (int i = (int)s.size()-1; i >= 0; --i) {
        char c = s[i];
        if      (c=='A') r+='T';
        else if (c=='T') r+='A';
        else if (c=='C') r+='G';
        else             r+='C';
    }
    return r;
}

/* =========================================================
   ROLLING HASH — Rabin-Karp  O(1) per k-mer
   Mersenne prime MOD = 2^61-1 avoids overflow via __uint128_t.
   ========================================================= */
class RollingHash {
    static constexpr uint64_t BASE = 5;
    static constexpr uint64_t MOD  = (1ULL << 61) - 1;
    uint64_t hv = 0, bp = 1;
    int k;

    static uint64_t cv(char c) {
        switch(c) {
            case 'A': return 1; case 'C': return 2;
            case 'G': return 3; case 'T': return 4;
        }
        return 1;
    }
    static uint64_t mm(uint64_t a, uint64_t b) {
        return (__uint128_t)a * b % MOD;
    }
public:
    explicit RollingHash(int k) : k(k) {}

    void init(const string& s, int pos) {
        hv = 0; bp = 1;
        for (int i = 0; i < k; ++i) {
            hv = (mm(hv, BASE) + cv(s[pos+i])) % MOD;
            if (i < k-1) bp = mm(bp, BASE);
        }
    }
    void slide(char old_c, char new_c) {
        hv = (hv + MOD - mm(bp, cv(old_c))) % MOD;
        hv = (mm(hv, BASE) + cv(new_c)) % MOD;
    }
    uint64_t value() const { return hv; }
};

/* =========================================================
   BLOOM FILTER — O(1) insert/query, ~8 MB
   Two-filter scheme: only add to graph if seen at least twice.
   Kept as a reference implementation. The pipeline currently
   uses an exact frequency counter (see Pass 1 in main) which
   provides threshold=3 needed for ONT sequencing error rates.
   ========================================================= */
class BloomFilter {
    static constexpr size_t BITS = 1ULL << 26;
    vector<uint8_t> b;
    int nh;

    size_t bit(size_t h, int i) const {
        uint64_t salt = (uint64_t)i * 2654435761ULL;
        return ((h ^ salt) * 6364136223846793005ULL
                          + 1442695040888963407ULL) % BITS;
    }
public:
    explicit BloomFilter(int hashes = 3) : b(BITS/8, 0), nh(hashes) {}
    void insert(size_t h) {
        for (int i = 0; i < nh; ++i) {
            size_t x = bit(h,i); b[x/8] |= (1u << (x%8));
        }
    }
    bool query(size_t h) const {
        for (int i = 0; i < nh; ++i) {
            size_t x = bit(h,i);
            if (!(b[x/8] & (1u << (x%8)))) return false;
        }
        return true;
    }
};

/* =========================================================
   DE BRUIJN GRAPH — adjacency list
   Node  = (k-1)-mer
   Edge  = k-mer
   edge_freq[u->v] = number of reads supporting that transition.
   ========================================================= */
struct Graph {
    unordered_map<string, vector<string>> adj;
    unordered_map<string, int> indeg, outdeg;
    unordered_map<string, int> edge_freq;     // "u->v" -> count

    void add_edge(const string& u, const string& v) {
        string key = u + "->" + v;
        if (edge_freq[key]++ == 0) {          // first time: register edge
            adj[u].push_back(v);
            outdeg[u]++;
            indeg[v]++;
            if (!indeg.count(u))  indeg[u]  = 0;
            if (!outdeg.count(v)) outdeg[v] = 0;
            if (!adj.count(v))    adj[v];
        }
    }

    // Eulerian source: outdeg - indeg == 1, else any node with out-edges
    string start_node() const {
        string best;
        for (auto& [n, od] : outdeg) {
            int id = indeg.count(n) ? indeg.at(n) : 0;
            if (od - id == 1)          return n;
            if (best.empty() && od > 0) best = n;
        }
        return best;
    }

    // Eulerian sink: indeg - outdeg == 1, else any node with no out-edges
    string sink_node() const {
        string best;
        for (auto& [n, id] : indeg) {
            int od = outdeg.count(n) ? outdeg.at(n) : 0;
            if (id - od == 1)                   return n;
            if (best.empty() && od == 0 && id > 0) best = n;
        }
        return best;
    }

    size_t V() const { return adj.size(); }
    size_t E() const {
        size_t e = 0;
        for (auto& [n, nb] : adj) e += nb.size();
        return e;
    }
};

/* =========================================================
   HIERHOLZER — Eulerian path  O(E)
   Visits every edge exactly once → longest possible assembly.
   ========================================================= */
vector<string> hierholzer(Graph& g) {
    string s = g.start_node();
    if (s.empty()) return {};

    unordered_map<string, int> idx;
    for (auto& [n, _] : g.adj) idx[n] = 0;

    vector<string> circuit;
    stack<string> stk;
    stk.push(s);

    while (!stk.empty()) {
        string v = stk.top();
        auto& nb = g.adj[v];
        if (idx[v] < (int)nb.size())
            stk.push(nb[idx[v]++]);
        else {
            circuit.push_back(v); stk.pop();
        }
    }
    reverse(circuit.begin(), circuit.end());
    return circuit;
}

/* =========================================================
   DIJKSTRA — Coverage-Weighted Assembly  O((V+E) log V)
   Edge weight  = (max_freq + 1) - freq(u,v)
                → high-coverage edges get low cost → Dijkstra
                  naturally selects the best-supported path.
   Returns the assembled sequence along the minimum-cost
   (= maximum-coverage) path from source to sink.
   ========================================================= */
string dijkstra_assemble(const Graph& g) {
    if (g.adj.empty()) return "";

    // Normalise weights: find global max frequency
    int max_f = 1;
    for (auto& [key, f] : g.edge_freq) max_f = max(max_f, f);

    string src  = g.start_node();
    string sink = g.sink_node();
    if (src.empty()) return "";

    // min-heap: (accumulated_cost, node)
    using Entry = pair<double, string>;
    priority_queue<Entry, vector<Entry>, greater<Entry>> pq;

    unordered_map<string, double> dist;
    unordered_map<string, string> prev;
    unordered_map<string, int>    hops;   // edge-hops from source

    for (auto& [n, _] : g.adj) dist[n] = 1e18;
    dist[src] = 0.0;
    hops[src] = 0;
    pq.push({0.0, src});

    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u] + 1e-9) continue;   // stale entry

        for (auto& v : g.adj.at(u)) {
            string key  = u + "->" + v;
            int    freq = g.edge_freq.count(key) ? g.edge_freq.at(key) : 1;
            double w    = double(max_f + 1 - freq);  // inverse: low cost = high freq
            double nd   = d + w;
            if (nd < dist[v] - 1e-9) {
                dist[v] = nd;
                prev[v] = u;
                hops[v] = hops.count(u) ? hops[u] + 1 : 1;
                pq.push({nd, v});
            }
        }
    }

    // Choose endpoint: prefer declared sink if reachable;
    // otherwise the deepest (most hops) reachable node.
    string end_node;
    if (!sink.empty() && dist.count(sink) && dist[sink] < 1e17) {
        end_node = sink;
    } else {
        int best_h = 0;
        for (auto& [n, h] : hops) {
            if (n == src) continue;
            if (dist[n] < 1e17 && h > best_h) { best_h = h; end_node = n; }
        }
    }
    if (end_node.empty()) return src;

    // Reconstruct path source → end_node via prev[] pointers
    vector<string> path;
    for (string cur = end_node; ; cur = prev[cur]) {
        path.push_back(cur);
        if (!prev.count(cur)) break;
    }
    reverse(path.begin(), path.end());

    // Stitch node labels into sequence (each subsequent node adds 1 char)
    string seq = path[0];
    for (size_t i = 1; i < path.size(); ++i) seq += path[i].back();
    return seq;
}

/* =========================================================
   MULTI-START GREEDY — coverage-weighted longest-path assembly

   Why this works when Hierholzer/Dijkstra give 1 k-mer:
   Real De Bruijn graphs from sequencing data are NOT Eulerian
   (many unbalanced nodes = dead-ends = tips). The graph is a
   forest of small chains, not a single connected Euler graph.

   Strategy:
   1. Sort all nodes by total outgoing edge frequency (most-covered
      nodes first — these are most likely to be in the real gene).
   2. For each candidate start, greedily follow the highest-frequency
      unused edge until stuck. Record sequence length.
   3. Return the longest sequence found across all starts.

   Using a GLOBAL used-edge set means independent runs each explore
   fresh territory — together they cover the whole graph.
   ========================================================= */
string multi_greedy_assemble(Graph& g, int top_n = 80) {
    if (g.adj.empty()) return "";

    // Score each node by sum of its outgoing edge frequencies
    vector<pair<int,string>> cands;
    cands.reserve(g.adj.size());
    for (auto& [n, neighbors] : g.adj) {
        int total_f = 0;
        for (auto& nb : neighbors) {
            string key = n + "->" + nb;
            total_f += g.edge_freq.count(key) ? g.edge_freq.at(key) : 1;
        }
        if (total_f > 0) cands.push_back({total_f, n});
    }
    sort(cands.rbegin(), cands.rend());   // highest coverage first

    string best;

    int tried = 0;
    for (auto& [tf, start] : cands) {
        if (tried++ >= top_n) break;

        string res = start;
        string cur = start;
        unordered_map<string,bool> used;

        while (true) {
            string best_nxt; int best_freq = -1;
            for (auto& nxt : g.adj[cur]) {
                string key = cur + "->" + nxt;
                if (!used[key]) {
                    int freq = g.edge_freq.count(key) ? g.edge_freq.at(key) : 1;
                    if (freq > best_freq) { best_freq = freq; best_nxt = nxt; }
                }
            }
            if (best_nxt.empty()) break;
            used[cur + "->" + best_nxt] = true;
            res += best_nxt.back();
            cur  = best_nxt;
        }

        if (res.size() > best.size()) best = res;
    }
    return best;
}

// Original single-start greedy kept for compatibility
string greedy_assemble(Graph& g) {
    return multi_greedy_assemble(g, 1);
}

/* =========================================================
   DP ERROR CORRECTION — sliding window majority vote  O(N)
   Fixes isolated single-base mismatches where both immediate
   neighbours agree (classic sequencing substitution error).
   ========================================================= */
string dp_correct(const string& seq) {
    string res = seq;
    for (int i = 1; i < (int)seq.size()-1; i++) {
        if (seq[i-1] == seq[i+1] && seq[i] != seq[i-1])
            res[i] = seq[i-1];
    }
    return res;
}

/* =========================================================
   GENOME METRICS
   ========================================================= */
double gc_content(const string& s) {
    if (s.empty()) return 0.0;
    int gc = 0;
    for (char c : s) if (c=='G' || c=='C') gc++;
    return 100.0 * gc / (double)s.size();
}

// N50: length L such that contigs >= L bp cover >= 50% of assembly.
// Splits on 'N' gap characters to handle multi-contig assemblies.
long long n50(const string& seq) {
    vector<long long> contigs;
    long long cur = 0;
    for (char c : seq) {
        if (c == 'N') { if (cur > 0) contigs.push_back(cur); cur = 0; }
        else            cur++;
    }
    if (cur > 0) contigs.push_back(cur);
    if (contigs.empty()) return 0;

    sort(contigs.rbegin(), contigs.rend());
    long long total = 0;
    for (auto l : contigs) total += l;
    long long half = total / 2, running = 0;
    for (auto l : contigs) { running += l; if (running >= half) return l; }
    return contigs.back();
}

/* =========================================================
   SUFFIX ARRAY — prefix-doubling  O(n log² n)

   Algorithm outline (Manber & Myers style):
   1. Initialise SA = [0..n-1], rank[i] = s[i].
   2. For gap = 1, 2, 4, ...:
      Sort SA by (rank[i], rank[i+gap]) — uses the O(gap)-correct
      ordering built in the previous round.
      Recompute rank[] from the new sorted order.
      Stop early once all ranks are distinct.
   ========================================================= */
vector<int> build_suffix_array(const string& s) {
    int n = (int)s.size();
    if (n == 0) return {};

    vector<int> sa(n), rank_(n), tmp(n);
    iota(sa.begin(), sa.end(), 0);
    for (int i = 0; i < n; i++) rank_[i] = (unsigned char)s[i];

    for (int gap = 1; gap < n; gap <<= 1) {
        // Comparator: primary key = rank_[a], secondary = rank_[a+gap] (-1 if OOB)
        auto cmp = [&](int a, int b) {
            if (rank_[a] != rank_[b]) return rank_[a] < rank_[b];
            int ra = (a + gap < n) ? rank_[a + gap] : -1;
            int rb = (b + gap < n) ? rank_[b + gap] : -1;
            return ra < rb;
        };
        sort(sa.begin(), sa.end(), cmp);

        // Rebuild ranks: sa[0] gets 0; sa[i] gets sa[i-1]'s rank + 1 if different
        tmp[sa[0]] = 0;
        for (int i = 1; i < n; i++)
            tmp[sa[i]] = tmp[sa[i-1]] + (cmp(sa[i-1], sa[i]) ? 1 : 0);
        rank_ = tmp;

        // Early exit: all ranks are now unique
        if (rank_[sa[n-1]] == n-1) break;
    }
    return sa;
}

/* =========================================================
   LCP ARRAY — Kasai's algorithm  O(n)

   Uses the inverse SA (rank array) to compute LCP in linear
   time by extending each match from the previous comparison.
   The key invariant: LCP drops by at most 1 per character,
   so h never decreases more than n times overall.
   ========================================================= */
vector<int> build_lcp(const string& s, const vector<int>& sa) {
    int n = (int)s.size();
    vector<int> rank_(n), lcp(n, 0);
    for (int i = 0; i < n; i++) rank_[sa[i]] = i;

    for (int i = 0, h = 0; i < n; i++) {
        if (rank_[i] > 0) {
            int j = sa[rank_[i] - 1];
            // Extend match character by character
            while (i+h < n && j+h < n && s[i+h] == s[j+h]) h++;
            lcp[rank_[i]] = h;
            if (h > 0) h--;    // Kasai's invariant: h decreases by at most 1
        }
    }
    return lcp;
}

/* =========================================================
   REPEAT FINDER — LCP array scan

   Two consecutive SA entries with LCP[i] >= min_len share a
   prefix of at least min_len characters — that prefix is a
   repeated substring.  We group maximal runs of high-LCP
   entries into "repeat regions", extract the shared pattern,
   and record all occurrence positions.
   ========================================================= */
struct RepeatRegion {
    string  pattern;      // the repeated substring
    int     length;       // length of shared prefix
    int     occurrences;  // how many times it appears
    vector<int> positions;// start positions in assembly
};

vector<RepeatRegion> find_repeats(const string& s,
                                   const vector<int>& sa,
                                   const vector<int>& lcp,
                                   int min_len = 20) {
    int n = (int)sa.size();
    vector<RepeatRegion> result;
    int orig_len = (int)s.size(); // s may have sentinel appended

    int i = 1;
    while (i < n) {
        if (lcp[i] >= min_len) {
            // Collect the maximal group of consecutive high-LCP entries
            int j         = i;
            int group_lcp = lcp[i];
            while (j < n && lcp[j] >= min_len) {
                group_lcp = min(group_lcp, lcp[j]);
                j++;
            }
            // Positions: sa[i-1], sa[i], ..., sa[j-1]
            vector<int> positions;
            auto valid = [&](int p) { return p + group_lcp <= orig_len; };
            if (valid(sa[i-1])) positions.push_back(sa[i-1]);
            for (int x = i; x < j; x++)
                if (valid(sa[x])) positions.push_back(sa[x]);

            if ((int)positions.size() >= 2) {
                string pat = s.substr(positions[0], group_lcp);
                result.push_back({pat, group_lcp,
                                  (int)positions.size(), positions});
            }
            i = j;
        } else {
            i++;
        }
    }

    // Sort: longest repeats first
    sort(result.begin(), result.end(),
         [](const RepeatRegion& a, const RepeatRegion& b) {
             return a.length > b.length;
         });
    return result;
}

/* =========================================================
   FASTQ STREAMING READER — O(N)
   Reads every 4th line (sequence line), normalises to ACGT.
   ========================================================= */
vector<string> read_fastq(const string& path, long long& nr) {
    ifstream f(path);
    if (!f) { cerr << "Cannot open " << path << "\n"; return {}; }
    vector<string> reads;
    string line; int ln = 0;
    while (getline(f, line)) {
        if (ln % 4 == 1) {
            string c; c.reserve(line.size());
            for (char ch : line) {
                char u = toupper(ch);
                c += (u=='A'||u=='C'||u=='G'||u=='T') ? u : 'A';
            }
            if (!c.empty()) reads.push_back(move(c));
        }
        ++ln;
    }
    nr = reads.size();
    return reads;
}

/* =========================================================
   JSON WRITER (manual, no external deps)
   Writes graph nodes, edges, and edge frequencies for
   the front-end visualiser.
   ========================================================= */
string escape_json(const string& s) {
    string o; o.reserve(s.size()+2);
    o += '"';
    for (char c : s) {
        if      (c == '"')  o += "\\\"";
        else if (c == '\\') o += "\\\\";
        else                o += c;
    }
    o += '"';
    return o;
}

// Smart graph export: selects structurally meaningful nodes instead of random ones.
// Priority: (1) assembly path nodes, (2) branch points, (3) path neighbours, (4) high-cov fill.
// Each node is tagged with a role; each edge is tagged on_path so the visualiser
// can highlight the assembly backbone in a distinct colour.
void write_graph_json(const string& path, const Graph& g,
                      const string& seq, int k) {
    const int MAX_NODES = 400, MAX_EDGES = 3000;

    // ── 1. Extract assembly path nodes in sequence order ──────────────────
    vector<string> path_ordered;
    unordered_set<string> path_set;
    if (!seq.empty() && (int)seq.size() >= k - 1) {
        for (int i = 0; i + k - 1 <= (int)seq.size(); ++i) {
            string nd = seq.substr(i, k - 1);
            if (g.adj.count(nd) && !path_set.count(nd)) {
                path_set.insert(nd);
                path_ordered.push_back(nd);
            }
        }
    }

    // ── 2. Identify branch nodes (combined degree > 2) ────────────────────
    unordered_set<string> branch_set;
    for (auto& [n, od] : g.outdeg) {
        int id = g.indeg.count(n) ? g.indeg.at(n) : 0;
        if (od + id > 2) branch_set.insert(n);
    }

    // ── 3. Build selected node map with roles ─────────────────────────────
    unordered_map<string, string> node_role;  // id → "path"|"branch"|"tip"|"other"
    unordered_map<string, int>    node_cov;

    auto add_node = [&](const string& n, const string& role) {
        if (node_role.count(n) || !g.adj.count(n)) return;
        node_role[n] = role;
        int cov = 0;
        for (auto& nb : g.adj.at(n)) {
            string key = n + "->" + nb;
            cov += g.edge_freq.count(key) ? g.edge_freq.at(key) : 1;
        }
        node_cov[n] = cov;
    };

    for (auto& n : path_ordered)  { if ((int)node_role.size() >= MAX_NODES) break; add_node(n, "path"); }
    for (auto& n : branch_set)    { if ((int)node_role.size() >= MAX_NODES) break; add_node(n, "branch"); }

    // Immediate neighbours of path nodes (shows branching off the backbone)
    for (auto& n : path_ordered) {
        if ((int)node_role.size() >= MAX_NODES) break;
        if (!g.adj.count(n)) continue;
        for (auto& nb : g.adj.at(n)) {
            if ((int)node_role.size() >= MAX_NODES) break;
            int od = g.outdeg.count(nb) ? g.outdeg.at(nb) : 0;
            add_node(nb, od == 0 ? "tip" : "other");
        }
    }

    // Fill remaining slots with highest-coverage nodes
    vector<pair<int,string>> cov_fill;
    for (auto& [n, _] : g.adj) {
        if (node_role.count(n)) continue;
        int cov = 0;
        for (auto& nb : g.adj.at(n)) {
            string key = n + "->" + nb;
            cov += g.edge_freq.count(key) ? g.edge_freq.at(key) : 1;
        }
        cov_fill.push_back({cov, n});
    }
    sort(cov_fill.rbegin(), cov_fill.rend());
    for (auto& [cov, n] : cov_fill) {
        if ((int)node_role.size() >= MAX_NODES) break;
        add_node(n, "other");
    }

    // ── 4. Collect edges between selected nodes ───────────────────────────
    // Build set of consecutive path-node pairs for on_path tagging
    unordered_set<string> path_edge_set;
    for (int i = 0; i + 1 < (int)path_ordered.size(); ++i)
        path_edge_set.insert(path_ordered[i] + "->" + path_ordered[i+1]);

    struct SelEdge { string u, v; int w; bool on_path; };
    vector<SelEdge> sel_edges;
    for (auto& [n, role] : node_role) {
        if (!g.adj.count(n)) continue;
        for (auto& nb : g.adj.at(n)) {
            if (!node_role.count(nb)) continue;
            if ((int)sel_edges.size() >= MAX_EDGES) goto done_collecting;
            string key = n + "->" + nb;
            int freq = g.edge_freq.count(key) ? g.edge_freq.at(key) : 1;
            sel_edges.push_back({n, nb, freq, path_edge_set.count(key) > 0});
        }
    }
    done_collecting:;

    // ── 5. Write JSON ─────────────────────────────────────────────────────
    ofstream f(path);
    f << "{\n  \"nodes\": [\n";
    bool first = true;
    for (auto& [n, role] : node_role) {
        if (!first) f << ",\n";
        f << "    {\"id\":"       << escape_json(n)
          << ",\"role\":"         << escape_json(role)
          << ",\"coverage\":"     << node_cov[n] << "}";
        first = false;
    }
    f << "\n  ],\n  \"edges\": [\n";
    first = true;
    for (auto& e : sel_edges) {
        if (!first) f << ",\n";
        f << "    {\"from\":"   << escape_json(e.u)
          << ",\"to\":"         << escape_json(e.v)
          << ",\"weight\":"     << e.w
          << ",\"on_path\":"    << (e.on_path ? "true" : "false") << "}";
        first = false;
    }
    f << "\n  ],\n";
    f << "  \"total_nodes\": "      << g.V()                 << ",\n";
    f << "  \"total_edges\": "      << g.E()                 << ",\n";
    f << "  \"path_node_count\": "  << path_ordered.size()   << "\n}\n";
}

/* =========================================================
   REPEAT REPORT WRITER
   ========================================================= */
void write_repeat_report(const string& path,
                          const vector<RepeatRegion>& repeats,
                          const string& seq,
                          int min_len) {
    ofstream f(path);
    f << "HelixForge — Repeat Analysis Report\n";
    f << string(60, '=') << "\n";
    f << "Assembly length : " << seq.size() << " bp\n";
    f << "Min repeat len  : " << min_len    << " bp\n";
    f << "Repeats found   : " << repeats.size() << "\n\n";

    int shown = min((int)repeats.size(), 25);
    if (shown == 0) { f << "No repeats found at this threshold.\n"; return; }

    f << "Top " << shown << " repeats (sorted by length):\n";
    f << string(60, '-') << "\n";
    f << left
      << setw(6)  << "Rank"
      << setw(10) << "Length"
      << setw(14) << "Occurrences"
      << "Positions (first 6)\n";
    f << string(60, '-') << "\n";

    for (int i = 0; i < shown; i++) {
        const auto& r = repeats[i];
        f << left
          << setw(6)  << (i+1)
          << setw(10) << r.length
          << setw(14) << r.occurrences;

        int pshow = min((int)r.positions.size(), 6);
        f << "[ ";
        for (int j = 0; j < pshow; j++) {
            f << r.positions[j];
            if (j+1 < pshow) f << ", ";
        }
        if ((int)r.positions.size() > 6) f << " ...";
        f << " ]\n";

        // Print up to 60 chars of the pattern
        f << "  Pattern: " << r.pattern.substr(0, 60);
        if (r.length > 60) f << "...";
        f << "\n\n";
    }
}

/* =========================================================
   MAIN — 7-stage pipeline
   ========================================================= */
int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << "Usage: assembler <input.fastq> <k> <output_dir>\n";
        return 1;
    }
    string fastq_path = argv[1];
    int    k          = stoi(argv[2]);
    string out_dir    = argv[3];

    if (k < 3 || k > 63) { cerr << "k must be in [3,63]\n"; return 1; }
    filesystem::create_directories(out_dir);

    Timer T; T.start();

    /* ── Stage 1 : Read FASTQ ─────────────────────────────── */
    cerr << "[1/7] Reading FASTQ...\n";
    long long nr = 0;
    auto reads = read_fastq(fastq_path, nr);
    if (reads.empty()) { cerr << "No reads.\n"; return 1; }
    cerr << "      " << nr << " reads loaded.\n";

    /* ── Stage 2 : Hash k-mers + build De Bruijn graph ────── */
    cerr << "[2/7] Hashing k-mers and building De Bruijn graph (k="
         << k << ")...\n";
    Timer ht; ht.start();
    RollingHash rh(k);

    // ── Pass 1: Exact k-mer frequency counting (forward + RC for accurate counts) ──
    // RC reads are counted here so k-mers from both strands meet the solid threshold,
    // but we do NOT add RC edges to the graph (that creates artificial cycles).
    unordered_map<uint64_t, int> kmer_freq;
    kmer_freq.reserve(8000000);
    for (auto& r : reads) {
        if ((int)r.size() < k) continue;
        // Forward
        rh.init(r, 0);
        kmer_freq[rh.value()]++;
        for (int i = 1; i+k <= (int)r.size(); ++i) {
            rh.slide(r[i-1], r[i+k-1]);
            kmer_freq[rh.value()]++;
        }
        // Reverse complement — improves solid-kmer detection on low-coverage strands
        string rc = rev_comp(r);
        rh.init(rc, 0);
        kmer_freq[rh.value()]++;
        for (int i = 1; i+k <= (int)rc.size(); ++i) {
            rh.slide(rc[i-1], rc[i+k-1]);
            kmer_freq[rh.value()]++;
        }
    }
    double hash_ms = ht.ms();

    // ── Pass 2: Build De Bruijn graph — forward edges only, solid k-mers only ──
    // Adding RC edges here would make the graph symmetric (every edge has a back-edge),
    // creating Eulerian circuits that traverse the entire graph instead of the target gene.
    // Forward-only edges form near-linear chains that assemble correctly.
    const int SOLID_MIN = 3;   // threshold = 3 filters ONT sequencing noise (~10-15% error)
    Timer gt; gt.start();
    Graph g;
    for (auto& r : reads) {
        if ((int)r.size() < k) continue;
        rh.init(r, 0);
        if (kmer_freq[rh.value()] >= SOLID_MIN) {
            string km = r.substr(0, k);
            g.add_edge(km.substr(0, k-1), km.substr(1, k-1));
        }
        for (int i = 1; i+k <= (int)r.size(); ++i) {
            rh.slide(r[i-1], r[i+k-1]);
            if (kmer_freq[rh.value()] >= SOLID_MIN) {
                string km = r.substr(i, k);
                g.add_edge(km.substr(0, k-1), km.substr(1, k-1));
            }
        }
    }
    double graph_ms = gt.ms();
    cerr << "      " << g.V() << " nodes, " << g.E() << " edges.\n";

    /* ── Stage 3 : Dijkstra coverage-weighted assembly ──────── */
    cerr << "[3/7] Dijkstra coverage-weighted assembly...\n";
    Timer djt; djt.start();
    string dijk_seq = dijkstra_assemble(g);
    double dijk_ms  = djt.ms();
    cerr << "      Dijkstra path  : " << dijk_seq.size() << " bp.\n";

    /* ── Stage 4 : Hierholzer Eulerian traversal ─────────────── */
    cerr << "[4/7] Hierholzer Eulerian traversal...\n";
    Timer tt; tt.start();
    auto   eul_path = hierholzer(g);
    string hier_seq;
    if (eul_path.size() > 1) {
        hier_seq = eul_path[0];
        for (size_t i = 1; i < eul_path.size(); ++i)
            hier_seq += eul_path[i].back();
    }
    double trav_ms = tt.ms();
    cerr << "      Hierholzer path: " << hier_seq.size() << " bp.\n";

    /* ── Multi-start greedy (always runs — handles fragmented graphs) ── */
    cerr << "      Running multi-start greedy (top-150 seed nodes)...\n";
    string greedy_seq = multi_greedy_assemble(g, 150);
    cerr << "      Greedy path    : " << greedy_seq.size() << " bp.\n";

    // Sanity-check Hierholzer: if its output is more than 10× the greedy output,
    // the graph has large cycles (non-Eulerian data) and Hierholzer traversed the
    // entire connected component instead of a single gene path — discard it.
    // A single barcode gene is typically ≤ 3000 bp; Dijkstra/greedy are better here.
    const size_t HIER_MAX_RATIO = 10;
    if (!greedy_seq.empty() && hier_seq.size() > greedy_seq.size() * HIER_MAX_RATIO) {
        cerr << "      Hierholzer discarded (full-graph circuit detected: "
             << hier_seq.size() << " bp >> greedy " << greedy_seq.size() << " bp).\n";
        hier_seq = "";
    }

    // Select the longest sane result across all three methods
    string seq;
    string method_used;
    struct Candidate { size_t len; string name; string seq; };
    vector<Candidate> results = {
        { hier_seq.size(),   "Hierholzer (Eulerian)",        hier_seq   },
        { dijk_seq.size(),   "Dijkstra (coverage-weighted)", dijk_seq   },
        { greedy_seq.size(), "Greedy (multi-start)",         greedy_seq },
    };
    auto best = max_element(results.begin(), results.end(),
        [](const Candidate& a, const Candidate& b){ return a.len < b.len; });
    seq         = best->seq;
    method_used = best->name;
    if (seq.empty()) { seq = greedy_seq; method_used = "Greedy (fallback)"; }

    cerr << "      Selected       : " << method_used
         << " → " << seq.size() << " bp.\n";

    /* ── Stage 5 : DP error correction ──────────────────────── */
    cerr << "[5/7] DP error correction...\n";
    Timer dt; dt.start();
    string seq_for_graph = seq;   // pre-correction: k-mers still exist in graph
    seq = dp_correct(seq);
    double dp_ms = dt.ms();

    /* ── Stage 6 : Suffix Array + LCP + repeat analysis ──────── */
    cerr << "[6/7] Building Suffix Array and LCP Array...\n";
    Timer sat; sat.start();

    // Minimum repeat length = k (repeats shorter than a k-mer are noise)
    const int REPEAT_MIN_LEN = max(k, 15);
    vector<int>          sa_arr, lcp_arr;
    vector<RepeatRegion> repeats;

    if (!seq.empty()) {
        // Append sentinel '$' (ASCII 36, smaller than A/C/G/T)
        // so no two suffixes are identical — required by SA algorithms.
        string sa_input = seq + "$";
        sa_arr  = build_suffix_array(sa_input);
        lcp_arr = build_lcp(sa_input, sa_arr);
        repeats = find_repeats(sa_input, sa_arr, lcp_arr, REPEAT_MIN_LEN);
    }
    double sa_ms = sat.ms();
    cerr << "      SA + LCP built. "
         << repeats.size() << " repeat region(s) found"
         << " (min_len=" << REPEAT_MIN_LEN << ").\n";

    /* ── Compute final metrics ───────────────────────────────── */
    double    gc      = gc_content(seq);
    long long n50_val = n50(seq);
    double    total_ms = T.ms();

    /* ── Stage 7 : Write outputs ─────────────────────────────── */
    cerr << "[7/7] Writing outputs...\n";

    // genome.fasta — 60-char wrapped FASTA
    {
        ofstream f(out_dir + "/genome.fasta");
        f << ">HelixForge_assembled_sequence method=" << method_used << "\n";
        for (size_t i = 0; i < seq.size(); i += 60)
            f << seq.substr(i, 60) << "\n";
    }

    // stats.txt — full pipeline summary
    {
        ofstream f(out_dir + "/stats.txt");
        f << "=== HelixForge Assembly Statistics ===\n\n"
          << "Input\n"
          << "  Reads processed   : " << nr          << "\n"
          << "  Unique k-mer hashes: " << kmer_freq.size() << "\n"
          << "  k-mer size (k)    : " << k           << "\n"
          << "  Solid threshold   : " << SOLID_MIN << "x\n\n"
          << "Graph\n"
          << "  Nodes (V)         : " << g.V()  << "\n"
          << "  Edges (E)         : " << g.E()  << "\n\n"
          << "Assembly\n"
          << "  Method selected   : " << method_used       << "\n"
          << "  Dijkstra length   : " << dijk_seq.size()   << " bp\n"
          << "  Hierholzer length : " << hier_seq.size()   << " bp\n"
          << "  Final length      : " << seq.size()        << " bp\n"
          << "  GC content        : " << fixed << setprecision(2)
                                      << gc            << " %\n"
          << "  N50               : " << n50_val       << " bp\n\n"
          << "Repeat Analysis (Suffix Array + LCP)\n"
          << "  Min repeat length : " << REPEAT_MIN_LEN    << " bp\n"
          << "  Repeat regions    : " << repeats.size()    << "\n";
        if (!repeats.empty()) {
            f << "  Longest repeat    : " << repeats[0].length
              << " bp (" << repeats[0].occurrences << "x)\n";
        }
        f << "\nTime Complexity (Theoretical)\n"
          << "  Rolling Hash      : O(N)           one slide per character\n"
          << "  Freq Counter      : O(N)           exact unordered_map count\n"
          << "  Graph Build       : O(V + E)       adjacency list\n"
          << "  Dijkstra          : O((V+E) log V) min-heap relaxation\n"
          << "  Hierholzer        : O(E)            each edge visited once\n"
          << "  DP Correction     : O(N)           single pass\n"
          << "  Suffix Array      : O(n log² n)    prefix doubling\n"
          << "  LCP Array         : O(n)           Kasai's algorithm\n"
          << "  Repeat Finding    : O(n)           LCP scan\n"
          << "  Overall           : O(N + (V+E) log V + n log² n)\n\n"
          << "Execution Time (Measured)\n"
          << "  Total             : " << total_ms << " ms\n"
          << "  Hashing           : " << hash_ms  << " ms\n"
          << "  Graph Build       : " << graph_ms << " ms\n"
          << "  Dijkstra          : " << dijk_ms  << " ms\n"
          << "  Hierholzer        : " << trav_ms  << " ms\n"
          << "  DP Correction     : " << dp_ms    << " ms\n"
          << "  Suffix Array+LCP  : " << sa_ms    << " ms\n";
    }

    // graph_data.json — nodes + weighted edges for visualiser
    write_graph_json(out_dir + "/graph_data.json", g, seq_for_graph, k);

    // repeats.txt — detailed repeat report
    write_repeat_report(out_dir + "/repeats.txt", repeats, seq, REPEAT_MIN_LEN);

    cerr << "Done in " << total_ms << " ms.\n";
    cerr << "Outputs: genome.fasta  stats.txt  graph_data.json  repeats.txt\n";
    return 0;
}
