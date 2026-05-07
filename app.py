"""
HelixForge — Genome Assembler  |  Cyberpunk Biopunk UI
Run:  streamlit run app.py
"""

import json, subprocess, tempfile, shutil
from pathlib import Path

import streamlit as st
import networkx as nx
import plotly.graph_objects as go
import pandas as pd

st.set_page_config(
    page_title="HelixForge — Genome Assembler",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

BASE_DIR   = Path(__file__).parent
ASSEMBLER  = BASE_DIR / ("assembler.exe" if (BASE_DIR / "assembler.exe").exists() else "assembler")
OUTPUT_DIR = BASE_DIR / "output"

# ══════════════════════════════════════════════════════════════════════════════
#  GLOBAL CSS — Cyberpunk Biopunk theme
# ══════════════════════════════════════════════════════════════════════════════
STYLE = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Orbitron:wght@400;700;900&family=Space+Mono:ital,wght@0,400;0,700;1,400&family=Inter:wght@300;400;600;700&display=swap');

/* ── CSS Variables ── */
:root {
  --bg:     #030010;
  --bg2:    #08001e;
  --card:   rgba(12, 0, 40, 0.75);
  --cyan:   #00f0ff;
  --purple: #cc00ff;
  --pink:   #ff0088;
  --green:  #00ff99;
  --orange: #ff6d00;
  --yellow: #ffe600;
  --teal:   #00c8a0;
  --text:   #d0dcff;
  --muted:  rgba(160,180,255,0.5);
  --border: rgba(150,100,255,0.18);
}

/* ── Keyframes ── */
@keyframes gradshift {
  0%,100% { background-position: 0% 50%; }
  50%      { background-position: 100% 50%; }
}
@keyframes pulseC  { 0%,100%{opacity:.7} 50%{opacity:1} }
@keyframes scanbar { 0%{transform:translateY(-100%)} 100%{transform:translateY(600%)} }
@keyframes float   { 0%,100%{transform:translateY(0)} 50%{transform:translateY(-6px)} }
@keyframes blink   { 0%,100%{opacity:1} 50%{opacity:0} }
@keyframes glow-in { from{opacity:0;transform:translateY(12px)} to{opacity:1;transform:translateY(0)} }
@keyframes spin    { from{transform:rotate(0deg)} to{transform:rotate(360deg)} }

/* ── Base ── */
html, [class*="css"] { font-family:'Inter',sans-serif; }
.stApp {
  background: var(--bg);
  background-image:
    radial-gradient(ellipse 90% 60% at 15%  8%, rgba(0,240,255,.07) 0%, transparent 55%),
    radial-gradient(ellipse 70% 50% at 85% 90%, rgba(204,0,255,.08) 0%, transparent 55%),
    radial-gradient(ellipse 50% 40% at 50% 50%, rgba(255,0,136,.04) 0%, transparent 60%);
}
#MainMenu, footer, header { visibility:hidden; }
.block-container { padding:1.8rem 2.2rem 4rem; max-width:1500px; }

/* ── Sidebar ── */
[data-testid="stSidebar"] {
  background: linear-gradient(180deg, #06001a 0%, #030010 100%) !important;
  border-right: 1px solid rgba(204,0,255,.2) !important;
}
[data-testid="stSidebar"] * { color: var(--text) !important; }

/* ── Inputs ── */
[data-testid="stFileUploader"] {
  background: rgba(0,240,255,.03) !important;
  border: 1.5px dashed rgba(0,240,255,.3) !important;
  border-radius: 10px !important;
  transition: border-color .3s !important;
}
[data-testid="stFileUploader"]:hover { border-color: rgba(0,240,255,.6) !important; }

/* Slider track */
[data-testid="stSlider"] [class*="track"]  { background: rgba(204,0,255,.2) !important; }
[data-testid="stSlider"] [class*="track--filled"] { background: linear-gradient(90deg,var(--cyan),var(--purple)) !important; }
[data-testid="stSlider"] [class*="thumb"]  {
  background: var(--purple) !important;
  border: 2px solid var(--cyan) !important;
  box-shadow: 0 0 12px var(--cyan) !important;
}

/* ── Buttons ── */
.stButton > button {
  background: linear-gradient(135deg, #00f0ff 0%, #9900ff 60%, #ff0088 100%) !important;
  background-size: 200% 200% !important;
  animation: gradshift 3s ease infinite !important;
  color: #fff !important;
  font-family: 'Orbitron', monospace !important;
  font-weight: 700 !important;
  font-size: .82rem !important;
  letter-spacing: .15em !important;
  text-transform: uppercase !important;
  border: none !important;
  border-radius: 8px !important;
  padding: .7rem 1.5rem !important;
  box-shadow: 0 0 24px rgba(0,240,255,.4), 0 0 48px rgba(204,0,255,.2) !important;
  transition: all .3s !important;
}
.stButton > button:hover {
  transform: translateY(-2px) scale(1.02) !important;
  box-shadow: 0 0 40px rgba(0,240,255,.6), 0 0 80px rgba(204,0,255,.35) !important;
}

/* ── Download buttons ── */
[data-testid="stDownloadButton"] > button {
  background: transparent !important;
  color: var(--cyan) !important;
  border: 1px solid rgba(0,240,255,.35) !important;
  font-family: 'Space Mono', monospace !important;
  font-size: .72rem !important;
  border-radius: 6px !important;
  letter-spacing: .08em !important;
  transition: all .25s !important;
}
[data-testid="stDownloadButton"] > button:hover {
  background: rgba(0,240,255,.08) !important;
  border-color: var(--cyan) !important;
  box-shadow: 0 0 16px rgba(0,240,255,.3) !important;
}

/* ── Expanders ── */
[data-testid="stExpander"] {
  background: rgba(12,0,40,.6) !important;
  border: 1px solid var(--border) !important;
  border-radius: 8px !important;
}
[data-testid="stExpander"] summary {
  color: var(--cyan) !important;
  font-family: 'Space Mono', monospace !important;
  font-size: .78rem !important;
}
[data-testid="stExpander"] > div { background: rgba(8,0,28,.8) !important; }

/* ── Code blocks ── */
pre, code { background: #050018 !important; color: var(--green) !important; font-family:'Space Mono',monospace !important; font-size:.75rem !important; }
.stCodeBlock * { background: #050018 !important; color: var(--green) !important; }

/* ── DataFrames ── */
.stDataFrame { border: 1px solid var(--border) !important; border-radius: 8px !important; }
.stDataFrame thead th { background:#100030 !important; color:var(--cyan) !important; font-family:'Space Mono',monospace !important; font-size:.75rem !important; }
.stDataFrame tbody td { background:#070020 !important; color:var(--text) !important; font-family:'Space Mono',monospace !important; font-size:.78rem !important; border-color: var(--border) !important; }

/* ── Spinner ── */
.stSpinner > div { border-top-color: var(--cyan) !important; }

/* ─────────────────────────────────────────────
   CUSTOM COMPONENTS
───────────────────────────────────────────── */

/* Hero */
.hf-hero {
  padding: 2rem 0 1rem;
  animation: glow-in .6s ease both;
}
.hf-logo {
  font-family: 'Orbitron', monospace;
  font-weight: 900;
  font-size: 3.4rem;
  letter-spacing: -.01em;
  background: linear-gradient(120deg, #00f0ff 0%, #cc00ff 40%, #ff0088 70%, #ff6d00 100%);
  background-size: 300% 300%;
  animation: gradshift 4s ease infinite;
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
  background-clip: text;
  line-height: 1.1;
  text-shadow: none;
}
.hf-sub {
  font-family: 'Space Mono', monospace;
  font-size: .72rem;
  color: var(--muted);
  letter-spacing: .22em;
  text-transform: uppercase;
  margin-top: .5rem;
}
.hf-divider {
  height: 1px;
  background: linear-gradient(90deg, var(--cyan) 0%, var(--purple) 40%, var(--pink) 70%, transparent 100%);
  margin: 1.2rem 0;
  opacity: .5;
}

/* Section labels */
.hf-section {
  font-family: 'Orbitron', monospace;
  font-size: .65rem;
  font-weight: 700;
  letter-spacing: .28em;
  text-transform: uppercase;
  color: var(--cyan);
  margin: 2rem 0 1rem;
  display: flex;
  align-items: center;
  gap: .7rem;
}
.hf-section::after {
  content: '';
  flex: 1;
  height: 1px;
  background: linear-gradient(90deg, rgba(0,240,255,.4), rgba(204,0,255,.2), transparent);
}

/* KPI cards */
.kpi-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(140px,1fr)); gap:.85rem; margin:.5rem 0 1.2rem; }
.kpi-card {
  background: var(--card);
  border: 1px solid rgba(255,255,255,.07);
  border-radius: 10px;
  padding: 1rem .9rem .85rem;
  position: relative;
  overflow: hidden;
  transition: transform .2s, box-shadow .2s;
  backdrop-filter: blur(8px);
}
.kpi-card::before {
  content:'';
  position:absolute; top:0; left:0; right:0; height:2px;
  background: var(--accent, var(--cyan));
  box-shadow: 0 0 12px var(--accent, var(--cyan));
}
.kpi-card::after {
  content:'';
  position:absolute; inset:0;
  background: radial-gradient(ellipse 80% 60% at 50% 0%, rgba(255,255,255,.04), transparent 70%);
  pointer-events:none;
}
.kpi-card:hover {
  transform: translateY(-3px);
  box-shadow: 0 8px 32px rgba(0,0,0,.5), 0 0 20px var(--accent, var(--cyan)) 40;
}
.kpi-label { font-family:'Space Mono',monospace; font-size:.62rem; letter-spacing:.14em; text-transform:uppercase; color:var(--muted); margin-bottom:.35rem; }
.kpi-value { font-family:'Orbitron',monospace; font-weight:700; font-size:1.5rem; color:#fff; line-height:1.1; }
.kpi-unit  { font-family:'Space Mono',monospace; font-size:.62rem; color:var(--muted); margin-top:.1rem; }

/* Method badge */
.method-badge {
  display:inline-flex; align-items:center; gap:.6rem;
  background: rgba(0,240,255,.06);
  border: 1px solid rgba(0,240,255,.3);
  border-radius: 8px;
  padding: .55rem 1.1rem;
  margin: .6rem 0 1rem;
}
.method-label { font-family:'Space Mono',monospace; font-size:.68rem; color:var(--muted); letter-spacing:.1em; text-transform:uppercase; }
.method-value { font-family:'Orbitron',monospace; font-size:.85rem; font-weight:700; color:var(--cyan); }

/* Sequence terminal */
.seq-terminal { background:#040018; border:1px solid rgba(0,240,255,.2); border-radius:10px; overflow:hidden; }
.seq-bar { background:#0a0030; padding:.5rem 1rem; display:flex; align-items:center; gap:.5rem; border-bottom:1px solid rgba(0,240,255,.12); }
.seq-dot { width:10px; height:10px; border-radius:50%; }
.seq-body { font-family:'Space Mono',monospace; font-size:.73rem; line-height:1.9; color:#00ffcc; padding:1rem 1.1rem; max-height:185px; overflow-y:auto; word-break:break-all; white-space:pre-wrap; }
.seq-body::-webkit-scrollbar { width:3px; }
.seq-body::-webkit-scrollbar-thumb { background:rgba(0,240,255,.3); border-radius:2px; }

/* Pipeline steps */
.pipeline { display:grid; grid-template-columns:repeat(7,1fr); gap:0; margin:1.2rem 0 1.8rem; }
.pipe-step {
  background: rgba(12,0,40,.5);
  border: 1px solid var(--border);
  border-right: none;
  padding: 1.2rem .7rem 1rem;
  text-align: center;
  position: relative;
  transition: background .25s;
  backdrop-filter: blur(4px);
}
.pipe-step:first-child { border-radius:10px 0 0 10px; }
.pipe-step:last-child  { border-right:1px solid var(--border); border-radius:0 10px 10px 0; }
.pipe-step:hover { background:rgba(0,240,255,.06); }
.pipe-step::after {
  content:'›'; position:absolute; right:-10px; top:50%; transform:translateY(-50%);
  color:rgba(0,240,255,.35); font-size:1.4rem; z-index:2;
}
.pipe-step:last-child::after { display:none; }
.pipe-num   { font-family:'Space Mono',monospace; font-size:.58rem; color:var(--muted); letter-spacing:.1em; }
.pipe-icon  { font-size:1.4rem; margin:.25rem 0; line-height:1; }
.pipe-title { font-family:'Orbitron',monospace; font-weight:700; font-size:.68rem; color:#fff; margin-bottom:.15rem; }
.pipe-sub   { font-family:'Space Mono',monospace; font-size:.58rem; color:var(--muted); }

/* Problem statement */
.problem-card {
  background: linear-gradient(135deg, rgba(0,240,255,.04) 0%, rgba(204,0,255,.05) 50%, rgba(255,0,136,.04) 100%);
  border: 1px solid rgba(204,0,255,.25);
  border-radius: 14px;
  padding: 1.8rem 2rem;
  margin: 1.2rem 0 1.5rem;
  position: relative;
  overflow: hidden;
}
.problem-card::before {
  content:'';
  position:absolute; top:0; left:0; right:0; height:2px;
  background: linear-gradient(90deg, var(--cyan), var(--purple), var(--pink));
}
.problem-title { font-family:'Orbitron',monospace; font-weight:700; font-size:1.1rem; color:var(--cyan); margin-bottom:1rem; letter-spacing:.05em; }
.problem-grid  { display:grid; grid-template-columns:repeat(3,1fr); gap:1rem; margin-top:1rem; }
.prob-item     { background:rgba(255,255,255,.03); border:1px solid var(--border); border-radius:8px; padding:.9rem 1rem; }
.prob-item-title { font-family:'Orbitron',monospace; font-size:.65rem; font-weight:700; letter-spacing:.12em; text-transform:uppercase; margin-bottom:.4rem; }
.prob-item-text  { font-family:'Inter',sans-serif; font-size:.8rem; color:var(--text); line-height:1.6; }

/* Kmer preview */
.kmer-preview {
  background: rgba(4,0,20,.9);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: .7rem 1rem;
  font-family: 'Space Mono', monospace;
  font-size: .72rem;
  color: var(--muted);
  margin: .3rem 0 .8rem;
  overflow-x: auto;
}
.kmer-seq  { font-size: .78rem; letter-spacing: .06em; line-height:1.8; }
.kmer-high { color: var(--cyan); background: rgba(0,240,255,.12); border-radius:3px; padding:0 2px; font-weight:700; text-shadow:0 0 8px var(--cyan); }
.kmer-info { margin-top:.35rem; font-size:.66rem; display:flex; gap:1.2rem; flex-wrap:wrap; }
.kmer-info-item { color:var(--muted); }
.kmer-info-item span { color:var(--purple); font-weight:700; }

/* Repeat card */
.repeat-card {
  background: rgba(255,109,0,.04);
  border: 1px solid rgba(255,109,0,.25);
  border-radius: 10px;
  padding: 1rem 1.2rem;
  margin-top:.6rem;
}
.repeat-title { font-family:'Orbitron',monospace; font-size:.68rem; font-weight:700; color:var(--orange); letter-spacing:.1em; text-transform:uppercase; margin-bottom:.6rem; }
.repeat-zero { font-family:'Space Mono',monospace; font-size:.78rem; color:var(--muted); }

/* Info banner */
.info-banner {
  background: rgba(0,240,255,.05);
  border: 1px solid rgba(0,240,255,.2);
  border-left: 3px solid var(--cyan);
  border-radius: 0 8px 8px 0;
  padding: .9rem 1.2rem;
  font-family: 'Space Mono', monospace;
  font-size: .8rem;
  color: var(--text);
  margin: .8rem 0;
}

/* Sidebar logo */
.sb-logo {
  font-family: 'Orbitron', monospace;
  font-weight: 900;
  font-size: 1.25rem;
  background: linear-gradient(120deg, var(--cyan), var(--purple), var(--pink));
  background-size: 200% 200%;
  animation: gradshift 4s ease infinite;
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
  background-clip: text;
  padding: .4rem 0 1rem;
}
.sb-section {
  font-family: 'Orbitron', monospace;
  font-size: .58rem;
  font-weight: 700;
  letter-spacing: .22em;
  text-transform: uppercase;
  color: var(--muted) !important;
  margin: 1rem 0 .4rem;
  padding-bottom: .3rem;
  border-bottom: 1px solid var(--border);
}

/* Overall badge */
.overall-badge {
  display:inline-flex; align-items:center; gap:.6rem;
  background: rgba(204,0,255,.06);
  border: 1px solid rgba(204,0,255,.3);
  border-radius: 8px;
  padding: .55rem 1.1rem;
  margin-top: .8rem;
}
.overall-label { font-family:'Space Mono',monospace; font-size:.68rem; color:var(--muted); letter-spacing:.1em; text-transform:uppercase; }
.overall-value { font-family:'Orbitron',monospace; font-size:.85rem; font-weight:700; color:var(--purple); }

/* Comparison row */
.cmp-row { display:flex; gap:1rem; margin:.5rem 0 1rem; }
.cmp-box { flex:1; background:rgba(255,255,255,.03); border:1px solid var(--border); border-radius:8px; padding:.7rem 1rem; }
.cmp-box-label { font-family:'Space Mono',monospace; font-size:.62rem; color:var(--muted); text-transform:uppercase; letter-spacing:.1em; }
.cmp-box-value { font-family:'Orbitron',monospace; font-size:1rem; font-weight:700; margin-top:.2rem; }

/* Tabs */
[data-testid="stTabs"] [role="tablist"] { gap:.4rem; border-bottom:1px solid var(--border); }
[data-testid="stTabs"] [role="tab"] {
  font-family:'Orbitron',monospace !important;
  font-size:.65rem !important;
  font-weight:700 !important;
  letter-spacing:.12em !important;
  color:var(--muted) !important;
  border-radius:6px 6px 0 0 !important;
  border: 1px solid transparent !important;
  border-bottom: none !important;
  padding:.5rem 1.1rem !important;
  transition: all .2s !important;
}
[data-testid="stTabs"] [role="tab"][aria-selected="true"] {
  color:var(--cyan) !important;
  border-color: var(--border) !important;
  background: rgba(0,240,255,.05) !important;
  border-bottom-color: var(--bg) !important;
}
</style>
"""
st.markdown(STYLE, unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
#  HELPERS
# ══════════════════════════════════════════════════════════════════════════════
def parse_stats(path: Path) -> dict:
    """Parse the new multi-section stats.txt format."""
    text  = path.read_text(errors="replace")
    lines = text.splitlines()
    data  = {
        "reads":"—","total_kmers":"—","kmer_size":"—",
        "nodes":"—","edges":"—",
        "method":"—","dijk_len":"—","hier_len":"—",
        "final_len":"—","gc":"—","n50":"—",
        "min_repeat":"—","repeat_regions":"0",
        "theoretical":[],"measured":{},"overall":"O(N + V + E)",
    }
    in_theo = False
    in_meas = False

    for line in lines:
        s = line.strip()
        if not s or s.startswith("="):
            continue
        if "Time Complexity" in s:
            in_theo = True;  in_meas = False;  continue
        if "Execution Time"  in s:
            in_theo = False; in_meas = True;   continue
        if ":" not in s:
            continue

        raw_key, _, raw_val = s.partition(":")
        k_str = raw_key.strip()
        v_str = raw_val.strip()
        kl    = k_str.lower()

        if in_theo:
            if "overall" in kl:
                data["overall"] = v_str
            elif v_str:
                tokens     = v_str.split()
                complexity = tokens[0] if tokens else ""
                desc_tail  = v_str[len(complexity):].strip().lstrip("-").strip()
                data["theoretical"].append({"step":k_str,"complexity":complexity,"reason":desc_tail})
        elif in_meas:
            data["measured"][k_str] = v_str
        else:
            if   "reads processed"   in kl: data["reads"]          = v_str
            elif "total k-mers"      in kl: data["total_kmers"]    = v_str
            elif "k-mer size"        in kl: data["kmer_size"]      = v_str
            elif "nodes (v)"         in kl: data["nodes"]          = v_str
            elif "edges (e)"         in kl: data["edges"]          = v_str
            elif "method selected"   in kl: data["method"]         = v_str
            elif "dijkstra length"   in kl: data["dijk_len"]       = v_str
            elif "hierholzer length" in kl: data["hier_len"]       = v_str
            elif "final length"      in kl: data["final_len"]      = v_str.replace(" bp","").strip()
            elif "gc content"        in kl: data["gc"]             = v_str.replace(" %","").strip()
            elif "n50"               in kl: data["n50"]            = v_str.replace(" bp","").strip()
            elif "min repeat length" in kl: data["min_repeat"]     = v_str
            elif "repeat regions"    in kl: data["repeat_regions"] = v_str
    return data


def parse_repeats(path: Path) -> list[dict]:
    """Parse repeats.txt — returns list of {rank, length, occ, positions, pattern}."""
    if not path.exists():
        return []
    rows = []
    current = {}
    for line in path.read_text(errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith("=") or s.startswith("-") or s.startswith("HelixForge") \
                or s.startswith("Assembly") or s.startswith("Min") or s.startswith("Repeats") \
                or s.startswith("Top") or s.startswith("Rank"):
            continue
        if s[0].isdigit() and len(s.split()) >= 2:
            if current:
                rows.append(current)
            parts = s.split()
            current = {"rank":parts[0],"length":parts[1],"occurrences":parts[2] if len(parts)>2 else "—","positions":"","pattern":""}
            pos_start = s.find("[")
            if pos_start != -1:
                current["positions"] = s[pos_start:]
        elif s.startswith("Pattern:"):
            current["pattern"] = s.replace("Pattern:","").strip()
    if current:
        rows.append(current)
    return rows


def build_3d_graph(graph_path: Path, max_nodes: int):
    gdata = json.loads(graph_path.read_text())
    nodes = gdata["nodes"][:max_nodes]
    node_set = set(nodes)

    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    edge_weights = {}
    for edge in gdata["edges"]:
        # Handle both old [u,v] and new {"from":u,"to":v,"weight":w} formats
        if isinstance(edge, list):
            u, v, w = edge[0], edge[1], 1
        else:
            u = edge.get("from",""); v = edge.get("to",""); w = edge.get("weight",1)
        if u in node_set and v in node_set:
            G.add_edge(u, v)
            edge_weights[(u, v)] = w

    pos  = nx.spring_layout(G, dim=3, seed=42, k=0.55)
    degs = [G.degree(n) for n in nodes]

    # Build edge traces with weight-based coloring
    max_w = max(edge_weights.values(), default=1)
    ex, ey, ez, ec = [], [], [], []
    for (u, v), w in edge_weights.items():
        if u not in pos or v not in pos:
            continue
        x0,y0,z0 = pos[u]; x1,y1,z1 = pos[v]
        ex += [x0,x1,None]; ey += [y0,y1,None]; ez += [z0,z1,None]

    edge_trace = go.Scatter3d(
        x=ex, y=ey, z=ez, mode="lines",
        line=dict(
            color=[w/max_w for (u,v),w in edge_weights.items() for _ in range(3)],
            colorscale=[[0,"rgba(0,30,50,.25)"],[0.4,"rgba(0,240,255,.5)"],[1,"rgba(204,0,255,.9)"]],
            width=1.5),
        hoverinfo="none", name="edges")

    node_trace = go.Scatter3d(
        x=[pos[n][0] for n in nodes],
        y=[pos[n][1] for n in nodes],
        z=[pos[n][2] for n in nodes],
        mode="markers",
        marker=dict(
            size=[3.5 + d*.8 for d in degs],
            color=degs,
            colorscale=[
                [0,   "#0a0030"],
                [0.2, "#00c8a0"],
                [0.5, "#00f0ff"],
                [0.8, "#cc00ff"],
                [1,   "#ff0088"]],
            colorbar=dict(
                title=dict(text="Degree", font=dict(color="#8090cc", size=9)),
                thickness=8, len=.45,
                tickfont=dict(color="#8090cc", size=8)),
            opacity=.92, line=dict(width=0)),
        text=nodes,
        hovertemplate="<b>%{text}</b><br>degree: %{marker.color}<extra></extra>",
        name="nodes")

    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        scene=dict(
            xaxis=dict(showgrid=False,zeroline=False,showticklabels=False,showbackground=False),
            yaxis=dict(showgrid=False,zeroline=False,showticklabels=False,showbackground=False),
            zaxis=dict(showgrid=False,zeroline=False,showticklabels=False,showbackground=False),
            bgcolor="rgba(0,0,0,0)"),
        legend=dict(x=.01, y=.99, font=dict(size=9,color="#8090cc"), bgcolor="rgba(0,0,0,0)"),
        margin=dict(l=0,r=0,t=0,b=0), height=500,
        paper_bgcolor="rgba(0,0,0,0)",
        hoverlabel=dict(bgcolor="#10003a", font_size=11, font_family="Space Mono"))
    return fig


def kpi_html(label, value, unit="", accent="var(--cyan)"):
    return (f'<div class="kpi-card" style="--accent:{accent}">'
            f'<div class="kpi-label">{label}</div>'
            f'<div class="kpi-value">{value}</div>'
            f'<div class="kpi-unit">{unit}</div></div>')


def kmer_preview_html(k: int) -> str:
    dna    = "ATCGATCGATCGATCGATCGATCG"
    start  = 3
    end    = min(start + k, len(dna))
    left   = dna[:start]
    middle = dna[start:end]
    right  = dna[end:]
    return f"""
    <div class="kmer-preview">
      <div style="font-size:.6rem;color:var(--muted);text-transform:uppercase;letter-spacing:.12em;margin-bottom:.35rem">
        Live k-mer preview (k={k})
      </div>
      <div class="kmer-seq"
           style="color:rgba(160,180,255,.5)">{left}<span class="kmer-high">{middle}</span>{right}</div>
      <div class="kmer-info">
        <div class="kmer-info-item">k-mer size: <span>{k} bp</span></div>
        <div class="kmer-info-item">node size: <span>{k-1} bp</span></div>
        <div class="kmer-info-item">sliding window → <span>O(1) / step</span></div>
      </div>
    </div>"""


# ══════════════════════════════════════════════════════════════════════════════
#  SIDEBAR
# ══════════════════════════════════════════════════════════════════════════════
with st.sidebar:
    st.markdown('<div class="sb-logo">⬡ HelixForge</div>', unsafe_allow_html=True)

    st.markdown('<div class="sb-section">Input File</div>', unsafe_allow_html=True)
    uploaded = st.file_uploader(
        "FASTQ file", type=["fastq","fq","txt"],
        help="Standard 4-line FASTQ. Supports large files.",
        label_visibility="collapsed")
    if not uploaded:
        st.markdown(
            '<div style="font-family:Space Mono,monospace;font-size:.68rem;'
            'color:rgba(160,180,255,.4);margin-top:.2rem">Drop .fastq / .fq file here</div>',
            unsafe_allow_html=True)

    st.markdown('<div class="sb-section">k-mer Parameters</div>', unsafe_allow_html=True)
    k = st.slider("k-mer size (k)", min_value=7, max_value=63, value=21, step=2,
                  help="Odd values avoid palindromes. Larger k = more specific but needs more coverage.")
    st.markdown(kmer_preview_html(k), unsafe_allow_html=True)

    max_vis = st.slider("Max graph nodes (visualiser)", 50, 600, 200, 50,
                        help="Capped for performance. Full graph saved to graph_data.json.")

    st.markdown("")
    run_btn = st.button("▶  RUN ASSEMBLY", type="primary", use_container_width=True)

    if uploaded:
        fsize = len(uploaded.getvalue()) / 1024
        st.markdown(
            f'<div style="font-family:Space Mono,monospace;font-size:.66rem;'
            f'color:rgba(0,240,255,.5);margin-top:.3rem;text-align:center">'
            f'📁 {uploaded.name} · {fsize:.1f} KB</div>',
            unsafe_allow_html=True)

    st.markdown('<div class="sb-section">Algorithm Notes</div>', unsafe_allow_html=True)
    with st.expander("Why De Bruijn graphs?"):
        st.write("Hamiltonian path (find a path visiting every **node** once) is NP-hard. "
                 "De Bruijn graphs reframe the problem as Eulerian path (visit every **edge** once) "
                 "which Hierholzer solves in O(E).")
    with st.expander("Rolling hash O(1)"):
        st.write("Rabin-Karp polynomial hash slides in O(1): remove the leftmost character, "
                 "add the new rightmost. Without this, hashing every k-mer is O(k·N) total.")
    with st.expander("Bloom filter memory"):
        st.write("Fixed 8 MB bit-array vs O(N·k) exact hash set. Two-filter scheme: "
                 "insert to 'once' then 'twice' — only k-mers in 'twice' are solid (real genome sequence).")
    with st.expander("Dijkstra + coverage"):
        st.write("Edge weight = (max_freq + 1) − observed_freq. "
                 "High coverage → low weight → Dijkstra finds the most-supported assembly path.")
    with st.expander("Suffix Array O(n log²n)"):
        st.write("Prefix-doubling (Manber & Myers): each round doubles the sorted prefix length. "
                 "O(log n) rounds × O(n log n) sort. Kasai's LCP is then O(n) by amortised analysis.")


# ══════════════════════════════════════════════════════════════════════════════
#  HERO
# ══════════════════════════════════════════════════════════════════════════════
st.markdown("""
<div class="hf-hero">
  <div class="hf-logo">HelixForge</div>
  <div class="hf-sub">
    De-Bruijn · Rolling Hash · Bloom Filter · Dijkstra · Hierholzer · Suffix Array · DP Correction
  </div>
</div>
<div class="hf-divider"></div>
""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
#  ASSEMBLE (on button click)
# ══════════════════════════════════════════════════════════════════════════════
if run_btn:
    if uploaded is None:
        st.error("⚠ Please upload a FASTQ file first."); st.stop()
    if not ASSEMBLER.exists():
        st.error(f"Assembler binary not found at `{ASSEMBLER}`.\n"
                 "Compile with:\n```\ng++ -O2 -std=c++17 -static -o assembler assembler_standalone.cpp\n```")
        st.stop()

    tmp = Path(tempfile.mkdtemp())
    fq  = tmp / "input.fastq"
    fq.write_bytes(uploaded.getvalue())
    out = tmp / "out"; out.mkdir()

    with st.spinner("🔬  Running C++ assembly engine — 7-stage pipeline …"):
        res = subprocess.run(
            [str(ASSEMBLER), str(fq), str(k), str(out)],
            capture_output=True, text=True, timeout=600)

    if res.returncode != 0:
        st.error("Assembly engine returned an error:")
        st.code(res.stderr, language="text")
        shutil.rmtree(tmp, ignore_errors=True); st.stop()

    OUTPUT_DIR.mkdir(exist_ok=True)
    for fname in ["genome.fasta","stats.txt","graph_data.json","repeats.txt"]:
        src = out / fname
        if src.exists(): shutil.copy(src, OUTPUT_DIR / fname)

    st.session_state["assembled"] = True
    shutil.rmtree(tmp, ignore_errors=True)
    st.rerun()


# ══════════════════════════════════════════════════════════════════════════════
#  TABS
# ══════════════════════════════════════════════════════════════════════════════
tab_results, tab_problem, tab_complexity = st.tabs(
    ["  ASSEMBLY RESULTS  ", "  PROBLEM STATEMENT  ", "  COMPLEXITY & ALGORITHMS  "])


# ─────────────────────────────────────────────────────────────────────────────
#  TAB 1 — RESULTS
# ─────────────────────────────────────────────────────────────────────────────
with tab_results:
    if not st.session_state.get("assembled"):
        # Landing state
        st.markdown(
            '<div class="info-banner">Upload a <code>.fastq</code> file in the sidebar, '
            'tune the k-mer slider, then click <strong>▶ RUN ASSEMBLY</strong>.</div>',
            unsafe_allow_html=True)

        st.markdown('<div class="hf-section">7-Stage Pipeline</div>', unsafe_allow_html=True)
        st.markdown("""
        <div class="pipeline">
          <div class="pipe-step"><div class="pipe-num">01</div><div class="pipe-icon">📂</div>
            <div class="pipe-title">FASTQ</div><div class="pipe-sub">Stream · O(N)</div></div>
          <div class="pipe-step"><div class="pipe-num">02</div><div class="pipe-icon">⚡</div>
            <div class="pipe-title">Hash</div><div class="pipe-sub">Rabin-Karp · O(1)</div></div>
          <div class="pipe-step"><div class="pipe-num">03</div><div class="pipe-icon">🌸</div>
            <div class="pipe-title">Bloom</div><div class="pipe-sub">8 MB · O(N)</div></div>
          <div class="pipe-step"><div class="pipe-num">04</div><div class="pipe-icon">🧬</div>
            <div class="pipe-title">De Bruijn</div><div class="pipe-sub">O(V+E)</div></div>
          <div class="pipe-step"><div class="pipe-num">05</div><div class="pipe-icon">🎯</div>
            <div class="pipe-title">Dijkstra</div><div class="pipe-sub">O((V+E)logV)</div></div>
          <div class="pipe-step"><div class="pipe-num">06</div><div class="pipe-icon">🛤</div>
            <div class="pipe-title">Hierholzer</div><div class="pipe-sub">Euler · O(E)</div></div>
          <div class="pipe-step"><div class="pipe-num">07</div><div class="pipe-icon">🔬</div>
            <div class="pipe-title">SA + LCP</div><div class="pipe-sub">O(n log²n)</div></div>
        </div>
        """, unsafe_allow_html=True)
        st.stop()

    # ── Load outputs ──────────────────────────────────────────────────────────
    fasta_p   = OUTPUT_DIR / "genome.fasta"
    stats_p   = OUTPUT_DIR / "stats.txt"
    graph_p   = OUTPUT_DIR / "graph_data.json"
    repeats_p = OUTPUT_DIR / "repeats.txt"

    if not stats_p.exists():
        st.error("Output files not found. Re-run the assembly."); st.stop()

    stats   = parse_stats(stats_p)
    repeats = parse_repeats(repeats_p)
    seq_txt = fasta_p.read_text() if fasta_p.exists() else ""
    sequence = "".join(l for l in seq_txt.splitlines() if not l.startswith(">"))

    # ── KPI grid ─────────────────────────────────────────────────────────────
    st.markdown('<div class="hf-section">Assembly Metrics</div>', unsafe_allow_html=True)
    kpis = [
        kpi_html("Reads",     stats["reads"],        "reads",  "var(--cyan)"),
        kpi_html("k-mers",    stats["total_kmers"],  "k-mers", "var(--purple)"),
        kpi_html("Nodes (V)", stats["nodes"],        "(k-1)-mers", "var(--pink)"),
        kpi_html("Edges (E)", stats["edges"],        "k-mers", "var(--green)"),
        kpi_html("Assembly",  stats["final_len"],    "bp",     "var(--orange)"),
        kpi_html("GC Content",stats["gc"],           "%",      "var(--yellow)"),
        kpi_html("N50",       stats["n50"],          "bp",     "var(--teal)"),
    ]
    st.markdown(f'<div class="kpi-grid">{"".join(kpis)}</div>', unsafe_allow_html=True)

    # ── Method + Dijkstra vs Hierholzer comparison ────────────────────────────
    method_color = "var(--cyan)" if "Hierholzer" in stats["method"] else "var(--purple)"
    st.markdown(
        f'<div class="method-badge">'
        f'<span class="method-label">Method selected</span>'
        f'<span class="method-value" style="color:{method_color}">⬡ {stats["method"]}</span>'
        f'</div>', unsafe_allow_html=True)

    if stats["dijk_len"] != "—" or stats["hier_len"] != "—":
        st.markdown(
            f'<div class="cmp-row">'
            f'<div class="cmp-box"><div class="cmp-box-label">Dijkstra (coverage-weighted)</div>'
            f'<div class="cmp-box-value" style="color:var(--purple)">{stats["dijk_len"]}</div></div>'
            f'<div class="cmp-box"><div class="cmp-box-label">Hierholzer (Eulerian path)</div>'
            f'<div class="cmp-box-value" style="color:var(--cyan)">{stats["hier_len"]}</div></div>'
            f'</div>', unsafe_allow_html=True)

    st.markdown('<div class="hf-divider"></div>', unsafe_allow_html=True)

    # ── Main columns ──────────────────────────────────────────────────────────
    col_left, col_right = st.columns([1, 1.15], gap="large")

    with col_left:
        # Sequence terminal
        st.markdown('<div class="hf-section">Assembled Sequence</div>', unsafe_allow_html=True)
        preview = sequence[:3600] + ("…" if len(sequence) > 3600 else "")
        st.markdown(f"""
        <div class="seq-terminal">
          <div class="seq-bar">
            <div class="seq-dot" style="background:#ff5f57"></div>
            <div class="seq-dot" style="background:#febc2e"></div>
            <div class="seq-dot" style="background:#28c840"></div>
            <span style="font-family:Space Mono,monospace;font-size:.68rem;
                         color:rgba(0,240,255,.5);margin-left:.5rem">
              genome.fasta — {len(sequence):,} bp
            </span>
          </div>
          <div class="seq-body">{preview}</div>
        </div>""", unsafe_allow_html=True)

        dl1, dl2 = st.columns(2)
        with dl1:
            st.download_button("⬇ genome.fasta", seq_txt,
                               "genome.fasta","text/plain", use_container_width=True)
        with dl2:
            st.download_button("⬇ stats.txt", stats_p.read_text(),
                               "stats.txt","text/plain", use_container_width=True)

        # Timing chart — all 7 stages
        st.markdown('<div class="hf-section">Execution Timing</div>', unsafe_allow_html=True)
        m = stats.get("measured", {})
        timing_stages = [
            ("Hashing",      m.get("Hashing","0 ms")),
            ("Graph Build",  m.get("Graph Build","0 ms")),
            ("Dijkstra",     m.get("Dijkstra","0 ms")),
            ("Hierholzer",   m.get("Hierholzer","0 ms")),
            ("DP Correct",   m.get("DP Correction","0 ms")),
            ("SA + LCP",     m.get("Suffix Array+LCP","0 ms")),
        ]
        t_labels, t_vals = [], []
        for label, raw in timing_stages:
            t_labels.append(label)
            try:   t_vals.append(float(raw.replace("ms","").strip()))
            except: t_vals.append(0.0)

        colors_bar = ["#00f0ff","#cc00ff","#ff0088","#00ff99","#ff6d00","#ffe600"]
        fig_bar = go.Figure(go.Bar(
            x=t_labels, y=t_vals,
            marker=dict(color=colors_bar, line=dict(width=0),
                        opacity=.85),
            text=[f"{v:.3f}" for v in t_vals], textposition="outside",
            textfont=dict(family="Space Mono", size=9, color="#8090cc")))
        fig_bar.update_layout(
            yaxis=dict(title="ms", color="#8090cc",
                       gridcolor="rgba(100,100,255,.08)",
                       tickfont=dict(family="Space Mono",size=8)),
            xaxis=dict(color="#8090cc",
                       tickfont=dict(family="Space Mono",size=8)),
            height=240, margin=dict(l=10,r=10,t=20,b=5),
            showlegend=False,
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(3,0,18,.6)")
        st.plotly_chart(fig_bar, use_container_width=True)

        total_t = m.get("Total","—")
        st.markdown(
            f'<div class="overall-badge">'
            f'<span class="overall-label">Total wall time</span>'
            f'<span class="overall-value">{total_t}</span></div>',
            unsafe_allow_html=True)

        # Repeat analysis
        st.markdown('<div class="hf-section">Repeat Analysis (SA + LCP)</div>',
                    unsafe_allow_html=True)
        n_rep = stats["repeat_regions"]
        if n_rep == "0" or n_rep == "—":
            st.markdown(
                f'<div class="repeat-card">'
                f'<div class="repeat-title">🔁 Repeat Regions</div>'
                f'<div class="repeat-zero">No repeats found at threshold {stats["min_repeat"]}. '
                f'Assembly appears non-repetitive.</div></div>',
                unsafe_allow_html=True)
        else:
            st.markdown(
                f'<div class="repeat-card">'
                f'<div class="repeat-title">🔁 {n_rep} Repeat Region(s) Detected</div></div>',
                unsafe_allow_html=True)
            if repeats:
                df_rep = pd.DataFrame(repeats)[["rank","length","occurrences","positions","pattern"]]
                df_rep.columns = ["Rank","Length (bp)","Occurrences","Positions","Pattern (preview)"]
                st.dataframe(df_rep, use_container_width=True, hide_index=True)

        if repeats_p.exists():
            st.download_button("⬇ repeats.txt", repeats_p.read_text(),
                               "repeats.txt","text/plain", use_container_width=True)

        with st.expander("Raw stats.txt"):
            st.code(stats_p.read_text(), language="text")

    with col_right:
        # 3D graph
        st.markdown('<div class="hf-section">3D De Bruijn Graph</div>', unsafe_allow_html=True)
        with st.spinner("Rendering graph …"):
            fig3d = build_3d_graph(graph_p, max_vis)
        st.plotly_chart(fig3d, use_container_width=True)

        gj = json.loads(graph_p.read_text())
        g_kpis = [
            kpi_html("Total Nodes", str(gj.get("total_nodes","—")), "(k-1)-mers","var(--cyan)"),
            kpi_html("Total Edges", str(gj.get("total_edges","—")), "overlaps",  "var(--purple)"),
        ]
        st.markdown(f'<div class="kpi-grid">{"".join(g_kpis)}</div>', unsafe_allow_html=True)
        st.markdown(
            '<div style="font-family:Space Mono,monospace;font-size:.64rem;'
            'color:rgba(160,180,255,.3);text-align:center;margin-top:.4rem">'
            'node = (k-1)-mer · edge = k-mer overlap · colour = degree · drag to rotate · '
            'edge brightness = coverage weight</div>',
            unsafe_allow_html=True)

        st.download_button("⬇ graph_data.json", graph_p.read_text(),
                           "graph_data.json","application/json", use_container_width=True)

        # Complexity table in results
        st.markdown('<div class="hf-section">Complexity Summary</div>', unsafe_allow_html=True)
        theo = stats.get("theoretical", [])
        if theo:
            df_th = pd.DataFrame(theo)
            df_th.columns = ["Algorithm","Complexity","Description"]
            st.dataframe(df_th, use_container_width=True, hide_index=True, height=295)
            st.markdown(
                f'<div class="overall-badge">'
                f'<span class="overall-label">Overall</span>'
                f'<span class="overall-value">{stats["overall"]}</span></div>',
                unsafe_allow_html=True)


# ─────────────────────────────────────────────────────────────────────────────
#  TAB 2 — PROBLEM STATEMENT
# ─────────────────────────────────────────────────────────────────────────────
with tab_problem:
    st.markdown('<div class="hf-section">The Genome Assembly Problem</div>',
                unsafe_allow_html=True)
    st.markdown("""
    <div class="problem-card">
      <div class="problem-title">🧬 What Is Genome Assembly?</div>
      <p style="font-family:Inter,sans-serif;font-size:.88rem;color:var(--text);line-height:1.75;margin:0 0 1rem">
        DNA sequencing machines cannot read an entire genome in one pass.
        Instead they produce millions of short overlapping fragments called
        <strong style="color:var(--cyan)">reads</strong> (75 – 300 base pairs each).
        Genome assembly is the computational challenge of reconstructing the original,
        complete genomic sequence from these millions of short, noisy, overlapping fragments —
        like reassembling a shredded book when you have millions of copies
        of random overlapping excerpts, some of which contain typos.
      </p>
      <div class="problem-grid">
        <div class="prob-item">
          <div class="prob-item-title" style="color:var(--cyan)">⚡ Why It's Hard</div>
          <div class="prob-item-text">
            • Reads are short (75–300 bp) vs genomes of billions of bp<br>
            • Sequencing errors: ~1% substitutions per base<br>
            • Repeat regions appear identically in multiple genomic locations<br>
            • Both DNA strands are sequenced simultaneously<br>
            • Volume: human genome needs ~1 TB of raw read data
          </div>
        </div>
        <div class="prob-item">
          <div class="prob-item-title" style="color:var(--purple)">🧩 Our Approach</div>
          <div class="prob-item-text">
            Break each read into k-mers (substrings of length k).
            Build a <strong>De Bruijn graph</strong> where each (k-1)-mer is a node
            and each k-mer is a directed edge. The assembled genome
            corresponds to an <strong>Eulerian path</strong> through this graph —
            a path that visits every edge exactly once. Euler's theorem
            guarantees this is solvable in O(E) time.
          </div>
        </div>
        <div class="prob-item">
          <div class="prob-item-title" style="color:var(--pink)">🔬 Why Not Naive?</div>
          <div class="prob-item-text">
            The naive overlap-layout-consensus (OLC) approach computes
            pairwise overlaps between all reads: <strong>O(N²)</strong> — infeasible
            for millions of reads. Hamiltonian path (visit every node once)
            is <strong>NP-hard</strong>. De Bruijn + Eulerian path reduces the
            problem to <strong>O(N)</strong> preprocessing + <strong>O(E)</strong> traversal.
          </div>
        </div>
      </div>
    </div>
    """, unsafe_allow_html=True)

    st.markdown('<div class="hf-section">Algorithm Design Decisions</div>',
                unsafe_allow_html=True)

    decisions = [
        ("Rolling Hash (Rabin-Karp)",   "var(--cyan)",
         "Naively hashing every k-mer requires O(k) operations each — giving O(Nk) total. "
         "Rolling hash maintains a polynomial hash and updates it in O(1) per slide by removing "
         "the outgoing character and adding the incoming one using a precomputed base power. "
         "Uses Mersenne prime MOD=2⁶¹−1 with __uint128_t to avoid overflow without modular arithmetic overhead.",
         "O(N·k) → O(N)"),

        ("Bloom Filter (Dual-Pass)",   "var(--purple)",
         "Sequencing errors create singleton k-mers that don't appear in the real genome. "
         "A naïve hash set storing all k-mers would require O(N·k) bytes of RAM — gigabytes for real data. "
         "Two Bloom filters (once, twice) at fixed 8 MB total classify k-mers as solid (seen ≥2×) "
         "with ~0.02% false-positive rate. Only solid k-mers enter the De Bruijn graph.",
         "O(N·k) RAM → 8 MB fixed"),

        ("Canonical K-mers",            "var(--green)",
         "DNA is double-stranded: a read ATCG and its reverse complement CGAT represent the same genomic region. "
         "Without canonicalization, these create duplicate conflicting graph nodes. "
         "canonical = min(kmer, rev_comp(kmer)) maps both strands to the same node, "
         "halving graph size and producing biologically correct assemblies.",
         "Graph size ÷ 2"),

        ("Dijkstra + Coverage Weights", "var(--orange)",
         "Edge weight = (max_freq + 1) − observed_freq. High-coverage edges get low cost, "
         "so Dijkstra's min-heap SSSP naturally finds the path supported by the most reads. "
         "This is the 'most confident' assembly backbone — useful when the graph is fragmented "
         "and no Eulerian path exists. Lazy deletion (stale-entry skip) keeps the heap clean.",
         "O((V+E) log V)"),

        ("Suffix Array (Prefix-Doubling)", "var(--pink)",
         "Built on the assembled sequence to find all repeat regions. "
         "Prefix-doubling: start with rank[i]=char(s[i]), then sort by pairs (rank[i], rank[i+gap]) "
         "doubling gap each round. After O(log n) rounds all suffixes are uniquely ranked. "
         "Sentinel '$' ensures no two suffixes are equal. Kasai's LCP uses the key invariant "
         "that LCP decreases by at most 1 per suffix — giving O(n) total.",
         "O(n log²n) + O(n) LCP"),
    ]

    for title, color, desc, complexity in decisions:
        with st.expander(f"◈  {title}  —  {complexity}"):
            st.markdown(
                f'<div style="font-family:Inter,sans-serif;font-size:.85rem;'
                f'color:var(--text);line-height:1.7;border-left:3px solid {color};'
                f'padding-left:.9rem">{desc}</div>',
                unsafe_allow_html=True)

    st.markdown('<div class="hf-section">Real-World Context</div>', unsafe_allow_html=True)
    st.markdown("""
    <div class="problem-card" style="border-color:rgba(0,255,153,.2)">
      <div class="problem-title" style="color:var(--green)">🌍 Where Is This Used?</div>
      <div class="problem-grid">
        <div class="prob-item">
          <div class="prob-item-title" style="color:var(--green)">Medical Genomics</div>
          <div class="prob-item-text">
            Cancer mutation calling, rare disease diagnosis, pharmacogenomics.
            GATK (Google/Broad Institute) and SPAdes (Saint Petersburg) use De Bruijn graphs at their core.
          </div>
        </div>
        <div class="prob-item">
          <div class="prob-item-title" style="color:var(--cyan)">Pandemic Response</div>
          <div class="prob-item-text">
            SARS-CoV-2 variant tracking required assembling thousands of viral genomes per day.
            Short-read assembly pipelines identical to HelixForge's were run at scale globally.
          </div>
        </div>
        <div class="prob-item">
          <div class="prob-item-title" style="color:var(--purple)">Biodiversity</div>
          <div class="prob-item-text">
            Earth BioGenome Project aims to sequence all ~1.5 million eukaryotic species.
            Efficient assembly algorithms like the ones implemented here are essential to that mission.
          </div>
        </div>
      </div>
    </div>
    """, unsafe_allow_html=True)


# ─────────────────────────────────────────────────────────────────────────────
#  TAB 3 — COMPLEXITY & ALGORITHMS
# ─────────────────────────────────────────────────────────────────────────────
with tab_complexity:
    st.markdown('<div class="hf-section">Full Complexity Table</div>', unsafe_allow_html=True)

    rows = [
        ("FASTQ Reader",          "O(N)",            "N = total characters across all reads",                        "Sequential I/O"),
        ("Reverse Complement",    "O(k)",            "Per k-mer; O(N) total across all reads",                      "String manipulation"),
        ("Rolling Hash (init)",   "O(k)",            "One-time per read; O(N/k) reads × O(k) = O(N)",               "Polynomial evaluation"),
        ("Rolling Hash (slide)",  "O(1)",            "Remove old char, add new char — amortised over all slides",    "Hash recurrence"),
        ("Bloom Filter insert",   "O(1)",            "h hash functions × O(1) bit-array access each",                "Bit manipulation"),
        ("Bloom Filter query",    "O(1)",            "Same as insert; stops early on first miss",                    "Bit manipulation"),
        ("Graph add_edge",        "O(1) avg",        "Hash map insertion; O(1) average for unordered_map",           "Hash table"),
        ("Graph Construction",    "O(V + E)",        "V = unique (k-1)-mers, E = unique solid k-mers",               "Adjacency list"),
        ("Dijkstra SSSP",         "O((V+E) log V)", "Binary min-heap; lazy deletion of stale entries",              "Priority queue"),
        ("Hierholzer Euler",      "O(E)",            "Each edge pushed/popped exactly once from stack",              "DFS + stack"),
        ("Greedy Fallback",       "O(E)",            "Scans each edge at most once; O(k) hash lookup per edge",      "Greedy + hash"),
        ("DP Error Correction",   "O(N)",            "Single pass; O(1) per position comparison",                    "DP / sliding window"),
        ("Suffix Array build",    "O(n log² n)",    "O(log n) doubling rounds × O(n log n) sort each round",        "Divide & conquer + sort"),
        ("LCP Array (Kasai)",     "O(n)",            "h increments ≤ n total; h decrements ≤ n total",               "Amortised analysis"),
        ("Repeat Finding",        "O(n)",            "Linear scan of LCP array; grouping is amortised O(n)",         "LCP scan"),
        ("GC Content",            "O(n)",            "Single pass counter over assembled sequence",                  "Linear scan"),
        ("N50 Metric",            "O(c log c)",      "c = number of contigs; dominated by sort",                    "Sort + prefix sum"),
        ("JSON Graph Writer",     "O(V + E)",        "Capped at 500 nodes / 2000 edges for visualiser performance", "Sequential I/O"),
        ("Overall Pipeline",      "O(N + (V+E)logV + n log²n)", "N dominates on large datasets; SA dominates on long assemblies", "All modules"),
    ]

    df_full = pd.DataFrame(rows, columns=["Module","Complexity","Explanation","Technique"])
    st.dataframe(df_full, use_container_width=True, hide_index=True, height=580)

    st.markdown('<div class="hf-section">Space Complexity</div>', unsafe_allow_html=True)
    space_rows = [
        ("Bloom Filters (×2)",  "16 MB fixed",  "2 × 8 MB bit-arrays — independent of N"),
        ("Rolling Hash",        "O(1)",         "Stores only current hash value and base power"),
        ("De Bruijn Graph",     "O(V + E)",     "Adjacency list + edge_freq map + indeg/outdeg"),
        ("Dijkstra structures", "O(V)",         "dist[], prev[], hops[] arrays + min-heap"),
        ("Suffix Array",        "O(n)",         "SA + rank + tmp vectors; 3n integers"),
        ("LCP Array",           "O(n)",         "One integer per suffix"),
        ("Overall",             "O(N + V + E + n)", "Bloom filter constant keeps large-N practical"),
    ]
    df_space = pd.DataFrame(space_rows, columns=["Structure","Space","Notes"])
    st.dataframe(df_space, use_container_width=True, hide_index=True)

    st.markdown('<div class="hf-section">APS Topics Applied</div>', unsafe_allow_html=True)
    aps_rows = [
        ("Graph Theory",         "De Bruijn graph — directed, weighted, (k-1)-mer nodes"),
        ("Eulerian Paths",       "Hierholzer O(E) — genome assembly = Eulerian path problem"),
        ("Shortest Path SSSP",   "Dijkstra with min-heap — coverage-weighted best-path assembly"),
        ("Priority Queue",       "std::priority_queue<pair<double,string>, greater<>> for Dijkstra"),
        ("Hashing",              "Rabin-Karp rolling hash — O(1) slide per k-mer position"),
        ("Probabilistic DS",     "Bloom Filter — O(1) insert/query at fixed 8 MB vs O(Nk) exact"),
        ("Divide & Conquer",     "SA prefix-doubling — each round doubles sorted prefix length"),
        ("String Algorithms",    "Suffix Array, LCP Array, Kasai — O(n log²n) + O(n)"),
        ("Dynamic Programming",  "DP error correction — optimal substructure per base position"),
        ("Amortised Analysis",   "Kasai's h-invariant — O(n) total via amortised decrement argument"),
        ("Greedy Algorithms",    "Frequency-weighted greedy walk — locally optimal edge selection"),
        ("Sliding Window",       "Rolling hash slide + DP correction window — amortised O(1)"),
        ("Two-Pointer",          "Kasai LCP extension — extend / decrement in linear total steps"),
        ("Space-Time Tradeoff",  "Bloom filter: accept ≈0.02% FP rate to save gigabytes of RAM"),
        ("Canonical Forms",      "min(kmer, rev_comp) — equivalence class for bidirectional k-mers"),
        ("Complexity Analysis",  "Every module profiled theoretically + measured; reported in stats.txt"),
    ]
    df_aps = pd.DataFrame(aps_rows, columns=["APS Topic","Application in HelixForge"])
    st.dataframe(df_aps, use_container_width=True, hide_index=True, height=520)
