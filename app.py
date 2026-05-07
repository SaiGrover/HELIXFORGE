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
@import url('https://fonts.googleapis.com/css2?family=Orbitron:wght@400;700;900&family=Space+Mono:wght@400;700&family=Inter:wght@300;400;500;600&display=swap');

/* ── CSS Variables — softer, slightly desaturated palette ── */
:root {
  --bg:      #07080f;
  --bg2:     #0d0f1c;
  --card:    rgba(16, 18, 36, 0.82);
  --cyan:    #4dd8e8;
  --purple:  #a855f7;
  --pink:    #e0508a;
  --green:   #3ecf8e;
  --orange:  #f5874f;
  --yellow:  #e8c84a;
  --teal:    #2db8a0;
  --text:    #c8d0e8;
  --muted:   rgba(170, 180, 210, 0.55);
  --border:  rgba(120, 100, 200, 0.16);
  --shadow:  rgba(0,0,0,0.45);
}

/* ── Keyframes ── */
@keyframes gradshift {
  0%,100% { background-position: 0% 50%; }
  50%      { background-position: 100% 50%; }
}
@keyframes glow-in { from{opacity:0;transform:translateY(10px)} to{opacity:1;transform:translateY(0)} }

/* ── Base — body text readable at 14px ── */
html, body, [class*="css"] {
  font-family: 'Inter', sans-serif;
  font-size: 14px;
  line-height: 1.6;
  color: var(--text);
}
.stApp {
  background: var(--bg);
  background-image:
    radial-gradient(ellipse 80% 50% at 10% 5%,  rgba(77,216,232,.05) 0%, transparent 60%),
    radial-gradient(ellipse 60% 40% at 90% 95%,  rgba(168,85,247,.06) 0%, transparent 60%);
}
#MainMenu, footer, header { visibility:hidden; }
.block-container { padding:1.6rem 2rem 4rem; max-width:1480px; }

/* ── Sidebar ── */
[data-testid="stSidebar"] {
  background: linear-gradient(180deg, #0b0d1e 0%, #07080f 100%) !important;
  border-right: 1px solid rgba(168,85,247,.15) !important;
}
[data-testid="stSidebar"] * { color: var(--text) !important; }
[data-testid="stSidebar"] p,
[data-testid="stSidebar"] li { font-size: 13px !important; line-height: 1.65 !important; }

/* ── File uploader ── */
[data-testid="stFileUploader"] {
  background: rgba(77,216,232,.03) !important;
  border: 1.5px dashed rgba(77,216,232,.25) !important;
  border-radius: 10px !important;
  transition: border-color .3s !important;
}
[data-testid="stFileUploader"]:hover { border-color: rgba(77,216,232,.5) !important; }

/* ── Sliders ── */
[data-testid="stSlider"] [class*="track"]        { background: rgba(168,85,247,.18) !important; }
[data-testid="stSlider"] [class*="track--filled"] { background: linear-gradient(90deg,var(--cyan),var(--purple)) !important; }
[data-testid="stSlider"] [class*="thumb"]         {
  background: var(--purple) !important;
  border: 2px solid var(--cyan) !important;
  box-shadow: 0 0 8px rgba(77,216,232,.4) !important;
}

/* ── Run button ── */
.stButton > button {
  background: linear-gradient(135deg, #2a9bb5 0%, #7c3aed 55%, #be4080 100%) !important;
  background-size: 200% 200% !important;
  animation: gradshift 4s ease infinite !important;
  color: #fff !important;
  font-family: 'Orbitron', monospace !important;
  font-weight: 700 !important;
  font-size: 13px !important;
  letter-spacing: .12em !important;
  text-transform: uppercase !important;
  border: none !important;
  border-radius: 8px !important;
  padding: .65rem 1.4rem !important;
  box-shadow: 0 0 18px rgba(77,216,232,.25), 0 0 36px rgba(168,85,247,.15) !important;
  transition: all .3s !important;
}
.stButton > button:hover {
  transform: translateY(-2px) !important;
  box-shadow: 0 4px 24px rgba(77,216,232,.35), 0 0 48px rgba(168,85,247,.22) !important;
}

/* ── Download buttons ── */
[data-testid="stDownloadButton"] > button {
  background: transparent !important;
  color: var(--cyan) !important;
  border: 1px solid rgba(77,216,232,.3) !important;
  font-family: 'Inter', sans-serif !important;
  font-size: 13px !important;
  border-radius: 6px !important;
  letter-spacing: .03em !important;
  transition: all .25s !important;
}
[data-testid="stDownloadButton"] > button:hover {
  background: rgba(77,216,232,.07) !important;
  border-color: var(--cyan) !important;
}

/* ── Expanders ── */
[data-testid="stExpander"] {
  background: rgba(16,18,36,.7) !important;
  border: 1px solid var(--border) !important;
  border-radius: 8px !important;
}
[data-testid="stExpander"] summary {
  color: var(--cyan) !important;
  font-family: 'Inter', sans-serif !important;
  font-size: 14px !important;
  font-weight: 600 !important;
}
[data-testid="stExpander"] p,
[data-testid="stExpander"] li { font-size: 13.5px !important; line-height: 1.7 !important; }
[data-testid="stExpander"] > div { background: rgba(10,12,26,.85) !important; }

/* ── Code / pre ── */
pre, code {
  background: #0c0e1e !important;
  color: #5de8b8 !important;
  font-family: 'Space Mono', monospace !important;
  font-size: 13px !important;
  line-height: 1.65 !important;
}

/* ── DataFrames ── */
.stDataFrame { border: 1px solid var(--border) !important; border-radius: 8px !important; }
.stDataFrame thead th {
  background: #0f1128 !important;
  color: var(--cyan) !important;
  font-family: 'Inter', sans-serif !important;
  font-size: 13px !important;
  font-weight: 600 !important;
}
.stDataFrame tbody td {
  background: #090b1a !important;
  color: var(--text) !important;
  font-family: 'Inter', sans-serif !important;
  font-size: 13px !important;
  line-height: 1.55 !important;
  border-color: var(--border) !important;
}

/* ── Spinner ── */
.stSpinner > div { border-top-color: var(--cyan) !important; }

/* ── Streamlit default text ── */
.stMarkdown p  { font-size: 14px; line-height: 1.7; color: var(--text); }
.stMarkdown li { font-size: 13.5px; line-height: 1.7; }

/* ═══════════════════════════════════════════
   CUSTOM COMPONENTS
═══════════════════════════════════════════ */

/* Hero */
.hf-hero { padding: 1.6rem 0 .8rem; animation: glow-in .5s ease both; }
.hf-logo {
  font-family: 'Orbitron', monospace;
  font-weight: 900;
  font-size: 2.6rem;
  background: linear-gradient(120deg, #4dd8e8 0%, #a855f7 45%, #e0508a 100%);
  background-size: 280% 280%;
  animation: gradshift 5s ease infinite;
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
  background-clip: text;
  line-height: 1.15;
}
.hf-sub {
  font-family: 'Inter', sans-serif;
  font-size: 12px;
  color: var(--muted);
  letter-spacing: .18em;
  text-transform: uppercase;
  margin-top: .4rem;
}
.hf-divider {
  height: 1px;
  background: linear-gradient(90deg, var(--cyan) 0%, var(--purple) 40%, var(--pink) 80%, transparent 100%);
  margin: 1rem 0;
  opacity: .35;
}

/* Section labels */
.hf-section {
  font-family: 'Inter', sans-serif;
  font-size: 11px;
  font-weight: 700;
  letter-spacing: .2em;
  text-transform: uppercase;
  color: var(--cyan);
  margin: 1.8rem 0 .9rem;
  display: flex;
  align-items: center;
  gap: .7rem;
  opacity: .85;
}
.hf-section::after {
  content: '';
  flex: 1;
  height: 1px;
  background: linear-gradient(90deg, rgba(77,216,232,.3), rgba(168,85,247,.15), transparent);
}

/* KPI cards */
.kpi-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(138px,1fr)); gap:.8rem; margin:.4rem 0 1.1rem; }
.kpi-card {
  background: var(--card);
  border: 1px solid rgba(255,255,255,.07);
  border-radius: 10px;
  padding: 1rem 1rem .9rem;
  position: relative;
  overflow: hidden;
  transition: transform .2s, box-shadow .2s;
  backdrop-filter: blur(6px);
}
.kpi-card::before {
  content:'';
  position:absolute; top:0; left:0; right:0; height:2px;
  background: var(--accent, var(--cyan));
  opacity: .75;
}
.kpi-card:hover { transform: translateY(-2px); box-shadow: 0 6px 24px var(--shadow); }
.kpi-label { font-family:'Inter',sans-serif; font-size:11px; font-weight:600; letter-spacing:.1em; text-transform:uppercase; color:var(--muted); margin-bottom:.4rem; }
.kpi-value { font-family:'Orbitron',monospace; font-weight:700; font-size:1.55rem; color:#e8ecff; line-height:1.1; }
.kpi-unit  { font-family:'Inter',sans-serif; font-size:11px; color:var(--muted); margin-top:.15rem; }

/* Method badge */
.method-badge {
  display:inline-flex; align-items:center; gap:.6rem;
  background: rgba(77,216,232,.05);
  border: 1px solid rgba(77,216,232,.22);
  border-radius: 8px;
  padding: .5rem 1rem;
  margin: .5rem 0 .9rem;
}
.method-label { font-family:'Inter',sans-serif; font-size:12px; color:var(--muted); letter-spacing:.08em; text-transform:uppercase; }
.method-value { font-family:'Inter',sans-serif; font-size:14px; font-weight:700; color:var(--cyan); }

/* Sequence terminal */
.seq-terminal { background:#0b0d1f; border:1px solid rgba(77,216,232,.18); border-radius:10px; overflow:hidden; }
.seq-bar { background:#0f1128; padding:.45rem 1rem; display:flex; align-items:center; gap:.5rem; border-bottom:1px solid rgba(77,216,232,.1); }
.seq-dot { width:9px; height:9px; border-radius:50%; }
.seq-body {
  font-family:'Space Mono',monospace; font-size:13px; line-height:1.85;
  color:#5de8b8; padding:.9rem 1.1rem;
  max-height:200px; overflow-y:auto; word-break:break-all; white-space:pre-wrap;
}
.seq-body::-webkit-scrollbar { width:3px; }
.seq-body::-webkit-scrollbar-thumb { background:rgba(77,216,232,.25); border-radius:2px; }

/* Pipeline steps */
.pipeline { display:grid; grid-template-columns:repeat(7,1fr); gap:0; margin:1rem 0 1.6rem; }
.pipe-step {
  background: rgba(16,18,36,.6);
  border: 1px solid var(--border);
  border-right: none;
  padding: 1rem .6rem .9rem;
  text-align: center;
  position: relative;
  transition: background .2s;
}
.pipe-step:first-child { border-radius:9px 0 0 9px; }
.pipe-step:last-child  { border-right:1px solid var(--border); border-radius:0 9px 9px 0; }
.pipe-step:hover { background:rgba(77,216,232,.05); }
.pipe-step::after {
  content:'›'; position:absolute; right:-10px; top:50%; transform:translateY(-50%);
  color:rgba(77,216,232,.3); font-size:1.3rem; z-index:2;
}
.pipe-step:last-child::after { display:none; }
.pipe-num   { font-family:'Inter',sans-serif; font-size:10px; font-weight:600; color:var(--muted); letter-spacing:.08em; }
.pipe-icon  { font-size:1.3rem; margin:.2rem 0; line-height:1; }
.pipe-title { font-family:'Inter',sans-serif; font-weight:700; font-size:12px; color:#d8dcf0; margin-bottom:.1rem; }
.pipe-sub   { font-family:'Space Mono',monospace; font-size:10px; color:var(--muted); }

/* Problem statement */
.problem-card {
  background: rgba(16,18,36,.7);
  border: 1px solid rgba(168,85,247,.2);
  border-radius: 12px;
  padding: 1.6rem 1.8rem;
  margin: 1rem 0 1.3rem;
  position: relative;
  overflow: hidden;
}
.problem-card::before {
  content:'';
  position:absolute; top:0; left:0; right:0; height:2px;
  background: linear-gradient(90deg, var(--cyan), var(--purple), var(--pink));
  opacity: .6;
}
.problem-title { font-family:'Inter',sans-serif; font-weight:700; font-size:17px; color:var(--cyan); margin-bottom:.9rem; }
.problem-grid  { display:grid; grid-template-columns:repeat(3,1fr); gap:.9rem; margin-top:.9rem; }
.prob-item     { background:rgba(255,255,255,.03); border:1px solid var(--border); border-radius:8px; padding:.85rem .95rem; }
.prob-item-title { font-family:'Inter',sans-serif; font-size:12px; font-weight:700; letter-spacing:.08em; text-transform:uppercase; margin-bottom:.35rem; }
.prob-item-text  { font-family:'Inter',sans-serif; font-size:13.5px; color:var(--text); line-height:1.65; }

/* K-mer preview */
.kmer-preview {
  background: rgba(10,12,26,.9);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: .65rem .9rem;
  font-family: 'Space Mono', monospace;
  font-size: 12px;
  color: var(--muted);
  margin: .3rem 0 .75rem;
  overflow-x: auto;
}
.kmer-seq  { font-size: 13px; letter-spacing: .05em; line-height:1.8; }
.kmer-high { color: var(--cyan); background: rgba(77,216,232,.1); border-radius:3px; padding:0 2px; font-weight:700; }
.kmer-info { margin-top:.3rem; font-size:11px; display:flex; gap:1rem; flex-wrap:wrap; }
.kmer-info-item { color:var(--muted); }
.kmer-info-item span { color:var(--purple); font-weight:700; }

/* Verification card */
.verify-card {
  background: rgba(16,18,36,.7);
  border: 1px solid rgba(62,207,142,.18);
  border-radius: 10px;
  padding: 1.1rem 1.3rem;
  margin: .5rem 0 .9rem;
}
.verify-title {
  font-family:'Inter',sans-serif; font-size:12px; font-weight:700;
  letter-spacing:.1em; text-transform:uppercase; margin-bottom:.8rem; color: var(--green);
}
.verify-row {
  display:flex; align-items:center; gap:.8rem;
  padding:.5rem 0; border-bottom:1px solid rgba(255,255,255,.04);
  font-family:'Inter',sans-serif; font-size:13.5px;
}
.verify-row:last-child { border-bottom:none; }
.verify-label { color:var(--muted); width:140px; flex-shrink:0; font-size:13px; }
.verify-bar-wrap { flex:1; height:5px; background:rgba(255,255,255,.06); border-radius:3px; overflow:hidden; }
.verify-bar { height:100%; border-radius:3px; transition:width .5s ease; }
.verify-status { font-size:12px; font-weight:700; letter-spacing:.05em; white-space:nowrap; }
.verify-val { color:#e0e4f8; width:85px; text-align:right; flex-shrink:0; font-size:13px; }
.blast-btn {
  display:inline-flex; align-items:center; gap:.45rem;
  background: rgba(77,216,232,.05);
  border: 1px solid rgba(77,216,232,.25);
  border-radius:7px; padding:.5rem 1rem;
  font-family:'Inter',sans-serif; font-size:13px; font-weight:600;
  color:var(--cyan); text-decoration:none; margin-top:.5rem;
  transition: all .2s;
}
.blast-btn:hover { background:rgba(77,216,232,.1); border-color:var(--cyan); }

/* Repeat card */
.repeat-card {
  background: rgba(16,18,36,.7);
  border: 1px solid rgba(245,135,79,.2);
  border-radius: 10px;
  padding: .95rem 1.1rem;
  margin-top:.5rem;
}
.repeat-title { font-family:'Inter',sans-serif; font-size:12px; font-weight:700; color:var(--orange); letter-spacing:.08em; text-transform:uppercase; margin-bottom:.5rem; }
.repeat-zero  { font-family:'Inter',sans-serif; font-size:13.5px; color:var(--muted); }

/* Info banner */
.info-banner {
  background: rgba(77,216,232,.04);
  border: 1px solid rgba(77,216,232,.18);
  border-left: 3px solid var(--cyan);
  border-radius: 0 8px 8px 0;
  padding: .85rem 1.1rem;
  font-family: 'Inter', sans-serif;
  font-size: 14px;
  color: var(--text);
  margin: .7rem 0;
}

/* Sidebar logo */
.sb-logo {
  font-family: 'Orbitron', monospace;
  font-weight: 900;
  font-size: 1.15rem;
  background: linear-gradient(120deg, var(--cyan), var(--purple), var(--pink));
  background-size: 200% 200%;
  animation: gradshift 5s ease infinite;
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
  background-clip: text;
  padding: .3rem 0 .9rem;
}
.sb-section {
  font-family: 'Inter', sans-serif;
  font-size: 10px;
  font-weight: 700;
  letter-spacing: .18em;
  text-transform: uppercase;
  color: var(--muted) !important;
  margin: .9rem 0 .35rem;
  padding-bottom: .3rem;
  border-bottom: 1px solid var(--border);
}

/* Overall badge */
.overall-badge {
  display:inline-flex; align-items:center; gap:.6rem;
  background: rgba(168,85,247,.05);
  border: 1px solid rgba(168,85,247,.22);
  border-radius: 8px;
  padding: .5rem 1rem;
  margin-top: .7rem;
}
.overall-label { font-family:'Inter',sans-serif; font-size:12px; color:var(--muted); letter-spacing:.08em; text-transform:uppercase; }
.overall-value { font-family:'Inter',sans-serif; font-size:14px; font-weight:700; color:var(--purple); }

/* Comparison row */
.cmp-row { display:flex; gap:.9rem; margin:.4rem 0 .9rem; }
.cmp-box { flex:1; background:rgba(255,255,255,.03); border:1px solid var(--border); border-radius:8px; padding:.65rem .95rem; }
.cmp-box-label { font-family:'Inter',sans-serif; font-size:11px; font-weight:600; color:var(--muted); text-transform:uppercase; letter-spacing:.08em; }
.cmp-box-value { font-family:'Orbitron',monospace; font-size:1.05rem; font-weight:700; margin-top:.2rem; }

/* Tabs */
[data-testid="stTabs"] [role="tablist"] { gap:.3rem; border-bottom:1px solid var(--border); }
[data-testid="stTabs"] [role="tab"] {
  font-family:'Inter',monospace !important;
  font-size:13px !important;
  font-weight:600 !important;
  letter-spacing:.06em !important;
  color:var(--muted) !important;
  border-radius:6px 6px 0 0 !important;
  border: 1px solid transparent !important;
  border-bottom: none !important;
  padding:.45rem 1rem !important;
  transition: all .2s !important;
}
[data-testid="stTabs"] [role="tab"][aria-selected="true"] {
  color:var(--cyan) !important;
  border-color: var(--border) !important;
  background: rgba(77,216,232,.05) !important;
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
                # Strip leading "- " dash bullet from key names
                step_name = k_str.lstrip("- ").strip()
                # Split on " - " to separate complexity expression from description
                if " - " in v_str:
                    complexity, desc_tail = v_str.split(" - ", 1)
                    complexity = complexity.strip(); desc_tail = desc_tail.strip()
                else:
                    parts = v_str.split(None, 1)
                    complexity = parts[0]; desc_tail = parts[1].strip() if len(parts) > 1 else ""
                data["theoretical"].append({"step":step_name,"complexity":complexity,"reason":desc_tail})
        elif in_meas:
            data["measured"][k_str] = v_str
        else:
            if   "reads processed"   in kl: data["reads"]          = v_str
            elif "total k-mer"       in kl or "unique k-mer" in kl: data["total_kmers"] = v_str
            elif "k-mer size"        in kl or "kmer size" in kl: data["kmer_size"] = v_str
            elif "nodes (v)"         in kl or "graph nodes" in kl: data["nodes"]   = v_str
            elif "edges (e)"         in kl or "graph edges" in kl: data["edges"]   = v_str
            elif "method selected"   in kl: data["method"]         = v_str
            elif "dijkstra length"   in kl: data["dijk_len"]       = v_str
            elif "hierholzer length" in kl: data["hier_len"]       = v_str
            elif "final length"      in kl or "assembly length" in kl:
                data["final_len"] = v_str.replace(" bp","").strip()
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
    raw_nodes = gdata["nodes"][:max_nodes]

    # Support old format (list of strings) and new format (list of dicts with role/coverage)
    if raw_nodes and isinstance(raw_nodes[0], dict):
        nodes    = [n["id"]               for n in raw_nodes]
        roles    = {n["id"]: n.get("role",     "other") for n in raw_nodes}
        coverages= {n["id"]: n.get("coverage", 1)       for n in raw_nodes}
    else:
        nodes    = raw_nodes
        roles    = {n: "other" for n in nodes}
        coverages= {n: 1       for n in nodes}

    node_set = set(nodes)

    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    edge_weights = {}
    path_edges   = set()

    for edge in gdata["edges"]:
        if isinstance(edge, list):
            u, v, w, on_path = edge[0], edge[1], 1, False
        else:
            u       = edge.get("from",    "")
            v       = edge.get("to",      "")
            w       = edge.get("weight",  1)
            on_path = edge.get("on_path", False)
        if u in node_set and v in node_set:
            G.add_edge(u, v, weight=float(w))
            edge_weights[(u, v)] = w
            if on_path:
                path_edges.add((u, v))

    pos = nx.spring_layout(G, dim=3, seed=42, k=0.65, weight="weight", iterations=90)

    ROLE_COLOR = {
        "path":   "#4dd8e8",
        "branch": "#f5874f",
        "tip":    "#a855f7",
        "other":  "rgba(100,110,150,0.4)",
    }
    ROLE_SIZE    = {"path": 6,   "branch": 3,    "tip": 2.5, "other": 2}
    ROLE_OPACITY = {"path": 1.0, "branch": 0.65, "tip": 0.5, "other": 0.3}
    ROLE_LABEL   = {
        "path":   "Assembly path",
        "branch": "Branch point",
        "tip":    "Dead-end tip",
        "other":  "Context node",
    }

    traces = []

    ox, oy, oz = [], [], []
    px, py, pz = [], [], []
    for (u, v) in edge_weights:
        if u not in pos or v not in pos:
            continue
        x0, y0, z0 = pos[u]; x1, y1, z1 = pos[v]
        if (u, v) in path_edges:
            px += [x0,x1,None]; py += [y0,y1,None]; pz += [z0,z1,None]
        else:
            ox += [x0,x1,None]; oy += [y0,y1,None]; oz += [z0,z1,None]

    if ox:
        traces.append(go.Scatter3d(
            x=ox, y=oy, z=oz, mode="lines",
            line=dict(color="rgba(70,80,130,0.18)", width=0.7),
            hoverinfo="none", name="Other edges", showlegend=True))
    if px:
        traces.append(go.Scatter3d(
            x=px, y=py, z=pz, mode="lines",
            line=dict(color="#4dd8e8", width=2.8),
            hoverinfo="none", name="Assembly path edges", showlegend=True))

    for role in ("other", "tip", "branch", "path"):
        rn = [n for n in nodes if roles.get(n) == role and n in pos]
        if not rn:
            continue
        traces.append(go.Scatter3d(
            x=[pos[n][0] for n in rn],
            y=[pos[n][1] for n in rn],
            z=[pos[n][2] for n in rn],
            mode="markers",
            marker=dict(
                size=ROLE_SIZE[role],
                color=ROLE_COLOR[role],
                opacity=ROLE_OPACITY[role],
                line=dict(width=0)),
            text=[f"{n}<br>cov: {coverages.get(n,0)}" for n in rn],
            hovertemplate=f"<b>%{{text}}</b><br>{ROLE_LABEL[role]}<extra></extra>",
            name=ROLE_LABEL[role]))

    fig = go.Figure(data=traces)
    fig.update_layout(
        scene=dict(
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, showbackground=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, showbackground=False),
            zaxis=dict(showgrid=False, zeroline=False, showticklabels=False, showbackground=False),
            bgcolor="rgba(0,0,0,0)"),
        legend=dict(
            x=0.01, y=0.99,
            font=dict(size=10, color="#a0b0cc"),
            bgcolor="rgba(10,12,28,0.6)",
            bordercolor="rgba(100,80,200,0.2)",
            borderwidth=1),
        margin=dict(l=0, r=0, t=0, b=0), height=520,
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

    st.markdown('<div class="sb-section">Why Each Algorithm?</div>', unsafe_allow_html=True)
    with st.expander("🧬 Why De Bruijn graphs?"):
        st.write("Think of the DNA reads like shredded sentence strips. "
                 "Instead of finding which strip connects to which (that's insanely hard — NP-hard to be exact), "
                 "we chop everything into small fixed-size words (k-mers) and ask: "
                 "which words flow into which? That gives us a graph where **finding the genome = tracing a path "
                 "that uses every edge once** — and that's easy to solve in O(E) time using Hierholzer's algorithm.")
    with st.expander("⚡ Why rolling hash?"):
        st.write("Imagine reading a 300-letter DNA read and computing a fingerprint (hash) for every "
                 "window of k letters. Doing it from scratch each time = very slow. "
                 "Rolling hash is a trick: **slide the window by one letter, subtract the old letter, "
                 "add the new one** — the hash updates in O(1) instead of O(k). "
                 "Across millions of reads, that's a massive speed-up.")
    with st.expander("🌸 Why Bloom filter?"):
        st.write("Sequencing machines make ~1% mistakes. Those errors create fake k-mers that appear "
                 "only once. We want to throw those away. "
                 "But storing all k-mers in memory = gigabytes of RAM. "
                 "A Bloom filter is a **tiny 8 MB checklist** — if a k-mer shows up twice, it's probably real. "
                 "If only once, it's probably a typo. We keep only the 'seen twice' ones.")
    with st.expander("🎯 Why Dijkstra for assembly?"):
        st.write("Sometimes the De Bruijn graph is too messy for a clean Eulerian path. "
                 "So we run Dijkstra as a backup — but instead of shortest distance, "
                 "we treat **high-coverage edges (seen many times in reads) as low-cost roads**. "
                 "Dijkstra naturally picks the path through the most-trusted data. "
                 "It's like choosing the highway with the most traffic — more reads = more confidence.")
    with st.expander("🔬 Why Suffix Array?"):
        st.write("After assembly, we want to know: are there repeated sections in the genome? "
                 "A Suffix Array sorts all possible suffixes of the assembled sequence alphabetically. "
                 "Then the LCP (Longest Common Prefix) array tells us **which neighbours share a long "
                 "starting stretch** — those are repeat regions. "
                 "Building it takes O(n log²n) time, but finding all repeats from the result is just O(n).")


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

    # ── Genome Verification panel ─────────────────────────────────────────────
    st.markdown('<div class="hf-section">Genome Verification</div>', unsafe_allow_html=True)

    def verify_gc(gc_str):
        """Return (value, bar_pct, color, status, note) for GC content."""
        try:
            v = float(gc_str)
        except:
            return None, 0, "#555", "NO DATA", "Run assembly first"
        # Broad range covers all common plant barcodes:
        # rbcL/MATK/trnL: 33-48%   ITS2/18S: 45-65%
        if 33 <= v <= 65:
            color, status = "#00ff99", "✓ NORMAL"
        elif 25 <= v < 33 or 65 < v <= 75:
            color, status = "#ffe600", "⚠ BORDERLINE"
        else:
            color, status = "#ff0088", "✗ UNUSUAL"
        note = ("GC ranges by marker: ITS2 ≈ 45–60% · 18S ≈ 50–60% · "
                "rbcL ≈ 33–45% · MATK ≈ 35–48% · trnL ≈ 33–42% · psbA3 ≈ 35–45%. "
                "Values below 30% or above 75% often indicate contamination or adapter sequence.")
        return v, min(v, 100), color, status, note

    def verify_n50(n50_str, asm_len_str, kmer_size_str):
        """Return (n50_val, asm_len, quality_str, color, note)."""
        try:
            n50 = int(n50_str)
        except:
            return None, None, "NO DATA", "#555", ""
        try:
            alen = int(asm_len_str)
        except:
            alen = n50
        try:
            k_used = int(kmer_size_str)
        except:
            k_used = 21
        ratio = n50 / alen if alen > 0 else 0
        # Flag as TRIVIAL only if assembly == 1 k-mer (alen <= k)
        if alen <= k_used:
            color   = "#ff0088"
            quality = "✗ TRIVIAL — assembly = 1 k-mer, lower k or increase coverage"
        elif ratio >= 0.5 and alen >= 200:
            color, quality = "#00ff99", "✓ GOOD"
        elif ratio >= 0.5:
            color, quality = "#ffe600", "⚠ SHORT but contiguous — try lower k for more length"
        elif ratio >= 0.2:
            color, quality = "#ffe600", "⚠ FRAGMENTED"
        else:
            color, quality = "#ff0088", "✗ HIGHLY FRAGMENTED"
        note = ("N50 = the length L where contigs ≥ L cover at least half the assembly. "
                "N50 = assembly length means one clean contig. "
                "Short assemblies on ONT data: try k=17–21. "
                "Very low N50 relative to length = many tiny fragments.")
        return n50, alen, quality, color, note

    gc_val, gc_bar, gc_color, gc_status, gc_note = verify_gc(stats["gc"])
    n50_v, alen_v, n50_quality, n50_color, n50_note = verify_n50(stats["n50"], stats["final_len"], stats["kmer_size"])

    # Build blast URL from assembled sequence
    blast_seq = sequence[:500] if sequence else ""
    import urllib.parse
    blast_url = ("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?"
                 "PAGE_TYPE=BlastSearch&PROGRAM=blastn&DATABASE=nt"
                 f"&QUERY={urllib.parse.quote(blast_seq)}&CMD=Put") if blast_seq else ""

    gc_display  = f"{gc_val:.2f} %" if gc_val is not None else "—"
    n50_display = f"{n50_v} bp"     if n50_v  is not None else "—"

    # Overall verdict
    try:
        kmer_size_int = int(stats["kmer_size"])
    except:
        kmer_size_int = 0
    asm_len_int = alen_v or 0
    gc_ok  = gc_val is not None and 25 <= gc_val <= 75
    asm_ok = asm_len_int >= 200

    # GC ranges: rbcL/trnL ~33-45%, ITS2 ~45-60%, 18S ~50-60%, MATK ~35-48%
    gc_borderline = gc_val is not None and (30 <= gc_val < 40 or 65 < gc_val <= 75)
    gc_unusual    = gc_val is not None and (gc_val < 30 or gc_val > 75)

    if asm_ok and gc_ok:
        verdict_color = "#00ff99"
        verdict_icon  = "✅"
        verdict_text  = "ASSEMBLY LOOKS VALID — GC content normal, length meaningful. BLAST to confirm organism match."
    elif not asm_ok and asm_len_int <= kmer_size_int:
        verdict_color = "#ff0088"
        verdict_icon  = "❌"
        verdict_text  = (f"TRIVIAL ASSEMBLY — output is just one {kmer_size_int}-bp k-mer. "
                         f"Lower k to 17–19 in the sidebar and re-run.")
    elif not asm_ok and gc_borderline:
        verdict_color = "#ffe600"
        verdict_icon  = "⚠️"
        verdict_text  = (f"PARTIAL ASSEMBLY ({asm_len_int} bp) with borderline GC ({gc_val:.1f}%). "
                         f"GC 33–45% is normal for rbcL/MATK/trnL markers. "
                         f"Try k=17–19 for a longer assembly.")
    elif not asm_ok:
        verdict_color = "#ffe600"
        verdict_icon  = "⚠️"
        verdict_text  = (f"PARTIAL ASSEMBLY ({asm_len_int} bp) — try k=17–19 for more length. "
                         f"Assembly is real but shorter than the full gene.")
    elif gc_unusual:
        verdict_color = "#ff0088"
        verdict_icon  = "❌"
        verdict_text  = "UNUSUAL GC CONTENT — likely adapter contamination or wrong organism. Check your FASTQ source."
    else:
        verdict_color = "#ffe600"
        verdict_icon  = "⚠️"
        verdict_text  = (f"BORDERLINE GC ({gc_val:.1f}%). "
                         f"Normal for rbcL (33–45%) and MATK (35–48%). BLAST to confirm.")

    st.markdown(
        f'<div style="background:rgba(255,255,255,.03);border:1px solid {verdict_color}40;'
        f'border-left:4px solid {verdict_color};border-radius:8px;padding:.8rem 1.1rem;'
        f'margin-bottom:.8rem;font-family:Space Mono,monospace;font-size:.8rem;'
        f'color:{verdict_color};line-height:1.5">'
        f'{verdict_icon} <strong>VERDICT:</strong> {verdict_text}</div>',
        unsafe_allow_html=True)

    st.markdown(f"""
    <div class="verify-card">
      <div class="verify-title" style="color:var(--green)">🔍 Assembly Quality Checks</div>

      <div class="verify-row">
        <div class="verify-label">GC Content</div>
        <div class="verify-bar-wrap">
          <div class="verify-bar" style="width:{gc_bar}%;background:{gc_color};box-shadow:0 0 6px {gc_color}60"></div>
        </div>
        <div class="verify-val">{gc_display}</div>
        <div class="verify-status" style="color:{gc_color}">{gc_status}</div>
      </div>
      <div style="font-family:Space Mono,monospace;font-size:.64rem;color:var(--muted);
                  padding:.15rem 0 .55rem 130px;line-height:1.5">{gc_note}</div>

      <div class="verify-row">
        <div class="verify-label">N50 Score</div>
        <div class="verify-bar-wrap">
          <div class="verify-bar" style="width:{min((n50_v or 0)/(alen_v or 1)*100,100):.1f}%;
               background:{n50_color};box-shadow:0 0 6px {n50_color}60"></div>
        </div>
        <div class="verify-val">{n50_display}</div>
        <div class="verify-status" style="color:{n50_color}">{n50_quality}</div>
      </div>
      <div style="font-family:Space Mono,monospace;font-size:.64rem;color:var(--muted);
                  padding:.15rem 0 .55rem 130px;line-height:1.5">{n50_note}</div>

      <div class="verify-row">
        <div class="verify-label">Assembly Length</div>
        <div class="verify-bar-wrap" style="background:transparent"></div>
        <div class="verify-val">{stats["final_len"]} bp</div>
        <div class="verify-status" style="color:{'#00ff99' if (alen_v or 0) >= 200 else ('#ff0088' if (alen_v or 0) <= kmer_size_int else '#ffe600')}">
          {'✓ Meaningful assembly' if (alen_v or 0) >= 200
           else ('✗ Trivial — only 1 k-mer assembled, lower k' if (alen_v or 0) <= kmer_size_int
                 else f'⚠ Partial assembly — {alen_v} bp (try k=17–19 for longer result)')}
        </div>
      </div>
      <div style="font-family:Space Mono,monospace;font-size:.64rem;color:var(--muted);
                  padding:.15rem 0 .55rem 130px;line-height:1.5">
        Expected lengths: ITS2 ≈ 200–500 bp · rbcL ≈ 550 bp · 18S ≈ 1800 bp · MATK ≈ 900 bp · psbA3 ≈ 450 bp.<br>
        If length = k (e.g. 15 bp for k=15), no path found — lower k. &nbsp;
        If 100–300 bp partial, try k=17–19 for better coverage.
        Good values: <strong style="color:#ffe600">k=17–19</strong> for most ONT barcode files.
      </div>
    </div>
    """, unsafe_allow_html=True)

    # BLAST + tips row
    vcol1, vcol2 = st.columns([1, 1])
    with vcol1:
        if blast_url and blast_seq:
            st.markdown(
                f'<a class="blast-btn" href="{blast_url}" target="_blank">'
                f'🔗 BLAST this sequence on NCBI</a>'
                f'<div style="font-family:Space Mono,monospace;font-size:.62rem;'
                f'color:var(--muted);margin-top:.35rem">'
                f'Opens NCBI BLASTn with your assembly — shows the closest matching organism in the database.</div>',
                unsafe_allow_html=True)
    with vcol2:
        with st.expander("💡 How to improve assembly quality"):
            st.markdown("""
**If assembly is too short (< 200 bp):**
- Lower k (try 15–19 for short-read data)
- Check the FASTQ has enough reads (ideally 1000+)
- The Bloom filter needs each k-mer seen ≥2× — higher coverage helps

**If GC content is unusual (< 30% or > 75%):**
- May indicate adapter contamination in your FASTQ
- Try filtering reads with quality < Q10 before assembling

**How to verify your assembly is biologically correct:**
1. Click BLAST → see if top hit matches your expected organism
2. Compare assembly length to known gene lengths for your marker (ITS2, rbcL, etc.)
3. GC content should roughly match known values for your species group
4. If N50 ≈ assembly length → assembly is one clean contig ✓
            """)

    st.markdown('<div class="hf-divider"></div>', unsafe_allow_html=True)

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

        def _ms(m, *keys):
            """Return float ms from the first matching key in dict m."""
            for key in keys:
                if key in m:
                    try: return float(m[key].replace("ms","").strip())
                    except: pass
            return 0.0

        timing_stages = [
            ("Hashing",     _ms(m, "Hashing", "Hashing Time")),
            ("Graph Build", _ms(m, "Graph Build", "Graph Build Time")),
            ("Dijkstra",    _ms(m, "Dijkstra", "Dijkstra Time")),
            ("Hierholzer",  _ms(m, "Hierholzer", "Traversal Time", "Hierholzer Time", "Traversal")),
            ("DP Correct",  _ms(m, "DP Correction", "DP Time", "DP Correct")),
            ("SA + LCP",    _ms(m, "Suffix Array+LCP", "SA Time", "SA+LCP Time", "Suffix Array")),
        ]
        t_labels = [s[0] for s in timing_stages]
        t_vals   = [s[1] for s in timing_stages]

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

        total_t = m.get("Total", m.get("Total Time", "—"))
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
        A DNA sequencing machine can't read your entire genome in one shot — it's like trying to read
        a book by only ever seeing random 200-letter chunks of it, over and over again.
        What comes out are millions of short DNA snippets called
        <strong style="color:var(--cyan)">reads</strong> (each ~75–300 letters long).
        <strong>Genome assembly</strong> is the job of stitching all those tiny overlapping pieces
        back into the full original sequence — with some pieces containing typos from the machine,
        and huge repeated sections that look identical no matter where you are in the genome.
        HelixForge solves this with a 7-stage algorithmic pipeline.
      </p>
      <div class="problem-grid">
        <div class="prob-item">
          <div class="prob-item-title" style="color:var(--cyan)">⚡ Why Is It Hard?</div>
          <div class="prob-item-text">
            • Each read is only 75–300 letters, but a genome can be <em>billions</em> long<br>
            • The sequencing machine makes ~1% mistakes — so the data has typos<br>
            • Some DNA sections repeat identically 100s of times — hard to place correctly<br>
            • The machine reads both strands of DNA at once, so you get mirrored copies<br>
            • A human genome produces ~1 TB of raw reads to process
          </div>
        </div>
        <div class="prob-item">
          <div class="prob-item-title" style="color:var(--purple)">🧩 How HelixForge Solves It</div>
          <div class="prob-item-text">
            We chop every read into tiny fixed-size chunks called <strong>k-mers</strong>
            (e.g. if k=21, each chunk is 21 letters). Then we build a
            <strong>De Bruijn graph</strong> — a map where each chunk connects to
            the next one that overlaps it. The genome sequence is hiding inside this graph
            as a path that visits every connection exactly once.
            We find that path using <strong>Hierholzer's algorithm</strong> in O(E) time.
          </div>
        </div>
        <div class="prob-item">
          <div class="prob-item-title" style="color:var(--pink)">🔬 Why Not Just Compare Everything?</div>
          <div class="prob-item-text">
            The obvious approach — compare every read to every other read to find overlaps —
            means O(N²) comparisons. With 5 million reads, that's
            <strong>25 trillion comparisons</strong>. Way too slow.
            Our approach converts the problem into a graph traversal, bringing it down to
            <strong>O(N)</strong> to build the graph and <strong>O(E)</strong> to find the answer.
            That's the difference between hours and milliseconds.
          </div>
        </div>
      </div>
    </div>
    """, unsafe_allow_html=True)

    st.markdown('<div class="hf-section">Algorithm Design Decisions</div>',
                unsafe_allow_html=True)

    decisions = [
        ("Rolling Hash — used in: k-mer extraction from every read",   "var(--cyan)",
         "<strong>The problem it solves:</strong> For every position in a DNA read, we need a fingerprint "
         "(hash) of the next k letters. If we recompute it from scratch each time, "
         "that's k operations × millions of positions = extremely slow.<br><br>"
         "<strong>What rolling hash does:</strong> It keeps the current fingerprint in memory and "
         "<em>slides</em> it one letter at a time — drop the old leftmost letter, add the new rightmost one. "
         "Each slide is just 2 math operations, so O(1) per step instead of O(k).<br><br>"
         "<strong>Where it runs:</strong> Stage 2 — every single read, every single position.",
         "O(N·k) → O(N)"),

        ("Bloom Filter — used in: filtering out sequencing errors",   "var(--purple)",
         "<strong>The problem it solves:</strong> Sequencing machines make typos (~1 per 100 letters). "
         "Those typos create k-mers that appear only once — they're fake, not real genome. "
         "If we kept all of them, the graph would be full of junk and the assembly would be wrong.<br><br>"
         "<strong>What Bloom filter does:</strong> It's a tiny 8 MB checklist. "
         "We pass over the reads twice: first pass marks a k-mer as 'seen once', "
         "second pass upgrades it to 'seen twice'. Only 'seen twice' k-mers are trusted as real genome sequence.<br><br>"
         "<strong>Where it runs:</strong> Stage 3 — before building the graph, so junk never gets in.",
         "8 MB fixed vs gigabytes of RAM"),

        ("Canonical K-mers — used in: De Bruijn graph construction",  "var(--green)",
         "<strong>The problem it solves:</strong> DNA is two-sided. A read <code>ATCG</code> and its "
         "mirror image <code>CGAT</code> are actually the same piece of genome — just read from opposite ends. "
         "Without handling this, we'd build two separate conflicting nodes for the same location.<br><br>"
         "<strong>What canonical k-mers do:</strong> For every k-mer, we also compute its reverse complement "
         "and keep whichever comes first alphabetically. Both strands collapse to one node. "
         "This cuts the graph size in half and makes the assembly biologically accurate.<br><br>"
         "<strong>Where it runs:</strong> Every time a k-mer is added to the graph (Stage 4).",
         "Graph size ÷ 2"),

        ("Dijkstra — used in: finding the best assembly path when the graph is messy", "var(--orange)",
         "<strong>The problem it solves:</strong> Sometimes the De Bruijn graph is fragmented — "
         "there's no single clean Euler path through everything. We need a backup strategy "
         "that still picks the most trustworthy route.<br><br>"
         "<strong>What Dijkstra does here:</strong> We assign each edge a weight based on how often "
         "that k-mer was seen in the reads. High coverage = low cost, low coverage = high cost. "
         "Dijkstra's min-heap naturally follows the path with the most read support — "
         "essentially finding the most confident assembly backbone.<br><br>"
         "<strong>Where it runs:</strong> Stage 5 — runs in parallel with Hierholzer. "
         "The longer result wins.",
         "O((V+E) log V)"),

        ("Suffix Array + LCP — used in: finding repeat regions in the assembled genome", "var(--pink)",
         "<strong>The problem it solves:</strong> Some DNA sections appear many times — "
         "these repeats confuse assemblers and are medically important (e.g. tandem repeat diseases). "
         "After assembly, we want to report exactly where repeats are and how long they are.<br><br>"
         "<strong>What Suffix Array does:</strong> Take the assembled sequence and sort every possible "
         "suffix (every possible starting position to the end) alphabetically. "
         "Suffixes that start with the same letters will sit next to each other. "
         "The LCP (Longest Common Prefix) array then records how many letters each neighbour shares — "
         "a long match = a repeat region.<br><br>"
         "<strong>Where it runs:</strong> Stage 7 — after assembly is complete, purely for analysis.",
         "O(n log²n) build + O(n) repeat scan"),
    ]

    for title, color, desc, complexity in decisions:
        with st.expander(f"◈  {title}"):
            st.markdown(
                f'<div style="font-family:Inter,sans-serif;font-size:.84rem;'
                f'color:var(--text);line-height:1.75;border-left:3px solid {color};'
                f'padding-left:.9rem">{desc}</div>'
                f'<div style="margin-top:.65rem;font-family:Space Mono,monospace;font-size:.68rem;'
                f'color:{color};background:rgba(255,255,255,.04);border-radius:5px;padding:.3rem .7rem;'
                f'display:inline-block">complexity: {complexity}</div>',
                unsafe_allow_html=True)

    st.markdown('<div class="hf-section">Real-World Impact</div>', unsafe_allow_html=True)
    st.markdown("""
    <div class="problem-card" style="border-color:rgba(0,255,153,.2)">
      <div class="problem-title" style="color:var(--green)">🌍 The Same Algorithms Power Real Science</div>
      <div class="problem-grid">
        <div class="prob-item">
          <div class="prob-item-title" style="color:var(--green)">🏥 Fighting Diseases</div>
          <div class="prob-item-text">
            Doctors use genome assembly to find mutations that cause cancer or rare diseases.
            Tools like GATK and SPAdes — used in hospitals worldwide — are built on exactly the
            same De Bruijn graph idea that HelixForge implements.
          </div>
        </div>
        <div class="prob-item">
          <div class="prob-item-title" style="color:var(--cyan)">🦠 COVID-19 Tracking</div>
          <div class="prob-item-text">
            Every time a new COVID variant was detected, labs had to assemble the viral genome
            from scratch in hours. The rolling hash + De Bruijn pipeline used here is the same
            approach those labs ran thousands of times per day during the pandemic.
          </div>
        </div>
        <div class="prob-item">
          <div class="prob-item-title" style="color:var(--purple)">🌿 Mapping All Life</div>
          <div class="prob-item-text">
            The Earth BioGenome Project wants to sequence every species on Earth (~1.5 million).
            Without fast assembly algorithms like ours, that would take centuries.
            Efficient k-mer hashing and graph traversal make it possible in years.
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
        ("📂 FASTQ Reader",        "O(N)",                   "Read the input file once, line by line. N = total letters in all reads.",                              "Used here: loading raw sequencing data"),
        ("🔄 Reverse Complement",  "O(k)",                   "Flip a k-mer and swap A↔T, C↔G. Needed because both DNA strands get sequenced.",                      "Used here: making canonical k-mers"),
        ("⚡ Rolling Hash (setup)","O(k)",                   "Compute the first hash of a k-letter window — done once per read at the start.",                       "Used here: k-mer fingerprinting"),
        ("⚡ Rolling Hash (slide)","O(1)",                   "Update the hash by dropping one letter and adding another. Much faster than recomputing from scratch.", "Used here: sliding across every read position"),
        ("🌸 Bloom Filter insert", "O(1)",                   "Flip a few bits in an 8 MB checklist. Constant time regardless of how many k-mers we've seen.",        "Used here: tracking which k-mers appeared twice"),
        ("🌸 Bloom Filter query",  "O(1)",                   "Check if those bits are set. If any bit is 0, the k-mer is definitely not there.",                     "Used here: filtering error k-mers before graph build"),
        ("🧬 Graph add_edge",      "O(1) avg",               "Insert a k-mer as an edge into a hash map — average O(1), very fast in practice.",                     "Used here: building the De Bruijn graph"),
        ("🧬 Graph Construction",  "O(V + E)",               "Build the full graph — V unique nodes (k-1-mers) and E unique edges (k-mers).",                        "Used here: Stage 4, after Bloom filtering"),
        ("🎯 Dijkstra SSSP",       "O((V+E) log V)",         "Min-heap shortest path, but here 'shortest' means most-covered. Picks the most trusted assembly path.","Used here: Stage 5, backup assembly when no Euler path"),
        ("🛤 Hierholzer Euler",    "O(E)",                   "Walk every edge exactly once using a stack. If the graph allows it, this gives the full genome path.",  "Used here: Stage 6, primary assembly strategy"),
        ("🔀 Greedy Fallback",     "O(E)",                   "If neither method gives a complete path, greedily follow the highest-frequency edge at each step.",     "Used here: last resort when graph is fragmented"),
        ("🩹 DP Error Correction", "O(N)",                   "Scan the assembled sequence and fix unlikely letters using a sliding window majority vote.",            "Used here: cleaning up the final sequence"),
        ("🔬 Suffix Array build",  "O(n log² n)",            "Sort all suffixes of the assembly alphabetically. Doubling trick means we only need log(n) rounds.",    "Used here: Stage 7, repeat region detection"),
        ("🔬 LCP Array (Kasai)",   "O(n)",                   "Find how many letters each pair of adjacent sorted suffixes share. Uses a clever invariant to stay O(n).","Used here: tells us where repeats start and end"),
        ("🔍 Repeat Finding",      "O(n)",                   "Scan the LCP array — neighbouring entries with a long common prefix are repeat regions.",               "Used here: reporting tandem repeats in the output"),
        ("📊 GC Content",          "O(n)",                   "Count G and C letters in the assembly, divide by total length. One simple pass.",                       "Used here: genome quality metric"),
        ("📊 N50 Metric",          "O(c log c)",             "Sort contig lengths, sum from longest until you hit 50% of total — that length is N50.",               "Used here: standard assembly quality score"),
        ("💾 JSON Graph Writer",   "O(V + E)",               "Write graph nodes and edges to a file for the 3D visualiser. Capped at 500 nodes for browser speed.",  "Used here: feeding the interactive graph display"),
        ("🏁 Overall Pipeline",    "O(N + (V+E)logV + n log²n)", "Dominated by Dijkstra on graph-heavy data, or by the Suffix Array on long assemblies.",           "All 7 stages combined"),
    ]

    df_full = pd.DataFrame(rows, columns=["Step","Complexity","Plain English — what it does & why","Where it's used"])
    st.dataframe(df_full, use_container_width=True, hide_index=True, height=600)

    st.markdown('<div class="hf-section">Memory Usage (Space Complexity)</div>', unsafe_allow_html=True)
    space_rows = [
        ("Bloom Filters (×2)",  "16 MB fixed",      "Two tiny 8 MB checklists — same size no matter how many reads you process. This is the key trick that makes HelixForge memory-efficient."),
        ("Rolling Hash",        "O(1)",              "Just two numbers in memory: the current hash value and a precomputed power. Doesn't grow with input size at all."),
        ("De Bruijn Graph",     "O(V + E)",          "One entry per unique k-mer (edge) and per unique (k-1)-mer (node). Grows with the genome, not the raw read count."),
        ("Dijkstra structures", "O(V)",              "Arrays for distances and a min-heap — one slot per node. Manageable even for large graphs."),
        ("Suffix Array",        "O(n)",              "Three arrays of length n (the assembly length). Longer assembly = more memory, but it's proportional and predictable."),
        ("LCP Array",           "O(n)",              "One number per suffix — exactly n integers for an assembly of length n."),
        ("Overall",             "O(N + V + E + n)", "The Bloom filter constant is the hero here — without it, storing all k-mers would need gigabytes instead of 16 MB."),
    ]
    df_space = pd.DataFrame(space_rows, columns=["What","Memory","Why it's this size"])
    st.dataframe(df_space, use_container_width=True, hide_index=True)

    st.markdown('<div class="hf-section">Algorithm & Data Structure Topics Used</div>', unsafe_allow_html=True)
    aps_rows = [
        ("Graph Theory",         "We model the genome as a De Bruijn graph — directed edges are k-mers, nodes are (k-1)-mers. Assembly = path through graph."),
        ("Eulerian Paths",       "Hierholzer's algorithm finds a path that uses every edge exactly once. This is exactly what genome assembly needs. O(E) time."),
        ("Shortest Path (SSSP)", "Dijkstra finds the path through edges with highest read coverage — our backup assembly strategy when no Euler path exists."),
        ("Priority Queue",       "Min-heap powers Dijkstra — always processes the lowest-cost node next. O(log V) per operation."),
        ("Hashing",              "Rabin-Karp rolling hash gives O(1) k-mer fingerprinting per position instead of O(k). Used on every single read."),
        ("Probabilistic DS",     "Bloom Filter — a space-efficient way to check 'have I seen this k-mer twice?' Uses only 8 MB regardless of data size."),
        ("Divide & Conquer",     "Suffix Array is built by repeatedly halving the problem — rank by 1 letter, then 2, then 4, then 8... until all suffixes are sorted. Each round doubles how much we know."),
        ("String Algorithms",    "Suffix Array sorts all suffixes; LCP Array finds shared prefixes between them; Kasai's algorithm builds LCP in O(n) by cleverly reusing work from the previous suffix."),
        ("Dynamic Programming",  "Error correction scans the assembly and fixes each letter based on what the surrounding window says — each fix depends only on the current position, so it's O(1) per step."),
        ("Amortised Analysis",   "Kasai's LCP algorithm seems like it could be slow, but a counter can only go down n times total — so even though individual steps vary, the total is always O(n)."),
        ("Greedy Algorithms",    "When both Dijkstra and Hierholzer fail to give a full path, we fall back to greedily always picking the most-seen edge next. Simple but effective."),
        ("Sliding Window",       "Rolling hash slides a fixed-size window across the read — add one letter on the right, drop one on the left. Same idea used in DP correction to check local base quality."),
        ("Two-Pointer",          "Kasai's LCP extension works like two pointers — one extending forward, one resetting back. Their total movements stay within n, keeping the whole thing O(n)."),
        ("Space-Time Tradeoff",  "Bloom filter accepts a tiny ~0.02% chance of a false positive in exchange for using 8 MB instead of gigabytes. We chose speed and memory over perfect accuracy — and it works."),
        ("Canonical Forms",      "For every k-mer, we pick min(kmer, reverse_complement) as the 'official' version. This groups both DNA strands under one name, cutting the graph size in half."),
        ("Complexity Analysis",  "Every module in HelixForge is measured both theoretically (what it should take) and in practice (how long it actually ran). Both are recorded in stats.txt for comparison."),
    ]
    df_aps = pd.DataFrame(aps_rows, columns=["Topic","How HelixForge uses it — in plain English"])
    st.dataframe(df_aps, use_container_width=True, hide_index=True, height=530)
