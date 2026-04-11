"""
HelixForge — Genome Assembler  |  Redesigned UI
Run:  streamlit run app.py
"""

import json, subprocess, tempfile, shutil
from pathlib import Path

import streamlit as st
import networkx as nx
import plotly.graph_objects as go
import pandas as pd

st.set_page_config(
    page_title="HelixForge",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

BASE_DIR   = Path(__file__).parent
ASSEMBLER  = BASE_DIR / ("assembler.exe" if not (BASE_DIR / "assembler").exists() else "assembler")
OUTPUT_DIR = BASE_DIR / "output"

STYLE = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Mono:ital,wght@0,400;0,700;1,400&family=Syne:wght@400;600;700;800&display=swap');
html,[class*="css"]{font-family:'Syne',sans-serif}
.stApp{background:#050a0f;background-image:radial-gradient(ellipse 80% 50% at 20% 10%,rgba(0,200,160,0.06) 0%,transparent 60%),radial-gradient(ellipse 60% 40% at 80% 80%,rgba(120,60,255,0.07) 0%,transparent 60%)}
#MainMenu,footer,header{visibility:hidden}
.block-container{padding:2rem 2.5rem 4rem;max-width:1400px}
[data-testid="stSidebar"]{background:#080e14!important;border-right:1px solid rgba(0,200,160,0.15)!important}
[data-testid="stSidebar"] *{color:#a8c4bc!important}
[data-testid="stSidebar"] h1,[data-testid="stSidebar"] h2,[data-testid="stSidebar"] h3{color:#00c8a0!important}
h1,h2,h3,h4,p,label,span,div{color:#d4e8e0}
.stButton>button{background:linear-gradient(135deg,#00c8a0 0%,#00a882 100%)!important;color:#050a0f!important;font-family:'Syne',sans-serif!important;font-weight:700!important;font-size:0.95rem!important;letter-spacing:0.08em!important;text-transform:uppercase!important;border:none!important;border-radius:4px!important;padding:0.65rem 1.5rem!important;transition:all 0.2s!important;box-shadow:0 0 20px rgba(0,200,160,0.3)!important}
.stButton>button:hover{box-shadow:0 0 35px rgba(0,200,160,0.55)!important;transform:translateY(-1px)!important}
[data-testid="stFileUploader"]{background:rgba(0,200,160,0.04)!important;border:1px dashed rgba(0,200,160,0.3)!important;border-radius:6px!important}
.stSlider>div>div>div{background:#00c8a0!important}
.stSlider>div>div{background:rgba(0,200,160,0.15)!important}
.streamlit-expanderHeader{background:rgba(0,200,160,0.06)!important;border:1px solid rgba(0,200,160,0.2)!important;border-radius:4px!important;color:#00c8a0!important;font-family:'Space Mono',monospace!important;font-size:0.82rem!important}
.streamlit-expanderContent{background:#060e0b!important;border:1px solid rgba(0,200,160,0.1)!important;border-top:none!important;font-size:0.88rem!important;line-height:1.7!important}
[data-testid="stExpander"]{background:#060e0b!important;border:1px solid rgba(0,200,160,0.2)!important;border-radius:4px!important}
[data-testid="stExpander"] details,[data-testid="stExpander"] summary{background:#060e0b!important}
[data-testid="stExpander"] summary{color:#00c8a0!important}
div[data-testid="stExpander"]>div{background:#060e0b!important}
pre{background:#030a07!important;color:#00ffb3!important;font-family:'Space Mono',monospace!important;font-size:0.76rem!important}
.stCodeBlock,.stCodeBlock *{background:#030a07!important;color:#00ffb3!important}
code{background:#030a07!important;color:#00ffb3!important}
.stDataFrame{border:1px solid rgba(0,200,160,0.2)!important;border-radius:6px!important}
.stDataFrame thead th{background:#0d1f1a!important;color:#00c8a0!important;font-family:'Space Mono',monospace!important;font-size:0.78rem!important}
.stDataFrame tbody td{background:#08120e!important;color:#a8c4bc!important;font-family:'Space Mono',monospace!important;font-size:0.82rem!important}
.stSpinner>div{border-top-color:#00c8a0!important}
[data-testid="stDownloadButton"]>button{background:transparent!important;color:#00c8a0!important;border:1px solid rgba(0,200,160,0.4)!important;font-family:'Space Mono',monospace!important;font-size:0.78rem!important;border-radius:4px!important;transition:all 0.2s!important}
[data-testid="stDownloadButton"]>button:hover{background:rgba(0,200,160,0.1)!important;border-color:#00c8a0!important}
.hf-hero{padding:2.5rem 0 1.5rem;margin-bottom:0.5rem}
.hf-logo{font-family:'Syne',sans-serif;font-weight:800;font-size:3.2rem;letter-spacing:-0.02em;background:linear-gradient(120deg,#00c8a0 0%,#00ffcc 40%,#7c3aff 100%);-webkit-background-clip:text;-webkit-text-fill-color:transparent;background-clip:text;line-height:1}
.hf-tagline{font-family:'Space Mono',monospace;font-size:0.8rem;color:rgba(0,200,160,0.6);letter-spacing:0.18em;text-transform:uppercase;margin-top:0.5rem}
.hf-divider{height:1px;background:linear-gradient(90deg,#00c8a0 0%,rgba(120,60,255,0.4) 50%,transparent 100%);margin:1.5rem 0;opacity:0.4}
.hf-section{font-family:'Space Mono',monospace;font-size:0.72rem;letter-spacing:0.22em;text-transform:uppercase;color:#00c8a0;margin:2rem 0 1rem;display:flex;align-items:center;gap:0.7rem}
.hf-section::after{content:'';flex:1;height:1px;background:linear-gradient(90deg,rgba(0,200,160,0.3),transparent)}
.kpi-grid{display:grid;grid-template-columns:repeat(5,1fr);gap:1rem;margin:1rem 0}
.kpi-card{background:rgba(0,200,160,0.04);border:1px solid rgba(0,200,160,0.15);border-radius:6px;padding:1.1rem 1rem 0.9rem;position:relative;overflow:hidden;transition:border-color 0.2s}
.kpi-card::before{content:'';position:absolute;top:0;left:0;right:0;height:2px;background:linear-gradient(90deg,#00c8a0,#7c3aff);opacity:0.7}
.kpi-card:hover{border-color:rgba(0,200,160,0.35)}
.kpi-label{font-family:'Space Mono',monospace;font-size:0.68rem;letter-spacing:0.14em;text-transform:uppercase;color:rgba(0,200,160,0.55);margin-bottom:0.4rem}
.kpi-value{font-family:'Syne',sans-serif;font-weight:700;font-size:1.6rem;color:#e8f5f0;line-height:1.1}
.kpi-unit{font-family:'Space Mono',monospace;font-size:0.7rem;color:rgba(0,200,160,0.4);margin-top:0.15rem}
.seq-terminal{background:#030a07;border:1px solid rgba(0,200,160,0.2);border-radius:6px;overflow:hidden}
.seq-terminal-bar{background:#0d1f1a;padding:0.5rem 1rem;display:flex;align-items:center;gap:0.5rem;border-bottom:1px solid rgba(0,200,160,0.15)}
.seq-dot{width:10px;height:10px;border-radius:50%}
.seq-content{font-family:'Space Mono',monospace;font-size:0.75rem;line-height:1.8;color:#00ffb3;padding:1rem;max-height:180px;overflow-y:auto;word-break:break-all;white-space:pre-wrap;letter-spacing:0.05em}
.seq-content::-webkit-scrollbar{width:4px}
.seq-content::-webkit-scrollbar-thumb{background:rgba(0,200,160,0.3);border-radius:2px}
.overall-badge{display:inline-flex;align-items:center;gap:0.6rem;background:rgba(0,200,160,0.07);border:1px solid rgba(0,200,160,0.25);border-radius:4px;padding:0.6rem 1rem;margin-top:0.8rem}
.overall-label{font-family:'Space Mono',monospace;font-size:0.72rem;color:rgba(0,200,160,0.6);letter-spacing:0.1em;text-transform:uppercase}
.overall-value{font-family:'Space Mono',monospace;font-size:0.95rem;font-weight:700;color:#00ffcc}
.sidebar-logo{font-family:'Syne',sans-serif;font-weight:800;font-size:1.3rem;background:linear-gradient(120deg,#00c8a0,#7c3aff);-webkit-background-clip:text;-webkit-text-fill-color:transparent;background-clip:text;padding:0.5rem 0 1rem}
.sidebar-section{font-family:'Space Mono',monospace;font-size:0.65rem;letter-spacing:0.2em;text-transform:uppercase;color:rgba(0,200,160,0.5)!important;margin:1.2rem 0 0.5rem}
.info-banner{background:rgba(0,200,160,0.06);border:1px solid rgba(0,200,160,0.2);border-left:3px solid #00c8a0;border-radius:0 6px 6px 0;padding:1rem 1.2rem;font-family:'Space Mono',monospace;font-size:0.82rem;color:#a8c4bc;margin:1rem 0}
.pipeline-steps{display:grid;grid-template-columns:repeat(5,1fr);gap:0;margin:1.5rem 0}
.pipeline-step{background:rgba(0,200,160,0.04);border:1px solid rgba(0,200,160,0.15);border-right:none;padding:1.4rem 1rem 1.2rem;text-align:center;position:relative;transition:background 0.2s}
.pipeline-step:last-child{border-right:1px solid rgba(0,200,160,0.15);border-radius:0 6px 6px 0}
.pipeline-step:first-child{border-radius:6px 0 0 6px}
.pipeline-step:hover{background:rgba(0,200,160,0.08)}
.pipeline-step::after{content:'›';position:absolute;right:-10px;top:50%;transform:translateY(-50%);color:rgba(0,200,160,0.4);font-size:1.4rem;z-index:2}
.pipeline-step:last-child::after{display:none}
.step-num{font-family:'Space Mono',monospace;font-size:0.65rem;color:rgba(0,200,160,0.4);letter-spacing:0.1em;margin-bottom:0.3rem}
.step-icon{font-size:1.5rem;margin-bottom:0.4rem;line-height:1}
.step-title{font-weight:700;font-size:0.85rem;color:#d4e8e0;margin-bottom:0.2rem}
.step-sub{font-family:'Space Mono',monospace;font-size:0.65rem;color:rgba(0,200,160,0.5)}
</style>
"""
st.markdown(STYLE, unsafe_allow_html=True)

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown('<div class="sidebar-logo">⬡ HelixForge</div>', unsafe_allow_html=True)
    st.markdown('<div class="sidebar-section">Input</div>', unsafe_allow_html=True)
    uploaded = st.file_uploader("FASTQ file", type=["fastq","fq","txt"],
        help="Standard 4-line FASTQ. Supports 100 MB+.", label_visibility="collapsed")
    if not uploaded:
        st.markdown('<div style="font-family:Space Mono,monospace;font-size:0.72rem;color:rgba(0,200,160,0.4);margin-top:0.3rem">Drop .fastq / .fq / .txt</div>', unsafe_allow_html=True)
    st.markdown('<div class="sidebar-section">Parameters</div>', unsafe_allow_html=True)
    k = st.slider("k-mer size", 11, 63, 21, 2)
    st.markdown(f'<div style="font-family:Space Mono,monospace;font-size:0.7rem;color:rgba(0,200,160,0.45);margin-top:-0.8rem;margin-bottom:0.8rem">k={k} · nodes={k-1}bp</div>', unsafe_allow_html=True)
    max_vis = st.slider("Graph nodes (viz)", 50, 500, 200, 50)
    st.markdown("")
    run_btn = st.button("▶  RUN ASSEMBLY", type="primary", use_container_width=True)
    st.markdown('<div class="sidebar-section">DAA Notes</div>', unsafe_allow_html=True)
    with st.expander("Rolling hash → O(N)"):
        st.write("Rabin-Karp slides in O(1) per step vs O(k) naive. Total: O(N).")
    with st.expander("De Bruijn vs O(N²)"):
        st.write("Naive overlap-layout is O(N²). De Bruijn encodes overlaps as edges; assembly = Eulerian path O(E).")
    with st.expander("Bloom Filter memory"):
        st.write("Fixed 8 MB vs O(N·k) hash map. Two-pass counting removes error k-mers.")
    with st.expander("Why Eulerian path?"):
        st.write("Hamiltonian is NP-hard. Eulerian is O(E) via Hierholzer. De Bruijn graphs are Eulerian by construction.")

# ── Helpers ───────────────────────────────────────────────────────────────────
def parse_stats(p):
    data, section = {}, None
    for line in p.read_text().splitlines():
        if not line.strip(): continue
        if "Time Complexity (Theoretical):" in line:
            section="theoretical"; data["theoretical"]=[]; continue
        if "Execution Time (Measured):" in line:
            section="measured"; data["measured"]={}; continue
        if section=="theoretical":
            if line.startswith("  -"):
                parts=line.strip("- ").split(":",1)
                if len(parts)>=2:
                    step=parts[0].strip(); rest=parts[1].strip()
                    complexity=rest.split()[0]
                    reason=rest.split("-",1)[1].strip() if "-" in rest else ""
                    data["theoretical"].append({"step":step,"complexity":complexity,"reason":reason})
            elif "Overall" in line:
                data["overall"]=line.split(":",1)[1].strip()
        elif section=="measured":
            if ":" in line:
                k_,v=line.split(":",1); data["measured"][k_.strip()]=v.strip()
        else:
            if ":" in line:
                k_,v=line.split(":",1); data[k_.strip()]=v.strip()
    return data

def build_3d_graph(graph_path, max_nodes):
    gdata=json.loads(graph_path.read_text())
    nodes=gdata["nodes"][:max_nodes]; node_set=set(nodes)
    G=nx.DiGraph(); G.add_nodes_from(nodes)
    for u,v in gdata["edges"]:
        if u in node_set and v in node_set: G.add_edge(u,v)
    pos=nx.spring_layout(G,dim=3,seed=42,k=0.5)
    ex,ey,ez=[],[],[]
    for u,v in G.edges():
        x0,y0,z0=pos[u]; x1,y1,z1=pos[v]
        ex+=[x0,x1,None]; ey+=[y0,y1,None]; ez+=[z0,z1,None]
    degs=[G.degree(n) for n in nodes]
    fig=go.Figure(data=[
        go.Scatter3d(x=ex,y=ey,z=ez,mode="lines",
            line=dict(color="rgba(0,200,160,0.22)",width=1),hoverinfo="none",name="k-mer edges"),
        go.Scatter3d(
            x=[pos[n][0] for n in nodes],y=[pos[n][1] for n in nodes],z=[pos[n][2] for n in nodes],
            mode="markers",
            marker=dict(size=3.5,color=degs,
                colorscale=[[0,"#001a10"],[0.3,"#00c8a0"],[0.7,"#7c3aff"],[1,"#ff6af0"]],
                colorbar=dict(
                    title=dict(text="Degree",font=dict(color="#a8c4bc",size=9)),
                    thickness=10,len=0.5,
                    tickfont=dict(color="#a8c4bc",size=9)),
                opacity=0.9,line=dict(width=0)),
            text=nodes,
            hovertemplate="<b>%{text}</b><br>degree: %{marker.color}<extra></extra>",
            name="(k-1)-mer nodes")
    ])
    fig.update_layout(
        scene=dict(
            xaxis=dict(showgrid=False,zeroline=False,showticklabels=False,showbackground=False),
            yaxis=dict(showgrid=False,zeroline=False,showticklabels=False,showbackground=False),
            zaxis=dict(showgrid=False,zeroline=False,showticklabels=False,showbackground=False),
            bgcolor="rgba(0,0,0,0)"),
        legend=dict(x=0.01,y=0.99,font=dict(size=10,color="#a8c4bc"),bgcolor="rgba(0,0,0,0)"),
        margin=dict(l=0,r=0,t=10,b=0),height=520,
        paper_bgcolor="rgba(0,0,0,0)",
        hoverlabel=dict(bgcolor="#0d1f1a",font_size=11,font_family="Space Mono"))
    return fig

def kpi(label,value,unit=""):
    return f'<div class="kpi-card"><div class="kpi-label">{label}</div><div class="kpi-value">{value}</div><div class="kpi-unit">{unit}</div></div>'

# ── Header ────────────────────────────────────────────────────────────────────
st.markdown("""
<div class="hf-hero">
  <div class="hf-logo">HelixForge</div>
  <div class="hf-tagline">De-novo genome assembler · Rolling Hash · Bloom Filter · Hierholzer · DP Correction</div>
</div>
<div class="hf-divider"></div>
""", unsafe_allow_html=True)

# ── Run ───────────────────────────────────────────────────────────────────────
if run_btn:
    if uploaded is None:
        st.error("No FASTQ file uploaded."); st.stop()
    if not ASSEMBLER.exists():
        st.error("Assembler binary not found.\n```bash\ng++ -O2 -std=c++17 -o assembler assembler_standalone.cpp\n```"); st.stop()
    tmp=Path(tempfile.mkdtemp()); fq=tmp/"input.fastq"; fq.write_bytes(uploaded.read())
    out=tmp/"output"; out.mkdir()
    with st.spinner("Assembling genome — running C++ DAA engine…"):
        result=subprocess.run([str(ASSEMBLER),str(fq),str(k),str(out)],capture_output=True,text=True,timeout=300)
    if result.returncode!=0:
        st.error("Assembly failed"); st.code(result.stderr); shutil.rmtree(tmp,ignore_errors=True); st.stop()
    OUTPUT_DIR.mkdir(exist_ok=True)
    for f in ["genome.fasta","stats.txt","graph_data.json"]:
        src=out/f
        if src.exists(): shutil.copy(src,OUTPUT_DIR/f)
    st.session_state["assembled"]=True
    shutil.rmtree(tmp,ignore_errors=True)
    st.rerun()

# ── Results ───────────────────────────────────────────────────────────────────
if st.session_state.get("assembled"):
    fasta_p=OUTPUT_DIR/"genome.fasta"; stats_p=OUTPUT_DIR/"stats.txt"; graph_p=OUTPUT_DIR/"graph_data.json"
    stats=parse_stats(stats_p); measured=stats.get("measured",{})

    st.markdown('<div class="hf-section">Assembly metrics</div>',unsafe_allow_html=True)
    reads=stats.get("Total Reads Processed","—"); kmers=stats.get("Total k-mers","—")
    nodes_v=stats.get("Graph Nodes (V)","—"); edges_v=stats.get("Graph Edges (E)","—")
    alen=stats.get("Assembly Length","—").replace(" bp","")
    st.markdown(f'<div class="kpi-grid">{kpi("Reads",reads,"reads")}{kpi("k-mers",kmers,"k-mers")}{kpi("Graph nodes",nodes_v,"(k-1)-mers")}{kpi("Graph edges",edges_v,"k-mers")}{kpi("Assembly",alen,"bp")}</div>',unsafe_allow_html=True)

    st.markdown('<div class="hf-section">Assembled sequence</div>',unsafe_allow_html=True)
    seq_text=fasta_p.read_text()
    sequence="".join(l for l in seq_text.splitlines() if not l.startswith(">"))
    preview=sequence[:4000]+("…" if len(sequence)>4000 else "")
    col_seq,col_dl=st.columns([6,1])
    with col_seq:
        st.markdown(f'''<div class="seq-terminal">
          <div class="seq-terminal-bar">
            <div class="seq-dot" style="background:#ff5f57"></div>
            <div class="seq-dot" style="background:#febc2e"></div>
            <div class="seq-dot" style="background:#28c840"></div>
            <span style="font-family:Space Mono,monospace;font-size:0.72rem;color:rgba(0,200,160,0.5);margin-left:0.5rem">genome.fasta — {len(sequence):,} bp</span>
          </div>
          <div class="seq-content">{preview}</div>
        </div>''',unsafe_allow_html=True)
    with col_dl:
        st.markdown("<div style='height:2.6rem'></div>",unsafe_allow_html=True)
        st.download_button("⬇ FASTA",seq_text,"genome.fasta","text/plain",use_container_width=True)

    st.markdown('<div class="hf-divider"></div>',unsafe_allow_html=True)
    left,right=st.columns([1,1.1],gap="large")

    with left:
        st.markdown('<div class="hf-section">Complexity analysis</div>',unsafe_allow_html=True)
        theo=stats.get("theoretical",[])
        if theo:
            df=pd.DataFrame(theo); df.columns=["Algorithm","Complexity","Rationale"]
            st.dataframe(df,use_container_width=True,hide_index=True,height=220)
            overall=stats.get("overall","O(N + V + E)")
            st.markdown(f'<div class="overall-badge"><span class="overall-label">Overall</span><span class="overall-value">{overall}</span></div>',unsafe_allow_html=True)

        st.markdown('<div class="hf-section" style="margin-top:1.8rem">Execution timing</div>',unsafe_allow_html=True)
        tlabels=["Hashing","Graph build","Traversal","DP correction"]
        tkeys=["Hashing Time","Graph Build Time","Traversal Time","DP Time"]
        tvals=[]
        for tk in tkeys:
            v=measured.get(tk,"0 ms").replace(" ms","").strip()
            try: tvals.append(float(v))
            except: tvals.append(0.0)
        fig_bar=go.Figure(go.Bar(x=tlabels,y=tvals,
            marker=dict(color=tvals,colorscale=[[0,"#003d2a"],[0.4,"#00c8a0"],[1,"#7c3aff"]],line=dict(width=0)),
            text=[f"{v:.2f}" for v in tvals],textposition="outside",
            textfont=dict(family="Space Mono",size=10,color="#a8c4bc")))
        fig_bar.update_layout(
            yaxis=dict(title="ms",color="#a8c4bc",gridcolor="rgba(0,200,160,0.07)",tickfont=dict(family="Space Mono",size=9)),
            xaxis=dict(color="#a8c4bc",tickfont=dict(family="Space Mono",size=9)),
            height=250,margin=dict(l=10,r=10,t=20,b=10),showlegend=False,
            paper_bgcolor="rgba(0,0,0,0)",plot_bgcolor="rgba(3,15,10,0.5)")
        st.plotly_chart(fig_bar,use_container_width=True)
        c1,c2=st.columns(2)
        with c1: st.download_button("⬇ stats.txt",stats_p.read_text(),"stats.txt","text/plain",use_container_width=True)
        with c2: st.download_button("⬇ graph.json",graph_p.read_text(),"graph_data.json","application/json",use_container_width=True)
        with st.expander("Raw stats.txt"): st.code(stats_p.read_text(),language="text")

    with right:
        st.markdown('<div class="hf-section">3D De Bruijn graph</div>',unsafe_allow_html=True)
        with st.spinner("Rendering 3D graph…"):
            fig3d=build_3d_graph(graph_p,max_vis)
        st.plotly_chart(fig3d,use_container_width=True)
        gj=json.loads(graph_p.read_text())
        st.markdown(f'<div style="display:flex;gap:1rem;margin-top:-0.5rem">{kpi("Total nodes",str(gj.get("total_nodes","—")),"(k-1)-mers")}{kpi("Total edges",str(gj.get("total_edges","—")),"overlaps")}</div>',unsafe_allow_html=True)
        st.markdown('<div style="font-family:Space Mono,monospace;font-size:0.68rem;color:rgba(0,200,160,0.35);margin-top:1.2rem;text-align:center">nodes=(k-1)-mers · edges=k-mer overlaps · colour=degree · drag to rotate</div>',unsafe_allow_html=True)

# ── Landing ───────────────────────────────────────────────────────────────────
else:
    st.markdown('<div class="info-banner">Upload a <code>.fastq</code> file in the sidebar, set your k-mer size, then click <strong>RUN ASSEMBLY</strong>.</div>',unsafe_allow_html=True)
    st.markdown('<div class="hf-section">Pipeline</div>',unsafe_allow_html=True)
    st.markdown("""<div class="pipeline-steps">
      <div class="pipeline-step"><div class="step-num">01</div><div class="step-icon">📂</div><div class="step-title">FASTQ Read</div><div class="step-sub">Streaming · O(N)</div></div>
      <div class="pipeline-step"><div class="step-num">02</div><div class="step-icon">⚡</div><div class="step-title">Rolling Hash</div><div class="step-sub">Rabin-Karp · O(1)</div></div>
      <div class="pipeline-step"><div class="step-num">03</div><div class="step-icon">🌸</div><div class="step-title">Bloom Filter</div><div class="step-sub">8 MB · O(N)</div></div>
      <div class="pipeline-step"><div class="step-num">04</div><div class="step-icon">🕸</div><div class="step-title">De Bruijn</div><div class="step-sub">Adj. list · O(V+E)</div></div>
      <div class="pipeline-step"><div class="step-num">05</div><div class="step-icon">🛤</div><div class="step-title">Hierholzer</div><div class="step-sub">Eulerian · O(E)</div></div>
    </div>""",unsafe_allow_html=True)
    st.markdown('<div class="hf-section" style="margin-top:2rem">Complexity summary</div>',unsafe_allow_html=True)
    rows=[("Rolling Hash","O(N)","Rabin-Karp O(1) slide per k-mer"),("Bloom Filter","O(N)","Fixed bit array, constant insert/query"),("Graph Construction","O(V+E)","Adjacency list, one insertion per k-mer"),("Hierholzer","O(E)","Each edge visited exactly once"),("DP Correction","O(n·m)","LCS on sliding windows n,m≤50")]
    st.dataframe(pd.DataFrame(rows,columns=["Algorithm","Complexity","Why"]),use_container_width=True,hide_index=True)
    st.markdown('<div class="overall-badge"><span class="overall-label">Overall</span><span class="overall-value">O(N + V + E)</span></div>',unsafe_allow_html=True)
