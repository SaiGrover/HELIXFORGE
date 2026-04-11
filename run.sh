#!/usr/bin/env bash
# HelixForge — build and launch script
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  HelixForge Setup"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# 1. Check GCC
if ! command -v g++ &> /dev/null; then
    echo "❌ g++ not found. Install GCC 11+."
    exit 1
fi
GCC_VER=$(g++ --version | head -1)
echo "✅ Compiler: $GCC_VER"

# 2. Compile C++ assembler
echo ""
echo "▶ Compiling C++ assembler..."
g++ -O2 -std=c++17 -o assembler assembler_standalone.cpp
echo "✅ assembler binary ready."

# 3. Install Python deps
echo ""
echo "▶ Installing Python dependencies..."
pip install -r requirements.txt --quiet
echo "✅ Python packages installed."

# 4. Launch Streamlit
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Launching HelixForge..."
echo "  Open: http://localhost:8501"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
streamlit run app.py
