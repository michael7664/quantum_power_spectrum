#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════
#  setup.sh — One-time Python environment setup for tng_data_prep
#  Run once before first use: bash setup.sh
# ═══════════════════════════════════════════════════════════════════
set -e

RED='\033[0;31m'; GREEN='\033[0;32m'
YELLOW='\033[1;33m'; CYAN='\033[0;36m'; NC='\033[0m'
info()    { echo -e "${CYAN}[INFO]${NC}  $1"; }
success() { echo -e "${GREEN}[OK]${NC}    $1"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $1"; }
error()   { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

echo ""
echo "══════════════════════════════════════════════"
echo "  TNG Data Prep — Setup"
echo "══════════════════════════════════════════════"
echo ""

# ── Check .env exists ────────────────────────────────────────────────
if [ ! -f ".env" ]; then
    if [ -f ".env.example" ]; then
        cp .env.example .env
        warn "Created .env from .env.example"
        warn "Edit .env and set your TNG_API_KEY, then re-run: bash setup.sh"
        exit 0
    else
        error ".env not found. Create it manually."
    fi
fi
success ".env found"

# ── Python check ─────────────────────────────────────────────────────
if ! command -v python3 &>/dev/null; then
    error "python3 not found. Install Python 3.9+ first."
fi
PY_VER=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
info "Python version: $PY_VER"

# ── Create virtual environment ───────────────────────────────────────
if [ ! -d "venv" ]; then
    info "Creating Python virtual environment..."
    python3 -m venv venv
    success "Created venv/"
else
    success "venv/ already exists"
fi

# ── Install packages ─────────────────────────────────────────────────
info "Installing Python packages..."
source venv/bin/activate
pip install --upgrade pip --quiet
pip install requests h5py numpy astropy astroquery --quiet
success "Installed: requests, h5py, numpy, astropy, astroquery"

# ── Create output dir in C++ project ─────────────────────────────────
source .env 2>/dev/null || true
OUT_DIR="${OUTPUT_DIR:-../quantum_power_spectrum/data}"
mkdir -p "$OUT_DIR"
success "Output directory ready: $OUT_DIR"

echo ""
echo "══════════════════════════════════════════════"
echo -e "  ${GREEN}Setup complete!${NC}"
echo "  Next: edit .env with your TNG_API_KEY"
echo "  Then: bash run.sh"
echo "══════════════════════════════════════════════"
echo ""
