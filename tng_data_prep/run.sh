#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════
#  run.sh — Fetch TNG/SDSS data → write .bin to C++ project's data/
#
#  Usage:
#    bash run.sh                          use config from .env
#    bash run.sh --dry-run                print plan, no downloads
#    bash run.sh --skip-fetch             reuse existing .bin
#    bash run.sh --inspect <file.bin>     inspect a .bin file
# ═══════════════════════════════════════════════════════════════════
set -e

RED='\033[0;31m'; GREEN='\033[0;32m'
YELLOW='\033[1;33m'; CYAN='\033[0;36m'
BOLD='\033[1m'; NC='\033[0m'
info()    { echo -e "${CYAN}[INFO]${NC}  $1"; }
success() { echo -e "${GREEN}[OK]${NC}    $1"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $1"; }
error()   { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step()    { echo -e "\n${BOLD}── $1 ──${NC}"; }

# ── Parse arguments ───────────────────────────────────────────────────
DRY_RUN=false
SKIP_FETCH=false
INSPECT_FILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --skip-fetch)
            SKIP_FETCH=true
            shift
            ;;
        --inspect)
            if [[ -z "$2" || "$2" == --* ]]; then
                error "--inspect requires a file path argument"
            fi
            INSPECT_FILE="$2"
            shift 2
            ;;
        *)
            warn "Unknown argument: $1"
            shift
            ;;
    esac
done

echo ""
echo "══════════════════════════════════════════════"
echo "  TNG Data Prep — Run"
echo "══════════════════════════════════════════════"

# ── Load .env safely ─────────────────────────────────────────────────
[ ! -f ".env" ] && error ".env not found. Run: bash setup.sh"

while IFS= read -r line; do
    # Skip blank lines and comment lines
    [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue

    # Split only on FIRST '=' to preserve '=' chars inside values
    key="${line%%=*}"
    value="${line#*=}"

    # Strip inline comments from value (everything after ' #' or '#')
    value="${value%%#*}"

    # Strip leading/trailing whitespace from both key and value
    key="${key#"${key%%[![:space:]]*}"}"
    key="${key%"${key##*[![:space:]]}"}"
    value="${value#"${value%%[![:space:]]*}"}"
    value="${value%"${value##*[![:space:]]}"}"

    # Skip lines with invalid key names
    [[ -z "$key" || "$key" =~ [^a-zA-Z0-9_] ]] && continue

    export "$key=$value"
done < .env

success "Loaded .env"

# ── Debug: confirm key loaded (shows first 6 chars only) ─────────────
if [ -z "$TNG_API_KEY" ]; then
    error "TNG_API_KEY failed to load from .env. Check the file format."
fi
info "TNG_API_KEY loaded: ${TNG_API_KEY:0:6}... (${#TNG_API_KEY} chars total)"


# ── Activate venv ────────────────────────────────────────────────────
[ ! -d "venv" ] && error "venv not found. Run: bash setup.sh"
source venv/bin/activate
success "Python venv activated"

# ── Inspect mode (works independently of fetch) ──────────────────────
if [ -n "$INSPECT_FILE" ]; then
    step "Inspecting Catalog"
    if [ ! -f "$INSPECT_FILE" ]; then
        error "File not found: $INSPECT_FILE"
    fi
    python fetch_tng_cluster.py --inspect "$INSPECT_FILE"
    exit 0
fi

# ── Validate API key ─────────────────────────────────────────────────
if [[ "$TNG_API_KEY" == "your_api_key_here"             || \
      "$TNG_API_KEY" == "your_regenerated_api_key_here" || \
      "$TNG_API_KEY" == "your_new_key_here"             || \
      ${#TNG_API_KEY} -lt 8                             || \
      -z "$TNG_API_KEY" ]]; then
    if [[ "$DATA_SOURCE" != "custom" ]]; then
        error "TNG_API_KEY not set in .env. Edit .env and add your real key."
    fi
fi


# ── Setup output dir ─────────────────────────────────────────────────
OUT_DIR="${OUTPUT_DIR:-../quantum_power_spectrum/data}"
mkdir -p "$OUT_DIR"

# ── Print plan ───────────────────────────────────────────────────────
step "Configuration"
echo "  DATA_SOURCE : $DATA_SOURCE"
echo "  OUTPUT_DIR  : $OUT_DIR"
echo "  API KEY     : ${TNG_API_KEY:0:6}... (${#TNG_API_KEY} chars)"
[ "$DATA_SOURCE" = "tng_cluster" ] && {
    echo "  SNAP        : $TNG_CLUSTER_SNAP"
    echo "  MODE        : $TNG_CLUSTER_MODE"
    echo "  MASS_MIN    : $TNG_CLUSTER_MASS_MIN Msun"
    echo "  N_CLUSTERS  : $TNG_CLUSTER_N_CLUSTERS"
    echo "  N_SUBS_MAX  : $TNG_CLUSTER_N_SUBS_MAX"
}
[ "$DRY_RUN"    = true ] && warn "DRY RUN — no data will be downloaded"
[ "$SKIP_FETCH" = true ] && warn "SKIP FETCH — will reuse existing .bin"

BIN_FILE=""

# ════════════════════════════════════════════════════════════════════
#  TNG-Cluster
# ════════════════════════════════════════════════════════════════════
if [ "$DATA_SOURCE" = "tng_cluster" ]; then
    step "Fetching TNG-Cluster Data"
    TAG="TNG-Cluster_snap${TNG_CLUSTER_SNAP}_${TNG_CLUSTER_MODE}"
    BIN_FILE="${OUT_DIR}/${TAG}.bin"

    if [ "$SKIP_FETCH" = true ] && [ -f "$BIN_FILE" ]; then
        success "Reusing existing: $BIN_FILE"
    else
        PY_ARGS=(
            --api-key       "$TNG_API_KEY"
            --snap          "$TNG_CLUSTER_SNAP"
            --mode          "$TNG_CLUSTER_MODE"
            --mass-min      "$TNG_CLUSTER_MASS_MIN"
            --sub-mass-min  "$TNG_CLUSTER_SUB_MASS_MIN"
            --N-clusters    "$TNG_CLUSTER_N_CLUSTERS"
            --N-subs-max    "$TNG_CLUSTER_N_SUBS_MAX"
            --out           "$BIN_FILE"
        )
        [ "$DRY_RUN" = true ] && PY_ARGS+=(--dry-run)

        # Show command with key partially hidden
        SAFE_CMD="python fetch_tng_cluster.py --api-key ${TNG_API_KEY:0:6}... ${PY_ARGS[*]:2}"
        info "Running: $SAFE_CMD"

        python fetch_tng_cluster.py "${PY_ARGS[@]}"
        PY_EXIT=$?
        [ $PY_EXIT -ne 0 ] && error "fetch_tng_cluster.py failed (exit $PY_EXIT)"
    fi

# ════════════════════════════════════════════════════════════════════
#  TNG100 / TNG300
# ════════════════════════════════════════════════════════════════════
elif [ "$DATA_SOURCE" = "tng100" ] || [ "$DATA_SOURCE" = "tng300" ]; then
    step "Fetching ${TNG_SIM} Data"
    TAG="${TNG_SIM}_snap${TNG_SNAP}"
    BIN_FILE="${OUT_DIR}/${TAG}.bin"

    if [ "$SKIP_FETCH" = true ] && [ -f "$BIN_FILE" ]; then
        success "Reusing existing: $BIN_FILE"
    else
        PY_ARGS=(
            --api-key   "$TNG_API_KEY"
            --sim       "$TNG_SIM"
            --snap      "$TNG_SNAP"
            --mass-min  "$TNG_MASS_MIN"
            --N-max     "$TNG_N_MAX"
            --out       "$BIN_FILE"
        )
        [ "$DRY_RUN" = true ] && PY_ARGS+=(--dry-run)

        info "Running: python fetch_tng.py --sim $TNG_SIM --snap $TNG_SNAP ..."
        python fetch_tng.py "${PY_ARGS[@]}"
        PY_EXIT=$?
        [ $PY_EXIT -ne 0 ] && error "fetch_tng.py failed (exit $PY_EXIT)"
    fi

# ════════════════════════════════════════════════════════════════════
#  SDSS
# ════════════════════════════════════════════════════════════════════
elif [ "$DATA_SOURCE" = "sdss" ]; then
    step "Fetching SDSS DR17 Data"
    TAG="SDSS_DR17_z${SDSS_Z_MIN}-${SDSS_Z_MAX}"
    BIN_FILE="${OUT_DIR}/${TAG}.bin"

    if [ "$SKIP_FETCH" = true ] && [ -f "$BIN_FILE" ]; then
        success "Reusing existing: $BIN_FILE"
    else
        PY_ARGS=(
            --z-min "$SDSS_Z_MIN"
            --z-max "$SDSS_Z_MAX"
            --N-max "$SDSS_N_MAX"
            --out   "$BIN_FILE"
        )
        info "Running: python fetch_sdss.py --z-min $SDSS_Z_MIN --z-max $SDSS_Z_MAX ..."
        python fetch_sdss.py "${PY_ARGS[@]}"
        PY_EXIT=$?
        [ $PY_EXIT -ne 0 ] && error "fetch_sdss.py failed (exit $PY_EXIT)"
    fi

# ════════════════════════════════════════════════════════════════════
#  Custom
# ════════════════════════════════════════════════════════════════════
elif [ "$DATA_SOURCE" = "custom" ]; then
    BIN_FILE="${CUSTOM_BIN_FILE:-${OUT_DIR}/custom_catalog.bin}"
    [ ! -f "$BIN_FILE" ] && \
        error "Custom .bin not found: $BIN_FILE\nRun convert_custom.py manually first."
    success "Using custom catalog: $BIN_FILE"

else
    error "Unknown DATA_SOURCE='$DATA_SOURCE' in .env"
fi

# ════════════════════════════════════════════════════════════════════
#  Final summary
# ════════════════════════════════════════════════════════════════════
step "Done"

if [ "$DRY_RUN" = false ]; then
    if [ -n "$BIN_FILE" ] && [ -f "$BIN_FILE" ]; then
        SIZE=$(du -h "$BIN_FILE" | cut -f1)
        success "Output : $BIN_FILE  ($SIZE)"
        echo ""
        info "Inspect the catalog:"
        echo "    bash run.sh --inspect $BIN_FILE"
        echo ""
        info "Run the C++ project:"
        echo "    cd ../quantum_power_spectrum"
        echo "    ./build/quantum_ps tng $BIN_FILE"
    else
        warn "Expected output not found: $BIN_FILE"
        warn "Check the Python fetch output above for errors."
    fi
fi

echo ""
echo "══════════════════════════════════════════════"
echo -e "  ${GREEN}Data preparation complete!${NC}"
echo "══════════════════════════════════════════════"
echo ""
