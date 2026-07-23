#!/usr/bin/env bash
#
# installOpenFOAM.sh
# ------------------
# Downloads and compiles ESI OpenFOAM (both DP and SP precision) for one or
# more user-specified versions, downloads a ParaView binary release, and
# rewrites every setupTemplates/*/defaultCluster/slurm/clusterDict so its
# paths point at THIS repository, the freshly built OpenFOAM and the ParaView
# install.
#
# ESI OpenFOAM compilation is only supported on Linux (native or WSL Ubuntu).
# The ParaView download and the clusterDict rewrite work on any POSIX shell.
#
# Usage:
#   ./installOpenFOAM.sh --versions "v2306 v2406" [options]
#
# See ./installOpenFOAM.sh --help for the full option list.

set -euo pipefail

# --------------------------------------------------------------------------
# Defaults
# --------------------------------------------------------------------------
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PREFIX="$HOME"                       # install root (OpenFOAM + ParaView live here)
VERSIONS=()                          # e.g. (v2306 v2406) - required for the OpenFOAM step
PRECISIONS=(DP SP)                   # which precisions to build
LABEL_SIZE=32                        # WM_LABEL_SIZE (Int32/Int64)
JOBS="$( { command -v nproc >/dev/null 2>&1 && nproc; } || echo 4)"

PARAVIEW_VERSION="5.11.1"
PARAVIEW_SERIES="v5.11"
PARAVIEW_FLAVOR="osmesa-MPI-Linux-Python3.9-x86_64"

OF_SOURCE_BASE="https://dl.openfoam.com/source"
PV_FILES_BASE="https://www.paraview.org/files"

CLUSTERDICT_PRECISION="DP"           # which build the clusterDicts source (DP is standard)
INJECT_FOAM_SOURCE=1                 # prepend a `source .../etc/bashrc` into clusterDict commands

INSTALL_DEPS=0
SKIP_OPENFOAM=0
SKIP_PARAVIEW=0
SKIP_CLUSTERDICT=0
DRY_RUN=0

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
log()  { printf '\033[1;34m[installOpenFOAM]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[installOpenFOAM] WARNING:\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31m[installOpenFOAM] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

run() {
    # Echo then execute, honouring --dry-run.
    log "\$ $*"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        "$@"
    fi
}

usage() {
    cat <<EOF
installOpenFOAM.sh - build ESI OpenFOAM (DP+SP) + ParaView and wire up clusterDicts

Required (unless --skip-openfoam):
  --versions "v2306 v2406"   Space-separated ESI OpenFOAM versions to build.
                             Each version is verified to exist before download.

Options:
  --prefix DIR               Install root for OpenFOAM and ParaView.
                             (default: \$HOME -> $HOME)
  --precision "DP SP"        Precisions to build. (default: DP SP)
  --label-size 32|64         WM_LABEL_SIZE integer width. (default: 32)
  --jobs N                   Parallel build jobs. (default: $JOBS)
  --paraview-version X.Y.Z   ParaView release. (default: $PARAVIEW_VERSION)
  --paraview-series vX.Y     ParaView series dir on paraview.org. (default: $PARAVIEW_SERIES)
  --paraview-flavor STR      ParaView binary flavor. (default: $PARAVIEW_FLAVOR)
  --clusterdict-precision DP|SP
                             Which OpenFOAM build the clusterDicts source. (default: DP)
  --no-foam-source           Do not inject a 'source .../etc/bashrc' into clusterDicts.
  --install-deps             apt-get install build prerequisites (Ubuntu/Debian, needs sudo).
  --skip-openfoam            Skip the OpenFOAM build stage.
  --skip-paraview            Skip the ParaView download stage.
  --skip-clusterdict         Skip rewriting the clusterDicts.
  --dry-run                  Print actions without executing them.
  -h, --help                 Show this help.

Examples:
  ./installOpenFOAM.sh --versions "v2406" --install-deps
  ./installOpenFOAM.sh --versions "v2306 v2606" --prefix \$HOME/openFoam --jobs 16
  ./installOpenFOAM.sh --skip-openfoam --skip-paraview   # only rewrite clusterDicts
EOF
}

# --------------------------------------------------------------------------
# Argument parsing
# --------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --versions)              read -r -a VERSIONS <<< "${2:-}"; shift 2 ;;
        --prefix)                PREFIX="${2:?}"; shift 2 ;;
        --precision)             read -r -a PRECISIONS <<< "${2:-}"; shift 2 ;;
        --label-size)            LABEL_SIZE="${2:?}"; shift 2 ;;
        --jobs)                  JOBS="${2:?}"; shift 2 ;;
        --paraview-version)      PARAVIEW_VERSION="${2:?}"; shift 2 ;;
        --paraview-series)       PARAVIEW_SERIES="${2:?}"; shift 2 ;;
        --paraview-flavor)       PARAVIEW_FLAVOR="${2:?}"; shift 2 ;;
        --clusterdict-precision) CLUSTERDICT_PRECISION="${2:?}"; shift 2 ;;
        --no-foam-source)        INJECT_FOAM_SOURCE=0; shift ;;
        --install-deps)          INSTALL_DEPS=1; shift ;;
        --skip-openfoam)         SKIP_OPENFOAM=1; shift ;;
        --skip-paraview)         SKIP_PARAVIEW=1; shift ;;
        --skip-clusterdict)      SKIP_CLUSTERDICT=1; shift ;;
        --dry-run)               DRY_RUN=1; shift ;;
        -h|--help)               usage; exit 0 ;;
        *)                       die "Unknown option: $1 (see --help)" ;;
    esac
done

PREFIX="${PREFIX%/}"
OF_ROOT="$PREFIX/OpenFOAM"
PV_ROOT="$PREFIX/ParaView"
PV_DIRNAME="ParaView-${PARAVIEW_VERSION}-${PARAVIEW_FLAVOR}"
PV_BIN="$PV_ROOT/$PV_DIRNAME/bin"

# --------------------------------------------------------------------------
# Prerequisite checks
# --------------------------------------------------------------------------
require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command '$1' not found. $2"
}

url_exists() {
    # HTTP existence check without downloading the body.
    curl -fsIL --retry 2 -o /dev/null "$1"
}

fetch() {
    # fetch URL DEST
    local url="$1" dest="$2"
    if [[ -f "$dest" ]]; then
        log "Already downloaded: $dest"
        return 0
    fi
    run curl -fL --retry 3 --progress-bar -o "$dest.part" "$url"
    [[ "$DRY_RUN" -eq 0 ]] && mv "$dest.part" "$dest"
}

check_prereqs() {
    require_cmd curl "Install curl and re-run."
    require_cmd tar  "Install tar and re-run."
    if [[ "$SKIP_OPENFOAM" -eq 0 ]]; then
        [[ "$(uname -s)" == "Linux" ]] || warn "ESI OpenFOAM compilation is only supported on Linux (native or WSL). Detected: $(uname -s)."
        for c in git make gcc g++ flex bison cmake mpirun; do
            command -v "$c" >/dev/null 2>&1 || warn "Build tool '$c' not found. Consider --install-deps (Ubuntu) or install it manually."
        done
    fi
}

install_deps() {
    [[ "$INSTALL_DEPS" -eq 1 ]] || return 0
    command -v apt-get >/dev/null 2>&1 || die "--install-deps only supports apt-based systems (Ubuntu/Debian)."
    log "Installing build prerequisites via apt-get (sudo required)..."
    run sudo apt-get update
    run sudo apt-get install -y \
        build-essential cmake git ca-certificates curl \
        flex bison zlib1g-dev libboost-system-dev libboost-thread-dev \
        libopenmpi-dev openmpi-bin gnuplot libreadline-dev libncurses-dev \
        libxt-dev libscotch-dev libptscotch-dev libfftw3-dev libcgal-dev
}

# --------------------------------------------------------------------------
# OpenFOAM
# --------------------------------------------------------------------------
normalise_version() {
    # Accept 2406 or v2406 -> v2406
    local v="$1"
    [[ "$v" == v* ]] || v="v$v"
    printf '%s' "$v"
}

verify_versions() {
    local ok=1 v url
    for v in "${VERSIONS[@]}"; do
        v="$(normalise_version "$v")"
        url="$OF_SOURCE_BASE/$v/OpenFOAM-$v.tgz"
        if url_exists "$url"; then
            log "Verified OpenFOAM source exists: $v"
        else
            warn "OpenFOAM source NOT found for '$v' at $url"
            ok=0
        fi
    done
    [[ "$ok" -eq 1 ]] || die "One or more requested OpenFOAM versions do not exist. Fix --versions and re-run."
}

build_precision() {
    # build_precision <foamDir> <precision>
    local foamDir="$1" prec="$2"
    local prefs="$foamDir/etc/prefs.sh"
    log "Configuring precision $prec (WM_LABEL_SIZE=$LABEL_SIZE) in $prefs"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        cat > "$prefs" <<EOF
# Written by installOpenFOAM.sh
export WM_PRECISION_OPTION=$prec
export WM_LABEL_SIZE=$LABEL_SIZE
export WM_COMPILE_OPTION=Opt
EOF
    fi
    local logFile="$foamDir/log.Allwmake.${prec}Int${LABEL_SIZE}"
    log "Building OpenFOAM ($prec) - this can take a long time. Log: $logFile"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        (
            set -e
            # shellcheck disable=SC1091
            source "$foamDir/etc/bashrc"
            log "  Platform: \${WM_OPTIONS:-unknown}"
            if [[ -x "$WM_THIRD_PARTY_DIR/Allwmake" ]]; then
                ( cd "$WM_THIRD_PARTY_DIR" && ./Allwmake -j "$JOBS" -q ) 2>&1 | tee "$foamDir/log.ThirdParty.${prec}Int${LABEL_SIZE}"
            fi
            ( cd "$WM_PROJECT_DIR" && ./Allwmake -j "$JOBS" -q -s ) 2>&1 | tee "$logFile"
        )
    fi
}

install_openfoam_version() {
    local v; v="$(normalise_version "$1")"
    local foamTgz="OpenFOAM-$v.tgz"
    local tpTgz="ThirdParty-$v.tgz"
    local foamUrl="$OF_SOURCE_BASE/$v/$foamTgz"
    local tpUrl="$OF_SOURCE_BASE/$v/$tpTgz"
    local foamDir="$OF_ROOT/OpenFOAM-$v"

    run mkdir -p "$OF_ROOT"

    if [[ -d "$foamDir" ]]; then
        log "OpenFOAM source already extracted: $foamDir"
    else
        log "Downloading $v ..."
        fetch "$foamUrl" "$OF_ROOT/$foamTgz"
        if url_exists "$tpUrl"; then
            fetch "$tpUrl" "$OF_ROOT/$tpTgz"
        else
            warn "No ThirdParty tarball for $v; relying on system libraries."
        fi
        log "Extracting $v ..."
        run tar -xzf "$OF_ROOT/$foamTgz" -C "$OF_ROOT"
        [[ -f "$OF_ROOT/$tpTgz" ]] && run tar -xzf "$OF_ROOT/$tpTgz" -C "$OF_ROOT"
    fi

    local prec
    for prec in "${PRECISIONS[@]}"; do
        build_precision "$foamDir" "$prec"
    done
    log "Finished OpenFOAM $v (${PRECISIONS[*]}). bashrc: $foamDir/etc/bashrc"
}

install_openfoam() {
    [[ "$SKIP_OPENFOAM" -eq 0 ]] || { log "Skipping OpenFOAM stage."; return 0; }
    [[ "${#VERSIONS[@]}" -gt 0 ]] || die "No versions given. Pass --versions \"v2406 ...\" or use --skip-openfoam."
    verify_versions
    local v
    for v in "${VERSIONS[@]}"; do
        install_openfoam_version "$v"
    done
}

# --------------------------------------------------------------------------
# ParaView
# --------------------------------------------------------------------------
install_paraview() {
    [[ "$SKIP_PARAVIEW" -eq 0 ]] || { log "Skipping ParaView stage."; return 0; }
    local tgz="$PV_DIRNAME.tar.gz"
    local url="$PV_FILES_BASE/$PARAVIEW_SERIES/$tgz"
    local dest="$PV_ROOT/$PV_DIRNAME"

    if [[ -x "$PV_BIN/pvbatch" ]]; then
        log "ParaView already installed: $PV_BIN/pvbatch"
        return 0
    fi

    log "Verifying ParaView download exists: $url"
    url_exists "$url" || die "ParaView release not found at $url (check --paraview-version/--paraview-series/--paraview-flavor)."

    run mkdir -p "$PV_ROOT"
    fetch "$url" "$PV_ROOT/$tgz"
    log "Extracting ParaView ..."
    run tar -xzf "$PV_ROOT/$tgz" -C "$PV_ROOT"

    if [[ "$DRY_RUN" -eq 0 && ! -x "$PV_BIN/pvbatch" ]]; then
        # Some tarballs extract to a differently cased/versioned dir; try to locate pvbatch.
        local found
        found="$(find "$PV_ROOT" -maxdepth 3 -type f -name pvbatch 2>/dev/null | head -n1 || true)"
        [[ -n "$found" ]] || die "pvbatch not found after extracting $tgz."
        PV_BIN="$(dirname "$found")"
        warn "pvbatch located at $PV_BIN (differs from expected layout)."
    fi
    log "ParaView ready: $PV_BIN/pvbatch"
}

# --------------------------------------------------------------------------
# clusterDict rewrite
# --------------------------------------------------------------------------
rewrite_clusterdicts() {
    [[ "$SKIP_CLUSTERDICT" -eq 0 ]] || { log "Skipping clusterDict rewrite."; return 0; }

    local pvbatch="$PV_BIN/pvbatch"
    local dicts=()
    local d
    while IFS= read -r d; do dicts+=("$d"); done < <(find "$REPO_DIR/setupTemplates" -path '*/defaultCluster/*/clusterDict' -type f | sort)
    [[ "${#dicts[@]}" -gt 0 ]] || { warn "No clusterDict files found under $REPO_DIR/setupTemplates."; return 0; }

    log "Rewriting ${#dicts[@]} clusterDict file(s):"
    printf '  %s\n' "${dicts[@]}"

    if [[ "$DRY_RUN" -eq 1 ]]; then
        log "(dry-run) would set:"
        log "  zeroTemplates -> $REPO_DIR/zeroTemplates"
        log "  postUtilities -> $REPO_DIR/postUtilities"
        log "  pvbatch       -> $pvbatch"
        [[ "$INJECT_FOAM_SOURCE" -eq 1 ]] && log "  foam source    -> per-template etc/bashrc ($CLUSTERDICT_PRECISION build)"
        return 0
    fi

    REPO_DIR="$REPO_DIR" PVBATCH="$pvbatch" OF_ROOT="$OF_ROOT" \
    INJECT_FOAM_SOURCE="$INJECT_FOAM_SOURCE" \
    CLUSTERDICT_PRECISION="$CLUSTERDICT_PRECISION" \
    python3 - "${dicts[@]}" <<'PYEOF'
import os, re, sys

repo   = os.environ["REPO_DIR"]
pvbatch = os.environ["PVBATCH"]
of_root = os.environ["OF_ROOT"]
inject  = os.environ.get("INJECT_FOAM_SOURCE", "0") == "1"

zero_tpl = f"{repo}/zeroTemplates"
post_dir = f"{repo}/postUtilities"

# Pick a default etc/bashrc for foam-source injection (primary = first version dir found).
def bashrc_for(path):
    # Choose a version-appropriate OpenFOAM build for this template.
    # Templates whose directory name contains a 4-digit ESI tag (e.g. 2606)
    # source that version; everything else uses the first available build.
    import glob
    builds = sorted(glob.glob(os.path.join(of_root, "OpenFOAM-v*")))
    if not builds:
        return None
    m = re.search(r'(\d{4})', os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(path)))))
    if m:
        tag = m.group(1)
        for b in builds:
            if tag in os.path.basename(b):
                return os.path.join(b, "etc", "bashrc")
    return os.path.join(builds[0], "etc", "bashrc")

# Token = run of chars that are not space/quote/semicolon.
TOKEN = r"[^ \"';]*"

def rewrite(text, path):
    # caseTemplates='...zeroTemplates'  (single-quoted value)
    text = re.sub(r"caseTemplates='[^']*'", f"caseTemplates='{zero_tpl}'", text)
    # any .../postUtilities prefix -> repo postUtilities (keeps trailing /script.py)
    text = re.sub(TOKEN + "postUtilities", post_dir, text)
    # any .../pvbatch -> installed pvbatch
    text = re.sub(TOKEN + "pvbatch", pvbatch, text)

    if inject:
        br = bashrc_for(path)
        if br:
            # Idempotent: strip a prior injected source, then prepend the new one.
            text = re.sub(
                r"(\. [^;]*etc/bashrc ; )?\. \$WM_PROJECT_DIR/bin/tools/RunFunctions",
                f". {br} ; . $WM_PROJECT_DIR/bin/tools/RunFunctions",
                text,
            )
    return text

for path in sys.argv[1:]:
    with open(path, "r") as fh:
        original = fh.read()
    updated = rewrite(original, path)
    if updated != original:
        with open(path, "w") as fh:
            fh.write(updated)
        print(f"  updated {path}")
    else:
        print(f"  unchanged {path}")
PYEOF
    log "clusterDict rewrite complete."
}

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
print_summary() {
    log "Done."
    echo
    echo "  Repository   : $REPO_DIR"
    echo "  Install root : $PREFIX"
    if [[ "$SKIP_OPENFOAM" -eq 0 ]]; then
        echo "  OpenFOAM     : $OF_ROOT/OpenFOAM-<version> (built: ${PRECISIONS[*]})"
        echo "  Source it    : source $OF_ROOT/OpenFOAM-<version>/etc/bashrc"
    fi
    [[ "$SKIP_PARAVIEW" -eq 0 ]] && echo "  ParaView     : $PV_BIN/pvbatch"
    [[ "$SKIP_CLUSTERDICT" -eq 0 ]] && echo "  clusterDicts : rewritten under $REPO_DIR/setupTemplates"
    echo
}

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
main() {
    log "Repository: $REPO_DIR"
    log "Install prefix: $PREFIX"
    [[ "$SKIP_OPENFOAM" -eq 0 ]] && log "OpenFOAM versions: ${VERSIONS[*]:-<none>} | precisions: ${PRECISIONS[*]}"
    check_prereqs
    install_deps
    install_openfoam
    install_paraview
    rewrite_clusterdicts
    print_summary
}

main "$@"
