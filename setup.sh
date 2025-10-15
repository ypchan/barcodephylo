#!/usr/bin/env bash
set -euo pipefail

# ========== Configuration ==========
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$REPO_DIR/scripts"
LINK_DIR="${LINK_DIR:-$HOME/bin}"     # Override with LINK_DIR=/some/dir
IQTREE_BIN="${IQTREE_BIN:-}"          # Optionally point to a specific iqtree binary
# ===================================

# ---------- Logging helpers ----------
c_green='\033[0;32m'; c_red='\033[0;31m'; c_yel='\033[0;33m'; c_blue='\033[0;34m'; c_end='\033[0m'
log()  { printf "${c_blue}[setup]${c_end} %s\n" "$*"; }
ok()   { printf "${c_green}[ok]${c_end} %s\n" "$*"; }
warn() { printf "${c_yel}[warn]${c_end} %s\n" "$*"; }
err()  { printf "${c_red}[err]${c_end} %s\n" "$*" >&2; }

require_dir() { [ -d "$1" ] || { err "Missing directory: $1"; exit 1; }; }

# ---------- 0) Sanity checks ----------
require_dir "$REPO_DIR"
require_dir "$SCRIPTS_DIR"
log "Repository: $REPO_DIR"
log "Scripts dir: $SCRIPTS_DIR"

if [[ "$(uname -s)" != "Linux" ]]; then
  warn "This setup is intended for Linux; proceeding anyway."
fi

# ---------- 1) Ensure LINK_DIR exists ----------
mkdir -p "$LINK_DIR"
ok "Link dir ready: $LINK_DIR"

# ---------- 2) Locate IQ-TREE (only dependency) ----------
resolve_iqtree() {
  local bin="$IQTREE_BIN"
  if [[ -n "$bin" ]]; then
    [[ -x "$bin" ]] || { err "IQTREE_BIN specified but not executable: $bin"; return 1; }
    echo "$bin"; return 0
  fi
  if command -v iqtree >/dev/null 2>&1; then
    command -v iqtree; return 0
  fi
  if command -v iqtree2 >/dev/null 2>&1; then
    command -v iqtree2; return 0
  fi
  return 1
}

IQBIN="$(resolve_iqtree || true)"
if [[ -z "${IQBIN:-}" ]]; then
  err "IQ-TREE not found (neither 'iqtree' nor 'iqtree2'). Please install IQ-TREE and re-run."
  exit 1
fi
ok "Found IQ-TREE binary: $IQBIN"

# Create a compatibility symlink 'iqtree' if only 'iqtree2' exists and 'iqtree' is missing.
if ! command -v iqtree >/dev/null 2>&1; then
  target="$LINK_DIR/iqtree"
  rm -f "$target"
  ln -s "$IQBIN" "$target"
  ok "Created compatibility symlink: $target -> $IQBIN"
fi

# ---------- 3) Symlink repo scripts into LINK_DIR ----------
log "Linking scripts into $LINK_DIR ..."
link_count=0
while IFS= read -r -d '' f; do
  base="$(basename "$f")"
  dest="$LINK_DIR/$base"

  # Ensure the script is executable (in repo)
  chmod +x "$f" || true

  # Remove existing target (file or symlink), then link
  if [[ -e "$dest" || -L "$dest" ]]; then
    rm -f "$dest"
  fi
  ln -s "$f" "$dest"
  ok "Linked: $dest -> $f"
  link_count=$((link_count+1))
done < <(find "$SCRIPTS_DIR" -maxdepth 1 -type f -print0)

log "Total linked scripts: $link_count"

# ---------- 4) PATH hint ----------
case ":$PATH:" in
  *":$LINK_DIR:"*) ok "\$PATH already includes $LINK_DIR" ;;
  *)
    warn "\$PATH does not include $LINK_DIR"
    echo "Add it for bash:  echo 'export PATH=\"\$HOME/bin:\$PATH\"' >> \$HOME/.bashrc && source \$HOME/.bashrc"
    echo "Add it for zsh:   echo 'export PATH=\"\$HOME/bin:\$PATH\"' >> \$HOME/.zshrc  && source \$HOME/.zshrc"
    ;;
esac

ok "Setup complete. You can now run scripts from anywhere."
