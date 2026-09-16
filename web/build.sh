#!/usr/bin/env bash
# Build the WASM binary and drop it next to web/index.html + mq_js_bundle.js,
# ready to serve.
#
# Usage:
#   web/build.sh          # build, then print how to serve it
#   web/build.sh --serve  # build, then serve on http://localhost:8000
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if ! rustup target list --installed 2>/dev/null | grep -q '^wasm32-unknown-unknown$'; then
  echo "Installing the wasm32-unknown-unknown target..."
  rustup target add wasm32-unknown-unknown
fi

cargo build --release --target wasm32-unknown-unknown
cp target/wasm32-unknown-unknown/release/billiards.wasm web/billiards.wasm

echo
echo "Built web/billiards.wasm"

if [ "${1:-}" = "--serve" ]; then
  echo "Serving http://localhost:8000/ (Ctrl-C to stop)"
  cd web && python3 -m http.server 8000
else
  echo "Serve it locally (wasm needs http://, not file://), e.g.:"
  echo "  cd web && python3 -m http.server 8000"
fi
