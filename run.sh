#!/usr/bin/env bash
# Launch the mRNA Designer app locally.
set -euo pipefail
cd "$(dirname "$0")"
PORT="${PORT:-8000}"
python3 -m uvicorn mrna_designer.api:app --host 127.0.0.1 --port "$PORT" --reload
