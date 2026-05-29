#!/usr/bin/env bash
# Native I/O only — no Boost required for kcollections 2.2+
set -euxo pipefail

if [ "$(uname)" = "Linux" ]; then
    if command -v yum &>/dev/null; then
        yum install -y zlib-devel
    elif command -v apt-get &>/dev/null; then
        apt-get update
        apt-get install -y zlib1g-dev
    elif command -v apk &>/dev/null; then
        apk add --no-cache zlib-dev
    fi
elif [ "$(uname)" = "Darwin" ]; then
    true
fi
