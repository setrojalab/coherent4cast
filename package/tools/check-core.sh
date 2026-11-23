#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cxx=${CXX:-g++}
extra=()
if [[ "${SANITIZE:-0}" == "1" ]]; then
  extra=(-fsanitize=address,undefined -fno-omit-frame-pointer -g)
fi
for src in "$repo"/tests/core/*.cpp; do
  exe="$work/$(basename "${src%.cpp}")"
  "$cxx" -std=c++17 -Wall -Wextra -pedantic -Wconversion -Werror \
    -I"$repo/inst/include" "${extra[@]}" "$src" -o "$exe"
  "$exe"
done
