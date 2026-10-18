#!/bin/bash
set -euo pipefail
package_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
docker build -f "${package_dir}/inst/docker/Dockerfile" "${package_dir}" -t vlemonidis/cment
docker push vlemonidis/cment
