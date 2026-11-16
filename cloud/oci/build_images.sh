#!/usr/bin/env bash
# Build every image the modules use from containers/<name>/Dockerfile, tagged
# exactly as the modules reference them. Skips images that already exist.
# Run on the VM:  bash ~/gencorr/cloud/oci/build_images.sh
set -euo pipefail
cd "$(dirname "$0")/../.."

grep -rhoE "ghcr\.io/mauricioroza/gencorr-[a-z]+:[A-Za-z0-9._-]+" modules/local | sort -u |
while read -r image; do
    name=${image#ghcr.io/mauricioroza/}; name=${name%%:*}
    if docker image inspect "$image" >/dev/null 2>&1; then
        echo "have  $image"
    else
        echo "build $image"
        docker build -q -t "$image" "containers/$name"
    fi
done
docker images --format '{{.Repository}}:{{.Tag}}  {{.Size}}' | grep gencorr
