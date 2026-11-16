#!/usr/bin/env bash
# =============================================================================
# Copy the pipeline code and the reference data from this machine (WSL) to
# the OCI VM, then build the container images there.
# Run in WSL:  cloud/oci/push_to_vm.sh ubuntu@<public-ip> [ssh-key] [refs-dir]
#
# Only files git knows about are copied (no results, biga_* folders or work/).
# The images are built on the VM from containers/*/Dockerfile, so the private
# gencorr-tools image (LAVA licence) never goes through a registry.
# =============================================================================
set -euo pipefail

HOST=${1:?usage: push_to_vm.sh ubuntu@<ip> [ssh-key] [refs-dir]}
KEY=${2:-$HOME/.ssh/oci_gencorr}
REFS=${3:-$HOME/refs}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
SSH="ssh -i $KEY -o StrictHostKeyChecking=accept-new"

echo "== code"
cd "$REPO"
git ls-files -co --exclude-standard -z \
    | rsync -az --from0 --files-from=- -e "$SSH" ./ "$HOST:gencorr/"
$SSH "$HOST" 'chmod +x ~/gencorr/bin/* ~/gencorr/cloud/oci/*.sh'

echo "== references ($REFS)"
rsync -az --info=progress2 -e "$SSH" \
    "$REFS/eur_w_ld_chr" "$REFS/w_hm3.snplist" "$REFS/ldak" "$REFS/g1000_eur.zip" \
    "$HOST:/data/refs/"
$SSH "$HOST" 'cd /data/refs && [ -f g1000_eur.bed ] || unzip -oq g1000_eur.zip'

echo "== images (built on the VM)"
$SSH "$HOST" 'bash ~/gencorr/cloud/oci/build_images.sh'
