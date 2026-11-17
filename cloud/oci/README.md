# gencorr on an Oracle Cloud (OCI) VM

Scripts to run gencorr on one x86 OCI VM (Ubuntu 24.04 or Oracle Linux 8/9)
with Nextflow's local executor and Docker. The full tutorial, with the console
steps, the reasons for each choice and the errors met on a new trial account,
is [docs/cloud_oci.md](../../docs/cloud_oci.md).

| File              | Runs on                           | What it does                                                                                                                        |
| ----------------- | --------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `bootstrap_vm.sh` | VM, `sudo`                        | `/data` (block volume, or boot disk with `none`), Docker, Java 21, Nextflow 26.10.5, `/data/gencorr_vm.config` with the VM's limits |
| `push_to_vm.sh`   | laptop (WSL)                      | git-tracked code → `~/gencorr`, `~/refs` → `/data/refs`, then `build_images.sh`                                                     |
| `build_images.sh` | VM                                | builds every image the modules reference from `containers/`                                                                         |
| `oci.config`      | (included by `gencorr_vm.config`) | LAVA memory 5 GB, one thread per task, `workDir = /data/work`                                                                       |

Quick start (VM already created, with a public IP):

```bash
scp -i ~/.ssh/oci_gencorr cloud/oci/bootstrap_vm.sh opc@<ip>:
ssh -i ~/.ssh/oci_gencorr opc@<ip> 'sudo bash bootstrap_vm.sh none'   # omit 'none' with a block volume
cloud/oci/push_to_vm.sh opc@<ip>
ssh -i ~/.ssh/oci_gencorr opc@<ip>
cd /data && nextflow run ~/gencorr/main.nf -profile test,docker -c /data/gencorr_vm.config --outdir results/test
```

Use `ubuntu@` instead of `opc@` on Ubuntu images.
