# Running gencorr on an Oracle Cloud (OCI) VM

One x86 VM runs everything with Nextflow's local executor and Docker. Data,
images and the Nextflow work directory live on a separate block volume, so the
VM can be stopped (no compute charge) or resized without losing anything.

## Why this layout

- **LAVA dominates the cost** (~1.5–3.5 CPU-hours per genome-wide pair). Its
  chunks are independent 1-CPU tasks that peak at ~3.7 GB RAM, so the VM
  should have many vCPUs with about 5 GB per vCPU.
- **x86, not Ampere (ARM)**: the images are amd64, and LDAK ships x86 binaries.
  The Always Free Ampere A1 shape would need arm64 rebuilds.
- **One VM, local executor**: no cluster or batch service to set up; OCI's
  Flex shapes can be resized (reboot) when you need more cores.
- **Images built on the VM** from `containers/`: the private `gencorr-tools`
  image (LAVA licence) never goes through a registry, and no tokens are needed.

## 1. Create the VM (OCI console)

1. **Compute → Instances → Create instance**
   - Image: **Canonical Ubuntu 24.04** (x86_64), or **Oracle Linux 9** if no
     shapes are offered for Ubuntu. The scripts handle both; on Oracle Linux the
     login user is `opc` instead of `ubuntu` in every command below.
   - Shape: **VM.Standard.E5.Flex** (or E4.Flex). Start with
     **8 OCPUs / 96 GB** (= 16 vCPUs; on x86 one OCPU is 2 vCPUs). Raise it
     for many pairs.
   - Networking: create a new VCN with a public subnet, assign a public IPv4.
   - SSH key: upload the public key of a new pair, e.g. made in WSL with
     `ssh-keygen -t ed25519 -f ~/.ssh/oci_gencorr`.
   - Boot volume: default (50 GB) is enough.
2. **Block Storage → Block volumes → Create**: 200 GB, same availability
   domain, Balanced performance. Then on the instance:
   **Attached block volumes → Attach**, type **Paravirtualized**, device path
   **/dev/oracleoci/oraclevdb**.
3. Optional but recommended: **Billing → Budgets** with an alert at e.g. $50,
   and in the VCN's security list restrict SSH (port 22) to your own IP.

If the shape is greyed out or the launch fails with a limit error, check
**Governance → Limits, Quotas and Usage** (trial accounts have low per-shape
limits) and request an increase, or use fewer OCPUs (4 OCPUs / 48 GB worked on
a new trial account where 8 did not). "No shapes are available for this image"
means the same limit, not a bad image.

Pitfalls seen on a new trial account:

- **No public IP**: the subnet must be public. Create the network with
  **Virtual cloud networks → Start VCN Wizard → Create VCN with Internet
  Connectivity** and pick its *public* subnet when creating the instance.
- **No block volume possible**: enlarge the boot volume instead (boot volume
  page → **Edit** → 200 GB), then on the VM
  `echo 1 | sudo tee /sys/class/block/sda/device/rescan; sudo /usr/libexec/oci-growfs -y`
  (Oracle Linux), and run the bootstrap as `sudo bash bootstrap_vm.sh none`.
  Data then lives on the boot volume: keep it if you ever terminate the VM.

## 2. Set it up

From WSL, in the repository:

```bash
scp -i ~/.ssh/oci_gencorr cloud/oci/bootstrap_vm.sh ubuntu@<ip>:
ssh -i ~/.ssh/oci_gencorr ubuntu@<ip> 'sudo bash bootstrap_vm.sh'
cloud/oci/push_to_vm.sh ubuntu@<ip>          # code + ~/refs + build images
```

`push_to_vm.sh` copies only files tracked by git (plus untracked, non-ignored
ones), the references from `~/refs` (~560 MB), and then builds the four images
on the VM (a few minutes; the LDSC build also runs its upstream tests).

## 3. Run

```bash
ssh -i ~/.ssh/oci_gencorr ubuntu@<ip>
tmux new -s gencorr                # keeps the run alive if SSH drops
cd /data
nextflow run ~/gencorr/main.nf -profile docker -c /data/gencorr_vm.config \
    --input input/samplesheet.csv --outdir results/<run> \
    --methods ldsc,lava \
    --ldsc_ld refs/eur_w_ld_chr --ldsc_merge_alleles refs/w_hm3.snplist \
    --lava_ref refs/g1000_eur --lava_snp_list refs/w_hm3.snplist \
    --gwascatalog_cache /data/gwascatalog_cache \
    --lava_loci_per_chunk 150 -resume
```

Put sumstats in `/data/input` (`rsync`/`scp`, or GWAS Catalog accessions in the
samplesheet, which the VM downloads itself at datacentre speed).

**Chunk size**: a pair is split into ~2,495 / `lava_loci_per_chunk` LAVA
tasks. Aim for at least as many tasks in total as vCPUs: for one pair on
16 vCPUs use 150; with many pairs the default 250 is fine. Results are
identical for any chunk size (per-locus seeds).

Copy results home and stop the VM when done:

```bash
rsync -az -e "ssh -i ~/.ssh/oci_gencorr" ubuntu@<ip>:/data/results/ ./results_oci/
```

Then **Stop** the instance in the console. Stopped Standard shapes are not
charged for CPU/memory; the volumes are (200 GB of block storage is in the
Always Free allowance). After resizing, re-run `sudo bash bootstrap_vm.sh`
to update `/data/gencorr_vm.config`. Clean `/data/work` between projects.

## Cost and the free trial

Check current prices in OCI's cost estimator; as a rough order of magnitude an
E5.Flex VM with 8 OCPUs / 96 GB costs well under $1 per hour, so a genome-wide
pair costs cents. Before the trial ends, copy everything you want to keep:
resources outside the Always Free tier are reclaimed if the account is not
upgraded to Pay As You Go.
