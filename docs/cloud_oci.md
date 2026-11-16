# Running gencorr on Oracle Cloud (OCI): setup tutorial

This is the exact setup used to run gencorr on an Oracle Cloud free-trial
account (November 2026), with the reason behind each choice and the errors met
on the way. The scripts it uses are in [`cloud/oci/`](../cloud/oci/).

The end result:

```
your laptop (WSL)                         OCI VM (Oracle Linux 9, x86)
──────────────────                        ──────────────────────────────────────
~/gencorr  (code)      ── push_to_vm.sh ─▶ ~/gencorr          code (git-tracked files)
~/refs     (references) ─────────────────▶ /data/refs         LD scores, g1000_eur, tagging file
                                           /data/docker       4 images, built on the VM
ssh / rsync  ◀── results ───────────────── /data/results      pipeline output
                                           /data/work         Nextflow work dir
                                           /data/gencorr_vm.config  limits for this VM
```

One VM runs everything: Nextflow with the local executor, every task in Docker.

## Contents

1. [Decisions and reasons](#1-decisions-and-reasons)
2. [SSH key](#2-ssh-key-on-the-laptop)
3. [Network](#3-network-vcn-with-a-public-subnet)
4. [The VM](#4-the-vm)
5. [Disk space](#5-disk-space-boot-volume-to-200-gb)
6. [Software on the VM](#6-software-on-the-vm-bootstrap_vmsh)
7. [Code, references and images](#7-code-references-and-images-push_to_vmsh)
8. [Test run](#8-test-run)
9. [Real runs](#9-real-runs)
10. [Daily use and costs](#10-daily-use-and-costs)
11. [Errors met and their fixes](#11-errors-met-and-their-fixes)

## 1. Decisions and reasons

| Decision | Reason |
|---|---|
| **One VM, local executor** (no cluster, no batch service) | gencorr's heavy part, LAVA, is many independent 1-CPU tasks: one machine with many cores runs them in parallel with no scheduler to set up. A bigger machine is a resize away (Flex shapes). |
| **x86 (AMD E5.Flex), not Ampere/ARM** | The images are built for amd64 and LDAK is distributed as an x86 binary. The free Ampere A1 shape would need arm64 rebuilds of all four images. |
| **4 OCPUs / 48 GB** (8 vCPUs) | LAVA chunks peak at ~3.7 GB RAM each on one CPU (measured in earlier runs with the genome-wide `g1000_eur` reference), so ~5 GB per vCPU keeps every core busy. 8 OCPUs was the plan, but the new trial account only allowed 4. On x86, 1 OCPU = 2 vCPUs. |
| **Oracle Linux 9** instead of Ubuntu | Ubuntu 24.04 offered no usable shape on the new account; Oracle Linux (OCI's own image) did. The scripts support both (`apt` or `dnf`); the login user is `opc` instead of `ubuntu`. |
| **Boot volume enlarged to 200 GB**, no separate block volume | Attaching a block volume was not available on the trial account at first. A 200 GB boot volume holds images (~4 GB), references (~4 GB unpacked) and work files, and stays inside OCI's Always Free 200 GB of storage. Kept that way when volumes became available ([why](#5-disk-space-boot-volume-to-200-gb)). |
| **Everything under `/data`** | One place for references, images, work dir and results; easy to check space (`df -h /data`) and to move onto a separate volume later. |
| **Images built on the VM** from `containers/` | `gencorr-tools` contains LAVA ("All rights reserved") and must stay private: building on the VM means it never goes through a registry, and the VM needs no GitHub token. The Dockerfiles pin every version, so the build matches the local images. |
| **Code copied with rsync of git-tracked files** | The repository is private; copying avoids putting credentials on the VM. Only tracked (and non-ignored) files go, so user data (`biga_*`, `results*`, `work/`) never leaves the laptop. |
| **Nextflow 26.10.5 pinned** | Same version as the development machine, where the test suite passes. Ignore the "newer version available" notice. |
| **One thread per task** (`OMP_NUM_THREADS=1`, data.table, BLAS) | Nextflow already runs one task per CPU; without the limit, R's data.table and BLAS start threads for every core in every task and overload the machine. Results do not change (LAVA seeds are per locus). |
| **LAVA memory request 5 GB** (pipeline default 6 GB) | Matches the measured ~3.7 GB peak with headroom; at 6 GB, 8 tasks would need 48 GB and one core would sit idle. An out-of-memory kill retries with 10 GB. |
| **SSH key with a passphrase + ssh-agent** | The key gives full access to the VM and its data; the agent means typing the passphrase once per terminal. |
| **Runs inside `tmux`** | A genome-wide run takes longer than an SSH session may last; tmux keeps it running if the connection drops. |
| **Stop the VM when idle, never terminate** | A stopped VM is not charged for CPU/memory. Terminating deletes the VM (and, by default, the disk with all data and images). |

## 2. SSH key (on the laptop)

In WSL:

```bash
ssh-keygen -t ed25519 -f ~/.ssh/oci_gencorr -C oci-gencorr     # set a passphrase
cat ~/.ssh/oci_gencorr.pub                                       # the line to give OCI
```

The key lives in WSL (`\\wsl.localhost\Ubuntu\home\<user>\.ssh` from Windows),
so pasting the `cat` output is easier than browsing for the file. Upload only
the `.pub` file, never `oci_gencorr`.

Once per terminal session, load it so the passphrase is asked only once:

```bash
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/oci_gencorr
```

## 3. Network: VCN with a public subnet

**☰ → Networking → Virtual cloud networks → Start VCN Wizard → Create VCN with
Internet Connectivity** (the wizard, *not* the plain "Create VCN" form).

- Name `gencorr-vcn`, root compartment, default CIDRs (`10.0.0.0/16`,
  public subnet `10.0.0.0/24`, private subnet `10.0.1.0/24`).

It creates the public subnet, an Internet Gateway, the route `0.0.0.0/0 → IGW`
and a security list that allows SSH. Check under the VCN's **Gateways** tab that
*Internet gateway-gencorr-vcn* is **Available** (the DRG section below it is
expected to be empty).

*Why:* the network the instance wizard created on the first try had no public
IP option, so the VM could not be reached from home.

Optional: in the default security list, change the SSH rule's source from
`0.0.0.0/0` to `<your-ip>/32`.

## 4. The VM

**☰ → Compute → Instances → Create instance**

| Field | Value |
|---|---|
| Name | `gencorr` |
| Placement | AD-1 (Vinhedo has one availability domain) |
| Shape | **Change shape → Virtual machine → AMD → VM.Standard.E5.Flex**; expand the row (⌄) to set **4 OCPUs / 48 GB** |
| Image | **Oracle Linux 9**, a build without `aarch64`/`GPU` in the name |
| Networking | existing VCN `gencorr-vcn`, subnet **public subnet-gencorr-vcn**, **Automatically assign public IPv4 address: on** |
| SSH keys | paste the `.pub` line |
| Boot volume | default (enlarged in the next step) |

Pick the shape **before** the image: the console then only lists images that
fit the shape. When it is **Running**, copy the public IP and test:

```bash
ssh -i ~/.ssh/oci_gencorr opc@<public-ip> 'uname -m; nproc; free -g | head -2'
# x86_64, 8, ~47 GB
```

## 5. Disk space: boot volume to 200 GB

1. Instance page → **Boot volume** → click the volume → **Edit** → size
   **200 GB** → **Save changes** (the VM keeps running).
2. Make the VM see the new size, then grow `/`:

   ```bash
   ssh -i ~/.ssh/oci_gencorr opc@<public-ip> 'echo 1 | sudo tee /sys/class/block/sda/device/rescan; lsblk /dev/sda'
   ssh -i ~/.ssh/oci_gencorr opc@<public-ip> 'sudo /usr/libexec/oci-growfs -y; df -h /'
   ```

   `lsblk` should show `sda` at 200G and `df` about 183G for `/`
   (15G stays on `/var/oled`, Oracle Linux's log volume). If `lsblk` still
   shows the old size, wait a minute or `sudo reboot`, then run `oci-growfs`
   again.

*Why the rescan:* resizing in the console enlarges the disk, but the running
kernel keeps the old size until told to re-read it. `oci-growfs` then grows the
three layers in order: partition (`sda3`) → LVM volume (`ocivolume-root`) →
XFS filesystem.

With a separate block volume instead (if your account allows it): create it in
the same AD (200 GB, Balanced), attach it as **Paravirtualized** with device
path **`/dev/oracleoci/oraclevdb`**, and run the bootstrap without `none`;
it formats and mounts it at `/data`.

*Why the boot disk was kept* when block volumes later became available: the
boot volume cannot be shrunk, so a second 200 GB volume would double the
storage (200 GB above the Always Free allowance) for little gain. A Flex VM is
resized in place with its disk, and when terminating the VM the boot volume
can be kept. Move to a block volume only if space runs out or the data should
outlive the VM.

## 6. Software on the VM: `bootstrap_vm.sh`

```bash
scp -i ~/.ssh/oci_gencorr cloud/oci/bootstrap_vm.sh opc@<public-ip>:
ssh -i ~/.ssh/oci_gencorr opc@<public-ip> 'sudo bash bootstrap_vm.sh none'
```

(`none` = no separate volume, use the boot disk.) The script:

- creates `/data/{refs,work,results,input,nxf_home,gwascatalog_cache}`, owned
  by `opc` (Nextflow writes `.nextflow/` in the folder runs start from);
- installs Docker (docker-ce repo on Oracle Linux), Java 21, rsync, unzip, tmux;
  Docker stores images in `/data/docker`; adds `opc` to the `docker` group;
- downloads the Nextflow 26.10.5 launcher from its GitHub release and fetches
  the runtime into `/data/nxf_home`;
- writes `/data/gencorr_vm.config`: `resourceLimits` and executor limits =
  all vCPUs and total RAM − 4 GB, plus `cloud/oci/oci.config` (LAVA memory,
  threads, `workDir = /data/work`).

It ends with `Done: 8 vCPUs, 46 GB RAM, data on /data (… free)`. It can be
re-run at any time; **re-run it after resizing the VM** to update the limits.

## 7. Code, references and images: `push_to_vm.sh`

In WSL, from the repository:

```bash
cloud/oci/push_to_vm.sh opc@<public-ip>
```

1. **Code**: `git ls-files` → rsync to `~/gencorr` on the VM.
2. **References** from `~/refs`: `eur_w_ld_chr/`, `w_hm3.snplist`, `ldak/`,
   `g1000_eur.zip` (≈595 MB, unzipped on the VM into `/data/refs`).
3. **Images**: `build_images.sh` on the VM builds every image the modules
   reference, tagged exactly as in the modules:

   | Image | Size |
   |---|---|
   | `gencorr-tools:0.1.0` (Python helpers, R + LAVA) | 1.96 GB |
   | `gencorr-popcorn:1.0.0-20a29e8` | 918 MB |
   | `gencorr-ldsc:3.0.1-128f996` (runs upstream tests while building) | 786 MB |
   | `gencorr-ldak:6.1` | 192 MB |

Re-run it after changing the code; existing images are skipped
(`docker rmi` an image to force a rebuild after editing its Dockerfile).

## 8. Test run

```bash
ssh -i ~/.ssh/oci_gencorr opc@<public-ip>
tmux new -s gencorr
cd /data && nextflow run ~/gencorr/main.nf -profile test,docker -c /data/gencorr_vm.config --outdir results/test
```

Result on this VM: 4 traits, 6 pairs, LDSC + LAVA + SumHer, **56 tasks,
1 min 45 s, "Pipeline completed successfully"**.

tmux: **Ctrl+B then D** leaves it running; `tmux attach -t gencorr` returns.

## 9. Real runs

Upload sumstats (or use GWAS Catalog accessions in the samplesheet, which the
VM downloads itself):

```bash
rsync -az --info=progress2 -e "ssh -i ~/.ssh/oci_gencorr" <files> opc@<public-ip>:/data/input/
```

The samplesheet must use VM paths (`/data/input/...`). Then, in tmux on the VM:

```bash
cd /data
nextflow run ~/gencorr/main.nf -profile docker -c /data/gencorr_vm.config \
    --input input/samplesheet.csv --outdir results/<run> \
    --methods ldsc,lava \
    --ldsc_ld refs/eur_w_ld_chr --ldsc_merge_alleles refs/w_hm3.snplist \
    --lava_ref refs/g1000_eur --lava_snp_list refs/w_hm3.snplist \
    --gwascatalog_cache /data/gwascatalog_cache -resume
```

**LAVA chunk size**: one pair = ~2,495 loci / `lava_loci_per_chunk` tasks.
With the default 250 that is ~10 tasks per pair, enough for 8 vCPUs. On a
bigger VM with few pairs, lower it so tasks ≥ vCPUs (e.g. 150 for 16 vCPUs).
Results do not depend on the chunk size.

After the first real run, check `results/<run>/pipeline_info/execution_trace_*.txt`
(`peak_rss` of `LAVA_RUN`) to confirm the 5 GB request still fits.

Copy results home:

```bash
rsync -az -e "ssh -i ~/.ssh/oci_gencorr" opc@<public-ip>:/data/results/ ./results_oci/
```

Clean `/data/work` between projects (`rm -rf /data/work/*`; `-resume` then
starts from scratch).

## 10. Daily use and costs

| | Running | Stopped | Terminated |
|---|---|---|---|
| Compute (OCPUs + memory) | billed per second | **not billed** | gone |
| Boot volume (all data) | storage | storage (within Always Free 200 GB) | deleted unless kept |
| Public IP | free | kept | released |

- **Stop** (not Force stop) on the instance page when idle, after checking no
  run is active (`docker ps` empty). **Start** brings everything back, same IP.
- Rough order of magnitude: 4 OCPUs + 48 GB ≈ US$0.20–0.25/h. The trial credit
  (R$1,500 ≈ US$300, billed in BRL for a Brazilian account) is far more than a
  month of runs needs; the limit is the **trial end date** (30 days). Check
  **Billing → Cost analysis**; set a budget alert under **Billing → Budgets**.
- Before the trial ends, copy results home, then either upgrade to Pay As You
  Go (pay only while the VM runs) or let it lapse; resources outside Always
  Free are stopped and later deleted.
- Delete leftovers from failed attempts (extra boot volumes, extra VCNs): they
  count against limits and storage.

## 11. Errors met and their fixes

| Message | Cause | Fix |
|---|---|---|
| Shape/OCPU settings not visible ("1 core OCPU, 12 GB") | The settings are inside the shape row | Change shape → expand the E5.Flex row (⌄) |
| `Shape VM.Standard.E5.Flex is not valid for image …` | ARM (aarch64) or GPU image selected | Pick an x86 build; choose the shape first, then the image |
| `Your resource creation limit has been reached` | Trial account limits (not missing credit) | Fewer OCPUs (4 worked), remove leftovers, check **Limits, Quotas and Usage**, request an increase |
| Subscription shows R$0.00 | That column is *usage*; the credit is "Total commitment value" (R$1,500) | Nothing to fix |
| `No shapes are available for this image` | Same limits, for that image | Try Oracle Linux 9 |
| "The instance requires a public IP address" | VM in a private subnet / no public IP | VCN wizard network, public subnet, auto-assign public IPv4 |
| `Permission denied (publickey…)` from a script | Key has a passphrase, no agent in that shell | `ssh-add` in the terminal you run commands from |
| Attach block volume unavailable | No volume / trial restriction | Enlarge the boot volume instead (step 5) |
| `oci-growfs`: `NOCHANGE: partition 3 could only be grown by 2015` | Kernel still sees the old disk size | Rescan (`echo 1 > /sys/class/block/sda/device/rescan`) or reboot |
| Bootstrap: `install: cannot stat 'nextflow'` | The get.nextflow.io installer put the launcher elsewhere | Fixed: the script downloads the launcher from the GitHub release |
| `ERROR ~ .nextflow/history.lock (No such file or directory)` | `/data` owned by root | Fixed in the script (`/data` owned by `opc`); by hand: `sudo chown opc:opc /data` |
| `tput: unknown terminal "unknown"` | No terminal in a non-interactive `ssh` command | Harmless |
| `Nextflow 26.10.6 is available` | Notice only | Keep 26.10.5 |
