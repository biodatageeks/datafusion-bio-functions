# VEP Benchmark Infrastructure (GCP)

Terraform + Ansible (via `ansible/ansible` provider) to run VEP benchmarks on a GCP VM. Single `terraform.tfvars` for all configuration.

## VM Spec

- **Machine**: `c3-highmem-8` (8 vCPUs, 64 GB RAM)
- **Disk**: 200 GB Hyperdisk Balanced (tunable IOPS/throughput)
- **OS**: Ubuntu 24.04 LTS
- **Scheduling**: SPOT (preemptible)
- **Software**: Docker, Rust 1.91, samtools, tabix, byobu

## Prerequisites

```bash
brew install terraform
curl -LsSf https://astral.sh/uv/install.sh | sh
gcloud auth application-default login
```

Ansible is installed automatically in a venv by `setup.sh` (via `uv`).

## Quick Start

### 1. Configure

```bash
cd benchmark-infra/terraform
cp terraform.tfvars.example terraform.tfvars
# Edit terraform.tfvars — set billing_account_id, project_id, SSH keys, and local data paths
```

All variables (GCP, SSH, data paths) live in a single `terraform.tfvars`.

### 2. Deploy

```bash
cd benchmark-infra/scripts
./setup.sh
```

Terraform provisions the VM, then the `ansible/ansible` provider automatically runs the playbook to install Docker + Rust, rsync data, and build the project.

### 3. Run benchmarks

```bash
ssh benchmark@<VM_IP>
```

#### Quick benchmark (automated Ensembl Docker vs native)

```bash
./run_benchmark.sh
```

#### Golden benchmark — Fjall backend

First build the fjall cache (one-time):

```bash
REPO=~/datafusion-bio-functions
DATA=/data/vep

# 1/4: Variation cache
$REPO/target/release/examples/load_cache_full \
  $DATA/chr1-vep/115_GRCh38_variation_1_vep.parquet \
  $DATA/chr1-vep/fjall_variation_cache 4

# 2/4: SIFT/PolyPhen
$REPO/target/release/examples/load_sift_cache \
  $DATA/chr1-vep/115_GRCh38_translation_sift_1_vep.parquet \
  $DATA/chr1-vep/fjall_variation_cache

# 3/4 + 4/4: Exons + Translations
$REPO/target/release/examples/load_context_kv \
  $DATA/chr1-vep/fjall_variation_cache \
  $DATA/chr1-vep/115_GRCh38_exon_1_vep.parquet \
  $DATA/chr1-vep/115_GRCh38_translation_core_1_vep.parquet
```

Then run with fjall:

```bash
mkdir -p /tmp/annotate_vep_golden_bench_chr1
cp $DATA/vcf/output/everything/HG002_chr1_0_vep115_golden.vcf \
   /tmp/annotate_vep_golden_bench_chr1/

$REPO/target/release/examples/annotate_vep_golden_bench \
  $DATA/vcf/HG002_chr1.vcf.gz \
  $DATA/chr1-vep/fjall_variation_cache \
  fjall 0 \
  $DATA/homo_sapiens \
  $DATA/vcf/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_golden_bench_chr1 \
  $DATA/chr1-vep \
  --steps=ensembl,datafusion --extended-probes --everything \
  --reference-fasta-path=$DATA/Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

#### Golden benchmark — Parquet backend

```bash
mkdir -p /tmp/annotate_vep_golden_bench_chr1_parquet
cp $DATA/vcf/output/everything/HG002_chr1_0_vep115_golden.vcf \
   /tmp/annotate_vep_golden_bench_chr1_parquet/

$REPO/target/release/examples/annotate_vep_golden_bench \
  $DATA/vcf/HG002_chr1.vcf.gz \
  $DATA/chr1-vep/115_GRCh38_variation_1_vep.parquet \
  parquet 0 \
  $DATA/homo_sapiens \
  $DATA/vcf/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_golden_bench_chr1_parquet \
  $DATA/chr1-vep \
  --steps=ensembl,datafusion --extended-probes --everything \
  --reference-fasta-path=$DATA/Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

#### Profile annotation (timing breakdown)

```bash
VEP_PROFILE=1 $REPO/target/release/examples/profile_annotation \
  $DATA/vcf/HG002_chr1.vcf.gz \
  $DATA/chr1-vep \
  320000 \
  --everything \
  --reference-fasta-path=$DATA/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --output=/tmp/HG002_chr1_annotated.vcf
```

Results saved to `/data/vep/results/`.

### 4. Tear down

```bash
./setup.sh destroy
```

Note: the GCP project is protected from destroy (`lifecycle { prevent_destroy = true }`). Only the VM and firewall are removed.

## Data transferred to the VM

| Data | VM path | Size |
|------|---------|------|
| VCF (HG002_chr1.vcf.gz + .tbi) | `/data/vep/vcf/` | ~12 MB |
| Golden VEP output | `/data/vep/vcf/output/everything/` | ~1.3 GB |
| Reference FASTA (+.fai) | `/data/vep/` | ~3 GB |
| chr1-vep parquet cache | `/data/vep/chr1-vep/` | ~3 GB |
| Ensembl VEP cache (115_GRCh38) | `/data/vep/homo_sapiens/115_GRCh38/` | ~15 GB |

## Cost

SPOT `c3-highmem-8` in `europe-west1`: ~$0.15/hr. Budget ~$1 per benchmark session.
