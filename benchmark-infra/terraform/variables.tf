variable "billing_account_id" {
  description = "GCP billing account ID (e.g. 012345-6789AB-CDEF01)"
  type        = string
}

variable "project_id" {
  description = "GCP project ID to create (must be globally unique)"
  type        = string
  default     = "vep-benchmark"
}

variable "project_name" {
  description = "Human-readable project name"
  type        = string
  default     = "VEP Benchmark"
}

variable "org_id" {
  description = "GCP organization ID (leave empty for projects outside an org)"
  type        = string
  default     = ""
}

variable "region" {
  description = "GCP region"
  type        = string
  default     = "europe-west1"
}

variable "zone" {
  description = "GCP zone"
  type        = string
  default     = "europe-west1-b"
}

variable "machine_type" {
  description = "GCP machine type (8 vCPUs, 64 GB RAM)"
  type        = string
  default     = "c3-highmem-8"
}

variable "disk_size_gb" {
  description = "Boot disk size in GB (needs space for VEP cache + data)"
  type        = number
  default     = 200
}

variable "disk_provisioned_iops" {
  description = "Provisioned IOPS for hyperdisk-balanced (baseline 3000, max 160000)"
  type        = number
  default     = 3000
}

variable "disk_provisioned_throughput" {
  description = "Provisioned throughput in MiB/s for hyperdisk-balanced (baseline 140, max 2400)"
  type        = number
  default     = 140
}

variable "ssh_user" {
  description = "SSH username for the VM"
  type        = string
  default     = "benchmark"
}

variable "ssh_pub_key_path" {
  description = "Path to SSH public key"
  type        = string
  default     = "~/.ssh/id_rsa.pub"
}

variable "ssh_priv_key_path" {
  description = "Path to SSH private key (for Ansible)"
  type        = string
  default     = "~/.ssh/id_rsa"
}

# ── Data paths (local machine) ───────────────────────────────────────
variable "local_vcf_dir" {
  description = "Local path to VCF files directory (contains HG002_chr1.vcf.gz)"
  type        = string
}

variable "local_fasta_path" {
  description = "Local path to reference FASTA file"
  type        = string
}

variable "local_chr1_vep_dir" {
  description = "Local path to chr1-vep parquet cache directory"
  type        = string
}

variable "local_ensembl_cache_dir" {
  description = "Local path to Ensembl VEP cache (homo_sapiens/115_GRCh38)"
  type        = string
}

variable "local_repo_dir" {
  description = "Local path to datafusion-bio-functions repo root"
  type        = string
}
