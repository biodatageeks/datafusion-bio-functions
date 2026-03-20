terraform {
  required_version = ">= 1.5"

  required_providers {
    google = {
      source  = "hashicorp/google"
      version = "~> 5.0"
    }
    ansible = {
      source  = "ansible/ansible"
      version = "~> 1.3"
    }
    null = {
      source  = "hashicorp/null"
      version = "~> 3.0"
    }
  }
}

# ── Providers ────────────────────────────────────────────────────────
# "bootstrap" has no project — used only for creating the project itself.
# Default provider uses the created project for all other resources.

provider "google" {
  alias  = "bootstrap"
  region = var.region
  zone   = var.zone
}

provider "google" {
  project = var.project_id
  region  = var.region
  zone    = var.zone
}

# ── GCP Project (survives terraform destroy) ─────────────────────────

resource "google_project" "benchmark" {
  provider = google.bootstrap

  name            = var.project_name
  project_id      = var.project_id
  billing_account = var.billing_account_id
  org_id          = var.org_id != "" ? var.org_id : null

  lifecycle {
    prevent_destroy = true
  }
}

resource "google_project_service" "compute" {
  project = google_project.benchmark.project_id
  service = "compute.googleapis.com"

  disable_on_destroy = false
}

# ── GCP Resources ────────────────────────────────────────────────────

resource "google_compute_firewall" "allow_ssh" {
  name    = "vep-benchmark-allow-ssh"
  network = "default"
  project = google_project.benchmark.project_id

  allow {
    protocol = "tcp"
    ports    = ["22"]
  }

  source_ranges = ["0.0.0.0/0"]
  target_tags   = ["vep-benchmark"]

  depends_on = [google_project_service.compute]
}

resource "google_compute_instance" "vep_benchmark" {
  name         = "vep-benchmark"
  machine_type = var.machine_type
  zone         = var.zone
  project      = google_project.benchmark.project_id

  tags = ["vep-benchmark"]

  boot_disk {
    initialize_params {
      image                  = "ubuntu-os-cloud/ubuntu-2404-lts-amd64"
      size                   = var.disk_size_gb
      type                   = "hyperdisk-balanced"
      provisioned_iops       = var.disk_provisioned_iops
      provisioned_throughput = var.disk_provisioned_throughput
    }
  }

  network_interface {
    network = "default"
    access_config {
      # Ephemeral public IP
    }
  }

  scheduling {
    preemptible        = true
    automatic_restart  = false
    provisioning_model = "SPOT"
  }

  metadata = {
    ssh-keys = "${var.ssh_user}:${file(var.ssh_pub_key_path)}"
  }

  service_account {
    scopes = ["cloud-platform"]
  }

  depends_on = [google_project_service.compute]
}

# ── Wait for SSH ─────────────────────────────────────────────────────

resource "null_resource" "wait_for_ssh" {
  triggers = {
    instance_id = google_compute_instance.vep_benchmark.id
  }

  connection {
    type        = "ssh"
    user        = var.ssh_user
    private_key = file(var.ssh_priv_key_path)
    host        = google_compute_instance.vep_benchmark.network_interface[0].access_config[0].nat_ip
    timeout     = "10m"
  }

  provisioner "remote-exec" {
    inline = ["echo 'SSH is ready'"]
  }
}

# ── Ansible Integration ──────────────────────────────────────────────

resource "ansible_host" "vep_benchmark" {
  name   = google_compute_instance.vep_benchmark.network_interface[0].access_config[0].nat_ip
  groups = ["benchmark"]

  variables = {
    ansible_user                 = var.ssh_user
    ansible_ssh_private_key_file = var.ssh_priv_key_path
    ansible_ssh_common_args      = "-o StrictHostKeyChecking=no"

    # Data paths — passed through to playbook as host vars
    local_vcf_dir           = var.local_vcf_dir
    local_fasta_path        = var.local_fasta_path
    local_chr1_vep_dir      = var.local_chr1_vep_dir
    local_ensembl_cache_dir = var.local_ensembl_cache_dir
    local_repo_dir          = var.local_repo_dir
  }
}

resource "ansible_playbook" "provision" {
  playbook   = "${path.module}/../ansible/playbook.yml"
  name       = google_compute_instance.vep_benchmark.network_interface[0].access_config[0].nat_ip
  replayable = true

  extra_vars = {
    ansible_user                 = var.ssh_user
    ansible_ssh_private_key_file = var.ssh_priv_key_path
    ansible_ssh_common_args      = "-o StrictHostKeyChecking=no"
    local_vcf_dir                = var.local_vcf_dir
    local_fasta_path             = var.local_fasta_path
    local_chr1_vep_dir           = var.local_chr1_vep_dir
    local_ensembl_cache_dir      = var.local_ensembl_cache_dir
    local_repo_dir               = var.local_repo_dir
  }

  timeouts {
    create = "120m"
  }

  depends_on = [
    null_resource.wait_for_ssh,
    ansible_host.vep_benchmark,
  ]
}

# ── Outputs ──────────────────────────────────────────────────────────

output "project_id" {
  description = "GCP project ID (protected from destroy)"
  value       = google_project.benchmark.project_id
}

output "instance_ip" {
  description = "Public IP of the benchmark VM"
  value       = google_compute_instance.vep_benchmark.network_interface[0].access_config[0].nat_ip
}

output "ssh_command" {
  description = "SSH command to connect"
  value       = "ssh ${var.ssh_user}@${google_compute_instance.vep_benchmark.network_interface[0].access_config[0].nat_ip}"
}

output "instance_name" {
  value = google_compute_instance.vep_benchmark.name
}
