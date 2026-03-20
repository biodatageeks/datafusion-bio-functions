#!/usr/bin/env bash
set -euo pipefail

# ── Single-command setup: terraform does everything including Ansible ──
# Usage: ./setup.sh [apply|destroy]

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INFRA_DIR="$(dirname "$SCRIPT_DIR")"
TF_DIR="${INFRA_DIR}/terraform"
ANSIBLE_DIR="${INFRA_DIR}/ansible"
VENV_DIR="${INFRA_DIR}/.venv"

ACTION="${1:-apply}"

# ── Python venv via uv ───────────────────────────────────────────────
if ! command -v uv &>/dev/null; then
  echo "ERROR: uv not found. Install it: curl -LsSf https://astral.sh/uv/install.sh | sh"
  exit 1
fi

if [ ! -d "${VENV_DIR}" ]; then
  echo "==> Creating venv with uv..."
  uv venv "${VENV_DIR}"
  echo "==> Installing Python dependencies..."
  uv pip install --python "${VENV_DIR}/bin/python" -r "${INFRA_DIR}/requirements.txt"
fi

# Activate venv for this script (ansible-playbook, ansible-galaxy)
source "${VENV_DIR}/bin/activate"

# ── Check prerequisites ──────────────────────────────────────────────
if ! command -v terraform &>/dev/null; then
  echo "ERROR: terraform not found. Install it: brew install terraform"
  exit 1
fi

if [ ! -f "${TF_DIR}/terraform.tfvars" ]; then
  echo "ERROR: ${TF_DIR}/terraform.tfvars not found."
  echo "Copy terraform.tfvars.example and fill in your values."
  exit 1
fi

cd "${TF_DIR}"

if [ "$ACTION" = "destroy" ]; then
  echo "==> Destroying infrastructure..."
  terraform destroy -auto-approve
  exit 0
fi

echo "==> Installing Ansible Galaxy roles..."
ansible-galaxy install -r "${ANSIBLE_DIR}/requirements.yml" --force

echo "==> Initializing Terraform..."
terraform init

echo "==> Applying (VM + Ansible provisioning)..."
terraform apply -auto-approve

echo ""
echo "============================================"
echo "  Benchmark VM ready!"
echo "  $(terraform output -raw ssh_command)"
echo "  Run: ~/run_benchmark.sh"
echo "============================================"
