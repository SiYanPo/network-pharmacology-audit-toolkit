#!/usr/bin/env bash
# ============================================================================
# network-pharmacology-pro 一键环境构建（Ubuntu 22.04+ / WSL Ubuntu-24.04）
# 用法： bash setup_env.sh [--with-test-data]
# 静默安装：scanpy/anndata, AutoDock Vina, OpenBabel, GROMACS, BLAST+, isPcr, acpype
# 幂等：重复运行安全；所有安装加 -y，无需交互。
# ============================================================================
set -euo pipefail
LOG="np_env_install.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== 环境构建开始 $(date) ==="

WITH_DATA=0
for a in "$@"; do [ "$a" = "--with-test-data" ] && WITH_DATA=1; done

# ---- 0. 系统包 ----
echo "[0/6] 系统依赖"
export DEBIAN_FRONTEND=noninteractive
sudo apt-get update -qq
sudo apt-get install -y -qq --no-install-recommends \
    python3 python3-pip python3-venv wget curl git ca-certificates \
    build-essential cmake pkg-config libboost-all-dev libeigen3-dev \
    gromacs ncbi-blast+ 2>/dev/null || {
    echo "[warn] 部分 apt 包失败，尝试单独安装 GROMACS/BLAST"
    sudo apt-get install -y -qq gromacs || true
    sudo apt-get install -y -qq ncbi-blast+ || true
}

# ---- 1. Python 虚拟环境 ----
echo "[1/6] Python venv + 核心包"
VENV="${NP_VENV:-$PWD/.npvenv}"
python3 -m venv "$VENV"
source "$VENV/bin/activate"
pip install -q --upgrade pip setuptools wheel
pip install -q numpy pandas scipy matplotlib scanpy anndata leidenalg igraph \
    rdkit meeko openbabel-wheel vina

# ---- 2. Conda 工具（acpype / isPcr）----
echo "[2/6] acpype + isPcr（conda-forge / bioconda）"
if ! command -v conda >/dev/null && ! command -v mamba >/dev/null; then
    echo "[info] 未检测到 conda，安装 Miniforge"
    wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh
    bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
    export PATH="$HOME/miniforge3/bin:$PATH"
    echo 'export PATH="$HOME/miniforge3/bin:$PATH"' >> ~/.bashrc
fi
CONDA=$(command -v mamba || command -v conda)
$CONDA install -y -q -c conda-forge acpype ParmEd
$CONDA install -y -q -c bioconda ispcr

# ---- 3. AutoDock Vina / OpenBabel 校验 ----
echo "[3/6] 校验 Vina / OpenBabel"
vina --version
obabel -V

# ---- 4. GROMACS 校验 ----
echo "[4/6] 校验 GROMACS"
if command -v gmx >/dev/null; then gmx --version | head -3; else echo "[warn] gmx 未找到，请手动 apt install gromacs"; fi

# ---- 5. 测试数据 ----
if [ "$WITH_DATA" -eq 1 ]; then
    echo "[5/6] 下载测试小数据集"
    mkdir -p tests/data
    # 1CRN 小蛋白（对接/MD 演示）
    [ -f tests/data/1crn.pdb ] || wget -q https://files.rcsb.org/download/1CRN.pdb -O tests/data/1crn.pdb
    # 阿司匹林 SDF（配体演示）
    [ -f tests/data/aspirin.sdf ] || wget -q https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/SDF?record_type=3d -O tests/data/aspirin.sdf
    # 极简 GMT（富集演示）
    cat > tests/data/demo.gmt <<'EOF'
APOPTOSIS	TP53	CASP3	BAX	BCL2
INFLAMMATION	IL6	TNF	IL1B	PTGS2
EOF
fi

# ---- 6. 写环境标记 + 自检 ----
echo "[6/6] 自检"
python scripts/np_envcheck.py --json np_env_after_setup.json || true
echo "=== 环境构建完成 $(date) ==="
echo "激活环境： source $VENV/bin/activate"
echo "运行自检： python scripts/np_envcheck.py"
