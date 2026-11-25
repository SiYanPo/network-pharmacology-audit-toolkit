# ============================================================================
# network-pharmacology-pro 容器化环境（v3.0）
#
# 构建： docker build -t np-pro:latest .
# 跑示例（离线，不依赖外部 API）：
#   docker run --rm -v $PWD:/work -w /work np-pro:latest \
#     python scripts/np_pipeline.py --inputs examples/inputs.demo.json -o /tmp/out --offline
#
# 设计要点（评审 §7）
#  1) 构建期**真跑一遍流水线**，而不是只 --version 打印工具存在——
#     "工具装了但流程跑不通"是最常见的假绿。
#  2) 运行期自检用 np_preflight（七类体检），失败要能定位到具体哪一类。
#  3) 不写死任何作者机器路径：容器内全部走 PATH，符合 np_tools 的解析顺序。
# ============================================================================
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    NP_OFFLINE=1 \
    CONDA_DIR=/opt/conda \
    PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin

# 1. 系统依赖（含 GROMACS + BLAST+；vina 走 pip 的 Python 绑定）
RUN apt-get update -qq && apt-get install -y -qq --no-install-recommends \
        python3 python3-pip python3-venv wget curl git ca-certificates \
        build-essential cmake pkg-config libboost-all-dev libeigen3-dev \
        gromacs ncbi-blast+ && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# 2. Miniforge（acpype/isPcr 通过 conda 装最稳）
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/mf.sh && \
    bash /tmp/mf.sh -b -p ${CONDA_DIR} && \
    rm /tmp/mf.sh && \
    ${CONDA_DIR}/bin/conda clean -afy

# 3. Python 栈：锁 Python 3.11（scanpy/numba 生态完整；3.13 装不上，见 PORTABILITY.md）
RUN conda install -y -q -c conda-forge \
        python=3.11 numpy pandas scipy matplotlib rdkit meeko openbabel \
        acpype parmed pyyaml h5py && \
    conda install -y -q -c bioconda ispcr && \
    pip install -q --no-cache-dir scanpy anndata leidenalg igraph vina pytest coverage

# 4. 代码
WORKDIR /work
COPY . /work

# 5. 构建期冒烟：工具存在 + preflight + **真跑示例流水线**（离线）
#    注意：vina 以 pip Python 绑定形态提供（无 CLI），所以检查 import 而非 --version。
RUN set -e && \
    python -c "import vina; print('vina py-OK')" && obabel -V && gmx --version | head -2 && \
    python -c "import scanpy, anndata, rdkit, meeko; print('python core OK')" && \
    python scripts/np_tools.py list && \
    (python scripts/np_preflight.py || true) && \
    python scripts/np_pipeline.py --inputs examples/inputs.demo.json -o /tmp/smoke --offline && \
    python scripts/np_contracts.py check --dir /tmp/smoke && \
    python scripts/np_fdr.py run --dir /tmp/smoke && \
    python scripts/np_evidence_ledger.py auto --dir /tmp/smoke && \
    test -s /tmp/smoke/09_evidence/evidence_ledger.csv && \
    echo "BUILD SMOKE OK"

HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
    CMD python scripts/np_tools.py resolve vina || exit 1

CMD ["bash"]
