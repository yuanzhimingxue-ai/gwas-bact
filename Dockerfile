# ---------------------------------------------------------
# Dockerfile  (放在仓库根目录 alongside environment.yml/gwas.sh)
# ---------------------------------------------------------
FROM mambaorg/micromamba:1.5.8

# 让 micromamba 在 RUN/ENTRYPOINT 中自动激活环境
SHELL ["/bin/bash", "-lc"]
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# ===== env =====
COPY environment.yml /tmp/env.yml
RUN micromamba create -y -n gwas -f /tmp/env.yml && micromamba clean -a -y
ENV PATH="/opt/conda/envs/gwas/bin:${PATH}"

# ===== script =====
# 关键修正：用 --chmod 赋执行权限，顺便避免再单独 RUN chmod
COPY --chmod=0755 gwas.sh /usr/local/bin/gwas.sh

# 入口：只需传一个“输入目录”参数（容器内路径）
ENTRYPOINT ["/usr/local/bin/gwas.sh"]
