# ===== base =====
FROM mambaorg/micromamba:1.5.8

# 让后续 RUN/ENTRYPOINT 自动激活 micromamba 环境
SHELL ["/bin/bash", "-lc"]
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# ===== env =====
COPY environment.yml /tmp/env.yml
RUN micromamba create -y -n gwas -f /tmp/env.yml && micromamba clean -a -y

# 方便直接找到工具
ENV PATH=/opt/conda/envs/gwas/bin:$PATH

# ===== script =====
COPY gwas.sh /usr/local/bin/gwas.sh
RUN chmod +x /usr/local/bin/gwas.sh

# 入口：只需传一个“输入目录”参数（容器内路径）
ENTRYPOINT ["/usr/local/bin/gwas.sh"]

