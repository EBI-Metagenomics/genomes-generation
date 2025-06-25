FROM mambaorg/micromamba:2.1.1

# Install git (needed for pip git+ installs)
USER root
RUN apt-get update && apt-get install -y git procps && rm -rf /var/lib/apt/lists/*

# Switch back to micromamba user
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml

RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

WORKDIR /home/$MAMBA_USER
SHELL ["/bin/bash", "-c"]
ENV PATH="/opt/conda/bin:$PATH"
