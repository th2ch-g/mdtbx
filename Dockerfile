FROM --platform=linux/amd64 ubuntu

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    ca-certificates \
    git \
    util-linux \
    && rm -rf /var/lib/apt/lists/*

ENV PIXI_HOME="/root/.pixi"
ENV PATH="${PIXI_HOME}/bin:${PATH}"
RUN curl -fsSL https://pixi.sh/install.sh | setarch x86_64 bash

WORKDIR /app

COPY ./pixi.lock .
COPY ./pyproject.toml .
COPY ./src/ ./src/

RUN setarch x86_64 pixi install --frozen

ENV PYTHONPATH /app

ENTRYPOINT ["pixi", "run"]
CMD ["bash"]
