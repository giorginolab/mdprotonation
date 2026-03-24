FROM python:3.13-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    UV_LINK_MODE=copy \
    NODE_OPTIONS=--openssl-legacy-provider \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_SERVER_PORT=7860 \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0

ARG APP_REPO=https://github.com/giorginolab/pkaScope.git
ARG APP_REF=main

RUN apt-get update \
    && apt-get install --yes --no-install-recommends git nodejs npm \
    && rm -rf /var/lib/apt/lists/*

RUN python -m pip install --no-cache-dir uv

WORKDIR /opt
RUN git clone --depth 1 --branch "${APP_REF}" --recurse-submodules "${APP_REPO}" app

WORKDIR /opt/app
RUN cd mdprotonation/vendor/streamlit-molstar && bash build.sh
RUN uv sync --frozen --no-dev

ENV PATH="/opt/app/.venv/bin:${PATH}"

EXPOSE 7860

CMD ["streamlit", "run", "main.py", "--server.port=7860", "--server.address=0.0.0.0"]
