---
title: pkaScope
sdk: docker
app_port: 7860
pinned: false
---

# pkaScope

This Space builds the app directly from [giorginolab/pkaScope](https://github.com/giorginolab/pkaScope).
The image clones the repository with git submodules so the vendored `streamlit-molstar` checkout is used during `uv sync`.
