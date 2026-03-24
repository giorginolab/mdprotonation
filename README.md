# mdprotonation

`mdprotonation` is a Streamlit web app for exploring protein protonation states as a function of pH. It runs [PROPka](https://github.com/jensengroup/propka) through its Python API, computes continuous protonation fractions for titratable sites, and exposes the result through an interactive structure viewer and residue tables.

Live demo: [Hugging Face Space](https://huggingface.co/spaces/tonigi/pkaScope)

## Current Features

- Load an example PDB from [`examples/`](/Users/toni/work/mdprotonation/examples), upload a local `.pdb` file, or fetch a structure by RCSB PDB ID.
- Run PROPka directly in Python and cache the parsed result in Streamlit.
- Sweep pH with a slider from 0.0 to 14.0.
- Inspect per-site pKa, protonated fraction, current charge, and a simple state label.
- View folded and unfolded charge profiles across pH.
- Inspect raw PROPka summary and determinant output in a dedicated pane.

## Running Locally

This repository uses `uv` for Python workflows.

```bash
git clone --recurse-submodules git@github.com:giorginolab/pkaScope.git
cd pkaScope
export NODE_OPTIONS=--openssl-legacy-provider
cd mdprotonation/vendor/streamlit-molstar && bash build.sh && cd ../../..
uv sync
uv run streamlit run main.py
```

If you already cloned the repository without submodules, run `git submodule update --init --recursive` once before building the vendored `streamlit-molstar` frontend and running `uv sync`.

Then open the local Streamlit URL shown in the terminal.

## How the First Viewer Version Works

The current viewer uses `streamlit-molstar` to render the structure. That component does not yet expose a convenient Python-side API for per-residue coloring or selection, so the app currently rewrites the loaded PDB before display:

- occupancy stores the average protonated fraction for each titratable residue
- B-factor stores transition intensity scaled to `0-100`

This keeps the pH-dependent state available in the structure artifact now, while leaving room for richer visual styling later.

## Project Layout

- [`main.py`](/Users/toni/work/mdprotonation/main.py): Streamlit app entrypoint
- [`mdprotonation/propka_analysis.py`](/Users/toni/work/mdprotonation/mdprotonation/propka_analysis.py): PROPka execution and result extraction
- [`mdprotonation/protonation.py`](/Users/toni/work/mdprotonation/mdprotonation/protonation.py): pH-dependent protonation and viewer encoding logic
- [`examples/7bcq.pdb`](/Users/toni/work/mdprotonation/examples/7bcq.pdb): example structure for local testing

## Next Steps

Planned additions include richer Mol* residue highlighting, chain and residue filters, structured PROPka tables, and additional panes for deeper inspection of PROPka outputs.
