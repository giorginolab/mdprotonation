# pkaScope

`pkaScope` is a Streamlit web app for exploring protein protonation
states as a function of pH. It runs
[PROPka](https://github.com/jensengroup/propka) through its Python API,
computes continuous protonation fractions for titratable sites, and
exposes the result through an interactive structure viewer and residue
tables.

Live demo: [Hugging Face Space](https://huggingface.co/spaces/tonigi/pkaScope)

## Current Features

- Load an example PDB from [`examples/`](/Users/toni/work/pkaScope/examples), upload a local `.pdb` file, or fetch a structure by RCSB PDB ID.
- Run PROPka directly in Python and cache the parsed result in Streamlit.
- Sweep pH with a slider from 0.0 to 14.0.
- Inspect per-site pKa, protonated fraction, current charge, and a simple state label.
- View folded and unfolded charge profiles across pH.
- Run `pdb2pqr` hydrogen-bond placement/optimization from a dedicated pane.
- Compare before/after structures and download optimized `.pdb` and `.pqr` outputs.
- Inspect raw PROPka summary and determinant output in a dedicated pane.

## Running Locally

This repository uses `uv` for Python workflows.

```bash
git clone git@github.com:giorginolab/pkaScope.git
cd pkaScope
uv sync
uv run streamlit run main.py
```

Then open the local Streamlit URL shown in the terminal.

## H-Bond Optimization Pane

The `H-Bond Optimize` tab runs `pdb2pqr` on demand when you click
`Run pdb2pqr optimization`. It uses the same global pH slider as the
other panes and caches successful results per structure and pH.

After completion, the pane shows:

- run summary metrics
- before/after Mol* previews
- download buttons for optimized PDB and PQR files
- a compact diagnostics section

## How the First Viewer Version Works

The current viewer uses `molviewspec` with the Streamlit interface shown
in `vendor/mol-view-spec/test-data/streamlit`. The app rewrites the
loaded PDB before display so per-residue protonation state is encoded in
the structure:

- occupancy stores the average protonated fraction for each titratable residue
- B-factor stores transition intensity scaled to `0-100`

This keeps the pH-dependent state available in the structure artifact now, while leaving room for richer visual styling later.
