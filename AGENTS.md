# Repository Guidelines

## Project Structure & Module Organization
This repository is being shaped into a Streamlit app for exploring protein protonation states as a function of pH. Keep [`main.py`](/Users/toni/work/mdprotonation/main.py) as the Streamlit entrypoint, and move reusable logic into a package such as `src/mdprotonation/` as the app grows. Use `examples/` for sample structures like [`examples/7bcq.pdb`](/Users/toni/work/mdprotonation/examples/7bcq.pdb), and keep design notes in top-level Markdown files such as [`streamlist-molstar-howto.md`](/Users/toni/work/mdprotonation/streamlist-molstar-howto.md).

Recommended modules:
- `propka.py` for PROPka API calls and parsed residue pKa data.
- `protonation.py` for pH-to-state logic.
- `viewer.py` for Mol* rendering and residue highlighting.
- `panes/` for UI panels; the first pane is the pH slider and structure viewer, with additional PROPka output panes added later.

## Architecture Plan
Build the first version around one uploaded or example PDB, processed through the PROPka Python API. The main interaction is a pH slider that updates residue protonation continuously from residue pKa values and pushes those state changes into the structure viewer. Keep calculation code separate from Streamlit widgets so later panes can reuse the same parsed PROPka result, residue table, and state model.

## Build, Test, and Development Commands
Use `uv` for all Python and package workflows.

- `uv sync`: install or refresh dependencies from [`pyproject.toml`](/Users/toni/work/mdprotonation/pyproject.toml) and [`uv.lock`](/Users/toni/work/mdprotonation/uv.lock).
- `uv run streamlit run main.py`: launch the local web app.
- `uv run python main.py`: run the lightweight entrypoint when debugging non-UI code.
- `uv run python -m pytest`: run tests once `tests/` exists.

## Coding Style & Naming Conventions
Target Python 3.13+, use 4-space indentation, and follow PEP 8 naming: `snake_case` for functions/modules and `PascalCase` for classes. Use type hints on any function that transforms PROPka output or residue state. Prefer pure functions for protonation logic; the UI layer should orchestrate input, caching, and display, not chemical calculations.

## Testing Guidelines
Add tests under `tests/` with names like `test_propka.py` and `test_protonation.py`. Prioritize unit tests for residue parsing, pKa mapping, and pH transition behavior before adding UI tests. Include edge cases around terminal residues, missing pKa values, and repeated slider updates.

## Commit & Pull Request Guidelines
The repository has no established Git history yet, so start with short imperative commits such as `Add PROPka residue parser` or `Wire pH slider to viewer state`. Pull requests should describe user-visible behavior, list verification commands, and include screenshots for Streamlit UI changes.
