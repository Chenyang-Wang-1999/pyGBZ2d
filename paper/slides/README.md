# Typst research talk: 2D GBZ algorithms

`gbz-algorithms.typ` is the editable English deck for a 20–25 minute physics research talk. `gbz-algorithms.pdf` is the compiled 16:9 version. It assumes basic NHSE, 1D non-Bloch band theory, and amoeba knowledge. The existing `pyGBZ2d.pptx` is unchanged.

The deck uses **Touying 0.7.4, University theme**, following the supplied academic example. It has serif typography, a centered title page, numbered section headings, a progress bar, and three-part academic footers. The theme call and slide helpers live directly in `gbz-algorithms.typ`, so the deck no longer depends on a local `theme.typ` file. Equations, slide text, and speaker notes remain editable. The main PDF has 20 pages, including an unnumbered title page and 19 numbered slides.

Only Touying and its uniwarn dependency are needed. The supplied example's CeTZ, Fletcher, and theorem/numbering packages are not imported unnecessarily: the existing SVG schematics remain editable vector assets, and native Typst heading numbering supplies the section numbers. `#speaker-note` is embedded in the source. Handout mode currently keeps one final view per slide; you can add Touying `#pause`/`#uncover` animation later and disable handout mode.

## Build

Use either the Typst CLI or its Python binding:

```powershell
python -m pip install typst
python paper/slides/build_slides.py
# Optional title-page metadata:
python paper/slides/build_slides.py --author "Your name" --institution "Your institution"
```

Or compile the main deck directly:

```powershell
typst compile paper/slides/gbz-algorithms.typ paper/slides/gbz-algorithms.pdf
```

On the first build, Typst downloads the pinned Touying package and its dependency. The builder uses `.build/packages` as the project-local package cache. Subsequent builds can use that cache offline. The builder also regenerates each schematic PDF from its SVG. It prefers an installed Python binding and falls back to a `typst` executable on PATH. The authoring-session compiler in `.build/python` is optional and ignored by Git. Compilation was checked with Typst 0.15.0. Use a current Typst release if older syntax support differs.

## Vector assets

`Figures/` contains four editable SVG schematics and matching vector PDFs:

- `strip-geometry`: chosen major direction and a finite transverse strip.
- `solver-flow`: the shared root engine and the two modulus constraints.
- `root-events`: equal modulus, multiple roots, and continuum detection.
- `amoeba-search`: inner mu2 search and outer mu1 update.

These are conceptual diagrams. None depicts a computed trajectory, spectrum, measured boundary, or performance curve. Diagram text stays editable in the SVG source; PDF versions embed vector paths/text.

## Numerical plots

See `DATA_REQUIREMENTS.md` for required model arrays, diagnostic E/mu1, energy grid, and example commands. `collect_slide_data.py` performs real calculations. `plot_slide_data.py` exports SVG/PDF plots from the saved JSON, with provenance and failure labels. Plotting requires NumPy and Matplotlib; collection additionally uses the repository's SciPy-based solver.

After producing all three numerical figure pairs, run `build_slides.py --with-data` to append them. The main PDF is not overwritten by this option. The optional data deck is `gbz-algorithms-with-data.pdf`.

`speaker-notes.md` gives timing and the physical explanation behind the equations. `sources.md` maps slides to theory and implementation. The main narrative emphasizes physical selection and numerical design: the Jensen/root-count reduction of the Ronkin function, and the middle-root/loop-winding reduction of the strip problem. Code-level details remain in the source map rather than the main slides. Older manuscript pseudocode does not override current code behavior.

