"""Compile the editable Typst deck and export each schematic as vector PDF.

Requires typst-py (`python -m pip install typst`) or a `typst` CLI on PATH.
The deck uses Touying 0.7.4 (downloaded by Typst on first build).
"""
from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
SCHEMATICS = ("strip-geometry", "solver-flow", "root-events", "amoeba-search")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--with-data", action="store_true")
    parser.add_argument("--author", default="pyGBZ2d research project")
    parser.add_argument("--institution", default=None)
    args = parser.parse_args()
    if args.with_data:
        required = [ROOT / "Figures" / (name + ".svg")
                    for name in ("root-tracks", "subsets", "spectrum")]
        for path in required:
            if not path.is_file():
                parser.error(f"missing computed figure: {path}; run plot_slide_data.py")
    build = ROOT / ".build"
    build.mkdir(exist_ok=True)
    # Optional project-local compiler installed for this authoring session.
    if (build / "python").exists():
        sys.path.insert(0, str(build / "python"))
    try:
        import typst
    except ImportError:
        typst = None
    cli = shutil.which("typst")
    if typst is None and cli is None:
        parser.error("install typst-py or add the Typst compiler to PATH")

    def compile_file(source, output, inputs=None):
        if typst is not None:
            typst.compile(str(source), output=str(output), root=str(ROOT),
                          sys_inputs=inputs or {},
                          package_cache_path=str(build / "packages"))
        else:
            command = [cli, "compile", "--root", str(ROOT)]
            command += ["--package-cache-path", str(build / "packages")]
            for k, v in (inputs or {}).items():
                command += ["--input", f"{k}={v}"]
            subprocess.run(command + [str(source), str(output)], check=True)

    for name in SCHEMATICS:
        wrapper = build / f"{name}.typ"
        wrapper.write_text(
            '#set page(width: auto, height: auto, margin: 0pt)\n'
            f'#image("../Figures/{name}.svg")\n', encoding="utf-8")
        compile_file(wrapper, ROOT / "Figures" / f"{name}.pdf")
    output = ROOT / ("gbz-algorithms-with-data.pdf" if args.with_data else "gbz-algorithms.pdf")
    inputs = {"with-data": str(args.with_data).lower(), "author": args.author}
    if args.institution:
        inputs["institution"] = args.institution
    compile_file(ROOT / "gbz-algorithms.typ", output, inputs)
    print(output)


if __name__ == "__main__":
    main()
