# Plot an `odgi` pangenome

The `odgi draw` SVG output is fairly limited currently, so this is a small tool as a stopgap. It requires two files, given an `in.og`:

1. A layout TSV from e.g. `odgi layout -i in.og -T out.layout.tsv`
2. A GFA from the same pangenome, e.g. `odgi view -i in.og -g > out.gfa`

## Main features

- Highlights specific (optionally multiple) paths through the pangenome. These are the `P` lines in the GFA.
- Clips x axes so that extreme values are compressed.
- Can subsample both edges and paths so the output SVG is smaller.

## Usage

For my use case, I used the following command:

```bash
# 1741217 lines in the layout tsv
odgi_path_svg -L data/component.717.layout.tsv \
  -g data/component.717.gfa \
  -o Logfia_minima.svg \
  -H Logfia_minima_u19 \
  --edge-simplify-eps 3 \
  --x-clip-low 0.05
```

Full usage below:

```txt
Render ODGI layouts + GFA paths as SVG, with optional path highlights

Usage: odgi_path_svg [OPTIONS] --layout <LAYOUT_TSV> --gfa <GFA> --output <SVG>

Options:
  -L, --layout <LAYOUT_TSV>     Layout TSV from `odgi layout -T`
  -g, --gfa <GFA>               GFA exported with `odgi view -g`
  -o, --output <SVG>            Output SVG file
  -H, --highlight <PATH>...     Path name to highlight (can be given multiple times)
      --width <PX>              Target SVG width in pixels [default: 1200]
      --padding <PX>            Padding in pixels around drawing [default: 20]
      --node-radius <PX>        Radius of background node circles (px) [default: 0.6]
      --max-aspect <RATIO>      Maximum allowed width/height ratio for the graph area [default: 3.0]
      --no-edges                Do not draw edges (only highlighted paths)
      --x-clip-low <PCT>        Lower X clip percentile (0-100); trims extreme left outliers [default: 0.0]
      --x-clip-high <PCT>       Upper X clip percentile (0-100); trims extreme right outliers [default: 100.0]
      --path-simplify-eps <PX>  Simplify highlighted paths: minimum distance (px) between kept points; 0 disables simplification [default: 0.0]
      --edge-simplify-eps <PX>  Simplify skeleton edges: minimum distance (px) between edge endpoints to draw; 0 draws all edges [default: 0.0]
  -h, --help                    Print help
  -V, --version                 Print version
```
