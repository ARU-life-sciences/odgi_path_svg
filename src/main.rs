use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::{Arg, ArgAction, Command};

#[derive(Debug, Clone, Copy)]
struct Point {
    x: f64,
    y: f64,
}

#[derive(Debug, Clone, Copy)]
struct PathStep {
    node_id: u64,
    reversed: bool,
}

#[derive(Debug, Clone)]
struct PathRec {
    name: String,
    steps: Vec<PathStep>,
}

#[derive(Debug, Clone, Copy)]
struct Edge {
    from: u64,
    to: u64,
}

#[derive(Debug, Clone)]
struct GraphData {
    paths: Vec<PathRec>,
    edges: Vec<Edge>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("odgi_path_svg")
        .version(clap::crate_version!())
        .about("Render ODGI layouts + GFA paths as SVG, with optional path highlights")
        .arg(
            Arg::new("layout")
                .long("layout")
                .short('L')
                .help("Layout TSV from `odgi layout -T`")
                .required(true)
                .value_name("LAYOUT_TSV")
                .value_parser(clap::value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("gfa")
                .long("gfa")
                .short('g')
                .help("GFA exported with `odgi view -g`")
                .required(true)
                .value_name("GFA")
                .value_parser(clap::value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("output")
                .long("output")
                .short('o')
                .help("Output SVG file")
                .required(true)
                .value_name("SVG")
                .value_parser(clap::value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("highlight")
                .long("highlight")
                .short('H')
                .help("Path name to highlight (can be given multiple times)")
                .value_name("PATH")
                .num_args(1..)
                .action(ArgAction::Append),
        )
        .arg(
            Arg::new("width")
                .long("width")
                .help("Target SVG width in pixels")
                .value_name("PX")
                .value_parser(clap::value_parser!(f64))
                .default_value("1200"),
        )
        .arg(
            Arg::new("padding")
                .long("padding")
                .help("Padding in pixels around drawing")
                .value_name("PX")
                .value_parser(clap::value_parser!(f64))
                .default_value("20"),
        )
        .arg(
            Arg::new("node_radius")
                .long("node-radius")
                .help("Radius of background node circles (px)")
                .value_name("PX")
                .value_parser(clap::value_parser!(f64))
                .default_value("0.6"),
        )
        .arg(
            Arg::new("max_aspect")
                .long("max-aspect")
                .help("Maximum allowed width/height ratio for the graph area")
                .value_name("RATIO")
                .value_parser(clap::value_parser!(f64))
                .default_value("3.0"),
        )
        .arg(
            Arg::new("no_edges")
                .long("no-edges")
                .help("Do not draw edges (only highlighted paths)")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("x_clip_low")
                .long("x-clip-low")
                .help("Lower X clip percentile (0-100); trims extreme left outliers")
                .value_name("PCT")
                .value_parser(clap::value_parser!(f64))
                .default_value("0.0"),
        )
        .arg(
            Arg::new("x_clip_high")
                .long("x-clip-high")
                .help("Upper X clip percentile (0-100); trims extreme right outliers")
                .value_name("PCT")
                .value_parser(clap::value_parser!(f64))
                .default_value("100.0"),
        )
        .arg(
            Arg::new("path_simplify_eps")
                .long("path-simplify-eps")
                .help("Simplify highlighted paths: minimum distance (px) between kept points; 0 disables simplification")
                .value_name("PX")
                .value_parser(clap::value_parser!(f64))
                .default_value("3.0"),
        )
        .arg(
            Arg::new("edge_simplify_eps")
                .long("edge-simplify-eps")
                .help("Simplify skeleton edges: minimum distance (px) between edge endpoints to draw; 0 draws all edges")
                .value_name("PX")
                .value_parser(clap::value_parser!(f64))
                .default_value("0.0"),
        )
        .get_matches();

    let layout_path: PathBuf = matches
        .get_one::<PathBuf>("layout")
        .expect("layout is required")
        .clone();
    let gfa_path: PathBuf = matches
        .get_one::<PathBuf>("gfa")
        .expect("gfa is required")
        .clone();
    let output_path: PathBuf = matches
        .get_one::<PathBuf>("output")
        .expect("output is required")
        .clone();

    let width: f64 = *matches.get_one::<f64>("width").unwrap();
    let padding: f64 = *matches.get_one::<f64>("padding").unwrap();
    let node_radius: f64 = *matches.get_one::<f64>("node_radius").unwrap();
    let max_aspect: f64 = *matches.get_one::<f64>("max_aspect").unwrap();
    let no_edges: bool = matches.get_flag("no_edges");

    let x_clip_low: f64 = *matches.get_one::<f64>("x_clip_low").unwrap();
    let x_clip_high: f64 = *matches.get_one::<f64>("x_clip_high").unwrap();
    let path_simplify_eps: f64 = *matches.get_one::<f64>("path_simplify_eps").unwrap();
    let edge_simplify_eps: f64 = *matches.get_one::<f64>("edge_simplify_eps").unwrap();

    let highlight_names: Vec<String> = matches
        .get_many::<String>("highlight")
        .map(|vals| vals.cloned().collect())
        .unwrap_or_default();

    eprintln!("[info] reading layout from {}", layout_path.display());
    let node_coords = read_layout_tsv(&layout_path)?;

    eprintln!("[info] reading GFA from {}", gfa_path.display());
    let graph = read_gfa(&gfa_path)?;

    eprintln!(
        "[info] found {} nodes, {} paths, {} edges",
        node_coords.len(),
        graph.paths.len(),
        graph.edges.len()
    );

    // Decide which paths to highlight
    let mut highlight_paths: Vec<PathRec> = Vec::new();
    if highlight_names.is_empty() {
        eprintln!("[info] no --highlight given; will draw edges only (graph skeleton)");
    } else {
        let wanted: HashSet<_> = highlight_names.iter().cloned().collect();
        for p in &graph.paths {
            if wanted.contains(&p.name) {
                highlight_paths.push(p.clone());
            }
        }
        if highlight_paths.is_empty() {
            eprintln!("[warn] none of the requested paths were found; only skeleton will be drawn");
        } else {
            eprintln!(
                "[info] {} highlighted paths will be drawn",
                highlight_paths.len()
            );
        }
    }

    // Compute bounding box
    let (min_x, max_x, min_y, max_y) = bounding_box(&node_coords)?;
    let dx = (max_x - min_x).max(1.0);
    let dy = (max_y - min_y).max(1.0);
    let aspect = dx / dy;
    eprintln!("[info] layout bounding box: dx={dx:.2}, dy={dy:.2}, aspect={aspect:.3}");

    // Aspect control: first pass (full bbox)
    let effective_dx = if aspect > max_aspect {
        eprintln!("[info] compressing x-dimension to respect max_aspect={max_aspect}");
        dy * max_aspect
    } else {
        dx
    };

    let usable_w = width - 2.0 * padding;
    let base_scale = usable_w / effective_dx;
    let height = 2.0 * padding + dy * base_scale;

    eprintln!(
        "[info] writing SVG to {} (width={:.1}, height={:.1})",
        output_path.display(),
        width,
        height
    );

    // Full bounding box (for y, and for reference)
    let (min_x_full, max_x_full, min_y, max_y) = bounding_box(&node_coords)?;

    // Clipped X bounds for scaling
    let (min_x_clip, max_x_clip) = x_clip_bounds(&node_coords, x_clip_low, x_clip_high)?;
    let dx_clip = (max_x_clip - min_x_clip).max(1.0);
    let dy = (max_y - min_y).max(1.0);
    let aspect_clip = dx_clip / dy;

    eprintln!("[info] full bbox: x=[{min_x_full:.2},{max_x_full:.2}], y=[{min_y:.2},{max_y:.2}]");
    eprintln!(
        "[info] clipped x-range ({x_clip_low}–{x_clip_high} pct): [{min_x_clip:.2},{max_x_clip:.2}], dx_clip={dx_clip:.2}, dy={dy:.2}, aspect_clip={aspect_clip:.3}"
    );

    // Aspect control based on clipped range
    let effective_dx = if aspect_clip > max_aspect {
        eprintln!(
            "[info] compressing x-dimension to respect max_aspect={max_aspect} (on clipped range)"
        );
        dy * max_aspect
    } else {
        dx_clip
    };

    let usable_w = width - 2.0 * padding;
    let base_scale = usable_w / effective_dx;
    let height = 2.0 * padding + dy * base_scale;

    eprintln!(
        "[info] writing SVG to {} (width={:.1}, height={:.1})",
        output_path.display(),
        width,
        height
    );

    let file = File::create(&output_path)?;
    let mut out = BufWriter::new(file);

    render_svg(
        &mut out,
        width,
        height,
        padding,
        node_radius,
        max_aspect,
        &node_coords,
        if no_edges { &[] } else { &graph.edges },
        &highlight_paths,
        &highlight_names,
        (min_x_clip, max_x_clip, min_y, max_y),
        path_simplify_eps,
        edge_simplify_eps,
    )?;

    Ok(())
}

/// Read TSV from `odgi layout -T`.
/// Expects first 3 columns: node_id, x, y.
fn read_layout_tsv(path: &PathBuf) -> Result<HashMap<u64, Point>, Box<dyn Error>> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);

    let mut map = HashMap::new();

    for line_res in reader.lines() {
        let line = line_res?;
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let id = match fields[0].parse::<u64>() {
            Ok(v) => v,
            Err(_) => continue, // skip header / junk
        };
        let x = match fields[1].parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let y = match fields[2].parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };

        map.insert(id, Point { x, y });
    }

    if map.is_empty() {
        return Err("layout TSV appears to contain no parsable coordinates".into());
    }

    Ok(map)
}

/// Read paths and edges from a GFA 1.0 file.
/// Assumes segment names are numeric node ids (ODGI-style).
fn read_gfa(path: &PathBuf) -> Result<GraphData, Box<dyn Error>> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);

    let mut paths = Vec::new();
    let mut edges = Vec::new();

    for line_res in reader.lines() {
        let line = line_res?;
        if line.is_empty() {
            continue;
        }
        let tag = line.as_bytes()[0] as char;

        match tag {
            'P' => {
                // P <name> <segment_list> <opt_tags>
                let mut fields = line.split('\t');
                let _tag = fields.next(); // "P"
                let name = match fields.next() {
                    Some(n) => n.to_string(),
                    None => continue,
                };
                let seglist = match fields.next() {
                    Some(s) => s,
                    None => continue,
                };

                let mut steps = Vec::new();
                for seg in seglist.split(',') {
                    if seg.is_empty() {
                        continue;
                    }
                    let (id_str, rev) = match seg.chars().last() {
                        Some('+') => (&seg[..seg.len() - 1], false),
                        Some('-') => (&seg[..seg.len() - 1], true),
                        _ => (seg, false),
                    };
                    if let Ok(id) = id_str.parse::<u64>() {
                        steps.push(PathStep {
                            node_id: id,
                            reversed: rev,
                        });
                    }
                }

                if !steps.is_empty() {
                    paths.push(PathRec { name, steps });
                }
            }
            'L' => {
                // L <from> <from_orient> <to> <to_orient> <overlap>
                let mut fields = line.split('\t');
                let _tag = fields.next(); // "L"
                let from = match fields.next() {
                    Some(s) => s,
                    None => continue,
                };
                let _from_orient = fields.next();
                let to = match fields.next() {
                    Some(s) => s,
                    None => continue,
                };
                let _to_orient = fields.next();
                let _overlap = fields.next();

                if let (Ok(f_id), Ok(t_id)) = (from.parse::<u64>(), to.parse::<u64>()) {
                    edges.push(Edge {
                        from: f_id,
                        to: t_id,
                    });
                }
            }
            _ => {
                // ignore S, H, other lines
            }
        }
    }

    Ok(GraphData { paths, edges })
}

/// Compute clipped X bounds by percentile (0–100).
fn x_clip_bounds(
    coords: &HashMap<u64, Point>,
    pct_low: f64,
    pct_high: f64,
) -> Result<(f64, f64), Box<dyn Error>> {
    if coords.is_empty() {
        return Err("no coordinates to compute x-clip bounds".into());
    }

    // Clamp percentiles to [0, 100]
    let pl = pct_low.clamp(0.0, 100.0);
    let ph = pct_high.clamp(0.0, 100.0);
    if pl > ph {
        return Err("x-clip-low must be <= x-clip-high".into());
    }

    // Collect and sort all x-values
    let mut xs: Vec<f64> = coords.values().map(|p| p.x).collect();
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = xs.len();

    let idx_low = ((n - 1) as f64 * (pl / 100.0)).round() as usize;
    let idx_high = ((n - 1) as f64 * (ph / 100.0)).round() as usize;

    let x_min_clip = xs[idx_low.min(n - 1)];
    let x_max_clip = xs[idx_high.min(n - 1)];

    Ok((x_min_clip, x_max_clip))
}

/// Compute bounding box from node coordinates
fn bounding_box(coords: &HashMap<u64, Point>) -> Result<(f64, f64, f64, f64), Box<dyn Error>> {
    let mut iter = coords.values();
    let first = iter
        .next()
        .ok_or("no coordinates to compute bounding box")?;
    let mut min_x = first.x;
    let mut max_x = first.x;
    let mut min_y = first.y;
    let mut max_y = first.y;

    for p in iter {
        if p.x < min_x {
            min_x = p.x;
        }
        if p.x > max_x {
            max_x = p.x;
        }
        if p.y < min_y {
            min_y = p.y;
        }
        if p.y > max_y {
            max_y = p.y;
        }
    }

    Ok((min_x, max_x, min_y, max_y))
}

/// Simple epsilon-based polyline simplification in SVG space.
/// Keeps first and last point, and any intermediate point at least `eps` away
/// (Euclidean distance) from the last kept point.
fn downsample_polyline(points: &[(f64, f64)], eps: f64) -> Vec<(f64, f64)> {
    if eps <= 0.0 || points.len() <= 2 {
        return points.to_vec();
    }

    let eps2 = eps * eps;
    let mut kept = Vec::with_capacity(points.len());
    let mut last = points[0];
    kept.push(last);

    for &p in &points[1..] {
        let dx = p.0 - last.0;
        let dy = p.1 - last.1;
        if dx * dx + dy * dy >= eps2 {
            kept.push(p);
            last = p;
        }
    }

    // Always keep final point
    if let Some(&last_p) = points.last() {
        if kept.last().copied() != Some(last_p) {
            kept.push(last_p);
        }
    }

    kept
}

/// Simplify edges by only keeping those whose endpoints are at least `eps` apart
/// (in SVG pixel space). If `eps <= 0`, keeps all edges.
fn simplify_edges<F>(
    edges: &[Edge],
    coords: &HashMap<u64, Point>,
    transform: &F,
    eps: f64,
) -> Vec<(f64, f64, f64, f64)>
where
    F: Fn(&Point) -> (f64, f64),
{
    let mut segments = Vec::new();

    if eps <= 0.0 {
        // No simplification: just transform everything
        for e in edges {
            if let (Some(p1), Some(p2)) = (coords.get(&e.from), coords.get(&e.to)) {
                let (x1, y1) = transform(p1);
                let (x2, y2) = transform(p2);
                segments.push((x1, y1, x2, y2));
            }
        }
        return segments;
    }

    let eps2 = eps * eps;

    for e in edges {
        if let (Some(p1), Some(p2)) = (coords.get(&e.from), coords.get(&e.to)) {
            let (x1, y1) = transform(p1);
            let (x2, y2) = transform(p2);
            let dx = x2 - x1;
            let dy = y2 - y1;
            if dx * dx + dy * dy >= eps2 {
                segments.push((x1, y1, x2, y2));
            }
        }
    }

    segments
}

#[allow(clippy::too_many_arguments)]
fn render_svg<W: Write>(
    out: &mut W,
    width: f64,
    height: f64,
    padding: f64,
    _node_radius: f64,
    max_aspect: f64,
    coords: &HashMap<u64, Point>,
    edges: &[Edge],
    highlight_paths: &[PathRec],
    highlight_names: &[String],
    bbox: (f64, f64, f64, f64),
    path_simplify_eps: f64,
    edge_simplify_eps: f64,
) -> Result<(), Box<dyn Error>> {
    let (min_x_clip, max_x_clip, min_y, max_y) = bbox;
    let dx_clip = (max_x_clip - min_x_clip).max(1.0);
    let dy = (max_y - min_y).max(1.0);
    let aspect_clip = dx_clip / dy;

    let usable_w = width - 2.0 * padding;
    let _usable_h = height - 2.0 * padding;

    // If aspect is too large, compress x so effective span is dy * max_aspect
    let effective_dx = if aspect_clip > max_aspect {
        dy * max_aspect
    } else {
        dx_clip
    };

    let base_scale = usable_w / effective_dx;
    let compress_factor = effective_dx / dx_clip; // <= 1.0 if compression
    let sx = base_scale * compress_factor;
    let sy = base_scale;

    // Transform from layout space → SVG space, with X clamped
    let transform = |p: &Point| -> (f64, f64) {
        let x_clamped = p.x.clamp(min_x_clip, max_x_clip);
        let x = padding + (x_clamped - min_x_clip) * sx;
        let y = padding + (max_y - p.y) * sy; // keep full Y range
        (x, y)
    };

    // SVG header
    writeln!(
        out,
        r##"<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     width="{width}"
     height="{height}"
     viewBox="0 0 {width} {height}">
  <rect x="0" y="0" width="{width}" height="{height}" fill="white" />
"##
    )?;

    // Draw edges (graph skeleton) as light grey, possibly simplified
    if !edges.is_empty() {
        writeln!(
            out,
            r##"  <g id="edges" stroke="black" stroke-width="0.4" stroke-linecap="round" opacity="0.7">"##
        )?;

        let eps2 = if edge_simplify_eps > 0.0 {
            edge_simplify_eps * edge_simplify_eps
        } else {
            0.0
        };

        let mut kept = 0usize;

        for (x1, y1, x2, y2) in simplify_edges(edges, coords, &transform, edge_simplify_eps) {
            if eps2 > 0.0 {
                let dx = x2 - x1;
                let dy = y2 - y1;
                if dx * dx + dy * dy < eps2 {
                    continue; // too short, skip
                }
            }

            kept += 1;
            writeln!(
                out,
                r#"    <line x1="{x1:.3}" y1="{y1:.3}" x2="{x2:.3}" y2="{y2:.3}" />"#
            )?;
        }

        eprintln!("[info] drew {kept} skeleton edges (from {})", edges.len());
        writeln!(out, "  </g>")?;
    }

    // Background nodes as faint dots (currently disabled)
    // writeln!(out, r##"  <g id="nodes" fill="#bbbbbb" opacity="0.7">"##)?;
    // for p in coords.values() {
    //     let (x, y) = transform(p);
    //     writeln!(
    //         out,
    //         r#"    <circle cx="{x:.3}" cy="{y:.3}" r="{node_radius:.3}" />"#
    //     )?;
    // }
    // writeln!(out, "  </g>")?;

    // Colour palette for highlighted paths
    let palette = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
        "#bcbd22", "#17becf",
    ];

    // Draw highlighted paths (if any)
    if !highlight_paths.is_empty() {
        writeln!(out, r#"  <g id="highlight_paths" fill="none">"#)?;

        for (idx, path) in highlight_paths.iter().enumerate() {
            // Build transformed polyline for this path
            let mut points: Vec<(f64, f64)> = Vec::with_capacity(path.steps.len());
            for step in &path.steps {
                if let Some(p) = coords.get(&step.node_id) {
                    points.push(transform(p));
                }
            }
            if points.len() < 2 {
                continue;
            }

            // Single downsample in SVG space
            let mut points = downsample_polyline(&points, path_simplify_eps);

            // Optional hard cap: don't let a single path go totally wild
            // (this is in SVG points *after* simplification)
            const MAX_POINTS_PER_PATH: usize = 200_000;
            if points.len() > MAX_POINTS_PER_PATH {
                let stride = (points.len() as f64 / MAX_POINTS_PER_PATH as f64).ceil() as usize;
                eprintln!(
                    "[info] path {} had {} points after eps-simplify; thinning by stride {stride}",
                    path.name,
                    points.len()
                );
                points = points.into_iter().step_by(stride.max(1)).collect();
            }

            let color = palette[idx % palette.len()];
            let stroke_width = 2.0; // thicker than skeleton

            // Chunk into smaller SVG paths *without* cloning each chunk
            let chunk_len: usize = 400;
            let mut start = 0usize;
            while start + 1 < points.len() {
                let end = (start + chunk_len).min(points.len());
                let slice = &points[start..end];
                if slice.len() < 2 {
                    break;
                }

                // Build path string for this chunk
                let (x0, y0) = slice[0];
                let mut d = format!("M {x0:.3} {y0:.3}");
                for &(x, y) in &slice[1..] {
                    use std::fmt::Write as FmtWrite;
                    write!(&mut d, " L {x:.3} {y:.3}").unwrap();
                }

                writeln!(
                    out,
                    r#"    <path d="{d}" stroke="{color}" stroke-width="{stroke_width}" opacity="0.95" />"#
                )?;

                if end == points.len() {
                    break;
                }
                // Overlap by 1 point to keep continuity
                start = end - 1;
            }
        }

        writeln!(out, "  </g>")?;

        // Legend
        if !highlight_names.is_empty() {
            writeln!(
                out,
                r#"  <g id="legend" font-family="monospace" font-size="10">"#
            )?;
            let mut y = padding + 10.0;
            let x = padding + 10.0;

            for (idx, name) in highlight_names.iter().enumerate() {
                let color = palette[idx % palette.len()];
                writeln!(
                    out,
                    r#"    <line x1="{x}" y1="{y}" x2="{x2}" y2="{y}" stroke="{color}" stroke-width="0.8" />"#,
                    x = x,
                    x2 = x + 20.0,
                    y = y,
                    color = color
                )?;
                writeln!(
                    out,
                    r#"    <text x="{tx}" y="{ty}" fill="black">{name}</text>"#,
                    tx = x + 25.0,
                    ty = y + 3.0,
                    name = name
                )?;
                y += 14.0;
            }
            writeln!(out, "  </g>")?;
        }
    }

    writeln!(out, "</svg>")?;
    Ok(())
}
