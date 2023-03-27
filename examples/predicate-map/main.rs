use robust::Coord;

// Directly evaluate the orient2d determinant.
// Refer: https://www.cs.cmu.edu/~quake/robust.html
fn simple_orient2d(p: Coord<f64>, q: Coord<f64>, r: Coord<f64>) -> f64 {
    (q.x - p.x) * (r.y - q.y) - (q.y - p.y) * (r.x - q.x)
}

// Directly evaluate the incircle determinant.
// Refer: https://www.cs.cmu.edu/~quake/robust.html
fn simple_incircle(a: Coord<f64>, b: Coord<f64>, c: Coord<f64>, d: Coord<f64>) -> f64 {
    let m11 = a.x - d.x;
    let m12 = a.y - d.y;
    let m13 = m11.powi(2) + m12.powi(2);

    let m21 = b.x - d.x;
    let m22 = b.y - d.y;
    let m23 = m21.powi(2) + m22.powi(2);

    let m31 = c.x - d.x;
    let m32 = c.y - d.y;
    let m33 = m31.powi(2) + m32.powi(2);

    m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31)
}

use std::cmp::Ordering;
fn orient2d_tests<F>(predicate: F, start: Coord<f64>, width: usize, height: usize) -> Vec<Ordering>
where
    F: Fn(Coord<f64>) -> f64,
{
    use float_extras::f64::nextafter;
    let mut yd = start.y;
    let mut data = Vec::with_capacity(width * height);

    for _ in 0..height {
        let mut xd = start.x;
        for _ in 0..width {
            let p = Coord { x: xd, y: yd };
            data.push(predicate(p).partial_cmp(&0.).unwrap());
            xd = nextafter(xd, std::f64::INFINITY);
        }
        yd = nextafter(yd, std::f64::INFINITY);
    }

    data
}

use std::path::Path;
fn write_png(data: &[Ordering], path: &Path, width: usize, height: usize) {
    assert_eq!(data.len(), width * height);

    use std::fs::File;
    use std::io::BufWriter;

    let file = File::create(path).unwrap();
    let ref mut w = BufWriter::new(file);

    let mut encoder = png::Encoder::new(w, width as u32, height as u32);
    encoder.set_color(png::ColorType::Grayscale);
    encoder.set_depth(png::BitDepth::Eight);

    let mut writer = encoder.write_header().unwrap();
    let data = data
        .iter()
        .map(|w| match w {
            Ordering::Less => 0u8,
            Ordering::Equal => 127,
            Ordering::Greater => 255,
        })
        .collect::<Vec<_>>();
    writer.write_image_data(&data).unwrap();
}

fn usage(name: &str) -> ! {
    eprintln!(
        "Usage: {} {{naive | robust}} {{incircle | orient2d}} <output.png>",
        name
    );
    std::process::exit(1);
}

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    if args.len() != 4 {
        usage(&args[0])
    }

    let p1 = Coord { x: 12., y: 12. };
    let p2 = Coord { x: 24., y: 24. };
    let p3 = Coord { x: -12., y: -12. };
    let predicate: Box<dyn Fn(Coord<f64>) -> f64> = match (args[1].as_str(), args[2].as_str()) {
        ("naive", "incircle") => Box::new(|p| simple_incircle(p1, p3, p2, p)),
        ("naive", "orient2d") => Box::new(|p| simple_orient2d(p1, p, p2)),
        ("robust", "incircle") => Box::new(|p| robust::incircle(p1, p3, p2, p)),
        ("robust", "orient2d") => Box::new(|p| robust::orient2d(p1, p, p2)),
        _ => usage(&args[0]),
    };

    let data = orient2d_tests(predicate, Coord { x: 0.5, y: 0.5 }, 256, 256);
    write_png(&data, Path::new(&args[3]), 256, 256);
}
