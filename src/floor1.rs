//! Floor type 1 (specification section 7): a piecewise-linear spectral
//! envelope on a dB scale, coded as predicted Y values at a list of X
//! positions.

use crate::bits::{BitReader, BitWriter, EndOfPacket, ilog};
use crate::codebook::Codebook;
use crate::error::{Result, invalid};
use crate::tables::FLOOR1_INVERSE_DB;

/// `{ 256, 128, 86, 64 }` indexed by `floor1_multiplier - 1`.
pub const RANGES: [i32; 4] = [256, 128, 86, 64];

/// One partition class of a floor 1 configuration.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Floor1Class {
    /// `floor1_class_dimensions`: Y values per partition, 1 to 8.
    pub dimensions: u8,
    /// `floor1_class_subclasses`: log2 of the number of subclass books, 0 to 3.
    pub subclasses: u8,
    /// `floor1_class_masterbooks`, read only when `subclasses > 0`.
    pub masterbook: u8,
    /// `floor1_subclass_books`: `1 << subclasses` book numbers, -1 for none.
    pub subclass_books: Vec<i16>,
}

/// A floor 1 configuration (7.2.2).
#[derive(Clone, Debug)]
pub struct Floor1 {
    /// `floor1_partition_class_list`.
    pub partition_class: Vec<u8>,
    /// The classes the partitions name.
    pub classes: Vec<Floor1Class>,
    /// `floor1_multiplier`, 1 to 4.
    pub multiplier: u8,
    /// `rangebits`: the last X position is `1 << rangebits`.
    pub rangebits: u8,
    /// `floor1_X_list`, in list (coding) order, including the implicit 0
    /// and `1 << rangebits` at its head.
    pub x_list: Vec<u32>,
    /// Indices of `x_list` in ascending order of X.
    sorted: Vec<usize>,
    /// `low_neighbor` and `high_neighbor` of each index from 2.
    low: Vec<usize>,
    high: Vec<usize>,
}

impl PartialEq for Floor1 {
    fn eq(&self, other: &Self) -> bool {
        self.partition_class == other.partition_class
            && self.classes == other.classes
            && self.multiplier == other.multiplier
            && self.rangebits == other.rangebits
            && self.x_list == other.x_list
    }
}

/// `low_neighbor` (9.2.4): the position of the greatest earlier value below `v[x]`.
pub fn low_neighbor(v: &[u32], x: usize) -> usize {
    let mut best: Option<usize> = None;
    for n in 0..x {
        if v[n] < v[x] && best.is_none_or(|b| v[n] > v[b]) {
            best = Some(n);
        }
    }
    best.unwrap_or(0)
}

/// `high_neighbor` (9.2.5): the position of the lowest earlier value above `v[x]`.
pub fn high_neighbor(v: &[u32], x: usize) -> usize {
    let mut best: Option<usize> = None;
    for n in 0..x {
        if v[n] > v[x] && best.is_none_or(|b| v[n] < v[b]) {
            best = Some(n);
        }
    }
    best.unwrap_or(1)
}

/// `render_point` (9.2.6).
pub fn render_point(x0: i32, y0: i32, x1: i32, y1: i32, x: i32) -> i32 {
    let dy = y1 - y0;
    let adx = x1 - x0;
    let err = dy.abs() * (x - x0);
    let off = err / adx;
    if dy < 0 { y0 - off } else { y0 + off }
}

/// `render_line` (9.2.7): writes `v[x]` for `x` in `[x0, x1)`, dropping
/// positions past the end of `v`.
pub fn render_line(x0: i32, y0: i32, x1: i32, y1: i32, v: &mut [i32]) {
    let dy = y1 - y0;
    let adx = x1 - x0;
    let mut ady = dy.abs();
    let base = dy / adx;
    let sy = if dy < 0 { base - 1 } else { base + 1 };
    let mut y = y0;
    let mut err = 0;
    ady -= base.abs() * adx;
    let len = v.len() as i32;
    if x0 >= 0 && x0 < len {
        v[x0 as usize] = y;
    }
    for x in x0 + 1..x1.min(len) {
        err += ady;
        if err >= adx {
            err -= adx;
            y += sy;
        } else {
            y += base;
        }
        v[x as usize] = y;
    }
}

impl Floor1 {
    /// Build a configuration, checking it as header decode does.
    pub fn new(partition_class: Vec<u8>, classes: Vec<Floor1Class>, multiplier: u8, rangebits: u8, x_list: Vec<u32>, codebooks: usize) -> Result<Self> {
        if partition_class.len() > 31 {
            return Err(invalid("floor 1 with more than 31 partitions"));
        }
        if !(1..=4).contains(&multiplier) || rangebits > 15 {
            return Err(invalid("floor 1 multiplier or rangebits out of range"));
        }
        for &c in &partition_class {
            if c as usize >= classes.len() {
                return Err(invalid("floor 1 partition names a missing class"));
            }
        }
        for c in &classes {
            if !(1..=8).contains(&c.dimensions) || c.subclasses > 3 || c.subclass_books.len() != 1 << c.subclasses {
                return Err(invalid("floor 1 class out of range"));
            }
            if c.subclasses > 0 && c.masterbook as usize >= codebooks {
                return Err(invalid("floor 1 master book does not exist"));
            }
            if c.subclass_books.iter().any(|&b| b >= codebooks as i16 || b < -1) {
                return Err(invalid("floor 1 subclass book does not exist"));
            }
        }
        let values: usize = 2 + partition_class.iter().map(|&c| classes[c as usize].dimensions as usize).sum::<usize>();
        if x_list.len() != values || values > 65 {
            return Err(invalid("floor 1 X list longer than 65 values or the wrong length"));
        }
        if x_list[0] != 0 || x_list[1] != 1 << rangebits {
            return Err(invalid("floor 1 X list does not start with 0 and 2^rangebits"));
        }
        let mut sorted: Vec<usize> = (0..values).collect();
        sorted.sort_by_key(|&i| x_list[i]);
        if sorted.windows(2).any(|w| x_list[w[0]] == x_list[w[1]]) {
            return Err(invalid("floor 1 X list values are not unique"));
        }
        let low = (0..values).map(|i| if i < 2 { 0 } else { low_neighbor(&x_list, i) }).collect();
        let high = (0..values).map(|i| if i < 2 { 1 } else { high_neighbor(&x_list, i) }).collect();
        Ok(Floor1 { partition_class, classes, multiplier, rangebits, x_list, sorted, low, high })
    }

    /// Header decode (7.2.2).
    pub(crate) fn read(r: &mut BitReader, codebooks: &[Codebook]) -> Result<Self> {
        Self::read_fields(r, codebooks.len()).map_err(|_| invalid("setup header ends inside a floor 1 configuration"))?
    }

    fn read_fields(r: &mut BitReader, codebooks: usize) -> std::result::Result<Result<Self>, EndOfPacket> {
        let partitions = r.read(5)? as usize;
        let mut partition_class = Vec::with_capacity(partitions);
        for _ in 0..partitions {
            partition_class.push(r.read(4)? as u8);
        }
        let max_class = partition_class.iter().map(|&c| c as i32).max().unwrap_or(-1);
        let mut classes = Vec::new();
        for _ in 0..=max_class {
            let dimensions = r.read(3)? as u8 + 1;
            let subclasses = r.read(2)? as u8;
            let masterbook = if subclasses > 0 { r.read(8)? as u8 } else { 0 };
            let mut subclass_books = Vec::new();
            for _ in 0..1 << subclasses {
                subclass_books.push(r.read(8)? as i16 - 1);
            }
            classes.push(Floor1Class { dimensions, subclasses, masterbook, subclass_books });
        }
        let multiplier = r.read(2)? as u8 + 1;
        let rangebits = r.read(4)? as u8;
        let mut x_list = vec![0, 1u32 << rangebits];
        for &c in &partition_class {
            for _ in 0..classes[c as usize].dimensions {
                x_list.push(r.read(rangebits as u32)?);
                if x_list.len() > 65 {
                    return Ok(Err(invalid("floor 1 X list longer than 65 values")));
                }
            }
        }
        Ok(Floor1::new(partition_class, classes, multiplier, rangebits, x_list, codebooks))
    }

    pub(crate) fn write(&self, w: &mut BitWriter) {
        w.write(self.partition_class.len() as u32, 5);
        for &c in &self.partition_class {
            w.write(c as u32, 4);
        }
        let max_class = self.partition_class.iter().map(|&c| c as i32).max().unwrap_or(-1);
        for c in &self.classes[..(max_class + 1) as usize] {
            w.write(c.dimensions as u32 - 1, 3);
            w.write(c.subclasses as u32, 2);
            if c.subclasses > 0 {
                w.write(c.masterbook as u32, 8);
            }
            for &b in &c.subclass_books {
                w.write((b + 1) as u32, 8);
            }
        }
        w.write(self.multiplier as u32 - 1, 2);
        w.write(self.rangebits as u32, 4);
        for &x in &self.x_list[2..] {
            w.write(x, self.rangebits as u32);
        }
    }

    /// `range` for this configuration's multiplier.
    pub fn range(&self) -> i32 {
        RANGES[self.multiplier as usize - 1]
    }

    /// Packet decode (7.2.3): `None` when the floor is unused this frame.
    pub(crate) fn decode(&self, r: &mut BitReader, books: &[Codebook]) -> std::result::Result<Option<Vec<i32>>, EndOfPacket> {
        if !r.read_flag()? {
            return Ok(None);
        }
        let bits = ilog(self.range() as i64 - 1);
        let mut y = Vec::with_capacity(self.x_list.len());
        y.push(r.read(bits)? as i32);
        y.push(r.read(bits)? as i32);
        for &class in &self.partition_class {
            let c = &self.classes[class as usize];
            let cbits = c.subclasses as u32;
            let csub = (1u32 << cbits) - 1;
            let mut cval = if cbits > 0 { books[c.masterbook as usize].decode_scalar(r)? } else { 0 };
            for _ in 0..c.dimensions {
                let book = c.subclass_books[(cval & csub) as usize];
                cval >>= cbits;
                y.push(if book >= 0 { books[book as usize].decode_scalar(r)? as i32 } else { 0 });
            }
        }
        Ok(Some(y))
    }

    /// Step 1 of curve computation (7.2.4): the final Y values and which
    /// posts step 2 draws.
    pub fn amplitudes(&self, y: &[i32]) -> (Vec<i32>, Vec<bool>) {
        let range = self.range();
        let n = self.x_list.len();
        let mut fin = vec![0i32; n];
        let mut step2 = vec![false; n];
        fin[0] = y[0];
        fin[1] = y[1];
        step2[0] = true;
        step2[1] = true;
        for i in 2..n {
            let lo = self.low[i];
            let hi = self.high[i];
            let predicted = render_point(self.x_list[lo] as i32, fin[lo], self.x_list[hi] as i32, fin[hi], self.x_list[i] as i32);
            let val = y[i];
            let highroom = range - predicted;
            let lowroom = predicted;
            let room = if highroom < lowroom { highroom * 2 } else { lowroom * 2 };
            if val != 0 {
                step2[lo] = true;
                step2[hi] = true;
                step2[i] = true;
                fin[i] = if val >= room {
                    if highroom > lowroom { val - lowroom + predicted } else { predicted - val + highroom - 1 }
                } else if val & 1 == 1 {
                    predicted - (val + 1) / 2
                } else {
                    predicted + val / 2
                };
            } else {
                fin[i] = predicted;
            }
        }
        // The guard 7.2.4 recommends against setups that abuse the coding.
        for v in fin.iter_mut() {
            *v = (*v).clamp(0, range - 1);
        }
        (fin, step2)
    }

    /// Step 2 of curve computation (7.2.4): the integer curve over `n`
    /// positions, before the dB table lookup.
    pub fn curve(&self, fin: &[i32], step2: &[bool], n: usize) -> Vec<i32> {
        let mult = self.multiplier as i32;
        let mut v = vec![0i32; n];
        let mut hx = 0i32;
        let mut hy = 0i32;
        let mut lx = 0i32;
        let mut ly = fin[self.sorted[0]] * mult;
        for &i in &self.sorted[1..] {
            if step2[i] {
                hy = fin[i] * mult;
                hx = self.x_list[i] as i32;
                render_line(lx, ly, hx, hy, &mut v);
                lx = hx;
                ly = hy;
            }
        }
        if (hx as usize) < n {
            render_line(hx, hy, n as i32, hy, &mut v);
        }
        v
    }

    /// The linear floor curve over `n` positions for decoded `y` values.
    pub fn synthesize(&self, y: &[i32], n: usize) -> Vec<f32> {
        let (fin, step2) = self.amplitudes(y);
        self.curve(&fin, &step2, n).into_iter().map(|v| FLOOR1_INVERSE_DB[v.clamp(0, 255) as usize]).collect()
    }

    /// The value to code at post `i` (list order, from 2) so that its final
    /// Y is `target`, given the final Y values of the posts before it: the
    /// inverse of step 1. `fin` must hold those earlier values.
    pub fn encode_value(&self, fin: &[i32], i: usize, target: i32) -> i32 {
        let range = self.range();
        let lo = self.low[i];
        let hi = self.high[i];
        let predicted = render_point(self.x_list[lo] as i32, fin[lo], self.x_list[hi] as i32, fin[hi], self.x_list[i] as i32);
        let highroom = range - predicted;
        let lowroom = predicted;
        let room = if highroom < lowroom { highroom * 2 } else { lowroom * 2 };
        let d = target - predicted;
        if d == 0 {
            return 0;
        }
        let interleaved = if d > 0 { 2 * d } else { -2 * d - 1 };
        if interleaved < room {
            interleaved
        } else if highroom > lowroom {
            d + lowroom
        } else {
            highroom - 1 - d
        }
    }

    /// The post order of `x_list` sorted by X, for the encoder.
    pub fn sorted_posts(&self) -> &[usize] {
        &self.sorted
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn floor(x_list: Vec<u32>, mult: u8, rangebits: u8) -> Floor1 {
        let extra = x_list.len() - 2;
        let class = Floor1Class { dimensions: 1, subclasses: 0, masterbook: 0, subclass_books: vec![0] };
        Floor1::new(vec![0; extra], vec![class], mult, rangebits, x_list, 1).unwrap()
    }

    #[test]
    fn neighbors_follow_the_definitions() {
        // The interleaved list of the 7.2.1 example.
        let x = [0u32, 128, 64, 32, 96, 16, 48, 80, 112];
        assert_eq!((low_neighbor(&x, 2), high_neighbor(&x, 2)), (0, 1));
        assert_eq!((low_neighbor(&x, 3), high_neighbor(&x, 3)), (0, 2));
        assert_eq!((low_neighbor(&x, 4), high_neighbor(&x, 4)), (2, 1));
        assert_eq!((low_neighbor(&x, 5), high_neighbor(&x, 5)), (0, 3));
        assert_eq!((low_neighbor(&x, 6), high_neighbor(&x, 6)), (3, 2));
        assert_eq!((low_neighbor(&x, 7), high_neighbor(&x, 7)), (2, 4));
        assert_eq!((low_neighbor(&x, 8), high_neighbor(&x, 8)), (4, 1));
    }

    #[test]
    fn render_point_truncates_toward_the_first_point() {
        assert_eq!(render_point(0, 10, 10, 20, 5), 15);
        assert_eq!(render_point(0, 10, 3, 20, 1), 13); // 10/3 = 3
        assert_eq!(render_point(0, 20, 3, 10, 1), 17); // 20 - 3
        assert_eq!(render_point(0, 7, 8, 7, 3), 7);
    }

    /// `render_line` draws within one unit of the exact line, starts on the
    /// first point and is monotonic between its ends.
    #[test]
    fn render_line_tracks_the_exact_line() {
        for (x0, y0, x1, y1) in [(0, 0, 10, 3), (0, 100, 7, 3), (5, 40, 300, 41), (3, 200, 4, 0), (0, 0, 128, 255), (10, 50, 20, 50)] {
            let mut v = vec![-1; x1 as usize + 1];
            render_line(x0, y0, x1, y1, &mut v);
            assert_eq!(v[x0 as usize], y0);
            for x in x0..x1 {
                let exact = y0 as f64 + (y1 - y0) as f64 * (x - x0) as f64 / (x1 - x0) as f64;
                assert!((v[x as usize] as f64 - exact).abs() < 1.0, "({x0},{y0})-({x1},{y1}) at {x}: {} vs {exact}", v[x as usize]);
                if x > x0 {
                    let step = v[x as usize] - v[x as usize - 1];
                    assert!(step * (y1 - y0).signum() >= 0);
                }
            }
            assert_eq!(v[x1 as usize], -1, "the end point is not drawn");
        }
    }

    /// Every target is reachable from every prediction: `encode_value`
    /// inverts step 1 exactly, for each multiplier.
    #[test]
    fn value_coding_inverts_the_prediction_for_every_target() {
        for mult in 1..=4u8 {
            let fl = floor(vec![0, 128, 64], mult, 7);
            let range = fl.range();
            for y0 in [0, 1, range / 3, range / 2, range - 1] {
                for y1 in [0, 2, range / 2, range - 1] {
                    for target in 0..range {
                        let fin = [y0, y1, 0];
                        let val = fl.encode_value(&fin, 2, target);
                        assert!(val >= 0 && val < range, "val {val}");
                        let (got, step2) = fl.amplitudes(&[y0, y1, val]);
                        assert_eq!(got[2], target, "mult {mult} y0 {y0} y1 {y1} target {target} val {val}");
                        assert_eq!(step2[2], val != 0);
                    }
                }
            }
        }
    }

    /// A hand-worked floor: posts 0, 128, 64 (n = 128) with Y 10, 20 and
    /// a coded difference of +4 at 64 (val 8), multiplier 2.
    #[test]
    fn hand_worked_curve() {
        let fl = floor(vec![0, 128, 64], 2, 7);
        let (fin, step2) = fl.amplitudes(&[10, 20, 8]);
        // predicted = 10 + 10*64/128 = 15; val 8 is even: 15 + 4.
        assert_eq!(fin, vec![10, 20, 19]);
        assert_eq!(step2, vec![true, true, true]);
        let c = fl.curve(&fin, &step2, 128);
        // Lines (0,20)-(64,38) and (64,38)-(128,40), in units of Y * 2.
        assert_eq!(c[0], 20);
        assert_eq!(c[64], 38);
        assert_eq!(c[32], 29);
        assert_eq!(c[127], 39);
        let lin = fl.synthesize(&[10, 20, 8], 128);
        assert_eq!(lin[64], FLOOR1_INVERSE_DB[38]);
        // An unchanged post (val 0) is not drawn: the line runs straight
        // from 0 to 128.
        let (fin, step2) = fl.amplitudes(&[10, 20, 0]);
        assert!(!step2[2]);
        let c = fl.curve(&fin, &step2, 128);
        assert_eq!(c[64], 30);
    }

    /// A curve shorter than the X range is truncated; a longer one is
    /// extended flat from the last post.
    #[test]
    fn curve_is_truncated_or_extended_to_n() {
        let fl = floor(vec![0, 64, 32], 1, 6);
        let (fin, step2) = fl.amplitudes(&[100, 100, 0]);
        assert_eq!(fl.curve(&fin, &step2, 32).len(), 32);
        let c = fl.curve(&fin, &step2, 100);
        assert!(c[64..].iter().all(|&v| v == 100));
    }

    #[test]
    fn setups_are_checked() {
        let class = Floor1Class { dimensions: 1, subclasses: 0, masterbook: 0, subclass_books: vec![0] };
        // Duplicate X.
        assert!(Floor1::new(vec![0], vec![class.clone()], 1, 7, vec![0, 128, 128], 1).is_err());
        // Book out of range.
        assert!(Floor1::new(vec![0], vec![class.clone()], 1, 7, vec![0, 128, 5], 0).is_err());
        assert!(Floor1::new(vec![0], vec![class], 1, 7, vec![0, 128, 5], 1).is_ok());
    }
}
