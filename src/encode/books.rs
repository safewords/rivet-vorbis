//! The encoder's codebooks and setup: Huffman codeword lengths designed
//! from modelled symbol probabilities, the floor 1 and residue 2
//! configurations, and the mappings and modes that tie them together.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use crate::codebook::{Codebook, VqLookup, float32_pack};
use crate::floor1::{Floor1, Floor1Class};
use crate::header::{Floor, Mapping, Mode, Setup};
use crate::residue::Residue;

/// Codebook numbers in the setup header.
pub(crate) const BOOK_FLOOR_MASTER: usize = 0;
pub(crate) const BOOK_FLOOR_SUB: [usize; 3] = [1, 2, 3];
pub(crate) const BOOK_CLASS: usize = 4;
/// Residue VQ books A..G (see [`RESIDUE_BOOKS`]).
pub(crate) const BOOK_VQ0: usize = 5;

/// Floor values coded by each subclass book (0 needs no book).
pub(crate) const FLOOR_SUB_LIMITS: [i32; 3] = [4, 16, 128];

/// A lattice VQ book: dimension, values per scalar, the smallest value
/// and the step between values.
#[derive(Clone, Copy)]
pub(crate) struct Lattice {
    pub dim: usize,
    pub values: usize,
    pub min: i32,
    pub step: i32,
    /// Laplacian scale of the modelled scalar distribution, in steps
    /// (chosen by measuring this encoder's bit rate with each).
    pub scale: f64,
}

/// The residue books, A to G.
pub(crate) const RESIDUE_BOOKS: [Lattice; 7] = [
    Lattice { dim: 4, values: 3, min: -1, step: 1, scale: 0.75 },
    Lattice { dim: 4, values: 5, min: -2, step: 1, scale: 0.8 },
    Lattice { dim: 2, values: 9, min: -4, step: 1, scale: 1.4 },
    Lattice { dim: 2, values: 17, min: -8, step: 1, scale: 2.6 },
    Lattice { dim: 2, values: 33, min: -16, step: 1, scale: 5.0 },
    Lattice { dim: 2, values: 15, min: -231, step: 33, scale: 0.3 },
    Lattice { dim: 1, values: 63, min: -15345, step: 495, scale: 0.5 },
];

/// Residue partition classes: the largest magnitude each covers and its
/// books per pass (indices into [`RESIDUE_BOOKS`]).
pub(crate) const CLASSES: [(i32, &[usize]); 8] = [
    (0, &[]),
    (1, &[0]),
    (2, &[1]),
    (4, &[2]),
    (8, &[3]),
    (16, &[4]),
    (247, &[5, 4]),
    (15592, &[6, 5, 4]),
];

/// The largest residue magnitude the books can code.
pub(crate) const MAX_RESIDUE: i32 = 15592;

/// Residue partition size (of the channel-interleaved vector).
pub(crate) const PARTITION: usize = 32;

/// Codeword lengths for symbols of the given weights: a Huffman code,
/// flattened until no codeword exceeds `max_len` bits. Every symbol gets a
/// codeword, so the tree is complete.
pub(crate) fn huffman_lengths(weights: &[f64], max_len: u32) -> Vec<u8> {
    let n = weights.len();
    if n == 1 {
        return vec![1];
    }
    let total: f64 = weights.iter().sum::<f64>().max(1e-300);
    let mut flatten = 0.0;
    loop {
        let w: Vec<f64> = weights.iter().map(|&x| x / total + flatten / n as f64).collect();
        let lengths = huffman(&w);
        if lengths.iter().all(|&l| l as u32 <= max_len) {
            return lengths;
        }
        flatten = if flatten == 0.0 { 1e-9 } else { flatten * 4.0 };
    }
}

fn huffman(w: &[f64]) -> Vec<u8> {
    // Weights as fixed point for a total order in the heap.
    let n = w.len();
    let mut parent = vec![usize::MAX; 2 * n - 1];
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = w.iter().enumerate().map(|(i, &x)| Reverse(((x * 1e15) as u64 + 1, i))).collect();
    let mut next = n;
    while heap.len() > 1 {
        let Reverse((a, i)) = heap.pop().expect("two nodes");
        let Reverse((b, j)) = heap.pop().expect("two nodes");
        parent[i] = next;
        parent[j] = next;
        heap.push(Reverse((a + b, next)));
        next += 1;
    }
    (0..n)
        .map(|i| {
            let mut d = 0u8;
            let mut x = i;
            while parent[x] != usize::MAX {
                x = parent[x];
                d += 1;
            }
            d
        })
        .collect()
}

fn laplace(v: f64, scale: f64) -> f64 {
    (-v.abs() / scale).exp()
}

fn scalar_book(weights: &[f64]) -> Codebook {
    Codebook::new(1, huffman_lengths(weights, 24), None).expect("designed scalar book")
}

fn lattice_book(l: &Lattice) -> Codebook {
    let entries = l.values.pow(l.dim as u32);
    let weights: Vec<f64> = (0..entries)
        .map(|e| {
            let mut p = 1.0;
            let mut rest = e;
            for _ in 0..l.dim {
                let digit = rest % l.values;
                rest /= l.values;
                p *= laplace((l.min + digit as i32 * l.step) as f64 / l.step as f64, l.scale);
            }
            p
        })
        .collect();
    let bits = usize::BITS - (l.values - 1).leading_zeros();
    let lookup = VqLookup {
        lookup_type: 1,
        minimum: float32_pack(l.min as f64),
        delta: float32_pack(l.step as f64),
        value_bits: bits as u8,
        sequence_p: false,
        multiplicands: (0..l.values as u32).collect(),
    };
    Codebook::new(l.dim as u16, huffman_lengths(&weights, 24), Some(lookup)).expect("designed lattice book")
}

/// Every codebook the encoder uses, in setup order.
pub(crate) fn codebooks() -> Vec<Codebook> {
    let mut books = Vec::new();
    // Floor master: three subclass digits (2 bits each, element 0 lowest).
    let digit = [0.40, 0.30, 0.22, 0.08];
    let master: Vec<f64> = (0..64).map(|c| digit[c & 3] * digit[c >> 2 & 3] * digit[c >> 4 & 3]).collect();
    books.push(scalar_book(&master));
    books.push(scalar_book(&[0.05, 0.40, 0.33, 0.22]));
    books.push(scalar_book(&(0..16).map(|v| laplace(v as f64, 5.0)).collect::<Vec<_>>()));
    books.push(scalar_book(&(0..128).map(|v| laplace(v as f64, 25.0)).collect::<Vec<_>>()));
    // Classbook: two classifications per word.
    let class = [0.30, 0.20, 0.20, 0.12, 0.08, 0.06, 0.03, 0.01];
    let words: Vec<f64> = (0..64).map(|w| class[w / 8] * class[w % 8]).collect();
    books.push(scalar_book(&words));
    for l in &RESIDUE_BOOKS {
        books.push(lattice_book(l));
    }
    books
}

/// Floor 1 posts (besides 0 and the end), sorted, for a vector of `half`
/// values: denser at low frequencies.
fn posts(half: usize) -> Vec<u32> {
    let long = [1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256, 320, 384, 448, 512, 640, 768, 896];
    let short = [1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64, 80, 96, 112];
    let (list, base): (&[u32], u32) = if half >= 1024 { (&long, 1024) } else { (&short, 128) };
    // Scale to this vector length, keeping posts distinct and inside it.
    let mut out: Vec<u32> = Vec::new();
    for &p in list {
        let x = (p as u64 * half as u64 / base as u64) as u32;
        if x > 0 && x < half as u32 && out.last() != Some(&x) {
            out.push(x);
        }
    }
    // A whole number of 3-post partitions.
    while out.len() % 3 != 0 {
        out.pop();
    }
    out
}

/// The posts in coding order: breadth-first bisection of the sorted list,
/// each interval split at the post nearest the geometric mean of its ends,
/// so each post is predicted from coarser ones on either side and the
/// lines between them follow a spectrum's roughly logarithmic shape.
fn bisection_order(sorted: &[u32], end: u32) -> Vec<u32> {
    let mut full = vec![0u32];
    full.extend_from_slice(sorted);
    full.push(end);
    let mut order = Vec::new();
    let mut queue = std::collections::VecDeque::from([(0usize, full.len() - 1)]);
    while let Some((lo, hi)) = queue.pop_front() {
        if hi - lo < 2 {
            continue;
        }
        let target = ((full[lo].max(1) as f64) * (full[hi] as f64)).sqrt();
        let mid = (lo + 1..hi)
            .min_by(|&a, &b| ((full[a] as f64).ln() - target.ln()).abs().total_cmp(&((full[b] as f64).ln() - target.ln()).abs()))
            .expect("interval has an inside");
        order.push(full[mid]);
        queue.push_back((lo, mid));
        queue.push_back((mid, hi));
    }
    order
}

/// The floor 1 configuration for vectors of `half` values.
pub(crate) fn floor(half: usize, codebooks: usize) -> Floor1 {
    let rangebits = half.trailing_zeros() as u8;
    let sorted = posts(half);
    let mut x_list = vec![0, half as u32];
    x_list.extend(bisection_order(&sorted, half as u32));
    let class = Floor1Class {
        dimensions: 3,
        subclasses: 2,
        masterbook: BOOK_FLOOR_MASTER as u8,
        subclass_books: vec![-1, BOOK_FLOOR_SUB[0] as i16, BOOK_FLOOR_SUB[1] as i16, BOOK_FLOOR_SUB[2] as i16],
    };
    Floor1::new(vec![0; sorted.len() / 3], vec![class], 2, rangebits, x_list, codebooks).expect("designed floor")
}

/// The residue 2 configuration for `channels` vectors of `half` values,
/// coding the interleaved vector up to `end`.
pub(crate) fn residue(end: usize) -> Residue {
    let mut cascade = Vec::new();
    let mut books = Vec::new();
    for (_, passes) in CLASSES {
        let mut row = [-1i16; 8];
        let mut bits = 0u8;
        for (pass, &b) in passes.iter().enumerate() {
            row[pass] = (BOOK_VQ0 + b) as i16;
            bits |= 1 << pass;
        }
        cascade.push(bits);
        books.push(row);
    }
    Residue {
        residue_type: 2,
        begin: 0,
        end: end as u32,
        partition_size: PARTITION as u32,
        classifications: CLASSES.len() as u8,
        classbook: BOOK_CLASS as u8,
        cascade,
        books,
    }
}

/// Coupling steps for a channel count (Vorbis channel order, 4.3.9): the
/// front pair, the rear or side pairs.
pub(crate) fn coupling(channels: u8) -> Vec<(u8, u8)> {
    match channels {
        2 => vec![(0, 1)],
        3 => vec![(0, 2)],
        4 => vec![(0, 1), (2, 3)],
        5..=7 => vec![(0, 2), (3, 4)],
        8 => vec![(0, 2), (3, 4), (5, 6)],
        _ => Vec::new(),
    }
}

/// The LFE channel of a layout, if it has one.
pub(crate) fn lfe(channels: u8) -> Option<usize> {
    match channels {
        6 => Some(5),
        7 => Some(6),
        8 => Some(7),
        _ => None,
    }
}

/// The whole setup: floors 0 (short) and 1 (long), residues 0 and 1 with
/// the given coded ends, mappings 0 and 1, modes 0 (short) and 1 (long).
pub(crate) fn setup(channels: u8, blocksize: [usize; 2], residue_end: [usize; 2]) -> Setup {
    let codebooks = codebooks();
    let floors = vec![Floor::One(floor(blocksize[0] / 2, codebooks.len())), Floor::One(floor(blocksize[1] / 2, codebooks.len()))];
    let residues = vec![residue(residue_end[0]), residue(residue_end[1])];
    let mapping = |i: u8| Mapping { coupling: coupling(channels), mux: vec![0; channels as usize], submap_floor: vec![i], submap_residue: vec![i] };
    Setup {
        codebooks,
        floors,
        residues,
        mappings: vec![mapping(0), mapping(1)],
        modes: vec![Mode { blockflag: false, mapping: 0 }, Mode { blockflag: true, mapping: 1 }],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn designed_lengths_make_complete_trees() {
        for w in [vec![1.0, 1.0], vec![0.5, 0.25, 0.125, 0.125], (0..300).map(|i| (-(i as f64) / 3.0).exp()).collect()] {
            let l = huffman_lengths(&w, 24);
            let kraft: f64 = l.iter().map(|&x| 0.5f64.powi(x as i32)).sum();
            assert!((kraft - 1.0).abs() < 1e-12);
            assert!(l.iter().all(|&x| (1..=24).contains(&x)));
        }
        // Steep weights are flattened to the length limit.
        let w: Vec<f64> = (0..64).map(|i| 0.5f64.powi(i)).collect();
        assert!(huffman_lengths(&w, 20).iter().all(|&x| x <= 20));
    }

    #[test]
    fn every_designed_book_is_valid() {
        let books = codebooks();
        assert_eq!(books.len(), BOOK_VQ0 + RESIDUE_BOOKS.len());
        for (l, b) in RESIDUE_BOOKS.iter().zip(&books[BOOK_VQ0..]) {
            // Lattice entry e has digits e % L, e / L ...
            let v = b.vector(1);
            assert_eq!(v[0], (l.min + l.step) as f32);
        }
        let s = setup(2, [256, 2048], [256, 2048]);
        for r in &s.residues {
            r.validate(&s.codebooks).unwrap();
        }
    }

    #[test]
    fn floor_posts_are_ordered_coarse_to_fine() {
        let f = floor(1024, 12);
        assert_eq!(&f.x_list[..3], &[0, 1024, 32]);
        assert_eq!((f.x_list.len() - 2) % 3, 0);
        let f = floor(128, 12);
        assert!(f.x_list.len() > 10);
        assert_eq!(floor(32, 12).x_list[1], 32);
    }

    #[test]
    fn class_ranges_are_covered_by_their_books() {
        for (max, passes) in CLASSES.iter().skip(1) {
            let reach: i32 = passes.iter().map(|&b| {
                let l = &RESIDUE_BOOKS[b];
                (l.min + (l.values as i32 - 1) * l.step).abs()
            }).sum();
            assert!(reach >= *max, "class reaching {max}: {reach}");
        }
    }
}
